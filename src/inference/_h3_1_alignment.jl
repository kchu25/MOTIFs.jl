# function make_effective_pos_mat(ms)
#     len_max = maximum(ms.lens)
#     effective_positions_mat = fill(false, (ms.num_motifs, len_max))
#     @inbounds for i = 1:ms.num_motifs
#         effective_positions_mat[i, 1:ms.lens[i]] = ms.effective_positions[i]
#     end
#     return cu(effective_positions_mat)
# end
    
const record_t = NTuple{3, UInt32}
const found_record_t = Vector{record_t};
const batch_size_greedy = 5000;

map_cartesian_n_increment(c, n)::record_t = c[1], c[2] + n -1, c[3]

############## greedy alignment ###############
# greedy search kernel
function greedy_search!(pwms, data_dat_gpu, lens, pos_scores)
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # kth pair
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y; # nth sequence
    l = (blockIdx().z - 1) * blockDim().z + threadIdx().z; # lth position

    L, N = size(data_dat_gpu); L_div_4 = CUDA.Int(L/4);
    K, _, _ = size(pwms);
    if k ≤ K && n ≤ N && l ≤ L_div_4-lens[k]+1
        @inbounds for (ind,i) in enumerate(l:l+lens[k]-1)
            # if eff_pos_mat[k,ind]
                for a = 1:4
                    pos_scores[k,n,l] += pwms[k,a,ind]*data_dat_gpu[(i-1)*4+a,n];                
                end
            # end
        end
        pos_scores[k,n,l] = pos_scores[k,n,l] > 0f0 ? pos_scores[k,n,l] : 0f0;
    end
    return nothing
end

function modify_w_found!(found_record, score_record, positions, scores, use_comp; rc=false)
    comp = rc ? true : false;
    @inbounds for (f, s) in zip(found_record, score_record)
        m, n, l = f[1], f[2], f[3];
        if haskey(positions[m], n)
            push!(positions[m][n], l)
            push!(scores[m][n], s)
            push!(use_comp[m][n],  comp)
        else
            positions[m][n] = [l];
            scores[m][n]    = [s];
            use_comp[m][n]  = [comp];
        end
    end
end

data_(data; test=false)  = test ? data.data_matrix_test : data.data_matrix
data_bg(data; test=false) = test ? data.data_matrix_bg_test : data.data_matrix_bg

function get_pos_scores_arr(ms, data; rc=false, bg=false, test=false)
    data_matrix = bg ? data_bg(data; test=test) : data_(data; test=test);
    found_record = found_record_t()
    score_record = float_type_retrieval[]
    
    # TODO fix this ugly hack
    length(size(data_matrix)) == 2 && (data_matrix = reshape(data_matrix, (size(data_matrix,1), 1, size(data_matrix,2))))
    L, _, N = size(data_matrix)

    maxlen = maximum(ms.lens);
    pwms = zeros(float_type_retrieval, ms.num_motifs, 4, maxlen);
    @inbounds for i = 1:ms.num_motifs pwms[i,:,1:ms.lens[i]] = 
        rc ? reverse(ms.pwms[i]) : ms.pwms[i]; end
        
    for n = 1:batch_size_greedy:N
        nend = min(n+batch_size_greedy-1, N)
        this_batch_size = nend-n+1
        data_matrix_gpu = reshape(cu(float_type_retrieval.(data_matrix[:,1,n:nend])), (L, this_batch_size));
        pos_scores = CUDA.zeros(float_type_retrieval, ms.num_motifs, this_batch_size, L);
        @cuda threads=ker_3d blocks=b_size_3d(pos_scores) greedy_search!(cu(pwms), 
                                                              data_matrix_gpu, 
                                                              cu(ms.lens), 
                                                              pos_scores
                                                              );
        pos_scores_arr = Array(pos_scores);
        found = findall(pos_scores_arr .> 0f0);
        append!(found_record, map_cartesian_n_increment.(found, n))
        append!(score_record, pos_scores_arr[found])
    end
    return found_record, score_record
end

function gpu_scan(ms, data; bg=false, test=false)
    found_record, score_record = 
        get_pos_scores_arr(ms, data; rc=false, bg=bg, test=test);
    found_record_rc, score_record_rc = 
        get_pos_scores_arr(ms, data; rc=true, bg=bg, test=test);
    
    positions, scores, use_comp = motifs_prep(ms);
    modify_w_found!(found_record, score_record, positions, scores, use_comp; rc=false)
    modify_w_found!(found_record_rc, score_record_rc, positions, scores, use_comp; rc=true)
    return positions, scores, use_comp
end

function scan_w_gpu!(ms, data; bg=false)
    positions, scores, use_comp = gpu_scan(ms, data; bg=bg)
    if bg 
        ms.positions_bg = positions;
        ms.scores_bg = scores;
        ms.use_comp_bg = use_comp;
    else
        ms.positions = positions;
        ms.scores = scores;
        ms.use_comp = use_comp;
    end
end

############## greedy alignment with non-overlap ###############

function max_score_ind(scores_n, scores_n_mask)
    max_score = -Inf;
    maxscore_ind = nothing;
    maxscore_ind_m = nothing;
    @inbounds for (ind_1,s) in enumerate(scores_n)
        for (ind_2, s_) in enumerate(s)
            if s_ > max_score && scores_n_mask[ind_1][ind_2]
                max_score = s_
                maxscore_ind = ind_2; 
                maxscore_ind_m = ind_1;
            end
        end
    end
    return maxscore_ind, maxscore_ind_m
end

max_pos_ind(positions_n, max_score_ind, max_score_m) = 
    positions_n[max_score_m][max_score_ind]

_intersect_(p, p_end, p_max, p_max_end) = 
    (p ≤ p_max ≤ p_end) || (p ≤ p_max_end ≤ p_end) || (p_max ≤ p ≤ p_max_end) || (p_max ≤ p_end ≤ p_max_end)


function __non_overlap_scan__!(positions, scores, use_comp, lens, N) 
    @inbounds for n = 1:N
        spans_pos = Int[]; 
        spans_len = Int[];

        indices_to_keep_n = [
            haskey(positions[i], n) ? fill(false, length(positions[i][n])) : Bool[]
                for i in eachindex(positions)];
        scores_n  = [haskey(scores[i], n) ? scores[i][n] : float_type_retrieval[] for i in eachindex(positions)];
        scores_n_mask  = [fill(true, length(s)) for s in scores_n]; # so entries in ms.scores aren't modified
        positions_n = [haskey(positions[i], n) ? positions[i][n] : Int[] for i in eachindex(positions)];   
        maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)

        if !isnothing(maxscore_ind)            
            while !isnothing(maxscore_ind)
                maxpos_ind = max_pos_ind(positions_n, maxscore_ind, maxscore_ind_m)
                intersect_ = false;
                for (p,l) in zip(spans_pos, spans_len)
                    # check whether this "segment" intersect with any previous segments
                    p_end     = p+l-1;
                    p_max     = positions_n[maxscore_ind_m][maxscore_ind];
                    p_max_end = p_max+lens[maxscore_ind_m]-1;
                    if _intersect_(p, p_end, p_max, p_max_end) 
                        intersect_ = true;
                        break
                    end
                end
                if !intersect_
                    indices_to_keep_n[maxscore_ind_m][maxscore_ind] = true;
                    push!(spans_pos, maxpos_ind);
                    push!(spans_len, lens[maxscore_ind_m]);
                end
                scores_n_mask[maxscore_ind_m][maxscore_ind]=false;
                maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)
            end
            for j in eachindex(positions)
                if haskey(positions[j], n)
                    positions[j][n] = positions[j][n][indices_to_keep_n[j]];
                    scores[j][n]    = scores[j][n][indices_to_keep_n[j]];
                    use_comp[j][n]  = use_comp[j][n][indices_to_keep_n[j]];
                    # so that the set of sequences covered by pwms are disjoint
                end
            end
        end
    end
    return positions, scores, use_comp
end
    
function non_overlap_scan!(ms::motifs{T,S}, N, bg) where {T <: Int, S <: Real}
    ms.positions, ms.scores, ms.use_comp = 
        __non_overlap_scan__!(ms.positions, ms.scores, ms.use_comp, ms.lens, N)
    keep = get_activate_counts(ms) .> 0;
    return filter_motifs_w_filter_vec(ms, keep, bg);
end

