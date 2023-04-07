#################################### obtain the enlarged matrix ####################################
function submat_comlement(data_matrix, start_, end_, k, len)
    piece = @view data_matrix[start_:end_,1,k];
    reverse(reshape(piece,(4,len)))
end

function enlarged_mat_add!(enlarged_mat, data_matrix, _start_, _end_, seq_ind, use_comp)
    if use_comp
        enlarged_mat .+= submat_comlement(data_matrix, four_based_ind(_start_), four_based_ind(_end_)+3, seq_ind, size(enlarged_mat,2))
    else
        enlarged_mat .+= reshape(view(data_matrix, four_based_ind(_start_):four_based_ind(_end_)+3,1,seq_ind), (4,size(enlarged_mat,2)))
    end
end

function enlarged_mat_increment!(ms, pfm_ind, seq_ind, ind_k, diff_front, diff_end, data_matrix, L, 
                                 enlarged_mat, enlarged_dict, enlarged_dict_use_comp; comp_match=false)
    use_comp = ms.use_comp[pfm_ind][seq_ind][ind_k];
    front_diff = use_comp ? diff_end : diff_front;
    end_diff   = diff_front + diff_end
    _start_    = ms.positions[pfm_ind][seq_ind][ind_k] - front_diff;
    _end_      = _start_ + ms.lens[pfm_ind] + end_diff - 1;

    if (1 ≤ _start_) && (_end_ ≤ L)
        enlarged_mat_add!(enlarged_mat, data_matrix, _start_, _end_, seq_ind, ms.use_comp[pfm_ind][seq_ind][ind_k])
        if !haskey(enlarged_dict, seq_ind)
            enlarged_dict[seq_ind] = [_start_]
            enlarged_dict_use_comp[seq_ind] = [ifelse(comp_match, !use_comp, use_comp)]
        else
            push!(enlarged_dict[seq_ind], _start_)
            push!(enlarged_dict_use_comp[seq_ind], comp_match)
        end        
    end
end

function obtain_enlarged_matrix(ms::motifs{T,S}, pfm_ind, diff_front, diff_end, data_matrix, L; comp_match=false) where {T,S}
    enlarged_length = ms.lens[pfm_ind] + diff_front + diff_end;
    enlarged_mat    = zeros(S, (4, enlarged_length))
    enlarged_dict = Dict{T,Vector{T}}(i=>[] for i in keys(ms.positions[pfm_ind]))
    enlarged_dict_use_comp = Dict{T,Vector{Bool}}(i=>[] for i in keys(ms.positions[pfm_ind]))
    @inbounds for seq_ind in keys(ms.positions[pfm_ind])
        for ind_k in 1:length(ms.positions[pfm_ind][seq_ind])
            enlarged_mat_increment!(ms, pfm_ind, seq_ind, ind_k, diff_front, diff_end, data_matrix, L, 
                enlarged_mat, enlarged_dict, enlarged_dict_use_comp; comp_match=comp_match)
        end
    end
    return ifelse(comp_match, reverse(enlarged_mat), enlarged_mat), enlarged_dict, enlarged_dict_use_comp
end

####################################################################################################


function this_cmat_pair_of_different_len_has_high_allr(cmat1, cmat2, i, j, bg; max_allowed_diff = max_allowed_diff, allr_thresh=allr_thresh)
    len_diff = size(cmat1,2)-size(cmat2,2)
    (abs(len_diff) > max_allowed_diff) && return nothing, nothing, nothing, nothing, nothing, false
    len_cond = len_diff < 0;
    smaller_mat = len_cond ? cmat1 : cmat2
    larger_mat  = len_cond ? cmat2 : cmat1
    smaller_ind = len_cond ? i : j
    larger_ind  = len_cond ? j : i
    smaller_mat_r = reverse(smaller_mat); # TODO later make it without a copy
    larger_mat_size = size(larger_mat,2)
    
    comp_match = false; diff_front = nothing; diff_end = nothing

    numcol_smaller_mat = size(smaller_mat,2);
    for i = 1:abs(len_diff)+1
        allr = avg_allr(smaller_mat, view(larger_mat, :,i:i+numcol_smaller_mat-1), bg)
        if allr > allr_thresh
            diff_front = i-1
            diff_end   = larger_mat_size - (i+numcol_smaller_mat-1)            
            break
        end
        allr_c = avg_allr(smaller_mat_r, view(larger_mat, :,i:i+numcol_smaller_mat-1), bg)
        if allr_c > allr_thresh
            comp_match = true
            diff_front = larger_mat_size - (i+numcol_smaller_mat-1)
            diff_end   = i-1
            break
        end
    end
    return diff_front, diff_end, smaller_ind, larger_ind, larger_mat,  comp_match
end

function alignment_merge!(ms, data, bg; allr_thresh=allr_thresh)
    count_matrix_each = posdicts2countmats(ms, data.data_matrix);
    ms.cmats = count_matrix_each;
    merged = fill(false, ms.num_motifs)
    for i = 1:ms.num_motifs
        for j = i+1:ms.num_motifs   
            merged[i] && continue
            merged[j] && continue
            diff_front, diff_end, smaller_ind, larger_ind, larger_mat, comp_match = 
                this_cmat_pair_of_different_len_has_high_allr(ms.cmats[i], ms.cmats[j], i, j, bg)
            if !isnothing(diff_front) && !isnothing(diff_end)

                enlarged_matrix, enlarged_dict, enlarged_dict_use_comp = 
                    obtain_enlarged_matrix(ms, smaller_ind, diff_front, diff_end, data.data_matrix, data.L; comp_match=comp_match)
                if avg_allr(enlarged_matrix, larger_mat, bg) > allr_thresh 
                    merge_add_counts && (ms.cmats[larger_ind] .+= enlarged_matrix)
                    # merge the positions
                    ms.positions[larger_ind] = mergewith(vcat, ms.positions[larger_ind], enlarged_dict)
                    # merge the use_comp
                    ms.use_comp[larger_ind] = mergewith(vcat, 
                                                    ms.use_comp[larger_ind], 
                                                    enlarged_dict_use_comp
                                                    )
                    merged[smaller_ind]  = true;          
                end               
            end            
        end
    end

    unmerged_cmats      = ms.cmats[.!merged];
    unmerged_positions  = ms.positions[.!merged];
    unmerged_use_comp   = ms.use_comp[.!merged];
    return countmats2motifs(unmerged_cmats, unmerged_positions, unmerged_use_comp, bg)
end
