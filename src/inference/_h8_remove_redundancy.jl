############# take out and merge procedure #############


pfm2pwm(pfm, bg) = log2.(pfm ./ bg);

#=
Delete the motifs that are indicated by the indicator vector
=#
function delete_indicated_motif!(ms, indicator_vec)
    !isnothing(ms.cmats) && (ms.cmats = ms.cmats[.!indicator_vec])
    !isnothing(ms.pfms)  && (ms.pfms = ms.pfms[.!indicator_vec])
    !isnothing(ms.pwms)  && (ms.pwms = ms.pwms[.!indicator_vec])
    !isnothing(ms.effective_segments) && (ms.effective_segments = ms.effective_segments[.!indicator_vec])
    !isnothing(ms.max_effective_lens) && (ms.max_effective_lens = ms.max_effective_lens[.!indicator_vec])
    !isnothing(ms.max_scores) && (ms.max_scores = ms.max_scores[.!indicator_vec])
    !isnothing(ms.min_scores) && (ms.min_scores = ms.min_scores[.!indicator_vec])
    !isnothing(ms.score_thresh) && (ms.score_thresh = ms.score_thresh[.!indicator_vec])
    !isnothing(ms.lens) && (ms.lens = ms.lens[.!indicator_vec])
    !isnothing(ms.positions) && (ms.positions = ms.positions[.!indicator_vec])
    !isnothing(ms.scores) && (ms.scores = ms.scores[.!indicator_vec])
    !isnothing(ms.use_comp) && (ms.use_comp = ms.use_comp[.!indicator_vec])
    !isnothing(ms.positions_bg) && (ms.positions_bg = ms.positions_bg[.!indicator_vec])
    !isnothing(ms.scores_bg) && (ms.scores_bg = ms.scores_bg[.!indicator_vec])
    !isnothing(ms.use_comp_bg) && (ms.use_comp_bg = ms.use_comp_bg[.!indicator_vec])
    ms.num_motifs = length(ms.cmats)
end

function delete_by_indices!(ms, indices_vec)
    deleteat!(ms.cmats, indices_vec)
    deleteat!(ms.effective_segments, indices_vec)
    deleteat!(ms.positions, indices_vec)
    deleteat!(ms.use_comp, indices_vec)
    deleteat!(ms.lens, indices_vec)
    ms.num_motifs -= length(indices_vec)
end

#=
take out the sub-motifs that are within the range of the length tuple

two objectives:
# modify the ms object so that it no longer contains the sub-motifs 
  that are within the range of the length tuple

# return the sub-motifs that are within the range of the length tuple
=#
function take_out_sub_ms_by_range_indicator!(ms, bg, range_low_bdd, range_uppper_bdd)
    len_this_range =  range_low_bdd .≤ ms.lens .≤ range_uppper_bdd
    all(len_this_range .== false) && return
    this_len_range_cmats = ms.cmats[len_this_range]
    this_len_positions   = ms.positions[len_this_range]
    this_len_use_comp    = ms.use_comp[len_this_range]
    
    this_len_ms = countmats2motifs(this_len_range_cmats, this_len_positions, this_len_use_comp, bg)
    delete_indicated_motif!(ms, len_this_range)
    return this_len_ms
end

#=
merge the two ms
=#
function merge_ms(ms1, ms2, bg)
    cat_cmats = vcat(ms1.cmats, ms2.cmats)
    cat_positions = vcat(ms1.positions, ms2.positions)
    cat_use_comp = vcat(ms1.use_comp, ms2.use_comp)
    return countmats2motifs(cat_cmats, cat_positions, cat_use_comp, bg)
end

############# calculate jaccard similarity and BFS for connected components #############
# ms_this_range = take_out_sub_ms_by_range_indicator!(ms, bg, 30, 35)
# # if !isnothing(ms_this_range) ...
# olap_ratio = get_overlap_ratio(ms_this_range)

function return_connected_components(olap_ratio, jaccard_thresh=0.8)
    A = olap_ratio .> jaccard_thresh;
    l = size(A, 1);
    marked = fill(false, l);
    trees = Vector{Int}[];
    
    for i = 1:l
        start = i # if not discovered; add conditions later
        queue = Deque{Int}();        
        @inbounds if !marked[i]
            marked[i] = true;
            subtree = Int[];
            push!(subtree, start)
            pushfirst!(queue, start)
            while length(queue) != 0
                v = pop!(queue);
                for j in findall(A[v,:] .> 0)
                    if !marked[j]
                        marked[j] = true;
                        push!(subtree, j);
                        pushfirst!(queue, j);
                    end
                end
            end                        
            # if size(subtree)[1] > 1 
                push!(trees, subtree)
            # end
        end
    end
    return trees
end


function allr(p, q, p_count, q_count, bg)
    allr_score = Float64[];
    for i = 1:size(p,2)
        view_p_col = Base.view(p, :, i);
        view_q_col = Base.view(q, :, i);
        nb_p = p_count .* view_p_col;
        nb_q = q_count .* view_q_col;
        a1=sum(nb_p .* pfm2pwm(view_q_col, bg)); 
        a2=sum(nb_q .* pfm2pwm(view_p_col, bg));
        push!(allr_score, (a1+a2)/(sum(nb_p)+sum(nb_q)))
    end
    return sum(allr_score)
end

function convolve_allr(pfm_c2, pfm,
                       counts_pfm_c2,
                       counts_pfm,                        
                       len_pfm_c2,
                       len_pfm,
                       bg; len_diff_tol=4
                       )
    #= len_pfm_c2 will always be smaller since we've select the ones
        with minimal length
    =#
    min_col = len_pfm_c2 - len_diff_tol;
    allrs = Float64[];
    # start and end indices for pwms
    # s1e2 for pfm_c2, s2e2 for pfm
    s1e1s = UnitRange{Int}[];
    s2e2s = UnitRange{Int}[];
    l_dec_1 = Int[]; l_dec_2 = Int[]; 
    r_inc_1 = Int[]; r_inc_2 = Int[];

    for i = 1:(len_pfm_c2+len_pfm-1)
        s1 = max(1, len_pfm_c2-i+1); e1 = min(len_pfm_c2, len_pfm_c2-(i-len_pfm));
        s2 = max(1, i-len_pfm_c2+1); e2 = min(i, len_pfm);
        overlap_count = e1-s1+1;
        push!(s1e1s, s1:e1); push!(s2e2s, s2:e2);
        #=  
            Note that:
            1) no need to calculate if the number of columns of the 
            pfm is less than min_col as specified
            2) no need to calculate the placements for which
            maximal value of the score is below the threshold
        =#
        if overlap_count ≥ min_col
            push!(allrs, allr(Base.view(pfm_c2,:,s1:e1), Base.view(pfm,:,s2:e2), 
                            counts_pfm_c2, counts_pfm, bg));
        else
            push!(allrs, -Inf);
        end        
        push!(l_dec_1, max(s2-1,0)); push!(l_dec_2, max(s1-1,0));
        push!(r_inc_1, max(0,len_pfm-i)); push!(r_inc_2, max(i-e2,0));
    end
    argmax_ind = argmax(allrs);
    return allrs[argmax_ind], 
           l_dec_1[argmax_ind], 
           r_inc_1[argmax_ind], 
           l_dec_2[argmax_ind], 
           r_inc_2[argmax_ind]
end

reduce_merge_with_two_dict(d1, d2) = mergewith(vcat, d1, d2)

# modify_subtree_positions_lens!(subtree_pos, subtree_lens, subtree_use_comb, subtree_num_ms, reverse_comp,
# ld1_matches, ri1_matches, ld2_matches, ri2_matches, data.L)

function modify_subtree_positions_lens!(
                    subtree_pos, subtree_lens, subtree_use_comb, subtree_num_ms, reverse_comp, 
                    ld1_matches, ri1_matches, ld2_matches, ri2_matches, L)
    # modify the positions
    max_start_decrement_root   = maximum(ld1_matches)
    max_start_decrement_root_c = maximum(ri1_matches)
    merged_len = subtree_lens[1] + max_start_decrement_root + max_start_decrement_root_c
    mark2delete = [Dict{Int, Vector{Int}}() for _ = 1:subtree_num_ms] 
    # mark2delete = [Dict{Int, Vector{Int}}(k=>[] for k in keys(subtree_pos[i])) for i = 1:subtree_num_ms] 
    # mark2delete = [Dict{Int, BitVector}(k=>trues(length(subtree_pos[i][k])) for k in keys(subtree_pos[i])) for i = 1:subtree_num_ms] 

    println("reverse_comp: ", reverse_comp)

    for key in keys(subtree_pos[1])
        for (ind, pos) in enumerate(subtree_pos[1][key])
            modified_start = subtree_use_comb[1][key][ind] ? 
                pos - max_start_decrement_root_c : pos - max_start_decrement_root;
            if modified_start > 0 && modified_start + merged_len -1 <= L
                subtree_pos[1][key][ind] = modified_start
            else
                # mark2delete[1][key][ind] = false
                if haskey(mark2delete[1], key)
                    push!(mark2delete[1][key], ind)
                else
                    mark2delete[1][key] = [ind]
                end
            end       
        end
    end

    # modify the rest of the positions
    for i = 2:subtree_num_ms
        for key in keys(subtree_pos[i])
            for (ind, pos) in enumerate(subtree_pos[i][key])
                diff = ifelse(reverse_comp[i-1],
                    subtree_use_comb[i][key][ind] ? 
                        ld2_matches[i-1]+max_start_decrement_root-ld1_matches[i-1] :
                        ri2_matches[i-1]+max_start_decrement_root_c-ri1_matches[i-1]
                ,
                    subtree_use_comb[i][key][ind] ? 
                        ri2_matches[i-1]+max_start_decrement_root_c-ri1_matches[i-1] : 
                        ld2_matches[i-1]+max_start_decrement_root-ld1_matches[i-1]
                )
                subtree_use_comb[i][key][ind] = ifelse(reverse_comp[i-1],
                        !subtree_use_comb[i][key][ind],
                        subtree_use_comb[i][key][ind])
                modified_start = pos - diff
                if modified_start > 0 && modified_start + merged_len - 1  <= L
                    subtree_pos[i][key][ind] = modified_start
                else
                    # mark2delete[1][key][ind] = false

                    if haskey(mark2delete[i], key)
                        push!(mark2delete[i][key], ind)
                    else
                        mark2delete[i][key] = [ind]
                    end
                end       
            end
        end
    end

    # for i = 1:subtree_num_ms
    #     for key in keys(mark2delete[i])
    #         println("=====================================")    
    #         println("i: ", i)
    #         println("key: ", key)
    #         println("mark2delete[i][key]: ", mark2delete[i][key])
    #         println("subtree_pos[i][key]: ", subtree_pos[i][key])
    #         println("subtree_use_comb[i][key]: ", subtree_use_comb[i][key])
    #         println("=====================================")    
    #     end
    # end


    # if haskey(subtree_pos[2], 3821)
    #     println(subtree_pos[2][3821]);
    #     println(subtree_use_comb[2][3821]);
    #     println("=====================================")
    # end

    # println(mark2delete)


    for i = 1:subtree_num_ms
        println(i)
        for key in keys(mark2delete[i])
            keep_indices= filter(x->!(x in mark2delete[i][key]), eachindex(subtree_pos[i][key]))
            # println(subtree_pos[i][key])
            # println(keep_indices)
            subtree_pos[i][key] = subtree_pos[i][key][keep_indices]
            subtree_use_comb[i][key] = subtree_use_comb[i][key][keep_indices]
            # try 
            #     deleteat!(subtree_pos[i][key], mark2delete[i][key])
            #     deleteat!(subtree_use_comb[i][key], mark2delete[i][key])            
            # catch e
            #     if isa(e, BoundsError)
            #         println("i: ", i)
            #         println("key: ", key)
            #         println("mark2delete[i][key]: ", mark2delete[i][key])
            #         println("subtree_pos[i][key]: ", subtree_pos[i][key])
            #         println("subtree_use_comb[i][key]: ", subtree_use_comb[i][key])
            #     end
            # end
        end
        # println("i: ", i)
        # if haskey(subtree_pos[2], 3821)
        #     println(subtree_pos[2][3821]);
        #     println(subtree_use_comb[2][3821]);
        # end
    end

    merged_pos = reduce(reduce_merge_with_two_dict, subtree_pos)
    merged_use_comp = reduce(reduce_merge_with_two_dict, subtree_use_comb)

    return merged_pos, merged_use_comp, merged_len
end

function update_matches!(i, ld1, ri1, ld2, ri2, rc, 
    ld1_matches, ri1_matches, ld2_matches, ri2_matches, reverse_comp)
    ld1_matches[i-1] = ld1
    ri1_matches[i-1] = ri1
    ld2_matches[i-1] = ld2
    ri2_matches[i-1] = ri2
    reverse_comp[i-1] = rc
end

function process_subtree(subtree, ms_this_range, data, bg)
    println("subtree: $subtree")
    # sort subtree by length
    sort_inds        = subtree[sortperm(ms_this_range.lens[subtree])]
    println("sort_inds: $sort_inds")
    subtree_cmats    = ms_this_range.cmats[sort_inds]
    subtree_pfms     = ms_this_range.pfms[sort_inds]
    subtree_pos      = ms_this_range.positions[sort_inds]
    subtree_use_comb = ms_this_range.use_comp[sort_inds]
    subtree_lens     = ms_this_range.lens[sort_inds]
    subtree_counts = [sum((@view cmat[:,1])) for cmat in subtree_cmats]
    subtree_num_ms = length(subtree)
    ld1_matches = Vector{Int}(undef, subtree_num_ms-1); 
    ri1_matches = Vector{Int}(undef, subtree_num_ms-1); 
    ld2_matches = Vector{Int}(undef, subtree_num_ms-1); 
    ri2_matches = Vector{Int}(undef, subtree_num_ms-1); 
    reverse_comp = Vector{Bool}(undef, subtree_num_ms-1);

    for i in 2:subtree_num_ms
        allr_score, ld1, ri1, ld2, ri2 = convolve_allr(
            subtree_pfms[1],
            subtree_pfms[i],
            subtree_counts[1],
            subtree_counts[i],
            subtree_lens[1],
            subtree_lens[i],
            bg);
        allr_score_c, ld1_c, ri1_c, ld2_c, ri2_c = convolve_allr(
                subtree_pfms[1],
                reverse(subtree_pfms[i]),
                subtree_counts[1],
                subtree_counts[i],
                subtree_lens[1],
                subtree_lens[i],
                bg);
        if allr_score_c > allr_score
            update_matches!(i, ld1_c, ri1_c, ld2_c, ri2_c, true, 
                ld1_matches, ri1_matches, ld2_matches, ri2_matches, reverse_comp)
        else
            update_matches!(i, ld1, ri1, ld2, ri2, false, 
                ld1_matches, ri1_matches, ld2_matches, ri2_matches, reverse_comp)    
        end
    end

    merged_pos, merged_use_comp, merged_len = 
        modify_subtree_positions_lens!(subtree_pos, subtree_lens, subtree_use_comb, subtree_num_ms, reverse_comp,
                        ld1_matches, ri1_matches, ld2_matches, ri2_matches, data.L)
                        
    merged_cmat = posdicts2countmats(merged_pos, merged_use_comp, merged_len, data.data_matrix)
    return merged_cmat, merged_pos, merged_use_comp
end

function merge_to_remove_redundancy!(ms, data, bg; max_len_diff=5, olap_ratio_thresh=0.95)
    _ms_ = deepcopy(ms);
    while(true)
        total_num_ms = _ms_.num_motifs;
        min_len = minimum(_ms_.lens)
        max_len = maximum(_ms_.lens)
        for i = min_len:(max_len-max_len_diff+1)
            @info "merging motifs in range $i to $(i+max_len_diff-1)"
            ms_this_range = take_out_sub_ms_by_range_indicator!(_ms_, bg, i, i+max_len_diff-1);
            if !isnothing(ms_this_range)
                olap_ratio = get_overlap_ratio(ms_this_range)
                trees = return_connected_components(olap_ratio, olap_ratio_thresh)
                merged_cmats, merged_poses, merged_use_comps = 
                    Vector{eltype(ms_this_range.cmats)}(undef, length(trees)), 
                    Vector{eltype(ms_this_range.positions)}(undef, length(trees)),
                    Vector{eltype(ms_this_range.use_comp)}(undef, length(trees))
                for (j, subtree) in enumerate(trees)
                    merged_cmat, merged_pos, merged_use_comp = nothing, nothing, nothing
                    if length(subtree) == 1
                        merged_cmat, merged_pos, merged_use_comp = 
                            ms_this_range.cmats[j], ms_this_range.positions[j], ms_this_range.use_comp[j]
                    else
                        merged_cmat, merged_pos, merged_use_comp = process_subtree(subtree, ms_this_range, data, bg)
                    end
                    merged_cmats[j] = merged_cmat; merged_poses[j] = merged_pos; merged_use_comps[j] = merged_use_comp
                end
                ms_this_tree = countmats2motifs(merged_cmats, merged_poses, merged_use_comps, bg)
                _ms_ = merge_ms(_ms_, ms_this_tree, bg)
            end
        end
        @info "total number of motifs: $total_num_ms, after merging: $(_ms_.num_motifs)"
        total_num_ms == _ms_.num_motifs && break        
    end
    return _ms_
end