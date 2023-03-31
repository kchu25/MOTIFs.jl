function pvec_from_test_data(ms, data; no_olap=false)
    positions_test, scores_test, use_comp_test = gpu_scan(ms, data; bg=false, test=true)
    positions_test_bg, scores_test_bg, use_comp_bg_test = gpu_scan(ms, data; bg=true, test=true)
    filter_position_by_best_thresh!.(positions_test, scores_test, use_comp_test, ms.score_thresh)
    filter_position_by_best_thresh!.(positions_test_bg, scores_test_bg, use_comp_bg_test, ms.score_thresh)
    # make it non-overlap

    if no_olap 
        positions_test, scores_test, use_comp_test = 
            __non_overlap_scan__!(positions_test, scores_test, use_comp_test, ms.lens, data.N_test)
    end

    union_test_pos = get_union_ranges.(positions_test, ms.lens)
    union_test_pos_bg = get_union_ranges.(positions_test_bg, ms.lens)

    total_test_pos    = get_total_occupied_positions.(union_test_pos)
    total_test_pos_bg = get_total_occupied_positions.(union_test_pos_bg)
    
    pvec = fisher_pvec(total_test_pos, total_test_pos_bg, data; test=true)
    uniq_active_counts_test = active_counts_position(get_uniq_pos.(positions_test))
    return pvec, uniq_active_counts_test
end

function modified_sort_perm(pvec, ms, alpha_fisher)
    #=
        1. partition the motifs into significant and non-significant
        2. sort the significant motifs by length and pvalue
        3. sort the non-significant motifs by length and pvalue
        4. combine the two sorted lists
    =#
    len_p = [(l,p) for (l,p) in zip(ms.lens, pvec)]
    len_p_significant_inds = findall(x->x[2] < alpha_fisher, len_p)
    len_p_non_significant_inds = findall(x->x[2] >= alpha_fisher, len_p)

    len_p_significant = len_p[len_p_significant_inds]
    len_p_non_significant = len_p[len_p_non_significant_inds]

    len_p_significant_sorted_inds = reverse!(sortperm(len_p_significant, by=x->(x[1], 1.0-x[2])))
    len_p_non_significant_sorted_inds = reverse!(sortperm(len_p_non_significant, by=x->(x[1], 1.0-x[2])))

    return vcat(len_p_significant_inds[len_p_significant_sorted_inds], len_p_non_significant_inds[len_p_non_significant_sorted_inds])
end



function get_pvec_and_related_info(ms, data, alpha_fisher; no_olap=false)
    pvec, uniq_active_counts_test = pvec_from_test_data(ms, data; no_olap=no_olap)
    use_vec = pvec .< alpha_fisher
    # sort_perm = sortperm(pvec); # sort it according to pvalues (small to big)
    sort_perm = modified_sort_perm(pvec, ms, alpha_fisher)
    p_vec_sort_perm, use_vec_sort_perm = pvec[sort_perm], use_vec[sort_perm]
    pvalues = get_rounded_pval.(p_vec_sort_perm, use_vec_sort_perm);
    return pvalues, sort_perm, uniq_active_counts_test
end
