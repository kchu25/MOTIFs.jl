function obtain_filtered_stored_code_components(data, cdl, hp, len, projs;
                                                filter_fcn = filter_code_components_using_quantile!,
                                                quantile_v = code_percentile)
    @info "Retrieving the non-zero code components..."
    stored_code_component           = code_retrieval(data, cdl, hp, len, projs)
    @info "Filtering out the low-magnitude non-zero code components..."
    stored_code_component_filtered  = filter_code_components(stored_code_component; 
                                                    filter_fcn=filter_fcn, 
                                                    quantile_v=quantile_v)
    @info "Get the scanning ranges..."
    scanning_range = 
        get_scanning_range_of_filtered_code_components(stored_code_component_filtered)
    return stored_code_component_filtered, scanning_range
end

function obtain_enriched_word_combinations(stored_code_component_filtered, scanning_range, hp)
    H = enumerate_triplets(stored_code_component_filtered, scanning_range, hp)
    @info "Obtaining the enriched patterns..."
    enriched_keys = get_enriched_keys(H);
    H_w_enriched_keys = getindices(H, enriched_keys)
    return H_w_enriched_keys
end

function obtain_H_w_enriched_keys(data, cdl, hp, len, projs)
    stored_code_component_filtered, scanning_range = 
        obtain_filtered_stored_code_components(data, cdl, hp, len, projs)
    H_w_enriched_keys = 
        obtain_enriched_word_combinations(stored_code_component_filtered, scanning_range, hp)
    return H_w_enriched_keys
end

function merge_trim_merge_H_w_enriched_keys(data, hp, H_w_enriched_keys, bg)
    H_w_enriched_keys_merged = merge_H(data, H_w_enriched_keys, hp, bg)
    H_w_enriched_keys_merged_and_trimmed = trim_H(data, H_w_enriched_keys_merged, bg)
    H_w_enriched_keys_merged_and_trimmed_and_merged = 
        merge_H(data, H_w_enriched_keys_merged_and_trimmed, hp, bg)
    return H_w_enriched_keys_merged_and_trimmed_and_merged
end

function greedy_align(H_w_enriched_keys_mtm, data, bg)
    try 
        ms = enriched_keys2motifs(H_w_enriched_keys_mtm, data, bg);
        ms = alignment_merge!(ms, data, bg);
        for i = 1:indep_run
            println("indep run $i")
            scan_w_gpu!(ms, data);
            scan_w_gpu!(ms, data; bg=true);
            filter_positions_scores_usecomp!(ms, data, bg);
            ms = filter_insignificant_motifs(ms, data, bg);
            expansions_ms!(ms, data, bg);
            ms = alignment_merge!(ms, data, bg);
            new_cmats = posdicts2countmats(ms, data.data_matrix);
            new_cmats = trim_cmats(new_cmats, bg);
            new_cmats = merge_count_matrices(new_cmats, bg);
            ms = countmats2motifs(new_cmats, bg);
        end

        for i = 1:2
            println("indep run $i")
            scan_w_gpu!(ms, data);
            scan_w_gpu!(ms, data; bg=true);
            filter_positions_scores_usecomp!(ms, data, bg);
            ms = filter_insignificant_motifs(ms, data, bg);
            ms = alignment_merge!(ms, data, bg);
            new_cmats = posdicts2countmats(ms, data.data_matrix);
            new_cmats = trim_cmats(new_cmats, bg);
            new_cmats = merge_count_matrices(new_cmats, bg);
            ms = countmats2motifs(new_cmats, bg);
        end
        return ms
    catch e
        if isa(e, MethodError)
            @info "caught a method error"
        end
    end
    return nothing
    # for j = 1:dep_run
    #     println("dep run $j")
    #     scan_w_gpu!(ms, data);
    #     scan_w_gpu!(ms, data; bg=true);
    #     # filter_positions_scores_usecomp!(ms, data);
    #     filter_positions_scores_usecomp!(ms, data);
    #     # ms = non_overlap_scan!(ms, data.N);
    #     ms = filter_insignificant_motifs(ms, data);
    #     expansions_ms!(ms, data);
    #     ms = alignment_merge!(ms, data);
    #     new_cmats = posdicts2countmats(ms, data.data_matrix);
    #     new_cmats = trim_cmats(new_cmats);
    #     new_cmats = merge_count_matrices(new_cmats);
    #     ms = countmats2motifs(new_cmats);
    # end
end

# function run_thru(cdl, data, hp, len, projs)
#     H_w_enriched_keys = obtain_H_w_enriched_keys(data, cdl, hp, len, projs)
#     H_w_enriched_keys_mtm = merge_trim_merge_H_w_enriched_keys(data, hp, H_w_enriched_keys)
#     ms = greedy_align(H_w_enriched_keys_mtm, data, hp)
#     return ms
# end

function get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.35, get_how_many=500)
    @info "Filtering out the low-magnitude non-zero code components... ($quantile)" 
    stored_code_component_filtered  = filter_code_components(stored_code_component; 
                                                    filter_fcn=filter_code_components_using_quantile!, 
                                                    quantile_v=quantile)
    @info "Get the scanning ranges..."
    scanning_range          = get_scanning_range_of_filtered_code_components(stored_code_component_filtered)
    H                       = enumerate_triplets(stored_code_component_filtered, scanning_range, hp)
    @info "Obtaining the enriched patterns..."
    enriched_keys           = get_enriched_keys(H; max_word_combinations=get_how_many);
    H_w_enriched_keys       = getindices(H, enriched_keys)
    H_w_enriched_keys_mtm   = merge_trim_merge_H_w_enriched_keys(data, hp, H_w_enriched_keys, this_bg)
    return H_w_enriched_keys_mtm
end


function get_H_base(stored_code_component, hp; base_thresh=0.01, enriched_atleast=3)
    stored_code_component_filtered  = filter_code_components(stored_code_component; 
                                                    filter_fcn=filter_code_components_using_quantile!, 
                                                    quantile_v=base_thresh)
    scanning_range          = get_scanning_range_of_filtered_code_components(stored_code_component_filtered)
    H                       = enumerate_triplets(stored_code_component_filtered, scanning_range, hp)
    enriched_keys = findall(x->length(x)> enriched_atleast, H)
    @info "Got H_base"
    return getindices(H, enriched_keys)
end

get_backup_enriched_keys(H_base, data, hp, this_bg) = merge_trim_merge_H_w_enriched_keys(data, hp, H_base, this_bg)


function run_thru(data, cdl, hp, len, projs, this_bg)

    @info "Retrieving the non-zero code components..."
    stored_code_component       = code_retrieval(data, cdl, hp, len, projs)
    H_base = get_H_base(stored_code_component, hp)

    H_w_enriched_keys_mtm_p25   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.25, get_how_many=1000)
    H_w_enriched_keys_mtm_p35   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.35, get_how_many=1000)
    H_w_enriched_keys_mtm_p45   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.45, get_how_many=1000)
    H_w_enriched_keys_mtm_p5    = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.5,  get_how_many=1000)
    H_w_enriched_keys_mtm_p65   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.65, get_how_many=1000)
    H_w_enriched_keys_mtm_p75   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.75, get_how_many=1000)

    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm_p75, H_w_enriched_keys_mtm_p65)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p5)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p45)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p35)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p25)

    for k in keys(H_w_enriched_keys_mtm)        
        if haskey(H_base, k)
            H_w_enriched_keys_mtm[k] = union(H_w_enriched_keys_mtm[k], H_base[k])
        else
            H_w_enriched_keys_mtm[k] = union(H_w_enriched_keys_mtm[k])
        end
    end
    
    length(H_w_enriched_keys_mtm) == 0 && (H_w_enriched_keys_mtm = get_backup_enriched_keys(H_base, data, hp, this_bg))
    length(H_w_enriched_keys_mtm) == 0 && return nothing

    ms = enriched_keys2motifs(H_w_enriched_keys_mtm, data, this_bg);
    expansions_ms!(ms, data, this_bg);
    ms = alignment_merge!(ms, data, this_bg);
    new_cmats = posdicts2countmats(ms, data.data_matrix);
    new_cmats = trim_cmats(new_cmats, this_bg);
    new_cmats = merge_count_matrices(new_cmats, this_bg);
    return countmats2motifs(new_cmats, this_bg);

end


function run_thru2(data, cdl, hp, len, projs, this_bg)

    @info "Retrieving the non-zero code components..."
    stored_code_component       = code_retrieval(data, cdl, hp, len, projs)
    H_base = get_H_base(stored_code_component, hp)

    H_w_enriched_keys_mtm_p25   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.25, get_how_many=1000)
    H_w_enriched_keys_mtm_p35   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.35, get_how_many=1000)
    H_w_enriched_keys_mtm_p45   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.45, get_how_many=1000)
    H_w_enriched_keys_mtm_p5    = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.5,  get_how_many=1000)
    H_w_enriched_keys_mtm_p65   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.65, get_how_many=1000)
    H_w_enriched_keys_mtm_p75   = get_H_w_enriched_keys_mtm_quantile(stored_code_component, data, this_bg, hp; quantile=0.75, get_how_many=1000)

    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm_p75, H_w_enriched_keys_mtm_p65)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p5)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p45)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p35)
    H_w_enriched_keys_mtm = merge(H_w_enriched_keys_mtm, H_w_enriched_keys_mtm_p25)

    for k in keys(H_w_enriched_keys_mtm)        
        if haskey(H_base, k)
            H_w_enriched_keys_mtm[k] = union(H_w_enriched_keys_mtm[k], H_base[k])
        else
            H_w_enriched_keys_mtm[k] = union(H_w_enriched_keys_mtm[k])
        end
    end
    
    length(H_w_enriched_keys_mtm) == 0 && (H_w_enriched_keys_mtm = get_backup_enriched_keys(H_base, data, hp, this_bg))
    length(H_w_enriched_keys_mtm) == 0 && return nothing

    ms = enriched_keys2motifs(H_w_enriched_keys_mtm, data, this_bg);
    # expansions_ms!(ms, data, this_bg);
    return ms
end