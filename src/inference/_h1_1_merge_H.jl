function merging_h!(count_matrices, H_w_enriched_keys, enriched_keys_vec, bg)
    uniq_lens, length_indicator = 
        obtain_uniq_len_and_len_indicator(count_matrices)
        
    @inbounds for i = 1:length(uniq_lens)
        cmat_w_len_i, cmat_w_len_i_c, num_of_cmats_this_len, unmerged = 
            obtain_cmats_and_its_info(count_matrices, length_indicator, i)

        enriched_keys_here  = enriched_keys_vec[view(length_indicator,:, i)]

        for j = 1:num_of_cmats_this_len, k=j+1:num_of_cmats_this_len
            if unmerged[j] && unmerged[k]
                allr    = avg_allr(cmat_w_len_i[j], cmat_w_len_i[k], bg)   > allr_thresh
                allr_c  = avg_allr(cmat_w_len_i[j], cmat_w_len_i_c[k], bg) > allr_thresh
                if allr || allr_c
                    unmerged[k] = false 
                    if merge_add_counts
                        if allr > allr_c                          
                            append!(H_w_enriched_keys[enriched_keys_here[j]], 
                                    H_w_enriched_keys[enriched_keys_here[k]])
                        else
                            append!(H_w_enriched_keys[enriched_keys_here[j]], 
                                    values_comp.(H_w_enriched_keys[enriched_keys_here[k]]))
                        end
                    end
                    H_w_enriched_keys[enriched_keys_here[k]] = Vector{value_type}() # make empty
                end
            end
        end
    end
    merged_into_keys    = findall(x->length(x)>0, H_w_enriched_keys)
    enriched_keys_vec   = collect(merged_into_keys)
    H_w_enriched_keys   = getindices(H_w_enriched_keys, merged_into_keys)
    return H_w_enriched_keys, enriched_keys_vec
end

function merge_H(data, H_w_enriched_keys, hp, bg)
    count_matrices      = obtain_count_matrices(data, H_w_enriched_keys)
    enriched_keys_vec   = collect(keys(H_w_enriched_keys))

    len_cmats = count_matrices |> length;
    H_w_enriched_keys, enriched_keys_vec = merging_h!(count_matrices, H_w_enriched_keys, enriched_keys_vec, bg)
    len_merged_cmats = keys(H_w_enriched_keys) |> length;


    while len_cmats > len_merged_cmats
        println("Removed $(len_cmats-len_merged_cmats) redundant count matrices; 
            left with $(len_merged_cmats) count matrices.")
        count_matrices = obtain_count_matrices(data, H_w_enriched_keys)
        len_cmats = count_matrices |> length;
        H_w_enriched_keys, enriched_keys_vec = merging_h!(count_matrices, H_w_enriched_keys, enriched_keys_vec, bg)
        len_merged_cmats = keys(H_w_enriched_keys) |> length;
    end
    return H_w_enriched_keys
end