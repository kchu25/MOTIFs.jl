function merging(count_matrices, bg)
    uniq_lens, length_indicator = 
        obtain_uniq_len_and_len_indicator(count_matrices)

    # merge the count matrices
    merged_cmats = Vector{Vector{Matrix{float_type_retrieval}}}();

    @inbounds for i = 1:length(uniq_lens)
        cmat_w_len_i, cmat_w_len_i_c, num_of_cmats_this_len, unmerged = 
            obtain_cmats_and_its_info(count_matrices, length_indicator, i)
        for j = 1:num_of_cmats_this_len, k=j+1:num_of_cmats_this_len
            if unmerged[j] && unmerged[k]
                allr    = avg_allr(cmat_w_len_i[j], cmat_w_len_i[k], bg)   > allr_thresh
                allr_c  = avg_allr(cmat_w_len_i[j], cmat_w_len_i_c[k], bg) > allr_thresh
                if allr || allr_c
                    unmerged[k] = false 
                    if merge_add_counts
                        if allr > allr_c  
                            cmat_w_len_i[j] .+= cmat_w_len_i[k]
                        else
                            cmat_w_len_i[j] .+= cmat_w_len_i_c[k]
                        end
                    end
                end
            end
        end
        push!(merged_cmats, cmat_w_len_i[unmerged])
    end
    return Iterators.flatten(merged_cmats) |> collect
end

function merge_count_matrices(count_matrices, bg)
    len_cmats = count_matrices |> length;
    merged_cmats = merging(count_matrices, bg)
    len_merged_cmats = merged_cmats |> length;

    while len_cmats > len_merged_cmats
        println("Removed $(len_cmats-len_merged_cmats) redundant count matrices; 
            left with $(len_merged_cmats) count matrices.")
        len_cmats = len_merged_cmats
        merged_cmats = merging(merged_cmats, bg)
        len_merged_cmats = merged_cmats |> length;
    end
    return merged_cmats
end