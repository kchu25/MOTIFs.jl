##### From the word combinations H, obtain the count matrices from enriched keys ########

function get_words(H, enriched_keys, len_enriched_keys; max_word_combinations=num_pfms2process)
    q = getindices(H, enriched_keys);
    q_sorted = sort(length.(q), rev=true)
    if len_enriched_keys > max_word_combinations
        return Indices(collect(keys(q_sorted))[1:max_word_combinations])
    else
        return Indices(collect(keys(q_sorted))[1:len_enriched_keys])
    end
end

function get_enriched_keys(H; max_word_combinations = num_pfms2process, 
            dec        = -5,
            count_from = cover_more_than,
            count_to   = cover_at_least
            )
    enriched_keys = nothing
    for count = count_from:dec:count_to
        enriched_keys = findall(x->length(x)> count, H)
        len_enriched_keys = length(enriched_keys)
        len_enriched_keys > max_word_combinations && 
            (return get_words(H, enriched_keys, len_enriched_keys; max_word_combinations=num_pfms2process))
    end
    return enriched_keys
end

function obtain_count_matrices(data, H)
    enriched_keys = keys(H)
    count_matrices_lengths = Array{Int}(undef, length(enriched_keys));
    @inbounds for (ind,k) in enumerate(enriched_keys)
        count_matrices_lengths[ind] = k.len
    end

    count_matrices = [zeros(float_type, (4, count_matrices_lengths[i])) for i = 1:length(enriched_keys)];
    @inbounds for (ind,k) in enumerate(enriched_keys)
        for v in H[k]
            pos_start   = four_based_ind(v.pos)
            pos_end     = four_based_ind(v.pos+count_matrices_lengths[ind]-1)+3
            onehot_code = reshape((@view data.data_matrix[pos_start:pos_end,1,v.seq_num]), 
                            (4, count_matrices_lengths[ind])) 
            count_matrices[ind] .+= v.comp ? reverse(onehot_code) : onehot_code
        end
    end
    return count_matrices
end



