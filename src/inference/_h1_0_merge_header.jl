function get_matrices(countmat1, countmat2, pseudocount, bg)
    countmat1 = countmat1 .+ pseudocount
    countmat2 = countmat2 .+ pseudocount
    pfm1 = countmat1 ./ sum(countmat1, dims=1)
    pfm2 = countmat2 ./ sum(countmat2, dims=1)
    pwm1 = log2.(pfm1 ./ bg);
    pwm2 = log2.(pfm2 ./ bg);
    return countmat1, countmat2, pwm1, pwm2
end

function avg_allr(countmat1, countmat2, bg;
              pseudocount=float_type(0.01))
    countmat1, countmat2, pwm1, pwm2 = get_matrices(countmat1, countmat2, pseudocount, bg)
    allr_vec = sum(countmat2 .* pwm1 .+ countmat1 .* pwm2, dims=1) ./ 
            sum(countmat1 + countmat2, dims=1)
    return sum(allr_vec)/length(allr_vec)
end

function obtain_uniq_len_and_len_indicator(count_matrices)
    count_matrices_lengths  = size.(count_matrices, 2)
    uniq_lens               = count_matrices_lengths |> unique |> sort;
    length_indicator        = BitMatrix(undef, (length(count_matrices), length(uniq_lens)));
    @inbounds for (ind, uniq_len) in enumerate(uniq_lens)
        length_indicator[:, ind] = count_matrices_lengths .== uniq_len
    end
    return uniq_lens, length_indicator
end

function obtain_cmats_and_its_info(count_matrices, length_indicator, i)
    cmat_w_len_i                        = count_matrices[view(length_indicator,:, i)]
    cmat_w_len_i_c                      = reverse.(cmat_w_len_i)
    number_of_count_matrices_this_len   = length(cmat_w_len_i)
    unmerged        = fill(true, number_of_count_matrices_this_len)
    return cmat_w_len_i, cmat_w_len_i_c, number_of_count_matrices_this_len, unmerged
end
