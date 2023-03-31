function active_counts_position(positions)
    activate_count = Vector{Float64}(undef, length(positions))
    @inbounds for i in eachindex(positions)
        sum_p = 0
        for key in keys(positions[i])
            sum_p += length(positions[i][key])
        end
        activate_count[i] = sum_p
    end
    return activate_count
end

function get_activate_counts(ms::motifs; bg = false)
    positions = bg ? ms.positions_bg : ms.positions;
    return active_counts_position(positions)
end

get_both_activate_counts(ms::motifs) =
    get_activate_counts(ms),  get_activate_counts(ms; bg=true)

function fisher_pvec(activate_counts, activate_counts_bg, data; test=false) 
    activate_sum = (test ? data.N_test : data.N) * data.L # total number of components that can be activated
    pvalues = fill(0.0, length(activate_counts));  

    for i in eachindex(activate_counts)
        a = activate_counts[i]; b = activate_counts_bg[i];
        c = activate_sum - a;   d = activate_sum - b;
        if a == 0 && b == 0
            pvalues[i] = 1.0 # if there's no activation in both just throw it away
        else
            q = FisherExactTest(promote_i(a, c, b, d)...);
            pvalues[i] = HypothesisTests.pvalue(q, tail=:right);
        end
    end
    return pvalues
end
   
function get_fisher_p_values(ms::motifs, data; test=false)
    union_positions, union_positions_bg = 
        get_union_ranges.(ms.positions, ms.lens), get_union_ranges.(ms.positions_bg, ms.lens)
    total_positions, total_positions_bg = 
        get_total_occupied_positions.(union_positions), get_total_occupied_positions.(union_positions_bg)
    return fisher_pvec(total_positions, total_positions_bg, data; test=test) 
end

function filter_insignificant_motifs(ms::motifs, data, this_bg; test=false)
    pvec = get_fisher_p_values(ms, data; test=test);
    keep = pvec .< pvalue_fisher_thresh
    return filter_motifs_w_filter_vec(ms, keep, this_bg)
end