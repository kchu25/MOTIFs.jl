# filter position by best threshold 

function get_hits(scores, thresh)
    hits = 0
    @inbounds for k in keys(scores) 
        hits += sum(scores[k] .> thresh)
    end
    return hits
end

function get_max_score(scores, scores_bg)
    max_score = -float_type_retrieval(Inf)

    @inbounds for k in keys(scores) 
        max_ik = maximum(scores[k]) 
        max_score = max_ik > max_score ? max_ik : max_score;
    end
    @inbounds for k in keys(scores_bg) 
        max_ik = maximum(scores_bg[k]) 
        max_score = max_ik > max_score ? max_ik : max_score;
    end
    return max_score
end

function get_min_score(scores, scores_bg)
    min_score = float_type_retrieval(Inf)
    @inbounds for k in keys(scores) 
        min_ik = minimum(scores[k]) 
        min_score = min_ik < min_score ? min_ik : min_score;
    end
    @inbounds for k in keys(scores_bg) 
        min_ik = minimum(scores_bg[k]) 
        min_score = min_ik < min_score ? min_ik : min_score;
    end
    return min_score
end

function get_pvalue(pwm)
    size_pwm = size(pwm,2)
    9 < size_pwm ≤ 11 && (return pvalue_Touzet_mid)
    size_pwm ≤ 9 && (return pvalue_Touzet_small)
    return pvalue_Touzet_large
end

function get_high_scoring_pwm_segments(pwm; seg=max_pwm_length_Touzet2)
    pwm_len = size(pwm,2);
    pwm_len % seg 
    max_score_each = [maximum(view(pwm, :, i)) for i = 1:pwm_len]
    sort_perm_max_cols = sortperm(max_score_each, rev=true)
    pwm_seg_views = [view(sort_perm_max_cols, i:i+seg-1) for i = 1:seg:(pwm_len-pwm_len % seg)]
    return [pwm[:, pwm_seg_views[i]] for i in eachindex(pwm_seg_views)]
end

# function get_best_thresh(eff_pos, pwm)
#     best_thresh = 0f0
#     for pos in eff_pos
#         length(pos) ≤ 1 && (continue)
#         if length(pos) > max_pwm_length_Touzet2
#             pwm_segments = get_high_scoring_pwm_segments(pwm[:, pos])
#             for pwm_seg in pwm_segments
#                 best_thresh += pvalue2score(pwm_seg, get_pvalue(pwm_seg))
#             end
#         else
#             this_pwm = pwm[:, pos];
#             best_thresh += pvalue2score(this_pwm, get_pvalue(this_pwm))
#         end
#     end
#     return best_thresh
# end

# function filter_position_by_best_thresh!(positions, scores, use_comp, best_thresh)
#     if !isnothing(positions) && !isnothing(scores) && !isnothing(use_comp)
#         @inbounds for k in keys(positions)
#             mask = scores[k] .> best_thresh
#             positions[k] = positions[k][mask]
#             scores[k] = scores[k][mask]
#             use_comp[k] = use_comp[k][mask]        
#         end
#     end
# end

# function filter_positions_scores_usecomp!(ms)
#     ms.score_thresh = get_best_thresh.(ms.effective_segments, ms.pwms);
#     filter_position_by_best_thresh!.(ms.positions, ms.scores, ms.use_comp, ms.score_thresh);
#     filter_position_by_best_thresh!.(ms.positions_bg, ms.scores_bg, ms.use_comp_bg, ms.score_thresh);
# end



function get_best_thresh(scores, bg_scores, max_score, min_score, max_eff_len, eff_pos, pwm, asum, bg)
    if any(length.(eff_pos) .< max_pwm_length_Touzet2)
        best_thresh = 0f0
        for pos in eff_pos
            (length(pos) > max_pwm_length_Touzet2 || length(pos) ≤ 1) && continue
            this_pwm = pwm[:, pos];
            best_thresh += pvalue2score(this_pwm, get_pvalue(this_pwm); bg=bg)
        end
        return best_thresh
    end    
    best_thresh = min_score
    score_thresh = min_score
    best_p = 1f0
    while score_thresh < max_score
        a = get_hits(scores, score_thresh)
        b = get_hits(bg_scores, score_thresh)
        c = asum-a
        d = asum-b
        q = FisherExactTest(promote_i(a, c, b, d)...);
        p = HypothesisTests.pvalue(q, tail=:right)
        p < best_p && (best_p = p; best_thresh = score_thresh)
        score_thresh += score_thresh_increment
    end
    return best_thresh
end

function filter_position_by_best_thresh!(positions, scores, use_comp, best_thresh)
    if !isnothing(positions) && !isnothing(scores) && !isnothing(use_comp)
        @inbounds for k in keys(positions)
            mask = scores[k] .> best_thresh
            positions[k] = positions[k][mask]
            scores[k] = scores[k][mask]
            use_comp[k] = use_comp[k][mask]        
        end
    end
end

function filter_positions_scores_usecomp!(ms, data, bg)
    @info "number of motifs: $(ms.num_motifs)"    
    @info "filtering motifs positions..."
    ms.max_scores = get_max_score.(ms.scores, ms.scores_bg);
    ms.min_scores = get_min_score.(ms.scores, ms.scores_bg);
    ms.score_thresh = Vector{float_type_retrieval}(undef, ms.num_motifs)
    for i = 1:ms.num_motifs
        ms.score_thresh[i] = get_best_thresh(ms.scores[i], ms.scores_bg[i], ms.max_scores[i], ms.min_scores[i], ms.max_effective_lens[i], ms.effective_segments[i], ms.pwms[i], data.N*data.L, bg);
    end
    filter_position_by_best_thresh!.(ms.positions, ms.scores, ms.use_comp, ms.score_thresh);
    filter_position_by_best_thresh!.(ms.positions_bg, ms.scores_bg, ms.use_comp_bg, ms.score_thresh);
    @info "done filtering motifs positions..."
end