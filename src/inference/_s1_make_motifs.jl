mutable struct motifs{T <: Integer, S <: Real}
    cmats::Vector{Matrix{S}}
    pfms::Vector{Matrix{S}}
    pwms::Union{Nothing,Vector{Matrix{S}}}
    effective_segments::Vector{Vector{UnitRange{Int}}}
    max_effective_lens::Vector{T}
    max_scores::Union{Nothing, Vector{S}}
    min_scores::Union{Nothing, Vector{S}}
    score_thresh::Union{Nothing, Vector{S}} # of the effective length
    lens::Vector{T}
    num_motifs::T
    positions::Union{Nothing, Vector{Dict{T, Vector{T}}}}
    scores::Union{Nothing, Vector{Dict{T, Vector{S}}}}
    use_comp::Union{Nothing, Vector{Dict{T, Vector{Bool}}}}
    positions_bg::Union{Nothing, Vector{Dict{T, Vector{T}}}}
    scores_bg::Union{Nothing, Vector{Dict{T, Vector{S}}}}
    use_comp_bg::Union{Nothing, Vector{Dict{T, Vector{Bool}}}}
end

ic(pfm,pwm) = reshape(sum(pfm .* pwm, dims=1), size(pfm,2))

pad_bit_ic_vec(q) = [false; q; false]

function moving_average(A::AbstractArray, m::Int)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end

function get_high_ic_segments(this_ic_vec; m=mv_avg_window, ic_threshold=effective_pos_ic_thresh)
    mv_this_ic_vec = pad_bit_ic_vec(moving_average(this_ic_vec, m) .> ic_threshold)
    diff_this_ic_vec = diff(push!(mv_this_ic_vec, false)) # append a false to the end
    starts = findall(diff_this_ic_vec .> 0)
    ends = findall(diff_this_ic_vec .< 0)
    # println("starts: ", starts)
    # println("ends: ", ends)
    return [start_:end_-1 for (start_, end_) in zip(starts, ends)]
end

filter_zero_eff_segs(eff_segs::Vector{Vector{UnitRange{Int}}}) = length.(eff_segs) .> 0

function get_filter_vec(cmats)
    # println(cmats)
    eff_segs = get_high_ic_segments.(cmat2ic.(cmats))
    # println("this eff seg: ", eff_segs)
    filter_vec = filter_zero_eff_segs(eff_segs)
    return filter_vec
end

function return_effective_segments_and_filter_vec(cmats)
    eff_segs = get_high_ic_segments.(cmat2ic.(cmats))
    filter_vec = filter_zero_eff_segs(eff_segs)
    return eff_segs, filter_vec
end

max_eff_len(range_vec) = maximum(length.(range_vec))

countmat2pfm(count_matrix, 
             countmat_sum_col; 
             ps=float_type_retrieval(0.01)) = 
        float_type_retrieval.((count_matrix .+ ps) ./ (countmat_sum_col .+ (4*ps)));

countmat2pfm(count_matrix; ps=float_type_retrieval(0.01)) = 
    countmat2pfm(count_matrix, sum(count_matrix, dims=1); ps=ps)

freq2pwm(pfm, bg) = log2.(pfm ./ bg)

function return_pfm_pwm_lens_nummats(countmats, bg)
    pfms = countmat2pfm.(countmats)
    pwms = [freq2pwm(pfm, bg) for pfm in pfms]
    lens = size.(pwms,2)
    return pfms, pwms, lens, length(lens)
end

function return_ms_others(ms, filter_vec)
    max_scores      = isnothing(ms.max_scores) ? nothing : ms.max_scores[filter_vec]
    min_scores      = isnothing(ms.min_scores) ? nothing : ms.min_scores[filter_vec]
    score_thresh    = isnothing(ms.score_thresh) ? nothing : ms.score_thresh[filter_vec]
    positions       = isnothing(ms.positions) ? nothing : ms.positions[filter_vec]
    scores          = isnothing(ms.scores) ? nothing : ms.scores[filter_vec]
    use_comp        = isnothing(ms.use_comp) ? nothing : ms.use_comp[filter_vec]
    positions_bg    = isnothing(ms.positions_bg) ? nothing : ms.positions_bg[filter_vec]
    scores_bg       = isnothing(ms.scores_bg) ? nothing : ms.scores_bg[filter_vec]
    use_comp_bg     = isnothing(ms.use_comp_bg) ? nothing : ms.use_comp_bg[filter_vec]
    return max_scores, min_scores, score_thresh, positions, 
           scores, use_comp, positions_bg, scores_bg, use_comp_bg
end


function countmats2motifs(count_mats, bg) # pval_thresh = 0.00027
    count_mats = count_mats[get_filter_vec(count_mats)]
    pfms, pwms, lens, num_pfms = return_pfm_pwm_lens_nummats(count_mats, bg)
    eff_segs = get_high_ic_segments.(cmat2ic.(count_mats))
    max_eff_lens = max_eff_len.(eff_segs)

    return motifs{Int, float_type_retrieval}(
                    count_mats,
                    pfms,
                    pwms,
                    eff_segs,
                    max_eff_lens,
                    nothing,
                    nothing, 
                    nothing,
                    lens,
                    num_pfms,
                    nothing, 
                    nothing,
                    nothing,
                    nothing,
                    nothing,
                    nothing
                );
end

function countmats2motifs(count_mats, positions, use_comp, bg) 
    filter_vec = get_filter_vec(count_mats)
    count_mats = count_mats[filter_vec]
    pfms, pwms, lens, num_pfms = return_pfm_pwm_lens_nummats(count_mats, bg)
    eff_segs = get_high_ic_segments.(cmat2ic.(count_mats))
    max_eff_lens = max_eff_len.(eff_segs)

    return motifs{Int, float_type_retrieval}( 
                    count_mats,
                    pfms,
                    pwms,
                    eff_segs,
                    max_eff_lens,
                    nothing,
                    nothing,
                    nothing, 
                    lens,
                    num_pfms,
                    positions[filter_vec], 
                    nothing,
                    use_comp[filter_vec],
                    nothing,
                    nothing,
                    nothing
                );
end

function filter_motifs_w_filter_vec(ms::motifs, filter_vec, bg)
    cmats = ms.cmats[filter_vec]

    eff_segs = get_high_ic_segments.(cmat2ic.(cmats))
    max_eff_lens = max_eff_len.(eff_segs)

    pfms, pwms, lens, num_pfms = return_pfm_pwm_lens_nummats(cmats, bg)

    max_scores, min_scores, score_thresh,  
    positions, scores, use_comp, 
    positions_bg, scores_bg, use_comp_bg = return_ms_others(ms, filter_vec)

    return motifs{Int, float_type_retrieval}( 
                    cmats,
                    pfms,
                    pwms,
                    eff_segs,
                    max_eff_lens,
                    max_scores,
                    min_scores,
                    score_thresh,
                    lens, 
                    num_pfms,
                    positions,
                    scores,
                    use_comp,
                    positions_bg,
                    scores_bg,
                    use_comp_bg
                );
end

function motifs_prep(ms::motifs{T,S}) where {T<:Integer,S<:Real}
    positions = [Dict{T,Vector{T}}()    for _ = 1:ms.num_motifs];
    scores    = [Dict{T,Vector{S}}()    for _ = 1:ms.num_motifs];
    use_comp  = [Dict{T,Vector{Bool}}() for _ = 1:ms.num_motifs];
    return positions, scores, use_comp
end

function push_position!(positions_i, seq_num, pos)
    if haskey(positions_i, seq_num)
        push!(positions_i[seq_num], pos)
    else
        positions_i[seq_num] = [pos]
    end
end

function push_use_comp!(use_comp_i, seq_num, use_comp)
    if haskey(use_comp_i, seq_num)
        push!(use_comp_i[seq_num], use_comp)
    else
        use_comp_i[seq_num] = [use_comp]
    end
end

function obtain_count_mat_pos_use_comp(H, data)
    count_matrices_lengths = Array{Int}(undef, length(H));
    enriched_keys = keys(H);
    @inbounds for (ind,k) in enumerate(enriched_keys)
        count_matrices_lengths[ind] = k.len
    end

    count_matrices = [zeros(Float32, (4, count_matrices_lengths[i])) for i in axes(count_matrices_lengths, 1)];
    
    positions = [Dict{Int, Vector{Int}}()  for _ in eachindex(count_matrices)]
    use_comp  = [Dict{Int, Vector{Bool}}() for _ in eachindex(count_matrices)]

    for (ind, k) in enumerate(keys(H))
        for v in H[k]
            push_position!(positions[ind], v.seq_num, v.pos)
            push_use_comp!(use_comp[ind], v.seq_num, v.comp)
            pos_start   = four_based_ind(v.pos)
            pos_end     = four_based_ind(v.pos+count_matrices_lengths[ind]-1)+3
            onehot_code = reshape(data.data_matrix[pos_start:pos_end,1,v.seq_num], 
                            (4, count_matrices_lengths[ind])) 
            count_matrices[ind] .+= v.comp ? reverse(onehot_code) : onehot_code
        end
    end

    return count_matrices, positions, use_comp
end

function enriched_keys2motifs(H, data, bg)
    count_matrices_lengths = Array{Int}(undef, length(H));
    enriched_keys = keys(H);
    @inbounds for (ind,k) in enumerate(enriched_keys)
        count_matrices_lengths[ind] = k.len
    end

    count_matrices = [zeros(Float32, (4, count_matrices_lengths[i])) for i in axes(count_matrices_lengths, 1)];
    
    positions = [Dict{Int, Vector{Int}}()  for _ in eachindex(count_matrices)]
    use_comp  = [Dict{Int, Vector{Bool}}() for _ in eachindex(count_matrices)]

    for (ind, k) in enumerate(keys(H))
        for v in H[k]
            push_position!(positions[ind], v.seq_num, v.pos)
            push_use_comp!(use_comp[ind], v.seq_num, v.comp)
            pos_start   = four_based_ind(v.pos)
            pos_end     = four_based_ind(v.pos+count_matrices_lengths[ind]-1)+3
            onehot_code = reshape(data.data_matrix[pos_start:pos_end,1,v.seq_num], 
                            (4, count_matrices_lengths[ind])) 
            count_matrices[ind] .+= v.comp ? reverse(onehot_code) : onehot_code
        end
    end
    return countmats2motifs(count_matrices, positions, use_comp, bg)
 end

