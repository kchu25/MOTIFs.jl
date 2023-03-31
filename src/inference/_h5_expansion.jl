function countvec2ic(count_vec, pos_ki; ps=0.01)
    count_vec = count_vec .+ ps;
    count_vec_sum = sum(count_vec);
    freq_vec = count_vec ./ count_vec_sum
    count_vec_sum/length(pos_ki), 2 + sum(freq_vec .* log2.(freq_vec))
end

get_char(data, p, char_ind, comp) =
        comp ?
            atcg2dummy[ atcg_comp[data.raw_data[p[2]][char_ind]] ] :
            atcg2dummy[ data.raw_data[p[2]][char_ind] ]

# get_char(data, p, char_ind, comp) =
#     comp ?
#         atcg2dummy[ atcg_comp[data.raw_data[p[2]].str[char_ind]] ] :
#         atcg2dummy[ data.raw_data[p[2]].str[char_ind] ]
            
count_vec_at_pos_ms(p, data, char_ind, comp) = 
    char_ind < 1 || char_ind > data.L ? 
        atcg2dummy['z'] :
        get_char(data, p, char_ind, comp)

get_expand_ind_left(p, Δ, comp)  = comp ? p[1][end]+Δ : p[1][1]-Δ;
get_expand_ind_right(p, Δ, comp) = comp ? p[1][1]-Δ   : p[1][end]+Δ;
get_shrink_ind_left(p, Δ, comp)  = comp ? p[1][end]-Δ : p[1][1]+Δ;
get_shrink_ind_right(p, Δ, comp) = comp ? p[1][1]+Δ   : p[1][end]-Δ;

function left_expansion_ms(data, pos_ki, dec; ic_expand_thresh=0.5, pc=1.0, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos_ms(p, data, get_expand_ind_left(p, dec, p[3]), p[3])
    end
    percentage_used, ic = countvec2ic(count_vec .+ pc, pos_ki)
    return percentage_used > tol && ic > ic_expand_thresh
end

function right_expansion_ms(data, pos_ki, inc; ic_expand_thresh=0.5, pc=1.0, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos_ms(p, data, get_expand_ind_right(p, inc, p[3]), p[3])
    end
    percentage_used, ic = countvec2ic(count_vec .+ pc, pos_ki)
    return percentage_used > tol && ic > ic_expand_thresh
end

function left_shrinkage_ms(data, pos_ki, inc; ic_shrink_thresh=0.2, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos_ms(p, data, get_shrink_ind_left(p, inc, p[3]), p[3])
    end
    percentage_used, ic = countvec2ic(count_vec, pos_ki)

    return percentage_used > tol && ic < ic_shrink_thresh
end

function right_shrinkage_ms(data, pos_ki, dec; ic_shrink_thresh=0.2, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos_ms(p, data, get_shrink_ind_right(p, dec, p[3]), p[3])
    end
    percentage_used, ic = countvec2ic(count_vec, pos_ki)
    return percentage_used > tol && ic < ic_shrink_thresh
end

function expansion_left_right_ms(data, pos_ki, ic_expand_thresh)
    expand_left = true; expand_right = true; 
    left_dec = 1; right_inc = 1;
    while expand_left 
        expand_left = 
            left_expansion_ms(data, pos_ki, left_dec; ic_expand_thresh=ic_expand_thresh)
        expand_left ? left_dec+=1 : left_dec-=1;
    end
    while expand_right
        expand_right = 
            right_expansion_ms(data, pos_ki, right_inc; ic_expand_thresh=ic_expand_thresh)
        expand_right ? right_inc+=1 : right_inc-=1;
    end
    return left_dec, right_inc
end

function trimming_left_right_ms(expansion::Tuple{Int, Int}, data, pos_ki, ic_shrink_thresh)
    # returns how much increment (decrement) from the left (from the right)
    shrink_left_  = expansion[1] == 0 ? true : false; 
    shrink_right_ = expansion[2] == 0 ? true : false;
    left_inc = shrink_left_ ? 1 : 0; 
    right_dec = shrink_right_ ? 1 : 0;
    while shrink_left_
        shrink_left_ = 
            left_shrinkage_ms(data, pos_ki, left_inc-1; ic_shrink_thresh=ic_shrink_thresh)
        shrink_left_ ? left_inc+=1 : left_inc -= 1;
    end
    while shrink_right_
        shrink_right_ = 
            right_shrinkage_ms(data, pos_ki, right_dec-1; ic_shrink_thresh=ic_shrink_thresh)
        shrink_right_ ? right_dec+=1 : right_dec-=1;
    end
    return left_inc, right_dec
end

function msa_expansion_ms(pos, data, ic_expand_thresh, ic_shrink_thresh)
    expansions = Vector{Tuple{Int,Int}}();
    shrinkage = Vector{Tuple{Int,Int}}();
    pos_lens = length.(pos)
    for (ind, pos_ki) in enumerate(pos)
        if pos_lens[ind] > 0
            push!(expansions, expansion_left_right_ms(data, pos_ki, ic_expand_thresh))
        else
            push!(expansions, (0,0))
        end
    end

    for (ind, expansion) in enumerate(expansions)
        if pos_lens[ind] > 0
            push!(shrinkage, trimming_left_right_ms(expansion, data, pos[ind], ic_shrink_thresh))
        else
            push!(shrinkage, (0,0))
        end
    end
    return expansions, shrinkage
end

function posdict2pos(ms)
    pos = [Vector{Tuple{UnitRange{Int64}, Int64, Bool}}() for _ = 1:ms.num_motifs]
    for i = 1:ms.num_motifs
        for k in keys(ms.positions[i])
            for (ind,w) in enumerate(ms.positions[i][k])
                push!(pos[i], (w:w+ms.lens[i]-1, k, ms.use_comp[i][k][ind]))
            end
        end
    end
    return pos
end

function edit_posdict_i!(expand__left, expand__right, 
                         shrink__left, shrink__right, ms, i, L)
    @assert expand__left == 0 || shrink__left == 0 "one left shink/expand has to be 0"
    @assert expand__right == 0 || shrink__right == 0 "one right shink/expand has to be 0"  

    left_Δ = -expand__left + shrink__left
    right_Δ = expand__right - shrink__right
    ms.lens[i] += -left_Δ + right_Δ    
    # @info "motif: $i, Left-Δ: $left_Δ, right-Δ: $right_Δ"
    for k in keys(ms.positions[i])
        indices_to_keep = fill(true, length(ms.positions[i][k]))
        for ind in 1:length(ms.positions[i][k])
            Δ = ifelse(ms.use_comp[i][k][ind], 
                ms.positions[i][k][ind] - right_Δ,
                ms.positions[i][k][ind] + left_Δ)
            if 1 ≤ Δ ≤ L && Δ + ms.lens[i] - 1 ≤ L
                ms.positions[i][k][ind] = Δ
            else
                indices_to_keep[ind] = false
            end
        end
        ms.positions[i][k] = ms.positions[i][k][indices_to_keep]
        ms.use_comp[i][k]  = ms.use_comp[i][k][indices_to_keep]
        # TODO scores as well?
    end    
end

function expansions_ms!(ms, data, bg; ic_expand_t=ic_expand_thresh, ic_shrink_t=ic_trim_thresh)
    use_next = fill(true, ms.num_motifs)
    ms_pos = posdict2pos(ms)
    expansions, shrinkage = 
        msa_expansion_ms(ms_pos, data, ic_expand_t, ic_shrink_t)
    for (ind, ((el,er),(sl,sr))) in enumerate(zip(expansions, shrinkage))
        if sl + sr > ms.lens[ind]
            use_next[ind] = false
        else
            edit_posdict_i!(el, er, sl, sr, ms, ind, data.L)
        end    
    end
    ms = countmats2motifs(ms.cmats[use_next], 
                            ms.positions[use_next], 
                            ms.use_comp[use_next], bg)
end