unique_positions(positions_array) = unique(positions_array)


#TODO speed this up
function get_uniq_pos(positions_i)
    positions_copy = typeof(positions_i)()
    for k in keys(positions_i)
        positions_copy[k] = unique_positions(positions_i[k])
    end
    return positions_copy
end

get_uniq_counts(ms) = 
    active_counts_position(get_uniq_pos.(ms.positions)), 
    active_counts_position(get_uniq_pos.(ms.positions_bg))

###########################################################################

function num_overlap(r1::UnitRange, r2::UnitRange)
    @assert r1[1] ≤ r1[end] && r2[1] ≤ r2[end] "range is not valid"
    @inbounds if r1[1] ≤ r2[1] ≤ r1[end]
        return min(r1[end],r2[end])-r2[1]+1;
    elseif  r2[1] ≤ r1[1] ≤ r2[end]
        return min(r1[end],r2[end])-r1[1]+1;
    else
        return 0;
    end
end

function total_active_position(position)
    activate_count = 0
    for key in keys(position)
        for r in position[key]
            activate_count += length(r)
        end          
    end
    return activate_count
end

function push_ranges!(ranges, _ranges_, i)
    if _ranges_[end][end] ≥ ranges[i][1]
        _ranges_[end] = _ranges_[end][1]:ranges[i][end]
    else
        push!(_ranges_, ranges[i])
    end
end

function union_ranges(ranges)
    length(ranges) == 0 && return eltype(ranges)[]
    ranges = sort(ranges, by = x->x[1])
    _ranges_ = [ranges[1]]
    @inbounds for i in eachindex(@view ranges[2:end])
        push_ranges!(ranges, _ranges_, i)
    end
    return _ranges_
end

function union_pos(positions_arr_k, len)
    ranges = [positions_arr_k[i]:positions_arr_k[i]+len-1 for i in eachindex(positions_arr_k)]
    return union_ranges(ranges)
end

function get_union_ranges(positions_i, len_i)
    positions_copy = Dict{Int, Vector{UnitRange{Int}}}(k=>[] for k in keys(positions_i));
    for k in keys(positions_copy)
        positions_copy[k] = union_pos(positions_i[k], len_i)
    end
    return positions_copy
end

function get_total_occupied_positions(position_ranges)
    total_positions = 0
    for k in keys(position_ranges)
        for r in position_ranges[k]
            total_positions += length(r)
        end
    end
    return total_positions
end


# function overlap_inc!(pos_i, pos_j, overlap_ij, mutual_keys)

# end

function get_overlap_ratio(ms)
    # @info "Computing union ranges..."
    # union_poses = get_union_ranges.(ms.positions, ms.lens)
    union_poses = Vector{Dict{Int64, Vector{UnitRange{Int64}}}}(undef, ms.num_motifs)
    @floop for i = 1:ms.num_motifs
        union_poses[i] = get_union_ranges(ms.positions[i], ms.lens[i])
    end
    
    # @info "Computing total active positions..."
    acs = total_active_position.(union_poses)
    olap_ratio_ij = zeros(Float32, (ms.num_motifs, ms.num_motifs))
    # @info "Computing overlap ratio..."
    @floop for i = 1:ms.num_motifs
        # println(i)
        @inbounds for j = i+1:ms.num_motifs
            pos_i = union_poses[i]
            pos_j = union_poses[j]
            overlap_ij = 0f0
            mutual_keys = intersect(keys(pos_i), keys(pos_j))
            for k in mutual_keys
                length(pos_i[k]) == 0 || length(pos_j[k]) == 0 && continue
                for pos_i_range in pos_i[k]
                    for pos_j_range in pos_j[k]
                        overlap_ij += num_overlap(pos_i_range, pos_j_range)
                    end
                end
            end
            olap_ratio_ij[i,j] = olap_ratio_ij[j,i] = overlap_ij / (acs[i]+acs[j]-overlap_ij)
        end
    end
    return olap_ratio_ij
end
