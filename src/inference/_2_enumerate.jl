create_key(c1,c2,c3,d12,d13,len) = (f1=c1, f2=c2, f3=c3, d12=d12, d13=d13, len=len)
create_value(seq_num, pos) = (seq_num=UInt32(seq_num), pos=pos, comp=false)
values_comp(v) = (seq_num=v.seq_num, pos=v.pos, comp=true)

function filter_code_components_using_median!(stored_code_components)
    mag_median = median(get_magnitude.(stored_code_components))
    filter!(x -> x.mag > mag_median, stored_code_components)
end

function filter_code_components_using_quantile!(stored_code_components; p=code_percentile)
    mag_percentile = quantile(get_magnitude.(stored_code_components), p)
    filter(x -> x.mag > mag_percentile, stored_code_components)
end

function filter_code_components(stored_code_components; 
                                filter_fcn=filter_code_components_using_quantile!, 
                                quantile_v=code_percentile)
    if filter_fcn != filter_code_components_using_quantile!
        return filter_fcn(stored_code_components)
    else
        return filter_fcn(stored_code_components; p=quantile_v)
    end
end

function get_scanning_range_of_filtered_code_components(stored_code_components) 
    cur_seq = 1; cur_range_start = 1; ranges = UnitRange{unit_range_int_t}[];
    @inbounds for i in eachindex(stored_code_components)
        if stored_code_components[i].seq != cur_seq
            push!(ranges, cur_range_start:i-1)
            cur_range_start = i
            cur_seq += 1
        end
    end
    return ranges
end

function insert_H!(H::Dictionary{composition_key_type, Vector{value_type}}, 
                   cartesian_inds_sorted, i::Int, j::Int, k::Int, ind::Int, h)
    c1,c2,c3 = get_fils(cartesian_inds_sorted, i,j,k)
    c1_pos = get_f1_pos(cartesian_inds_sorted, i)
    d12,d13 = get_d12_d13(cartesian_inds_sorted, i,j,k)
    len = d13 + h
    key = create_key(c1,c2,c3,d12,d13,len)
    value = create_value(ind, c1_pos)
    haskey(H, key) ? push!(H[key], value) : insert!(H, key, [value])
end

sort_by_x1(x) = x[1]

function enumerate_triplets(stored_code_component_filtered, seq_ranges, hp)
    H = Dictionary{composition_key_type, Vector{value_type}}();
    @inbounds for ind in eachindex(seq_ranges)
        whats_in_store = @view stored_code_component_filtered[seq_ranges[ind]]
        cartesian_inds_sorted = sort(whats_in_store, by=sort_by_x1)
        car_len = length(cartesian_inds_sorted)
        for i = 1:car_len-2
            for j = i+1:car_len-1
                for k = j+1:car_len
                    insert_H!(H, cartesian_inds_sorted, i, j, k, ind, hp.h)
                end
            end
        end
    end
    return H
end
