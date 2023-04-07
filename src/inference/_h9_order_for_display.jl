const ind2dna_str = Dict{Int, String}(1=>"A", 2=>"C", 3=>"G", 4=>"T")

function get_relaxed_consensus_str(pfm; any_regex=".", prob_thresh=0.5)
    argmax_inds = reshape(argmax(pfm, dims=1), (size(pfm,2),));
    char_array = [ind2dna_str[i[1]] for i in argmax_inds]
    char_array[findall((@view pfm[argmax_inds]) .< prob_thresh)] .= any_regex
    return join(char_array)
end

function unitranges_to_bitarray(uranges::Vector{UnitRange{Int}}, len::Int)
    bits = falses(len)
    for urange in uranges
        bits[urange] .= true
    end
    return bits
end

function apply_bitarray(s::AbstractString, b::BitArray)
    @assert length(s) == length(b)
    res = ""
    for i in eachindex(s)
        if b[i]
            res *= s[i]
        else
            res *= "-"
        end
    end
    return res
end

get_dashed_strings(consensus_str, str_bitarray) = apply_bitarray(consensus_str, str_bitarray)

function get_dashed_strings(ms)
    consensus_strs = get_relaxed_consensus_str.(ms.pfms)
    str_bitarrays = unitranges_to_bitarray.(ms.effective_segments, ms.lens)
    return get_dashed_strings.(consensus_strs, str_bitarrays)
end

function clear_islands_of_len!(char_array, str, sep, len)
    for i = 1:length(str)-(len+2)+1
        if char_array[i] == sep && char_array[i+len+1] == sep
            for j = 1:len
                char_array[i+j] = sep
            end            
        end
    end
end

function clean_str(str; sep='-', islands2clear=[2,3])
    char_array = Vector{Char}(undef, length(str))
    char_array[1] = str[1]; char_array[end] = str[end];
    for i = 2:length(str)-1
        if str[i] == sep
            char_array[i] = sep
        elseif str[i-1] == sep && str[i+1] == sep
            char_array[i] = sep
        else
            char_array[i] = str[i]
        end
    end
    # clear out all the islands of length 2
    for len in islands2clear
        clear_islands_of_len!(char_array, str, sep, len)
    end
    join(char_array)
end

# clear all the dashes in the beginning and end of the string given a vector of strings
function clear_dash(str; sep1 = '-', sep2 = '.')
    strlen = length(str)
    front_break_ind = 0
    for i = 1:strlen
        (str[i] != sep1 && str[i] != sep2) && break        
        front_break_ind = i
    end
    back_break_ind = strlen+1
    for i = strlen:-1:1
        (str[i] != sep1 && str[i] != sep2) && break        
        back_break_ind = i
    end
    front_break_ind > back_break_ind && return ""    
    return str[front_break_ind+1:back_break_ind-1]
end

get_cleaned_dashed_strs(ms) = clear_dash.(clean_str.(get_dashed_strings(ms)))

function edit_distance(s, t)
    s, t = (length(t) > length(s)) ? (s, t) : (t, s);
    slen = length(s);
    tlen = length(t);
    @inbounds while slen > 0 && s[slen] == t[tlen]
        slen -= 1; tlen -= 1;
    end
    start = 0;
    @inbounds if s[1]==t[1] || slen == 0
        while start < slen && s[start+1] == t[start+1]
            start += 1;
        end
        slen -= start;
        tlen -= start;
        slen == 0 && return tlen;
        # t = t[start+1:start+tlen];
    end

    v_0 = [i for i = 1:tlen]; # preallocate this later; needed for GPU kernel
    cur = 0
    @inbounds for i = 1:slen
        schar = s[start+i]; # check this
        left = cur = i-1;
        for j = 1:tlen
            abv = cur;
            cur = left;
            left = v_0[j]
            if schar != t[start+j]
                cur += 1;
                insdel = abv + 1;
                cur = (insdel < cur) ? insdel : cur;
                insdel = left + 1;
                cur = (insdel < cur) ? insdel : cur;
            end
            v_0[j] = cur;            
        end
    end
    return cur
end

function fill_pseudo_normalized_edit_similarity_matrix(dashed_strs)
    len_dashed_strs = length(dashed_strs)
    len_each_dashed_strs = length.(dashed_strs)
    str_similarity_mat = zeros(Float64, (len_dashed_strs, len_dashed_strs))
    for i = 1:len_dashed_strs
        for j = i+1:len_dashed_strs
            str_similarity_mat[i,j] = str_similarity_mat[j,i] = 
                min(len_each_dashed_strs[i], len_each_dashed_strs[j]) / edit_distance(dashed_strs[i], dashed_strs[j])
        end
    end
    return str_similarity_mat
end

mean_length_of_connected_components(trees) = mean.(trees)

# function obtain_groupings_for_display(ms; component_cut_off=2.75)
#     dashed_strs = get_cleaned_dashed_strs(ms)
#     str_similarity_mat = fill_pseudo_normalized_edit_similarity_matrix(dashed_strs)
#     trees = return_connected_components(str_similarity_mat, component_cut_off)
#     mean_lengths = mean_length_of_connected_components(trees)
#     # display the groups that have more than 1 member first
#     # within these groups, display the ones with the longest mean length first
#     groups_with_more_than_1_member = findall(length.(trees) .> 1)
#     groups_with_more_than_1_member_order = sortperm(mean_lengths[groups_with_more_than_1_member], rev=true)
#     groups_with_just_1_member = findall(length.(trees) .== 1)    
#     groups_with_just_1_member_order = sortperm(mean_lengths[groups_with_just_1_member], rev=true)
#     order_for_display = Vector{Int}(undef, ms.num_motifs)

#     k = 1
#     for i in groups_with_more_than_1_member_order
#         group_indices = 
#             trees[groups_with_more_than_1_member[i]]
#         len_group_indices = length(group_indices)
#         order_for_display[k:k+len_group_indices-1] = group_indices
#         k += len_group_indices
#     end
#     for i in groups_with_just_1_member_order
#        group_index = trees[groups_with_just_1_member[i]][1]
#        order_for_display[k] = group_index
#        k += 1
#     end
#     return order_for_display
# end  

function get_length_tree(dashed_strs_here, tree)
    length_tree = Vector{Vector{Int}}(undef, length(tree))
    for (i, l) in enumerate(tree)
        length_tree[i] = length.(dashed_strs_here[l])
    end
    length_tree
end

function get_order_for_this_group(dash_strs_here, component_cut_off)
    str_similarity_mat_abv = fill_pseudo_normalized_edit_similarity_matrix(dash_strs_here)
    trees_here = return_connected_components(str_similarity_mat_abv, component_cut_off)
    mean_lengths_here = mean_length_of_connected_components(get_length_tree(dash_strs_here, trees_here))
    groups_with_more_than_1_member = findall(length.(trees_here) .> 1)
    groups_with_more_than_1_member_order = sortperm(mean_lengths_here[groups_with_more_than_1_member], rev=true)
    groups_with_just_1_member = findall(length.(trees_here) .== 1)    
    groups_with_just_1_member_order = sortperm(mean_lengths_here[groups_with_just_1_member], rev=true)
    order_for_display_here = Vector{Int}(undef, length(dash_strs_here))

    k = 1
    for i in groups_with_more_than_1_member_order
        group_indices = 
        trees_here[groups_with_more_than_1_member[i]]
        len_group_indices = length(group_indices)
        order_for_display_here[k:k+len_group_indices-1] = group_indices
        k += len_group_indices
    end
    for i in groups_with_just_1_member_order
        group_index = trees_here[groups_with_just_1_member[i]][1]
        order_for_display_here[k] = group_index
        k += 1
    end
    return order_for_display_here
end

# function obtain_groupings_for_display(ms, component_cut_off=2.75)
#     dashed_strs = get_cleaned_dashed_strs(ms)
#     more_than_one_island_indicator = length.(ms.effective_segments) .> 1
#     one_island_indicator = .!more_than_one_island_indicator
#     # organize all the inds of the motifs with more than one island to be display first
#     inds_abv            = findall(more_than_one_island_indicator)
#     inds_below          = findall(one_island_indicator)
#     dashed_strs_abv     = @view dashed_strs[more_than_one_island_indicator]
#     dashed_strs_below   = @view dashed_strs[one_island_indicator]
#     # get the sort perm
#     inds_abv_sort_perm = get_order_for_this_group(dashed_strs_abv, component_cut_off) 
#     inds_below_sort_perm = get_order_for_this_group(dashed_strs_below, component_cut_off)
#     display_order = vcat(inds_abv[inds_abv_sort_perm], inds_below[inds_below_sort_perm])
#     return display_order
# end


function obtain_groupings_for_display1(ms; component_cut_off=2.75, big_pattern_thresh=80)
    dashed_strs = get_cleaned_dashed_strs(ms)
    big_pattern_indicator = length.(dashed_strs) .> big_pattern_thresh
    not_big_pattern_indicator = .!big_pattern_indicator
    more_than_one_island_indicator = (length.(ms.effective_segments) .> 1) .& not_big_pattern_indicator
    one_island_indicator = (.!more_than_one_island_indicator) .& not_big_pattern_indicator
    # organize all the inds of the motifs with more than one island to be display first
    inds_abv            = findall(big_pattern_indicator)
    inds_middle         = findall(more_than_one_island_indicator)
    inds_below          = findall(one_island_indicator)
    dashed_strs_abv     = @view dashed_strs[big_pattern_indicator]
    dashed_strs_middle  = @view dashed_strs[more_than_one_island_indicator]
    dashed_strs_below   = @view dashed_strs[one_island_indicator]
    # get the sort perm
    inds_abv_sort_perm = get_order_for_this_group(dashed_strs_abv, component_cut_off) 
    inds_middle_sort_perm = get_order_for_this_group(dashed_strs_middle, component_cut_off)
    inds_below_sort_perm = get_order_for_this_group(dashed_strs_below, component_cut_off)
    display_order = vcat(inds_abv[inds_abv_sort_perm], inds_middle[inds_middle_sort_perm], inds_below[inds_below_sort_perm])
    return display_order
end
