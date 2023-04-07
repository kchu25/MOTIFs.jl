function get_start_end(positions, len4, i, k, ind_k)
    start_ = (positions[i][k][ind_k]-1)*4+1;
    end_ = start_+len4-1;
    return start_, end_
end

function msa_add!(i, len4, positions, use_complement, 
    lens, msa, data_matrix; ps=0.01, return_count_mat=false)
    @inbounds for k in keys(positions[i])        
        for ind_k in 1:length(positions[i][k])
            # println("i: $i, k: $k, ind_k: $ind_k")
            start_, end_ = get_start_end(positions, len4, i, k, ind_k)
            if use_complement[i][k][ind_k]
                msa .+= 
                    reshape(submat_comlement(data_matrix,start_,end_,k, lens[i]),(4, lens[i]));
            else
                msa .+= reshape((@view data_matrix[start_:end_,1,k]), (4, lens[i]));
            end
        end
    end
    return return_count_mat ? msa .+ ps : countmat2pfm(msa .+ ps)
end

f_retrieval_t(x) = float_type_retrieval.(x)

function posdicts2countmats(ms::motifs, data_matrix::Array{S,3}) where {S<:Real}
    count_mats = Vector{Matrix{S}}(undef, ms.num_motifs);
    @inbounds for i in 1:ms.num_motifs
        # println(i)
        msa = zeros(S, (4, ms.lens[i]));
        count_mat = msa_add!(i, 4*ms.lens[i], ms.positions, 
            ms.use_comp, ms.lens, msa, data_matrix; return_count_mat=true)
        count_mats[i] = count_mat
    end

    return f_retrieval_t.(count_mats)
end

function posdicts2countmats(pos, use_comp, len, data_matrix::Array{S,3}) where {S<:Real}
    len4 = 4*len;
    count_mat = zeros(eltype(data_matrix), (4, len) );
    for k in keys(pos)
        for ind_k in 1:length(pos[k])
            start_  = (pos[k][ind_k]-1)*4+1;
            end_    = start_+len4-1;
            if use_comp[k][ind_k]
                count_mat .+= 
                    reshape(submat_comlement(data_matrix, start_, end_, k, len),(4, len));
            else
                count_mat .+= reshape((@view data_matrix[start_:end_, 1, k]), (4, len));
            end
        end
    end
    return Float16.(count_mat)
end