#=
Couple things:
1. Get code components
2. Filter code components if necessary: 
    a. filter out low magnitude components by setting the threshold 
       to the median of the magnitude or the v percentile of the magnitude
3. Get the scanning range of the filtered code components
4. Enumerate the triplet (word combinations) of the filtered code components
    using both the scanning range and the filtered code components
5. Get the start and end range for each syntax filter
=#

get_code_component_this_batch(cartesian_ind, magnitude, i) = 
    (position=cartesian_ind[1], fil=cartesian_ind[3], seq=cartesian_ind[4]+i-1, mag=magnitude)

get_magnitude(x) = x.mag

get_fils(cartesian_inds_sorted, i,j,k) = 
    cartesian_inds_sorted[i][2], cartesian_inds_sorted[j][2], cartesian_inds_sorted[k][2]

get_d12_d13(cartesian_inds_sorted, i,j,k) = 
    cartesian_inds_sorted[j][1] - cartesian_inds_sorted[i][1], 
    cartesian_inds_sorted[k][1] - cartesian_inds_sorted[i][1]

get_f1_pos(cartesian_inds_sorted, i) = cartesian_inds_sorted[i][1]

function append_code_component!(X, stored_code_components, i)
    cartesian_inds = findall(X .> 0);
    append!(stored_code_components, 
        get_code_component_this_batch.(cartesian_inds |> cpu, float_type_retrieval.(X[cartesian_inds]) |> cpu, i))
end

function code_retrieval(data, cdl, hp, len, projs)
    lambda_sparsity_warmup, lambda_sparsity, _,
    lambda_stepsize_warmup, omega_stepsize_warmup, lambda_stepsize, omega_stepsize, _,
    penalty_xyz, _, D, F_orig = prep_params(cdl, hp, projs)

    data_load = Flux.DataLoader(data.data_matrix, batchsize=hp.batch_size, shuffle=false, partial=false)
    stored_code_components = stored_code_component_t[]

    i = 1;
    @inbounds for S in data_load
        S = S |> gpu;
        _, _, X = ADMM_XYZ(S, D, F_orig,  
                    lambda_stepsize_warmup, lambda_stepsize, 
                    lambda_sparsity_warmup, lambda_sparsity,
                    omega_stepsize_warmup, omega_stepsize, 
                    penalty_xyz, 
                    hp, len, projs
                    )
        append_code_component!(X, stored_code_components, i)
        i += hp.batch_size;
    end

    return stored_code_components
end
