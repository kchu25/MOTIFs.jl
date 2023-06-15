function setup_num_epochs(number_training_samples)
    if number_training_samples < 1000
        return 25
    elseif number_training_samples < 10000
        return 10
    elseif number_training_samples < 100000
        return 5
    else
        return 3
    end
end

function train_ucdl(data; 
                    num_epochs=nothing,
                    # filter_len::Int = 8,
                    # f_len::Int = filter_len*4,
                    # M::Int = 45,
                    # twoM::Int = 2*M,
                    # h::Int = 12,
                    # K::Int = 25,
                    # q::Int = 20,
                    # batch_size::Int = 16,
                    # num_pass_xyz::Int = 6,
                    # num_pass_df::Int = 3,
                    # magnifying_factor::float_type = 10,
                    # gamma::float_type = 0.1
                    l1_loss_thresh=float_type(95.0)
    )
    hp          = Hyperparam();
    len         = length_info(hp, data);
    projs       = projectors(hp, len);
    cdl         = ucdl(hp);
    data_load   = Flux.DataLoader(data.data_matrix, batchsize=hp.batch_size, shuffle=true, partial=false);
    ps          = Flux.params(cdl);
    opt         = Flux.AdaBelief();

    num_epochs = isnothing(num_epochs) ? setup_num_epochs(data.N) : num_epochs;
    break_condition = false;
    for i in 1:num_epochs
        for S in data_load
            S = S |> gpu;
            gs = gradient(ps) do
                forward_pass_return_loss(S, cdl, hp, len, projs)
            end

            Flux.Optimise.update!(opt, ps, gs) # update parameters
            l1_loss = sum(abs.(prep_syntax_filters(cdl.F))) 
            # "l1 loss: $l1_loss" |> println
            if l1_loss < l1_loss_thresh 
                break_condition = true
                break
            end
        end
        break_condition && break
        println("Epoch: $i completed")
    end
    return cdl, hp, len, projs
end