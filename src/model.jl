Base.@kwdef mutable struct Hyperparam
    filter_len::Int = 8                 # length of each filter
    f_len::Int = filter_len*4           # length of the filter adjusted for one hot DNA sequence
    M::Int = 50                         # number of filters
    twoM::Int = 2*M                     # 2*M
    h::Int = 12                         # height of the syntax filter
    K::Int = 24                         # number of syntax filters
    q::Int = 32                         # how sparse the syntax code should be for each sequence
    batch_size::Int = 6                 # batch size
    num_pass_xyz::Int = 6               # number of passes for the x,y,z parameters
    num_pass_df::Int = 3                # number of passes for the d,f parameters
    magnifying_factor::float_type = 10  # magnifying factor for the sparse code Z, Y
    gamma::float_type = 0.1             # regularization for filter incoherence
end

struct length_info
    L::Int                              # length of the sequence
    C::Int                              # code length
    c::Int                              # code length divided by 4
    l::Int                              # syntax code length 
    MB::Int                             # hp.M * hp.batch_size
    KB::Int                             # hp.K * hp.batch_size
    CS_vlen::Int                        # C + L - 1
    last_avail_ind::Int                 # last data point in the code that can be used for a full batch

    function length_info(hp, data)
        L               = 4*data.L
        C               = L-hp.f_len+1
        c               = data.L-hp.filter_len+1
        l               = c-hp.h+1
        MB              = hp.M * hp.batch_size
        KB              = hp.K * hp.batch_size
        CS_vlen         = C + L - 1;
        last_avail_ind  = (data.N - data.N % hp.batch_size);
        new(L, C, c, l, MB, KB, CS_vlen, last_avail_ind)
    end
end

struct projectors
    mapdrange::CuArray{float_type, 2}
    mapclarge::CuArray{float_type, 2}
    z_mask_n::CuArray{float_type, 3, CUDA.Mem.DeviceBuffer}
    pseudocount_matrix::CuArray{float_type, 3}

    function projectors(hp, len)
        mapdrange                               = zeros(float_type, (hp.f_len, len.C+len.L-1));
        mapdrange[:, len.C:len.C+hp.f_len-1]    = Matrix(I, hp.f_len, hp.f_len);
        mapdrange                               = cu(mapdrange);

        mapclarge                               = zeros(float_type, (len.C, len.c));    
        mapclarge[1:4:end,:]                    = Matrix(I, len.c, len.c);
        mapclarge                               = cu(mapclarge);

        z_mask_col  = 1:len.C .∈ [1:4:len.C]; 
        z_mask_n    = cu(float_type.(repeat(z_mask_col, outer=(1,hp.M,hp.batch_size))));
        pseudocount_matrix = fill(0.001f0, (4, hp.filter_len, hp.M))

        new(
            mapdrange,
            mapclarge,
            z_mask_n,
            pseudocount_matrix
        )
    end
end

struct ucdl
    lambda_sparsity_warmup::float_type                  # sparsity param for Z and Y when their initial values are zero
    lambda_sparsity::Array{float_type, 1}               # sparsity param for Z and Y during iterations
    kappa_sparsity::Array{float_type, 1}                # sparsity param for F during iterations

    lambda_stepsize_warmup::float_type                  # step size for Z and Y when their initial values are zero
    omega_stepsize_warmup::float_type                   # step size for X when its initial values is zero
    lambda_stepsize::Array{float_type, 1}               # step size for Z and Y during iterations
    omega_stepsize::Array{float_type, 1}                # step size for X during iterations
    kappa_stepsize::Array{float_type, 1}                # step size for F during iterations

    D::CuArray{float_type, 3}                           # filters (f_len, 1, M, K)    
    F::CuArray{float_type, 4}                           # syntax filters (hp.h, hp.twoM, 1, hp.K)
    
    penalty_xyz::Array{float_type, 1}                   # penalty for x,y,z
    mu::Array{float_type, 1}                            # step size for D during iterations

    function ucdl(hp; η₁=float_type(0.05))
        D                       = cu(randomly_initialize_filters(
                                   repeats=hp.filter_len, 
                                   how_many_filters=hp.M,
                                   float_type=float_type));
        D                       = sqrt.(D);
        F = cu(abs.(0.1 .* randn(float_type, (hp.h, hp.twoM, 1, hp.K))));
        lambda_sparsity_warmup  = η₁ * rand(float_type);
        lambda_sparsity         = η₁ * rand(float_type, hp.num_pass_xyz);
        kappa_sparsity          = η₁ * rand(float_type, hp.num_pass_df);
        lambda_stepsize_warmup  = η₁ * rand(float_type);
        omega_stepsize_warmup   = η₁ * rand(float_type);
        lambda_stepsize         = η₁ * rand(float_type, hp.num_pass_xyz);
        omega_stepsize          = η₁ * rand(float_type, hp.num_pass_xyz);
        kappa_stepsize          = η₁ * rand(float_type, hp.num_pass_df);
        penalty_xyz             = η₁ * rand(float_type, hp.num_pass_xyz);
        mu                      = η₁ * rand(float_type, hp.num_pass_df);
        new(
            lambda_sparsity_warmup, lambda_sparsity, kappa_sparsity,
            lambda_stepsize_warmup, omega_stepsize_warmup, lambda_stepsize, omega_stepsize, kappa_stepsize,
            D, F, penalty_xyz, mu
        )
    end
    function ucdl(lambda_sparity_warmup, 
                        lambda_sparsity, 
                        kappa_sparsity, 
                        lambda_stepsize_warmup, 
                        omega_stepsize_warmup, 
                        lambda_stepsize, 
                        omega_stepsize, 
                        kappa_stepsize, 
                        D, 
                        F, 
                        penalty_xyz, 
                        mu)
        return new(
            lambda_sparity_warmup, 
            lambda_sparsity, 
            kappa_sparsity, 
            lambda_stepsize_warmup, 
            omega_stepsize_warmup, 
            lambda_stepsize, 
            omega_stepsize, 
            kappa_stepsize, 
            D, 
            F, 
            penalty_xyz, 
            mu
        )
    end

end

Flux.@functor ucdl

function prep_filters(D, hp, projs) 
    D_init      = D .^ 2  
    D_init_r    = reshape(D_init, (4, hp.filter_len, hp.M))
    D_init_r    += projs.pseudocount_matrix               # pseudocounts to avoid division by zero
    D_init_r    = D_init_r ./ sum(D_init_r, dims=1)
    D_init      = reshape(D_init_r, (hp.f_len, 1, hp.M))
    return D_init
end

function prep_syntax_filters(F)
    F = F.^2
    return F ./ (sqrt.(sum(F.^2, dims=(1,2)))) # normalize F
end

function prep_params(ucdl, hp, projs)
    lambda_sparsity_warmup  = ucdl.lambda_sparsity_warmup^2
    lambda_sparsity         = ucdl.lambda_sparsity.^2
    kappa_sparsity          = ucdl.kappa_sparsity.^2
    lambda_stepsize_warmup  = ucdl.lambda_stepsize_warmup^2
    omega_stepsize_warmup   = ucdl.omega_stepsize_warmup^2
    lambda_stepsize         = ucdl.lambda_stepsize.^2
    omega_stepsize          = ucdl.omega_stepsize.^2
    kappa_stepsize          = ucdl.kappa_stepsize.^2
    penalty_xyz             = ucdl.penalty_xyz.^2
    mu                      = ucdl.mu.^2
    D                       = prep_filters(ucdl.D, hp, projs)
    F                       = prep_syntax_filters(ucdl.F)
    return lambda_sparsity_warmup, lambda_sparsity, kappa_sparsity,
           lambda_stepsize_warmup, omega_stepsize_warmup, lambda_stepsize, omega_stepsize, kappa_stepsize,
           penalty_xyz, mu, D, F
end

function warmup_ZY(S, D, lambda_stepsize_warmup, lambda_sparsity_warmup, projs)
    DᵀS             = convolution(S, D, pad=0, flipped=true)
    DS              = convolution(S, D, pad=0)
    Z_update        = (lambda_stepsize_warmup .* DᵀS) .- (lambda_sparsity_warmup * lambda_stepsize_warmup)
    Y_update        = (lambda_stepsize_warmup .* DS)  .- (lambda_sparsity_warmup * lambda_stepsize_warmup)
    Z               = Flux.NNlib.relu.(projs.z_mask_n .* Z_update)
    Y               = Flux.NNlib.relu.(projs.z_mask_n .* Y_update)
    return Z, Y
end

function generate_bitmat(X, hp)
    bitmat      = CUDA.zeros(eltype(X), (size(X, 1), 1, hp.K, hp.batch_size))
    X_reshape   = reshape(X, (size(X, 1)*hp.K, hp.batch_size))
    vals        = reshape(partialsort.(eachcol(X_reshape), hp.q, rev=true) |> cu, (1,1,1, hp.batch_size))
    bitmat[X .≥ vals] .= 1;
    return bitmat
end

function project_X(X, hp)
    bitmat = @ignore generate_bitmat(X, hp)
    return X .* bitmat
end

function create_ZY_mask(ZY)
    ZY_mask                     = CUDA.zeros(eltype(ZY), size(ZY))
    Z_nz                 = ZY[ZY .> 0] 
    if isempty(Z_nz)
        return nothing
    else
        Z_nz_median             = median(Z_nz)
        ZY_mask[ZY .≥ Z_nz_median] .= 1
        return ZY_mask
    end
end

function cat_ZY(Z, Y, hp, len) 
    ZY      = reshape(hcat(Z[1:4:end,:,:], Y[1:4:end,:,:]), (len.c, hp.twoM, 1, hp.batch_size))
    ZY_mask = @ignore create_ZY_mask(ZY);
    return isnothing(ZY_mask) ? hp.magnifying_factor .* ZY : hp.magnifying_factor .* (ZY_mask .* ZY)
end

function warmup_X(F, Z, Y, omega_stepsize_warmup, hp, len)
    ZY          = cat_ZY(Z, Y, hp, len) # (len.c, 2hp.M, 1, hp.batch_size
    X_updated   = omega_stepsize_warmup .* convolution(ZY, F, pad=0, flipped=true) # (len.l, 1, hp.K, hp.batch_size)
    return project_X(X_updated, hp)
end

function return_left_right_FX(FX, hp, len, projs)
    left_FX     = reshape(FX[:, 1:hp.M, :, :],     (len.c, hp.M, hp.batch_size))
    right_FX    = reshape(FX[:, hp.M+1:end, :, :], (len.c, hp.M, hp.batch_size))
    return left_FX, right_FX
end

function warmup_XYZ(S, D, F, lambda_stepsize_warmup, lambda_sparsity_warmup, 
                    omega_stepsize_warmup, hp, len, projs
                    )
    Z, Y                = warmup_ZY(S, D, lambda_stepsize_warmup, lambda_sparsity_warmup, projs)
    X                   = warmup_X(F, Z, Y, omega_stepsize_warmup, hp, len)
    FX                  = sum(convolution(X, F, pad=(hp.h-1, hp.twoM-1), groups=hp.K), dims=3)
    left_FX, right_FX   = return_left_right_FX(FX, hp, len, projs)
    return Z, Y, X, FX, left_FX, right_FX
end

#=
eta, mu: scaled dual variables 
=#
function update_ZY(S, Z, Y, D, left_FX, right_FX, alpha, beta, lambda_sparsity, lambda_stepsize, penalty_xyz, hp, projs, num_pass)    
    ZD, YD      = convolution(Z, D, pad=hp.f_len-1, groups=hp.M), convolution(Y, D, pad=hp.f_len-1, groups=hp.M, flipped=true)
    diff        = sum(ZD + YD, dims=2) - S
    z_grad      = convolution(diff, D, pad=0, flipped=true) + penalty_xyz[num_pass] .* (Z - batched_mul(projs.mapclarge, left_FX  + alpha))
    y_grad      = convolution(diff, D, pad=0)               + penalty_xyz[num_pass] .* (Y - batched_mul(projs.mapclarge, right_FX + beta))
    Z_updated   = Z - lambda_stepsize[num_pass] .* z_grad .- (lambda_sparsity[num_pass] .* lambda_stepsize[num_pass])
    Y_updated   = Y - lambda_stepsize[num_pass] .* y_grad .- (lambda_sparsity[num_pass] .* lambda_stepsize[num_pass])
    return Flux.NNlib.relu.(projs.z_mask_n .* Z_updated), Flux.NNlib.relu.(projs.z_mask_n .* Y_updated)
end

function update_X(FX, Z, Y, X, F, alpha, beta, omega_stepsize, hp, len, num_pass)
    alpha_beta  = reshape(hcat(alpha, beta), (len.c, hp.twoM, 1, hp.batch_size))
    ZY          = cat_ZY(Z, Y, hp, len)
    diff        = sum(FX, dims=3) - (ZY - alpha_beta)
    x_grad      = convolution(diff, F, pad=0, flipped=true)
    X_updated   = X - omega_stepsize[num_pass] .* x_grad
    return project_X(X_updated, hp)
end

function one_forward_step_XYZ(S, Z, Y, D, X, F, left_FX, right_FX, FX, alpha, beta, 
                              lambda_sparsity, lambda_stepsize, 
                              omega_stepsize, penalty_xyz, 
                              hp, len, projs, num_pass
                              )
    Z, Y                = update_ZY(S, Z, Y, D, left_FX, right_FX, alpha, beta, lambda_sparsity, lambda_stepsize, penalty_xyz, hp, projs, num_pass)
    X                   = update_X(FX, Z, Y, X, F, alpha, beta, omega_stepsize, hp, len, num_pass)
    FX                  = sum(convolution(X, F, pad=(hp.h-1, hp.twoM-1), groups=hp.K), dims=3)
    left_FX, right_FX   = return_left_right_FX(FX, hp, len, projs)
    alpha               = alpha + left_FX  - (@view Z[1:4:end,:,:])
    beta                = beta  + right_FX - (@view Y[1:4:end,:,:])
    return Z, Y, X, FX, left_FX, right_FX, alpha, beta
end

conv_code_diff(code, diff, hp, len) = reshape(
                                convolution(reshape(upsample_nearest(diff, (1,hp.M,1)), (len.L, len.MB, 1)), 
                                                reshape(code, (len.C, 1, len.MB)) , pad=len.C-1, groups=len.MB, flipped=true),
                                                (len.CS_vlen, hp.M, hp.batch_size));

function update_D(S, Z, Y, D, mu, hp, len, projs, num_pass)
    sumZD       = sum(convolution(Z, D, pad=hp.f_len-1, groups=hp.M), dims=2)
    sumYRD      = sum(convolution(Y, D, pad=hp.f_len-1, groups=hp.M, flipped=true), dims=2)
    ZᵀsumZD     = conv_code_diff(Z, sumZD, hp, len)
    YᵀsumZD     = conv_code_diff(Y, sumZD, hp, len)
    ZᵀsumYRD    = conv_code_diff(Z, sumYRD, hp, len)
    YᵀsumYRD    = conv_code_diff(Y, sumYRD, hp, len)
    ZᵀS         = conv_code_diff(Z, S, hp, len)
    YᵀS         = conv_code_diff(Y, S, hp, len)
    D_grad      = reshape(
                    projs.mapdrange*reshape(sum(ZᵀsumZD + ZᵀsumYRD + ZᵀS + reverse(YᵀsumZD + YᵀsumYRD +YᵀS, dims=1), dims=3), (len.C+len.L-1, hp.M)),
                    (hp.f_len, 1, hp.M));
    Breg_num    = reshape(D .* exp.(-mu[num_pass] .* D_grad), (4, hp.filter_len, 1, hp.M));
    D_updated   = reshape((Breg_num ./ sum(Breg_num,dims=1)), (hp.f_len, 1, hp.M));  
    return D_updated
end

function F_gradient(ZY, X, F, hp, len, theta)
    diff_X_upsampled = 
        upsample_nearest(sum(convolution(X, F, pad=(hp.h-1, hp.twoM-1), groups=hp.K), dims=3) - (ZY  + theta)
        , (1, 1, hp.K, 1))
    diff_r      = reshape(diff_X_upsampled, (len.c, hp.twoM, hp.K*hp.batch_size, 1))
    X_r         = reshape(X, (len.l, 1, 1, hp.K*hp.batch_size))
    conv_diff_X = convolution(diff_r, X_r, pad=0, flipped=true, groups=len.KB)
    F_conv      = reshape(conv_diff_X, (hp.h, hp.twoM, hp.K, hp.batch_size))
    F_grad      = reshape(sum(F_conv, dims=4), (hp.h, hp.twoM, 1, hp.K))
    return F_grad
end

function update_F(ZY, X, F, hp, len, theta, kappa_stepsize, kappa_sparsity, num_pass)
    F_grad      = F_gradient(ZY, X, F, hp, len, theta)
    F_updated   = Flux.NNlib.relu.(F - kappa_stepsize[num_pass] * F_grad .-(kappa_stepsize[num_pass] * kappa_sparsity[num_pass]))
    return F_updated ./  (sqrt.(sum(F_updated.^2, dims=(1,2))))
end

function loss(S, Z, Y, X, D, ZY, F, F_orig, hp)
    normalize_factor = (1.0f0/float_type(hp.batch_size));

    DZ                          = sum(convolution(Z, D, pad=hp.f_len-1, groups=hp.M), dims=2)
    DY                          = sum(convolution(Y, D, pad=hp.f_len-1, groups=hp.M, flipped=true), dims=2)
    reconstruction_loss         = normalize_factor*sum((DZ+DY-S).^2)
    FX                          = sum(convolution(X, F, pad=(hp.h-1, hp.twoM-1), groups=hp.K), dims=3)
    syntax_reconstruction_loss  = normalize_factor*sum((FX - ZY).^2)
    
    # abs_F_orig = abs.(F_orig)
    # l1l2_loss_orig = (sum(abs_F_orig) + sum(sqrt.(sum(abs_F_orig.^2, dims=(1,2)))))
    # l1l2_loss_orig = sum(abs_F_orig) 
    # println("reconstruction_loss: ", reconstruction_loss)
    # println("syntax_reconstruction_loss: ", syntax_reconstruction_loss)
    # println("l1l2_loss_orig: ", l1l2_loss_orig)
    return reconstruction_loss + syntax_reconstruction_loss
end

create_alpha_beta_as_zeros(S, hp, len) = 
    CUDA.zeros(eltype(S), (len.c, hp.M, hp.batch_size)), CUDA.zeros(eltype(S), (len.c, hp.M, hp.batch_size));

function ADMM_XYZ(S, D, F,  
                 lambda_stepsize_warmup, lambda_stepsize, 
                 lambda_sparsity_warmup, lambda_sparsity,
                 omega_stepsize_warmup, omega_stepsize, 
                 penalty_xyz, 
                 hp, len, projs
                 )
    # create initial scaled dual variables as zeros
    alpha, beta = @ignore create_alpha_beta_as_zeros(S, hp, len)

    # warm up
    Z, Y, X, FX, left_FX, right_FX = 
        warmup_XYZ(S, D, F, 
                   lambda_stepsize_warmup, lambda_sparsity_warmup, 
                   omega_stepsize_warmup, hp, len, projs)

    # iterations
    for num_pass = 1:hp.num_pass_xyz
        Z, Y, X, FX, left_FX, right_FX, alpha, beta = 
            one_forward_step_XYZ(S, Z, Y, D, X, F, 
                                 left_FX, right_FX, FX, 
                                 alpha, beta, 
                                 lambda_sparsity, lambda_stepsize, 
                                 omega_stepsize, penalty_xyz, 
                                 hp, len, projs, num_pass)
    end
    return Z, Y, X
end

create_theta_as_zeros(S, hp, len) = 
    CUDA.zeros(eltype(S), (len.c, hp.twoM, 1, hp.batch_size));

function ADMM_DF(S, Z, Y, X, D, F, mu,
                 kappa_sparsity, kappa_stepsize, hp, len, projs)
    # create the initial scaled dual variable as zeros
    theta   = @ignore create_theta_as_zeros(S, hp, len)
    ZY      = cat_ZY(Z, Y, hp, len) 
    for num_pass = 1:hp.num_pass_df
        D = update_D(S, Z, Y, D, mu, hp, len, projs, num_pass)
        F = update_F(ZY, X, F, hp, len, theta, kappa_stepsize, kappa_sparsity, num_pass)
        theta = theta + sum(convolution(X, F, pad=(hp.h-1, hp.twoM-1), groups=hp.K), dims=3) - ZY
    end
    return ZY, D, F
end

function forward_pass_return_loss(S, cdl, hp, len, projs)
    lambda_sparsity_warmup, lambda_sparsity, kappa_sparsity,
    lambda_stepsize_warmup, omega_stepsize_warmup, lambda_stepsize, omega_stepsize, kappa_stepsize,
    penalty_xyz, mu, D, F_orig = prep_params(cdl, hp, projs)

    Z, Y, X = ADMM_XYZ(S, D, F_orig,  
                 lambda_stepsize_warmup, lambda_stepsize, 
                 lambda_sparsity_warmup, lambda_sparsity,
                 omega_stepsize_warmup, omega_stepsize, 
                 penalty_xyz, hp, len, projs
                 )

    ZY, D, F = ADMM_DF(S, Z, Y, X, D, F_orig, mu,
                       kappa_sparsity, kappa_stepsize,
                       hp, len, projs
                       )

    return loss(S, Z, Y, X, D, ZY, F, F_orig, hp)
end


function retrieve_code(S, cdl, hp, len, projs)
    lambda_sparsity_warmup, lambda_sparsity, kappa_sparsity,
    lambda_stepsize_warmup, omega_stepsize_warmup, lambda_stepsize, omega_stepsize, kappa_stepsize,
    penalty_xyz, _, D, F_orig = prep_params(cdl, hp, projs)

    Z, Y, X = ADMM_XYZ(S, D, F_orig,  
                 lambda_stepsize_warmup, lambda_stepsize, 
                 lambda_sparsity_warmup, lambda_sparsity,
                 omega_stepsize_warmup, omega_stepsize, 
                 penalty_xyz, 
                 hp, len, projs
                 )
    return F_orig, Z, Y, X
end
