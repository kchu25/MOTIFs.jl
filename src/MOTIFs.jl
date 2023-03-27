module MOTIFs

# Write your package code here.
using Flux, CUDA, Zygote, StatsBase, LinearAlgebra,
    Dictionaries, Random

using Zygote: @ignore


const float_type = Float32
const convolution = Flux.NNlib.conv;

function randomly_initialize_filters(; 
                       dim=4,
                       rng=Random.GLOBAL_RNG, 
                       repeats=5, 
                       how_many_filters=10,
                       float_type=Float16)    
    arr = zeros(float_type,(dim+1,
                            repeats,
                            1,
                            how_many_filters));    
    for i = 1:repeats, j = 1:how_many_filters
        unif = rand(rng, dim-1);
        arr[2:dim,i,1,j] .= sort(unif);
    end
    arr[dim+1,:,:,:] .= 1;
    return reshape(diff(arr, dims=1), (dim*repeats,1,how_many_filters));
end

include("model.jl")


end
