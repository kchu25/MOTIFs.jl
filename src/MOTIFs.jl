module MOTIFs

# modules for motif discovery
using Flux, CUDA, Zygote, StatsBase, LinearAlgebra, Random, SeqShuffle
# modules for inferring the results
using Dictionaries, DataStructures, StaticArrays
# modules for rendering the results
using HypothesisTests, FLoops, Mustache, DataFrames

using Zygote: @ignore

export discover_motifs

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

function get_data_bg(data)
    this_bg = reshape(sum(reshape(data.data_matrix, (4, data.L, data.N)), dims=(2,3)), (4,))
    this_bg = float_type.(this_bg ./ sum(this_bg))
    return this_bg
end

include("loadfasta/helpers.jl")
include("loadfasta/fasta.jl")
include("model.jl")
include("train.jl")
include("inference/_0_const.jl")
include("inference/_1_code_retrieval.jl")
include("inference/_2_enumerate.jl")
include("inference/_3_make_pfms.jl")
include("inference/_s1_make_motifs.jl")
include("inference/_s2_filter_pos_w_scores.jl")
include("inference/_h0_trim.jl")
include("inference/_h1_0_merge_header.jl")
include("inference/_h1_1_merge_H.jl")
include("inference/_h1_2_merge_countmats.jl")
include("inference/_h2_Touzet.jl")
include("inference/_h3_1_alignment.jl")
include("inference/_h3_2_alignment_merge.jl")
include("inference/_h4_overlap_ratio.jl")
include("inference/_h5_expansion.jl")
include("inference/_h6_positions2countmat.jl")
include("inference/_h7_fisher.jl")
include("inference/_h8_remove_redundancy.jl")
include("inference/_h9_order_for_display.jl")
include("inference/_g1_obtain_coutmats.jl")
include("render/const.jl")
include("render/helpers.jl")
include("render/html_template_olap.jl")
include("render/html_template_no_olap.jl")
include("render/plotting.jl")
include("render/pvec_calculations.jl")
include("render/render.jl")
include("wrap.jl")

end
