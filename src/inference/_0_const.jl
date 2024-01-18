const float_type_retrieval = Float16

const stored_code_component_t = 
    NamedTuple{(:position, :fil, :seq, :mag), Tuple{UInt16, UInt16, UInt32, float_type_retrieval}}

const composition_key_type =
     NamedTuple{(:f1,:f2,:f3,:d12,:d13,:len), Tuple{Int8, Int8, Int8, UInt16, UInt16, UInt16}}

const value_type = 
    NamedTuple{(:seq_num, :pos, :comp), Tuple{UInt32, UInt16, Bool}}

# const se_type = 
#     NamedTuple{(:start_, :end_), Tuple{Int16, Int16}}

const unit_range_int_t = UInt32

const cover_more_than = 200
const cover_at_least = 10
const syntax_filter_thresh = 0.01
const max_pwm_length_Touzet = 15
# Touzet's threshold calulation pvalue_Touzet2
const _granularity_ = 1e-1; # initial granularity for score2pvalue and pval2score
const _k_ = 100; # decreasing factor for finer granularity in each iteration
const _bg_ = [.25,.25,.25,.25]; # default background
const pvalue_Touzet = 0.0002

# for finding the best socre threshold
const score_thresh_increment = float_type_retrieval(0.5)

const background = float_type_retrieval.([0.25 for _ = 1:4])

const threads_1d = 512;
const threads_2d = 32;
const threads_3d = 10;
const ker_1d = threads_1d;
const ker_2d = (threads_2d, threads_2d);
const ker_3d = (threads_3d, threads_3d, threads_3d);

b_size_1d(X) = ceil.(Int, size(X) ./ threads_1d)
b_size_2d(X) = ceil.(Int, size(X) ./ threads_2d)
b_size_3d(X) = ceil.(Int, size(X) ./ threads_3d)

const atcg2dummy = Dict{Char, Vector{Float32}}('a'=>[1,0,0,0], 'A'=>[1,0,0,0],
                                               'c'=>[0,1,0,0], 'C'=>[0,1,0,0],
                                               'g'=>[0,0,1,0], 'G'=>[0,0,1,0],
                                               't'=>[0,0,0,1], 'T'=>[0,0,0,1],
                                               'z'=>[0,0,0,0])

const atcg_comp = Dict{Char, Char}('a'=>'t', 'A'=>'T', 
                                   'c'=>'g', 'C'=>'G',
                                   'g'=>'c', 'G'=>'C',
                                   't'=>'a', 'T'=>'A')

four_based_ind(ind) = (ind-1)*4+1
promote_i(x...) = Int.(x);

const num_pfms2process = 500
const pvalue_fisher_thresh      = 1e-10 # cannot use float16 since it's pvalue ( use log scale instead?)
const _pseudocounts_            = float_type_retrieval(0.0001)
const allr_thresh               = float_type_retrieval(0.875)
const ic_trim_thresh            = float_type_retrieval(1.0)
const ic_expand_thresh          = float_type_retrieval(0.1)
const max_pwm_length_Touzet2    = 15
const pvalue_Touzet_large = 0.0001
const pvalue_Touzet_mid   = 0.0001
const pvalue_Touzet_small = 0.0003
const max_allowed_diff  = 80
const indep_run         = 13;
const dep_run           = 10;
const merge_add_counts  = true
const code_percentile   = 0.01
const effective_pos_ic_thresh = 0.5
const mv_avg_window = 3
const pfm_minimal_length = 8