######## fasta load settings: ##################
const max_num_read_fasta = 100000;
const dat_t = Float32;          # data_matrix_type
const int_t = Int32;            # integer type

#= each label must associated with at least 200 sequences 
   for classification 
=#
const label_count_thresh = 250; 
const DNA_complement = Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A');        

reverse_complement(s::String) = 
    join(islowercase(s[si]) ? s[si] : DNA_complement[s[si]] for si = length(s):-1:1)

get_count_map(v) = countmap(v)

"""
Assume all dna strings are of the same length and contains no masked strings
"""
function reading_for_DNA_regression(filepath::String; parse_float_type=Float32)
    reads = read(filepath, String);
    splits = split(reads, '>');
    dna_reads = Vector{String}();
    labels = Vector{parse_float_type}();
    for i in splits
        if !isempty(i)
            splits_i = split(i, "\n");
            if !occursin("N", splits_i[2]) && !occursin("n", splits_i[2])
                push!(labels, parse(parse_float_type, splits_i[1]));
                push!(dna_reads, splits_i[2]);
            end
        end
    end
    return labels, dna_reads
end

function permute_dataset(labels, seqs; train_test_ratio=0.8)
    shuffled_indices = randperm(labels |> length)
    n_rows = length(shuffled_indices);
    n_train = round(Int, n_rows*train_test_ratio)
    n_test = n_rows - n_train;    
    labels_train = @view labels[shuffled_indices][1:n_train]
    seqs_train   = @view seqs[shuffled_indices][1:n_train]
    labels_test  = @view labels[shuffled_indices][n_train+1:end]
    seqs_test    = @view seqs[shuffled_indices][n_train+1:end]
    return n_train, n_test, labels_train, seqs_train, labels_test, seqs_test
end

function read_and_permute(filepath::String; train_test_ratio=0.8, parse_float_type=Float32)
    labels, seqs = reading_for_DNA_regression(filepath; parse_float_type=parse_float_type)
    return permute_dataset(labels, seqs; train_test_ratio=train_test_ratio)
end


"""
Read the strings from the datasets provided by 
https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks
and processed by 
https://github.com/kchu25/genomic_benchmark_datasets

Assume all dna strings are of the same length and contains no masked strings
"""
function reading_for_DNA_classification(filepath::String)
    reads = read(filepath, String);
    splits = split(reads, '>');
    dna_reads = Vector{String}();
    labels = Vector{Int}();
    for i in splits
        if !isempty(i)
            splits_i = split(i, "\n");
            if !occursin("N", splits_i[2]) && !occursin("n", splits_i[2])
                push!(labels, parse(Int, splits_i[1]));
                push!(dna_reads, splits_i[2]);
            end
        end
    end

    uniq_labels = unique(labels)
    labels = reduce(hcat, [uniq_labels .== v for v in labels]); # bitarrays for class indicates

    return labels, dna_reads
end

function reading(filepath::String;
                 max_entries=max_num_read_fasta)
    # read the file; process the fasta file to get the DNA part; # rid of sequences that contains "n"
    reads = read(filepath, String);
    dna_reads = Vector{String}();
    for i in split(reads, '>')
        if !isempty(i)
            splits = split(i, "\n");
            this_read = join(splits[2:end]);        
            !occursin("N", this_read) && !occursin("n", this_read) && (push!(dna_reads, this_read);)            
        end
    end    
    dna_reads = length(dna_reads) > max_entries ? dna_reads[1:max_entries] : dna_reads;    
    # rid of all dna sequences that's not the same length as sequence 1
    # o/w markov mat assignment may report error
    dna_reads = [s for s in dna_reads if length(s) == length(dna_reads[1])];
    return dna_reads
end

function read_fasta(filepath::String; 
                    max_entries=max_num_read_fasta
                    )::Vector{String}
    #= read a fasta file =#
    dna_reads = reading(filepath; max_entries);   
    return [uppercase(i) for i in dna_reads]
end

function dna2dummy(dna_string::String, dummy::Dict; F=Float32)    
    v = Array{F,1}(undef, 4*length(dna_string));
    for (index, alphabet) in enumerate(dna_string)
        start = (index-1)*4+1;
        v[start:start+3] = dummy[uppercase(alphabet)];
    end
    return v
end

#=
get the set of dummy-vectors from a set of dna-strings
the dummy-vectors are all of same length (for now)
=#
function data_2_dummy(dna_strings; F=Float32)

    dummy = Dict('A'=>Array{F}([1, 0, 0, 0]), 
                   'C'=>Array{F}([0, 1, 0, 0]),
                   'G'=>Array{F}([0, 0, 1, 0]), 
                   'T'=>Array{F}([0, 0, 0, 1]));

    how_many_strings = length(dna_strings);
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    how_many_strings == 0  && return nothing;
    _len_ = length(dna_strings[1]); # length of each dna string in data    
    _S_ = Array{F, 2}(undef, (4*_len_, how_many_strings));
    for i = 1:how_many_strings
        length(dna_strings[i]) == _len_ && (@inbounds _S_[:, i] = dna2dummy(dna_strings[i], dummy))
    end
    return _S_
end

function get_train_test_inds(dna_read, train_test_split_ratio, shuffle)
    # println(shuffle)
    len_dna_read = length(dna_read)
    how_may_entries_in_test = Int(floor((1-train_test_split_ratio)*len_dna_read));
    test_set_inds = nothing;

    shuffled_inds = randperm(len_dna_read)
    # inds was 1:len_dna_read
    if shuffle 
        test_set_inds = sample(shuffled_inds, 
                      how_may_entries_in_test, 
                      replace=false)
    else
        test_set_inds = collect(len_dna_read-how_may_entries_in_test+1:len_dna_read)
    end
    # println("test_set_inds: ", test_set_inds)
    train_set_inds = setdiff(shuffled_inds, test_set_inds)
    return train_set_inds, test_set_inds
end

function get_data_matrices2(dna_read, labels; 
                           k_train=1, k_test=2, 
                           FloatType=dat_t, 
                           train_test_split_ratio=0.85,
                           shuffle=true)
    # set train_test_split_ratio = 1.0 if no test set is needed    
    train_set_inds, test_set_inds = get_train_test_inds(dna_read, train_test_split_ratio, shuffle)
    # println(train_set_inds)
    dna_read_train = @view dna_read[train_set_inds]
    dna_read_test = @view dna_read[test_set_inds]    
    labels_train = labels[train_set_inds]
    labels_test = labels[test_set_inds]
    
    shuffled_dna_read_train = seq_shuffle.(dna_read_train; k=k_train);
    data_matrix_train = data_2_dummy(dna_read_train; F=FloatType);
    data_matrix_bg_train = data_2_dummy(shuffled_dna_read_train; F=FloatType);

    shuffled_dna_read_test = seq_shuffle.(dna_read_test; k=k_test);
    data_matrix_test = data_2_dummy(dna_read_test; F=FloatType);
    data_matrix_bg_test = data_2_dummy(shuffled_dna_read_test; F=FloatType);
    
    # estimate the Markov background (order 1)
    acgt_freq_train, markov_mat_train = est_1st_order_markov_bg(shuffled_dna_read_train; F=FloatType);
    acgt_freq_test, markov_mat_test = est_1st_order_markov_bg(shuffled_dna_read_test; F=FloatType);
    data_bg_prob_train = SeqShuffle.assign_bg_prob(shuffled_dna_read_train, markov_mat_train, acgt_freq_train);
    data_bg_prob_test = SeqShuffle.assign_bg_prob(shuffled_dna_read_test, markov_mat_test, acgt_freq_test);
   
    return data_matrix_train, 
           data_matrix_bg_train, 
           data_bg_prob_train, 
           acgt_freq_train, 
           markov_mat_train, 
           data_matrix_test,
           data_matrix_bg_test,
           data_bg_prob_test,
           acgt_freq_test,
           markov_mat_test,
           length(dna_read_train),
           length(dna_read_test),
           train_set_inds,
           test_set_inds,
           labels_train,
           labels_test
end

function get_data_matrices(dna_read; 
                           k_train=1, k_test=2, 
                           FloatType=dat_t, 
                           train_test_split_ratio=0.85,
                           shuffle=true)
    # set train_test_split_ratio = 1.0 if no test set is needed    
    train_set_inds, test_set_inds = get_train_test_inds(dna_read, train_test_split_ratio, shuffle)
    # println(train_set_inds)
    dna_read_train = @view dna_read[train_set_inds]
    dna_read_test = @view dna_read[test_set_inds]    
    
    shuffled_dna_read_train = seq_shuffle.(dna_read_train; k=k_train);
    data_matrix_train = data_2_dummy(dna_read_train; F=FloatType);
    data_matrix_bg_train = data_2_dummy(shuffled_dna_read_train; F=FloatType);

    shuffled_dna_read_test = seq_shuffle.(dna_read_test; k=k_test);
    data_matrix_test = data_2_dummy(dna_read_test; F=FloatType);
    data_matrix_bg_test = data_2_dummy(shuffled_dna_read_test; F=FloatType);
    
    # estimate the Markov background (order 1)
    acgt_freq_train, markov_mat_train = est_1st_order_markov_bg(shuffled_dna_read_train; F=FloatType);
    acgt_freq_test, markov_mat_test = est_1st_order_markov_bg(shuffled_dna_read_test; F=FloatType);
    data_bg_prob_train = SeqShuffle.assign_bg_prob(shuffled_dna_read_train, markov_mat_train, acgt_freq_train);
    data_bg_prob_test = SeqShuffle.assign_bg_prob(shuffled_dna_read_test, markov_mat_test, acgt_freq_test);
   
    return data_matrix_train, 
           data_matrix_bg_train, 
           data_bg_prob_train, 
           acgt_freq_train, 
           markov_mat_train, 
           data_matrix_test,
           data_matrix_bg_test,
           data_bg_prob_test,
           acgt_freq_test,
           markov_mat_test,
           length(dna_read_train),
           length(dna_read_test),
           train_set_inds,
           test_set_inds
end