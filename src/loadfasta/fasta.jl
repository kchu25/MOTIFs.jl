
const dna_meta_data = Vector{NamedTuple{(:str, :motif_where, :mode), 
                            Tuple{String, UnitRange{Int64}, Int64}}}


mutable struct FASTA_DNA{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    raw_data_test::Vector{String}
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}, Nothing}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
    data_matrix_bg_gpu::Union{CuArray{S,3}, CuArray{S,2}, Nothing}
    labels::Union{Nothing, Vector{String}, Vector{Int}}
    meta_data::Union{Nothing, dna_meta_data}
    acgt_freq_test::Union{Nothing, Vector{S}}
    markov_bg_mat_test::Union{Nothing, Matrix{S}}
    data_matrix_test::Union{Nothing, Array{S,3}, Array{S,2}}
    data_matrix_bg_test::Union{Nothing, Array{S,3}, Array{S,2}}
    N_test::Int

    function FASTA_DNA{S}(dna_read::Vector{String};
                          k_train=1, k_test=2, # kmer frequency in the test set 
                          train_test_split_ratio=0.9,
                          shuffle=true) where {S <: Real}
        data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
            data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
                markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(dna_read; k_train=k_train, k_test=k_test, 
                                  train_test_split_ratio=train_test_split_ratio, 
                                  shuffle=shuffle, 
                                  FloatType=S);
        L = Int(size(data_matrix,1)/4);
        data_matrix = reshape(data_matrix, 4*L, 1, N_train);
        data_matrix_test = reshape(data_matrix_test, 4*L, 1, N_test)
        data_matrix_bg = reshape(data_matrix_bg, 4*L, 1, N_train)
        new(        
            N_train,
            L,
            acgt_freq,
            markov_bg_mat,
            dna_read[train_set_inds],
            dna_read[test_set_inds],
            data_matrix,
            nothing,
            data_matrix_bg,
            nothing,
            nothing,
            nothing,
            acgt_freq_test,
            markov_bg_mat_test,
            data_matrix_test,
            data_matrix_bg_test,
            N_test
            )
    end

    function FASTA_DNA{S}(fasta_location::String; 
                        max_entries=max_num_read_fasta,
                        k_train=1, k_test=2, # kmer frequency in the test set 
                        train_test_split_ratio=0.9,
                        shuffle=true
                        ) where {S <: Real}

        dna_read = nothing; labels = nothing;
        dna_read = read_fasta(fasta_location; max_entries);
        # dna_read[1] |> println
        data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
            data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
                markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(dna_read; k_train=k_train, k_test=k_test, 
                                  train_test_split_ratio=train_test_split_ratio, 
                                  shuffle=shuffle, 
                                  FloatType=S);
        L = Int(size(data_matrix,1)/4);
        data_matrix = reshape(data_matrix, 4*L, 1, N_train);
        data_matrix_test = reshape(data_matrix_test, 4*L, 1, N_test)
        data_matrix_bg = reshape(data_matrix_bg, 4*L, 1, N_train)
        new(        
            N_train,
            L,
            acgt_freq,
            markov_bg_mat,
            dna_read[train_set_inds],
            dna_read[test_set_inds],
            data_matrix,
            nothing,
            data_matrix_bg,
            nothing,
            labels,
            nothing,
            acgt_freq_test,
            markov_bg_mat_test,
            data_matrix_test,
            data_matrix_bg_test,
            N_test
            )
    end
end