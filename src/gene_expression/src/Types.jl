mutable struct TranscriptionalControlModel

    # data for gene model -
    transcriptional_weight_matrix::Array{Float64,2}
    transcriptional_parameter_matrix::Array{Float64,2}
    transcriptional_actor_index_array::Array{Int64,1}

    # length factors -
    nt_seq_length::Float64
    aa_seq_length::Float64

    # constructor -
    function TranscriptionalControlModel()
        this = new()
    end
end
