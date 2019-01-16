module VLTestModule

    # include -
    top_level_path = pwd()
    include("$(top_level_path)/test/src/Include.jl")

    # export -
    export load_default_data_dictionary_test, sample_flux_space_test, network_decomposition_test
end
