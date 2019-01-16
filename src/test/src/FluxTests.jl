function sample_flux_space_test()

    # load a data data_dictionary -
    dd = load_default_data_dictionary_test()

    # sample -
    return sample_flux_space(dd,10)
end

function network_decomposition_test()

    # load a data data_dictionary -
    dd = load_default_data_dictionary_test()

    # get the STM -
    stoichiometric_matrix = dd["stoichiometric_matrix"]

    return calculate_basis_null_space(stoichiometric_matrix)
end

function flux_sample_test()

    

end
