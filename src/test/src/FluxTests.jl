const path_to_measurements_file = "$(pwd())/experimental_data/test_data/Glucose.json"

function sample_flux_space_test(organism_id, number_of_runs, number_of_samples, sweep_index)

    # load a data data_dictionary -
    dd = load_default_data_dictionary_test(organism_id)

    # pass the default dictionary to a customization method -
    updated_data_dictionary = optimize_specific_growth_rate(dd)

    # update dictionary with experimental data?
    updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file)

    # sweep the growth rate -
    flux_ensemble = zeros(95,1)
    for run_index = 1:number_of_runs

        # flux_ensemble = objective_function_sweep([13], updated_data_dictionary, 3)
        soln_bounds_array = updated_data_dictionary["flux_bounds_array"]
        soln_bounds_array[13,1] = (0.0 + (run_index - 1)*0.1)
        soln_bounds_array[13,2] = 1.8

        soln_set = sample_flux_space(soln_bounds_array,updated_data_dictionary,number_of_samples)
        flux_ensemble = [flux_ensemble soln_set]
    end

    # sample -
    return flux_ensemble[:,2:end]
end
