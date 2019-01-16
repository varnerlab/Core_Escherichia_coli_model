# includes -
include("Include.jl")

# setup some global consts -
const path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
const model_file_name = "modelReg"
const organism_id = :eco    # we use the KEGG organism symbols
const path_to_measurements_file = "$(pwd())/experimental_data/test_data/Measurements.json"

# declare a return type -
mutable struct VLFluxResult


    objective_value::Float64
    flux_array::Array{Float64,1}
    dual_array::Array{Float64,1}
    uptake_array::Array{Float64,1}
    flux_variability_array::Array{Float64,1}

    exit_flag::Int64
    status_flag::Int64

    # constructor -
    function VLFluxResult()
        this = new()
    end
end


"""
TODO: Fill me in with some stuff ...
"""
function maximize_specific_growth_rate()

    # load the default data_dictionary -
    default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id)

    # pass the default dictionary to a customization method -
    updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary)

    # update dictionary with experimental data?
    updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file)

    # calculate the flux variability -
    return calculate_flux_variabilty(updated_data_dictionary,[])
end

function sample_flux_space_with_experimental_constraints(solution_bounds_array::Array{Float64,2}, number_of_samples::Int64)

    # load the default data_dictionary -
    default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id)

    # update dictionary with experimental data?
    updated_data_dictionary = constrain_measured_fluxes(default_data_dictionary, path_to_measurements_file)

    # call the sample method -
    return sample_flux_space(solution_bounds_array, updated_data_dictionary,number_of_samples)
end

# call an executable function -
# (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag) = maximize_specific_growth_rate()
(calculated_flux_array,dual_value_array) = maximize_specific_growth_rate()
#flux_ensemble = sample_flux_space_with_experimental_constraints(calculated_flux_array[:,2:end],1000)
