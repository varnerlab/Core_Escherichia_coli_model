# includes -
include("Include.jl")

# setup some global consts -
const path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
const model_file_name = "modelReg"
const organism_id = :eco    # we use the KEGG organism symbols
const path_to_measurements_file = "$(pwd())/experimental_data/test_data/Glucose.json"

# declare a return type for the flux variability calculation -
mutable struct VLFluxVariabilityResult

    optimal_flux_array::Array{Float64,1}
    flux_variability_array::Array{Float64,1}

    # constructor -
    function VLFluxVariabilityResult()
        this = new()
    end
end

# declare a return type -
mutable struct VLOptimalFluxResult

    objective_value::Float64
    flux_array::Array{Float64,1}
    dual_array::Array{Float64,1}
    uptake_array::Array{Float64,1}
    exit_flag::Int64
    status_flag::Int64
    flux_bounds_array::Array{Float64,2}

    # constructor -
    function VLOptimalFluxResult()
        this = new()
    end
end


"""
TODO: Fill me in with some stuff ...
"""
function maximize_specific_growth_rate(path_to_measurements_file::String; number_of_samples::Int64 = 100)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id);

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);

        # estimate the optimal flux distrubution -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # build a return type -
        fluxResult = VLOptimalFluxResult();
        fluxResult.objective_value = objective_value;
        fluxResult.flux_array = calculated_flux_array;
        fluxResult.dual_array = dual_value_array;
        fluxResult.uptake_array = uptake_array;
        fluxResult.exit_flag = exit_flag;
        fluxResult.status_flag = status_flag;

        # get some problem setup information -
        fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

        # cache -
        push!(results_array, fluxResult);

        # user message -
        msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,msg)]);
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end

# call an executable function -
number_of_samples = 1000;
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file; number_of_samples = 1000);
