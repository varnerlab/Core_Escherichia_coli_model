function calculate_residual(flux_vector,STM,W)

    # compute the species error -
    species_error_vector = STM*flux_vector

    # compute the residual -
    residual = transpose(species_error_vector)*W*species_error_vector

    # return -
    return residual
end

function calculate_gradient!(G,flux_vector,STM,W)

    # compute the residual -
    residual = STM*flux_vector

    # compute the grad -
    tmp_vector = transpose(residual)*(transpose(W)+W)*STM
    index = 1
    for value in tmp_vector
        G[index] = value
        index = index + 1
    end
end

function sample_flux_space(solution_bounds_array::Array{Float64,2}, data_dictionary::Dict{String,Any}, number_of_samples::Int64)

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # get some stuff from the data dictionary -
    STM = data_dictionary["stoichiometric_matrix"]

    # Compute the W matrix -
    (number_of_species,number_of_reactions) = size(STM)
    W = 10000*Matrix{Float64}(I,number_of_species,number_of_species)

    # setup the call to the optimizer --
    objective_function(x) = calculate_residual(x,STM,W)
    gradient_function!(G,x) = calculate_gradient!(G,x,STM,W)

    # get the bounds array -
    lower_bound_array = solution_bounds_array[:,1]
    upper_bound_array = solution_bounds_array[:,2]

    # initialize the ensemble -
    number_of_fluxes = length(lower_bound_array)
    flux_ensemble = zeros(number_of_fluxes,1)

    sobol_sequence = SobolSeq(number_of_fluxes)
    good_sample_count = 1
    q_vector = Float64[]
    flux = Float64[]
    for sample in sobol_sequence

        # generate a q_vector -
        q_vector = lower_bound_array.*sample + (1 .- sample).*upper_bound_array

        # make a call to the optim package -
        result = optimize(objective_function,gradient_function!,lower_bound_array, upper_bound_array, q_vector, Fminbox(LBFGS()))

        # get the flux -
        flux = Optim.minimizer(result)

        # check the error?
        converged_flag = Optim.converged(result)
        local_error = STM*flux

        if (converged_flag == true && norm(local_error)<1e-8)

            # user message -
            msg = "Completed $(good_sample_count) of $(number_of_samples) trials ..."

            # update the progress bar -
            ProgressMeter.next!(p; showvalues = [(:status,msg)])

            # package -
            flux_ensemble = [flux_ensemble flux]

            # update the sample count -
            good_sample_count = good_sample_count + 1
        end

        # should we stop?
        if (good_sample_count>number_of_samples)
            @info "Completed ...\r"
            break
        end
    end

    # return -
    return flux_ensemble[:,2:end]
end
