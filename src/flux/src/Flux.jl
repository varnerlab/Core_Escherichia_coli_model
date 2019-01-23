function objective_function_sweep(sweep_flux_index_array::Array{Int64,1},data_dictionary::Dict{String,Any}, number_of_divisions::Int64)

    # how many fluxes do we have?
    number_of_fluxes = data_dictionary["number_of_reactions"]
    number_of_fluxes_to_sample = length(sweep_flux_index_array)

    # make a copy of the flux_bounds_array -
    fba_original = copy(data_dictionary["flux_bounds_array"])

    # initialize -
    flux_ensemble = zeros(number_of_fluxes,1)

    # do max -
    for flux_index in sweep_flux_index_array

        # initialize -
        c_vector = zeros(number_of_fluxes)
        c_vector[flux_index] = -1.0

        # update the data d -
        data_dictionary["objective_coefficient_array"] = c_vector

        # evaluate lp max and min -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag) = calculate_optimal_flux_distribution(data_dictionary)

        # add this solution to the ensemble -
        flux_ensemble = [flux_ensemble calculated_flux_array]

        # what is the original max solution?
        flux_solution_max = calculated_flux_array

        # sample the between the bounds for this flux -
        theta = 0.0
        fba_new = fba_original
        for sample_index = 1:number_of_divisions

            theta = (number_of_divisions - sample_index)*(1/number_of_divisions)
            upper_bound = theta*flux_solution_max[flux_index]
            fba_new[flux_index,2] = upper_bound
            data_dictionary["flux_bounds_array"] = fba_new

            # call the lp again -
            (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag) = calculate_optimal_flux_distribution(data_dictionary)

            # cache -
            flux_ensemble = [flux_ensemble calculated_flux_array]
        end
    end

    return flux_ensemble[:,2:end]
end

function calculate_optimal_flux_distribution(data_dictionary)

    # set some constants -
    TIME_RESTART_LIM = 60

    # Get the stoichiometric_matrix from data_dictionary -
    tmp_stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
    additional_constraint_array = data_dictionary["additional_constraint_array"];
    stoichiometric_matrix = [tmp_stoichiometric_matrix ; additional_constraint_array];
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

    # # Setup the GLPK problem -
    lp_problem = GLPK.Prob();
    GLPK.set_prob_name(lp_problem, "sample");
    GLPK.set_obj_name(lp_problem, "objective")

    # Set solver parameters
    solver_parameters = GLPK.SimplexParam();
    GLPK.init_smcp(solver_parameters);
    solver_parameters.msg_lev = GLPK.MSG_OFF;

    # Are we doing min -or- max?
    min_flag = data_dictionary["is_minimum_flag"];
    if min_flag == true
    	GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    else
    	GLPK.set_obj_dir(lp_problem, GLPK.MAX);
    end

    # Set the number of constraints and fluxes -
    GLPK.add_rows(lp_problem, number_of_species);
    GLPK.add_cols(lp_problem, number_of_fluxes);

    # Setup flux bounds, and objective function -
    default_bounds_array = data_dictionary["flux_bounds_array"]
    (number_of_fluxes,number_of_bounds) = size(default_bounds_array)
    for flux_index = 1:number_of_fluxes

    	flux_lower_bound = default_bounds_array[flux_index,1]
    	flux_upper_bound = default_bounds_array[flux_index,2]

    	# Check bounds type ... default is DB -
    	if (flux_upper_bound == flux_lower_bound)
    		flux_constraint_type = GLPK.FX
    	else
    		flux_constraint_type = GLPK.DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(flux_index)

    	# Set the bounds in GLPK -
    	GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
    	GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end

    # Setup objective function -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)

    	# Set the objective function value in GLPK -
    	GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
    end

    # Setup problem constraints for the metabolites -
    species_bounds_array = data_dictionary["species_bounds_array"]
    for species_index = 1:number_of_species

    	species_lower_bound = species_bounds_array[species_index,1]
    	species_upper_bound = species_bounds_array[species_index,2]

    	# defualt
    	species_constraint_type = GLPK.FX
    	if (species_lower_bound != species_upper_bound)
    		species_constraint_type = GLPK.DB
    	end

    	# set the symbol -
    	species_symbol = "x_"*string(species_index)

    	# Set the species bounds in GLPK -
    	GLPK.set_row_name(lp_problem, species_index, species_symbol);
    	GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);

    end

    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Int,number_of_species*number_of_fluxes);
    col_index_array = zeros(Int,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
    	for flux_index in flux_index_vector
    		row_index_array[counter] = species_index;
    		col_index_array[counter] = flux_index;
    		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
    		counter+=1;
    	end
    end
    GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);

    # Call the solver -
    exit_flag = GLPK.simplex(lp_problem, solver_parameters);

    # Get the objective function value -
    objective_value = GLPK.get_obj_val(lp_problem);

    # Get the calculated flux values from GLPK -
    calculated_flux_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	calculated_flux_array[flux_index] = GLPK.get_col_prim(lp_problem, flux_index);
    end

    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
    end

    # is this solution optimal?
    status_flag = GLPK.get_status(lp_problem)

    # Calculate the uptake array -
    uptake_array = stoichiometric_matrix*calculated_flux_array;

    # Formulate the return tuple -
    return (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag,status_flag);
end

function calculate_flux_variabilty(data_dictionary,additional_flux_constraints_array)

    # --- Phase 0: GET SOME STUFF ---------------------------------------------------#
    # Get the stoichiometric_matrix from data_dictionary -
    tmp_stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
    stoichiometric_matrix = Float64[]

    # do we have additional constraints?
    if haskey(data_dictionary,"additional_constraint_array") == true
        additional_constraint_array = data_dictionary["additional_constraint_array"]
        stoichiometric_matrix = [tmp_stoichiometric_matrix ; additional_constraint_array]
    else
        stoichiometric_matrix = tmp_stoichiometric_matrix
    end

    # get the size -
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

    # initialize the output array -
    calculated_flux_array = zeros(number_of_fluxes,3)

    # set some constants -
    TIME_RESTART_LIM = 60
    EXIT_FLAG_SUCCESS = 0
    TOL = 1e-6
    # ------------------------------------------------------------------------------- #

    # --- Phase 1: SOLVE THE ORIGINAL PROBLEM --------------------------------------- #

    # Setup the GLPK problem -
    lp_problem = GLPK.Prob();
    GLPK.set_prob_name(lp_problem, "sample");
    GLPK.set_obj_name(lp_problem, "objective")

    # Are we doing min -or- max?
    min_flag = data_dictionary["is_minimum_flag"];
    if min_flag == true
    	GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    else
    	GLPK.set_obj_dir(lp_problem, GLPK.MAX);
    end

    # Set the number of constraints and fluxes -
    GLPK.add_rows(lp_problem, number_of_species);
    GLPK.add_cols(lp_problem, number_of_fluxes);

    # Set solver parameters
    solver_parameters = GLPK.SimplexParam();
    solver_parameters.msg_lev = GLPK.MSG_ERR;
    solver_parameters.presolve = GLPK.ON;

    # Setup flux bounds, and objective function -
    default_bounds_array = data_dictionary["flux_bounds_array"]
    (number_of_fluxes,number_of_bounds) = size(default_bounds_array)
    for flux_index = 1:number_of_fluxes

    	flux_lower_bound = default_bounds_array[flux_index,1]
    	flux_upper_bound = abs(default_bounds_array[flux_index,2])

    	# Check bounds type ... default is DB -
    	if (flux_upper_bound == flux_lower_bound)
    		flux_constraint_type = GLPK.FX
    	else
    		flux_constraint_type = GLPK.DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(flux_index)

    	# Set the bounds in GLPK -
    	GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
    	GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end


    # Setup objective function -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)

    	# Set the objective function value in GLPK -
    	GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
    end

    # Setup problem constraints for the metabolites -
    species_bounds_array = data_dictionary["species_bounds_array"]
    for species_index = 1:number_of_species

    	species_lower_bound = species_bounds_array[species_index,1]
    	species_upper_bound = species_bounds_array[species_index,2]

    	# defualt
    	species_constraint_type = GLPK.FX
    	if (species_lower_bound != species_upper_bound)
    		species_constraint_type = GLPK.DB
    	end

    	# set the symbol -
    	species_symbol = "x_"*string(species_index)

    	# Set the species bounds in GLPK -
    	GLPK.set_row_name(lp_problem, species_index, species_symbol);
    	GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);

    end

    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Int,number_of_species*number_of_fluxes);
    col_index_array = zeros(Int,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
    	for flux_index in flux_index_vector
    		row_index_array[counter] = species_index;
    		col_index_array[counter] = flux_index;
    		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
    		counter+=1;
    	end
    end
    GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);

    # solve --
    # Call the solver -
    exit_flag = GLPK.simplex(lp_problem, solver_parameters);

    # Get the objective function value -
    objective_value = GLPK.get_obj_val(lp_problem);

    # Get the calculated flux values for the optimal distrubution from GLPK -
    for flux_index in flux_index_vector
    	calculated_flux_array[flux_index,1] = GLPK.get_col_prim(lp_problem, flux_index);
    end

    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
    end
    # ------------------------------------------------------------------------------- #

    # --- Phase 2: SOLVE THE FVA PROBLEM -------------------------------------------- #
    target_value = 0.0
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        target_value = floor(objective_value/TOL)*(TOL);
    else
        target_value = ceil(objective_value/TOL)*(TOL);
    end

    # do we have any additional flux constraints?
    number_of_extra_flux_constraints = length(additional_flux_constraints_array)

    # Add an extra row for the objective -
    m = GLPK.add_rows(lp_problem,number_of_extra_flux_constraints+1)
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        GLPK.set_row_bnds(lp_problem, m, GLPK.LO,target_value,0.0);
    else
        GLPK.set_row_bnds(lp_problem, m, GLPK.UP,0.0,target_value);
    end

    # populate new row, add the object value -
    index_array = Int[]
    value_array = Float64[]
    for flux_index = 1:number_of_fluxes
        push!(index_array,flux_index)
        push!(value_array,GLPK.get_obj_coef(lp_problem, flux_index))
    end
    GLPK.set_mat_row(lp_problem,m,number_of_fluxes,index_array,value_array)

    # Add additional rows for the extra flux constraints -
    for flux_constraint_index = 1:number_of_extra_flux_constraints

        # what rate do I want to constrain?
        target_flux_index = additional_flux_constraints_array[flux_constraint_index]

        # What is my current row index?
        current_row_index = m + flux_constraint_index

        # What is my target value?
        target_value = calculated_flux_array[target_flux_index]

        # add the row bound -
        if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
            GLPK.set_row_bnds(lp_problem, current_row_index, GLPK.LO,target_value,0.0);
        else
            GLPK.set_row_bnds(lp_problem, current_row_index, GLPK.UP,0.0,target_value);
        end

        # Add the extra row to the st array -
        index_array = Int[]
        value_array = Float64[]
        for flux_index = 1:number_of_fluxes
            push!(index_array,flux_index)

            if (flux_index == target_flux_index)
                push!(value_array,1.0)
            else
                push!(value_array,0.0)
            end
        end
        GLPK.set_mat_row(lp_problem,current_row_index,number_of_fluxes,index_array,value_array)
    end

    # zero out the objective coefficients -
    for flux_index = 1:number_of_fluxes
        GLPK.set_obj_coef(lp_problem, flux_index, 0.0);
    end

    solver_parameters.presolve = GLPK.OFF;
    solver_parameters.msg_lev = GLPK.MSG_OFF;
    solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;

    # main FVA loop -
    for round_index = 1:2

        # are we min are maximizing the flux?
        if round_index == 1
            GLPK.set_obj_dir(lp_problem, GLPK.MIN);
        else
            GLPK.set_obj_dir(lp_problem, GLPK.MAX);
        end

        for flux_index = 1:number_of_fluxes

            # run the calc -
            GLPK.set_obj_coef(lp_problem, flux_index, 1.0);
            exit_flag = GLPK.simplex(lp_problem, solver_parameters);

            # reset the objective coefficient -
            GLPK.set_obj_coef(lp_problem,flux_index,0.0)

            if (exit_flag != EXIT_FLAG_SUCCESS)

                solver_parameters.tm_lim = 10000*TIME_RESTART_LIM;
                solver_parameters.presolve = GLPK.ON
                GLPK.adv_basis(lp_problem, 0);
                exit_flag = GLPK.simplex(lp_problem, solver_parameters);
                if (exit_flag != EXIT_FLAG_SUCCESS)
                    # ok, something muy malo is happening. throw an error ...
                    error("ERROR: Unable to solve modified problem. Exiting ... ")
                end
            end

            # reset the problem options -
            solver_parameters.presolve = GLPK.OFF;
            solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;

            if (GLPK.get_obj_dir(lp_problem) == GLPK.MIN)
                calculated_flux_array[flux_index,2] = GLPK.get_obj_val(lp_problem)
            else
                calculated_flux_array[flux_index,3] = GLPK.get_obj_val(lp_problem)
            end
        end
    end
    # ------------------------------------------------------------------------------- #

    # return the array -
    return (calculated_flux_array,dual_value_array)
end
