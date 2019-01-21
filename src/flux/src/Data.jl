"""
TODO: Fill me in with some stuff ...
"""
function optimize_specific_growth_rate(data_dictionary::Dict{String,Any})

    # get some stuff from dd -
    number_of_reactions = data_dictionary["number_of_reactions"]

    # update the c-vector -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]
    objective_coefficient_array[13] = -1

    # return -
    return data_dictionary
end

function constrain_measured_fluxes(data_dictionary::Dict{String,Any}, path_to_measurements_file::String)

    # TODO: is the path_to_measurements_file legit?
    measurements_dictionary = JSON.parsefile(path_to_measurements_file)

    # get some stuff from the dd -
    rxn_name_list = data_dictionary["list_of_reaction_name_strings"]
    flux_bounds_array = data_dictionary["flux_bounds_array"]
    (number_of_reactions,nc) = size(flux_bounds_array)

    # do we have any flux ratios?
    flux_ratio_dictionary_array = measurements_dictionary["flux_ratio_measurements"]

    # initialize -
    number_of_additional_constraints = length(flux_ratio_dictionary_array)
    constraint_counter = 1
    additional_constraint_array = zeros(number_of_additional_constraints,number_of_reactions)
    for (reaction_key, local_measurement_dict) in flux_ratio_dictionary_array

        # ok, so we have a reaction key - find the index of this reaction in the reaction list -
        idx_reaction_match = findall(rxn_name_list .== reaction_key)

        # if we have this reaction, then update the bounds array -
        if (isempty(idx_reaction_match) == false)

            # all fluxes are relative to glucose -
            glucose_uptake_index = 50   # GLCpts -

            # get ratio value -
            ratio_value = parse(Float64,local_measurement_dict["mean_value"])

            # get the actual index -
            idx_reaction = (getindex(idx_reaction_match))[1]

            # update the constraint -
            additional_constraint_array[constraint_counter,idx_reaction] = 1.0
            additional_constraint_array[constraint_counter,glucose_uptake_index] = -1*ratio_value

            # update -
            constraint_counter = constraint_counter + 1
        end
    end

    # cache the additional constraints -
    data_dictionary["additional_constraint_array"] = additional_constraint_array

    # we need to add some additional "species bounds" -
    tmp_species_bounds_array = data_dictionary["species_bounds_array"]
    additional_constraint_block = zeros(number_of_additional_constraints,2)
    species_bounds_array = [tmp_species_bounds_array ; additional_constraint_block]
    data_dictionary["species_bounds_array"] = species_bounds_array

    # iterate through the measurements_dictionary and set values -
    individual_measuerment_dictionaries = measurements_dictionary["exchange_flux_measurements"]
    for (reaction_key, local_measurement_dict) in individual_measuerment_dictionaries

        # ok, so we have a reaction key - find the index of this reaction in the reaction list -
        idx_reaction_match = findall(rxn_name_list .== reaction_key)

        # if we have this reaction, then update the bounds array -
        if (isempty(idx_reaction_match) == false)

            # get the measured value -
            measured_value = parse(Float64,local_measurement_dict["mean_value"])

            # get the actual index -
            idx_reaction = (getindex(idx_reaction_match))[1]

            # what is the directionality of this measurement?
            directionality = Symbol(local_measurement_dict["directionality"])
            if directionality == :input

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = 0.0

            elseif directionality == :output

                flux_bounds_array[idx_reaction,1] = 0
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :free

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :input_fixed

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = -1*measured_value

            elseif directionality == :output_fixed

                flux_bounds_array[idx_reaction,1] = measured_value
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :nil

                flux_bounds_array[idx_reaction,1] = -1e-6
                flux_bounds_array[idx_reaction,2] = 1e-6

            else
                # TODO: issue a warning ...
            end
        end
    end

    # update -
    data_dictionary["flux_bounds_array"] = flux_bounds_array

    # return modified dictionary -
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function generate_exchange_flux_index_array(reaction_name_array::Array{String,1},pattern::String)

    # initialize -
    exchange_flux_index_array = Int64[]

    # we are looking for a name starting w/EX_
    number_of_reactions = length(reaction_name_array)
    for reaction_index = 1:number_of_reactions

        # grab the name -
        reaction_name = reaction_name_array[reaction_index]

        # check -
        if (startswith(reaction_name,pattern) == true)
            push!(exchange_flux_index_array,reaction_index)
        end
    end

    # return -
    return exchange_flux_index_array
end


"""
TODO: Fill me in with some stuff ...
"""
function generate_default_data_dictionary(path_to_cobra_mat_file::String,model_name::String, organism_id::Symbol)

    # Wher are we?
    top_level_path = pwd()

    # load the biophysical_constants dictionary -
    default_biophysical_dictionary = load_default_biophysical_dictionary(organism_id)

    # TODO: check if string is a legit path to the cobra file, and the model name is ok
    # load the *.mat code from the cobra code folder -
    file = matopen(path_to_cobra_mat_file)
    cobra_dictionary = read(file,model_name)
    close(file)

    # Setup: the stoichiometric matrix -
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # get the species symbol list -
    list_of_metabolite_symbols = cobra_dictionary["mets"]

    # get list of gene symbols -
    list_of_gene_symbols = cobra_dictionary["genes"]

    # setup the objective_coefficient_array -
    objective_coefficient_array = cobra_dictionary["c"]

    # Setup: default flux bounds array, including the "exhange fluxes"
    default_flux_bounds_array = zeros(number_of_reactions,2)
    lb = cobra_dictionary["lb"] # lower bound -
    ub = cobra_dictionary["ub"] # upper bound -
    for reaction_index = 1:number_of_reactions
        default_flux_bounds_array[reaction_index,1] = lb[reaction_index]
        default_flux_bounds_array[reaction_index,2] = ub[reaction_index]
    end

    # calculate the default vamx (based upon biophysical_constants/literature) -
    (default_vmax, default_enzyme_concentration) = calculate_default_vmax(default_biophysical_dictionary)

    # load the ec_number file -
    path_to_ec_file::String = top_level_path*"/cobra/config/data/ec_numbers.dat"
    ec_number_dictionary =  load_ec_mapping_file(path_to_ec_file)

    # load the gene order array -
    path_to_gene_file = top_level_path*"/cobra/config/data/gene_list.dat"
    gene_order_array = load_gene_order_file(path_to_gene_file)

    # calculate the vamx array -
    path_to_brenda_data = top_level_path*"/cobra/config/data/ECN-EcoliCoreModel-Palsson-EcoSal-Plus-2013.json"
    total_vmax_array = generate_vmax_array(path_to_brenda_data, default_biophysical_dictionary, gene_order_array, ec_number_dictionary)

    # correct the total_vmax_array for this model, using the rules -
    local_data_dictionary = Dict{String,Any}()
    local_data_dictionary["default_vmax_value"] = default_vmax
    model_vmax_array = calculate_rules_vector(local_data_dictionary, total_vmax_array)

    # update the default bounds array w/our "default" biophysical_constants -
    flux_bounds_array = update_default_flux_bounds_array(default_flux_bounds_array, model_vmax_array, cobra_dictionary)

    # species bounds array - default, everything is bounded 0,0
    # species in the [e] (extracellular) compartment are unbounded
    # we have *updated* the STM w/exchange reactions [] -> a and a -> []
    # all species are bounded to 0
    species_bounds_array = zeros(number_of_species,2)

    # create list of reaction strings -
    list_of_reaction_name_strings = cobra_dictionary["rxns"]

    # What sense do we do? (by default we min)
    is_minimum_flag = true

    # construct the reaction string list -
    list_of_chemical_reaction_strings = reconstruct_reaction_string_list(cobra_dictionary)

    # grab the reversible reaction list -
    reversible_reaction_flag_array = cobra_dictionary["rev"]

    # calculate the reaction name -
    reaction_name_array_tmp = cobra_dictionary["rxns"]
    reaction_name_array = String[]
    for rxn_name in reaction_name_array_tmp
        push!(reaction_name_array,rxn_name)
    end
    exchange_flux_index_array = generate_exchange_flux_index_array(reaction_name_array,"EX_")
    exchange_flux_index_array = [exchange_flux_index_array ; 13]    # add the growth rate -



    # =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
    data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
    data_dictionary["flux_bounds_array"] = flux_bounds_array
    data_dictionary["species_bounds_array"] = species_bounds_array
    data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
    data_dictionary["list_of_gene_symbols"] = list_of_gene_symbols
    data_dictionary["list_of_reaction_name_strings"] = list_of_reaction_name_strings
    data_dictionary["list_of_chemical_reaction_strings"] = list_of_chemical_reaction_strings
	data_dictionary["is_minimum_flag"] = is_minimum_flag
    data_dictionary["reversible_reaction_flag_array"] = reversible_reaction_flag_array
    data_dictionary["exchange_flux_index_array"] = exchange_flux_index_array

    # stuff for rules -
    data_dictionary["default_vmax"] = default_vmax

    # in case we need something later -
    data_dictionary["cobra_dictionary"] = cobra_dictionary
    # ========================================================================================= #
    return data_dictionary
end

# - PRIVATE HELPER METHODS -------------------------------------------------------------------- #
function load_gene_order_file(path_to_gene_file)

    # gene array -
    gene_symbol_array = String[]

    # open the gene file -
    # TODO: We need to check to see if this is a legit path -
    open(path_to_gene_file) do f

        # load the file into the buffer -
        buffer = read(f, String)

        # split along new line -
        list_of_records = split(buffer,'\n')

        # how many records do we have?
        number_of_records = length(list_of_records)
        for record_index = 1:number_of_records

            # convert the record to a string -
            record = string(list_of_records[record_index])

            # is the record empty?
            if (isempty(record) == false)

                # split into fields -
                field_array = split(record,',')

                # grab the ec number -
                order_index = field_array[1]
                gene_symbol = field_array[2]

                # add a record -
                push!(gene_symbol_array, gene_symbol)
            end
        end
    end

    # return -
    return gene_symbol_array
end

function load_ec_mapping_file(path_to_ec_file::String)

    # initialize the ec number array -
    ec_number_dictionary = Dict{String,String}()

    # open the ec file -
    open(path_to_ec_file) do f

        # load the file into the buffer -
        buffer = read(f, String)

        # split along new line -
        list_of_records = split(buffer,'\n')

        # how many records do we have?
        number_of_records = length(list_of_records)
        for record_index = 1:number_of_records

            # convert the record to a string -
            record = string(list_of_records[record_index])

            # is the record empty?
            if (isempty(record) == false)

                # split into fields -
                field_array = split(record,'\t')

                # grab the ec number -
                bg_number = field_array[1]
                ec_number = field_array[2]

                # add a record -
                ec_number_dictionary[bg_number] = ec_number

            end
        end
    end

    # return -
    return ec_number_dictionary
end

function generate_vmax_array(path_to_brenda_data::String, biophysical_dictionary::Dict{String,Any}, gene_symbol_array::Array{String,1}, ec_number_dictionary::Dict{String,String})

    # TODO: is the path_to_brenda_data legit?
    kcat_measurements_dictionary = JSON.parsefile(path_to_brenda_data)

    # initialize -
    vmax_array = Float64[]

    # calculate a default value -
    (default_vmax, default_enzyme_concentration) = calculate_default_vmax(biophysical_dictionary)

    # for each gene symbol, lookup an ec number, then get a range of measured kcats -
    for gene_symbol in gene_symbol_array

        # do we have this gene key?
        if (haskey(ec_number_dictionary, gene_symbol) == true)

            # what is the ec number for this key?
            ec_number_key = ec_number_dictionary[gene_symbol]

            # ok, so we have an ec number. Do we have a kcat measurement for this ec number?
            if (haskey(kcat_measurements_dictionary, ec_number_key) == true)


                # ok, we have an ec_number w/data - get the data
                # this will come back as an array of dictionaries
                local_vmax_array = Float64[]
                array_data_dictionaries = kcat_measurements_dictionary[ec_number_key]
                for local_data_dictionary::Dict{String,Any} in array_data_dictionaries

                    # get the turnoverNumber -
                    turnoverNumber = parse(Float64,local_data_dictionary["turnoverNumber"])*(3600) # hr^-1

                    # there is a "strange" data issue w/the Brenda data. They use -999 to mean no value
                    if (turnoverNumber>0)
                        # cache -
                        push!(local_vmax_array, turnoverNumber*default_enzyme_concentration)
                    end
                end

                if isempty(local_vmax_array) == true
                    push!(local_vmax_array, default_vmax)
                end

                # what is the min/max -
                min_value = minimum(local_vmax_array)
                max_value = maximum(local_vmax_array)

                # calculate a value (uniform random between the measurements)
                sampled_vmax_value = min_value + (max_value - min_value)*rand()

                # yes - we have an experimental value -
                push!(vmax_array, sampled_vmax_value)
            else

                # nope - bummer. Use the default value -
                push!(vmax_array, default_vmax)
            end
        else

            # nope - bummer. Use the default value -
            push!(vmax_array, default_vmax)
        end
    end

    # return -
    return vmax_array
end

function calculate_default_vmax(biophysical_dictionary::Dict{String,Any})

    # ok, so we need to get some stuff from the dictionary -
    # TODO: Check these keys are contained in the dictionary
    default_turnover_number = parse(Float64,biophysical_dictionary["biophysical_constants"]["default_turnover_number"]["value"])              # convert to h^-1
    default_enzyme_concentration = parse(Float64,biophysical_dictionary["biophysical_constants"]["default_enzyme_concentation"]["value"])     # mumol/gDW

    # calculate the default VMax -
    default_vmax = (default_turnover_number)*(default_enzyme_concentration)*(3600)

    return (default_vmax, default_enzyme_concentration)
end

function update_default_flux_bounds_array(flux_bounds_array::Array{Float64,2}, vmax_array::Array{Float64,1}, cobra_dictionary::Dict{String,Any})

    # get the reverse flag array from the cobra dictionary -
    rev_flag_array = cobra_dictionary["rev"]

    # how many bounds do we have?
    (number_of_bounds, number_of_columns) = size(flux_bounds_array)

    # initialize and update the array -
    updated_flux_bounds_array = zeros(number_of_bounds,number_of_columns)
    for flux_index = 1:number_of_bounds

        # get vmax value -
        vmax_value = vmax_array[flux_index]

        # setup new bounds -
        updated_flux_bounds_array[flux_index,1] = 0.0
        updated_flux_bounds_array[flux_index,2] = vmax_value
        if (rev_flag_array[flux_index] == 1.0)
            updated_flux_bounds_array[flux_index,1] = -1*vmax_value
        end
    end

    # return -
    return updated_flux_bounds_array
end

function load_default_biophysical_dictionary(organism_id::Symbol)

    # load the organism_id mapping file -
    top_level_path = pwd()
    path_to_mapping_file = top_level_path*"/flux/config/Mapping.json"

    # TODO: Check do we have a mapping file?
    organism_id_map = JSON.parsefile(path_to_mapping_file)

    # TODO: Throw an error if we are missing this organism
    local_path_to_biophysical_dictionary = organism_id_map[String(organism_id)]["path_to_default_constants"]

    # TODO: check - does this path point to a file?
    global_path_to_biophysical_dictionary = "$(pwd())"*local_path_to_biophysical_dictionary
    default_biophysical_dictionary = JSON.parsefile(global_path_to_biophysical_dictionary)

    # return -
    return default_biophysical_dictionary
end
