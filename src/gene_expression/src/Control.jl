function calculate_transcription_transfer_function_vector(x,binding_parameter_matrix)

    # how many factors do we have?
    number_of_effectors = length(x)

    # initialize -
    f_vector = Float64[]

    # main loop -
    for effector_index = 1:number_of_effectors

        # get the parameters for this function -
        K = binding_parameter_matrix[effector_index,1]
        n = binding_parameter_matrix[effector_index,2]

        # get the effector value -
        x_value = x[effector_index]

        # compute -
        f_value = (x_value^n)/(K + x_value^n)

        # cache -
        push!(f_vector,f_value)
    end

    # return -
    return f_vector
end

function calculate_transcription_control_vector(x,data_dictionary::Dict{String,Any})

    # How many control terms are we going to have?
    gene_species_array = data_dictionary["gene_species_array"]
    number_of_genes = length(gene_species_array)
    control_array = Float64[]

    # get the gene_order dictionary -
    control_structure_dictionary = data_dictionary["control_structure_dictionary"]

    # b-vector -
    b = [0 ; 1]

    # iterate through the gene_data_dictionary -
    for gene_symbol in gene_species_array

        # get the gene model -
        gene_model = control_structure_dictionary[gene_symbol]

        # get the gain and parameter matrix -
        transcriptional_gain_matrix = gene_model.transcriptional_weight_matrix
        binding_parameter_matrix = gene_model.transcriptional_parameter_matrix

        # which range of effectors influence this gene?
        effector_range = gene_model.transcriptional_actor_index_array

        # grab the effector concentration -
        effector_concentration_vector = x[effector_range]

        # calculate the transfer function vector -
        f = calculate_transcription_transfer_function_vector(effector_concentration_vector,binding_parameter_matrix)
        f_vector = [1 ; f]

        # caclulate the control term -
        local_control_array = transcriptional_gain_matrix*f_vector + b
        control_value = local_control_array[1]/local_control_array[2]

        # cache -
        push!(control_array,control_value)
    end

    # return -
    return control_array
end

function calculate_translation_control_vector(x,data_dictionary::Dict{String,Any})

    # for now - assume no translational control -
    range_of_proteins = data_dictionary["protein_range_array"]

    # return a vectors of 1's -
    return ones(length(collect(range_of_proteins)))
end

function calculate_system_control_vector(x,data_dictionary::Dict{String,Any})

    # initialize control vector -
    control_vector = Float64[]

    # grab the external state vector, append to internal state vector -
    external_state_vector = data_dictionary["external_state_vector"]
    state_vector = [x ; external_state_vector]

    # 1: Calculate the u-variables (control the transcription rate) -
    u_variable_array = calculate_transcription_control_vector(state_vector,data_dictionary)
    for value in u_variable_array
        push!(control_vector, value)
    end

    # 2: Calculate the w-variables (control the translation rate) -
    w_variable_array = calculate_translation_control_vector(state_vector,data_dictionary)
    for value in w_variable_array
        push!(control_vector, value)
    end

    # return control vector -
    return control_vector
end
