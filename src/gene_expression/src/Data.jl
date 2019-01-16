
"""
TODO: Fill me in with some stuff ...
"""
function build_biophysical_dictionary(path_to_biophysical_dictionary::String)

    # TODO: Check to see if the path points to a legit JSON file -
    # ...

    # load the database file w/the species information in it -
    raw_biophysical_dict = JSON.parsefile(path_to_biophysical_dictionary)

    # remove the top level dictionary -
    biophysical_data_dictionary = raw_biophysical_dict["biophysical_constants"]

    # -------------------------------------------------------------------------------------------------------------------------------------- #
    # alias the constants, so we can do some computation on them -
    av_number = 6.02e23                                                                                                                # copies/mol
    cell_diameter = parse(Float64,raw_biophysical_dict["cell_diameter"]["value"])                                                      # mum
	mass_of_single_cell = parse(Float64,raw_biophysical_dict["mass_of_single_cell"]["value"])                                          # g/cell
	copy_number_of_rnapII = parse(Float64,raw_biophysical_dict["copies_of_rnapII_per_cell"]["value"])                                  # copies/cells
    copy_number_of_ribosome = parse(Float64,raw_biophysical_dict["copies_of_ribosome_per_cell"]["value"])                              # copies/cells
    mRNA_half_life_TF = parse(Float64,raw_biophysical_dict["average_mRNA_half_life"]["value"])                                         # hr
    protein_half_life = parse(Float64,raw_biophysical_dict["average_protein_half_life"]["value"])                                      # hr
    infrastructure_half_life = parse(Float64,raw_biophysical_dict["average_infrastructure_half_life"]["value"])                        # hr
    doubling_time_cell = parse(Float64,raw_biophysical_dict["doubling_time_cell"]["value"])                                            # hr
    max_translation_rate = parse(Float64,raw_biophysical_dict["translation_elongation_rate"]["value"])                                 # aa/sec
	max_transcription_rate = parse(Float64,raw_biophysical_dict["transcription_elongation_rate"]["value"])                             # nt/sec
	transcription_initiation_time_contstant = parse(Float64,raw_biophysical_dict["transcription_initiation_time_contstant"]["value"])  # sec
	average_transcript_length = parse(Float64,raw_biophysical_dict["average_transcript_length"]["value"])                              # nt
	average_protein_length = parse(Float64,raw_biophysical_dict["average_protein_length"]["value"])                                    # aa
	fraction_nucleus = parse(Float64,raw_biophysical_dict["fraction_nucleus"]["value"])                                                # dimensionless
	avg_gene_number = parse(Float64,raw_biophysical_dict["average_gene_copies_per_cell"]["value"])                                     # copies/cell
	avg_polysome_number = parse(Float64,raw_biophysical_dict["average_polysome_number"]["value"])                                      # ribsomoses/per transcript
    copy_number_of_basal_TF = parse(Float64,raw_biophysical_dict["copies_of_basal_TF_per_cell"]["value"])                              # copies/cell
    default_enzyme_concentation = parse(Float64,raw_biophysical_dict["default_enzyme_concentation"]["value"])                          # mumol/gDW
    default_turnover_number = parse(Float64,raw_biophysical_dict["default_turnover_number"]["value"])                                  # s^-1

    # do some calculations ...

    # Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

    # concentration scale factor -
    CSF = 1e6

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = copy_number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*CSF       # umol/gdw
	ribosome_concentration = copy_number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*CSF   # umol/gdw

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                                  # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                               # hr^-1
	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(0.5)			         # hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)                # hr^-1
	kcat_translation = avg_polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1       # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                              # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                                # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*CSF                  # umol/gdw

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                     # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 91*(1/av_number)*(1/mass_of_single_cell)*CSF                     # umol/gdw
	saturation_translation = 160000*(1/av_number)*(1/mass_of_single_cell)*CSF                   # umol/gdw

    # what is the sigma70 concentration?
    basal_TF_concentration = copy_number_of_basal_TF*(1/av_number)*(1/mass_of_single_cell)*CSF  # umol/gdw
	# -------------------------------------------------------------------------------------------------------------------------------------- #

	# Package the txtl parameters -
	txtl_parameter_dictionary = Dict{AbstractString,Any}()
	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	txtl_parameter_dictionary["degrdation_constant_infrastructure"] = degrdation_constant_infrastructure  # hr^-1
	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	txtl_parameter_dictionary["death_rate_constant"] = death_rate_constant
	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
	txtl_parameter_dictionary["average_transcript_length"] = average_transcript_length
	txtl_parameter_dictionary["average_protein_length"] = average_protein_length
	txtl_parameter_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation
	txtl_parameter_dictionary["kcat_translation_initiation"] = kcat_translation_initiation
    txtl_parameter_dictionary["basal_TF_concentration"] = basal_TF_concentration
    txtl_parameter_dictionary["default_enzyme_concentation"] = default_enzyme_concentation
    txtl_parameter_dictionary["default_turnover_number"] = default_turnover_number

    # return -
    return txtl_parameter_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function build_transcription_control_dictionary(configuration_dictionary::Dict{String,Any})




    # Initialize -
    control_structure_dictionary = Dict{AbstractString,TranscriptionalControlModel}()
    gene_species_array = String[]

    # setup the paths -
    model_directory_name = configuration_dictionary["model_directory_name"]
    model_files_path = "$(pwd())/models/$(model_directory_name)"

    # load the mapping array -
    gene_mapping_array_filename = configuration_dictionary["gene_mapping_array"]
    gene_mapping_array = readdlm("$(model_files_path)/$(gene_mapping_array_filename)")

    # get the gene_control_models -
    gene_model_array = configuration_dictionary["gene_control_models"]
    gene_index = 1
    for gene_model in gene_model_array

        # Build the control model -
        gene_control_model = TranscriptionalControlModel()

        # build the actor index array -
        idx_array = findall(gene_mapping_array[gene_index,:].==1)
        gene_control_model.transcriptional_actor_index_array = idx_array

        # Get the control and parameter arrays -
        gain_array_filename = gene_model["control_array"]
        actor_file_name = gene_model["actor_parameter_array"]

        # load the arrays -
        gain_array = readdlm("$(model_files_path)/$(gain_array_filename)")
        actor_parameter_array = readdlm("$(model_files_path)/$(actor_file_name)")

        # cache -
        gene_control_model.transcriptional_weight_matrix = gain_array
        gene_control_model.transcriptional_parameter_matrix = actor_parameter_array

        # grab the lengths -
        gene_control_model.aa_seq_length = parse(Float64,gene_model["aa_seq_length"])
        gene_control_model.nt_seq_length = parse(Float64,gene_model["nt_seq_length"])

        # get the id -
        gene_id = gene_model["id"]

        # cache -
        control_structure_dictionary[gene_id] = gene_control_model

        # grab the species -
        push!(gene_species_array, gene_id)

        # update index -
        gene_index = gene_index + 1
    end

    # return -
    return (control_structure_dictionary, gene_species_array)
end

function build_steady_state_data_dictionary(path_to_json_config_file::String)

    # = 0: Initialize  =========================================================================== %

    # TODO: Need to check the path to see if this is a legit file.
    # if not - then throw an error

    # load the database file w/the species information in it -
    raw_config_dict = JSON.parsefile(path_to_json_config_file)
    # ============================================================================================ %

    # Load the biophysical_dictionary -
    biophysical_dictionary_file_name = raw_config_dict["biophysical_model_filename"]
    path_to_biophysical_dictionary = "$(pwd())/gene_expression/config/constants/$(biophysical_dictionary_file_name)"
    biophysical_dictionary = build_biophysical_dictionary(path_to_biophysical_dictionary)

    # Construct the system control structure dictionary -
    (control_structure_dictionary, gene_species_array) = build_transcription_control_dictionary(raw_config_dict)


    # = 2: Setup the sytem size ================================================================== %
    number_of_internal_species = parse(Int,raw_config_dict["number_of_internal_species"])
    number_of_external_species = parse(Int,raw_config_dict["number_of_external_species"])
    transcript_range_array = 1:4
    protein_range_array = 5:8
    # ============================================================================================ %

    # = 3: Bounds ================================================================================ %
    lower_bound_array = zeros(number_of_internal_species)
    upper_bound_array = 100*ones(number_of_internal_species)
    # ============================================================================================ %

    # = 4: External vector values ================================================================ %
    external_state_vector = [
    ]
    # ============================================================================================ %

    # ======== DO NOT EDIT BELOW THIS LINE ======================================================= %
    data_dictionary = Dict()
    data_dictionary["control_structure_dictionary"] = control_structure_dictionary
    data_dictionary["gene_species_array"] = gene_species_array



    data_dictionary["number_of_internal_species"] = number_of_internal_species
    data_dictionary["number_of_external_species"] = number_of_external_species
    data_dictionary["transcript_range_array"] = transcript_range_array
    data_dictionary["protein_range_array"] = protein_range_array
    data_dictionary["lower_bound_array"] = lower_bound_array
    data_dictionary["upper_bound_array"] = upper_bound_array
    data_dictionary["external_state_vector"] = external_state_vector
    return data_dictionary
    # ============================================================================================ %
end
