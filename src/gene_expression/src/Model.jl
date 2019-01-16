
# - PUBLIC FUNCTIONS -------------------------------------------------------------------------------------------------------------- #
function solve_steady_state_system(ic_array::Array{Float64,1},data_dictionary::Dict{String,Any})

    # setup the call to the optimizer --
    objective_function(x) = evaluate_steady_state_balances(x,data_dictionary)

    # get the lower, and upper_bounds for the states -
    lower_bound_array = data_dictionary["lower_bound_array"]
    upper_bound_array = data_dictionary["upper_bound_array"]

    # how many species do we have?
    number_of_species = data_dictionary["number_of_internal_species"]


    # make a call to the optim package -
    result = optimize(objective_function,lower_bound_array, upper_bound_array, ic_array,Fminbox(LBFGS()))

    # get the extent -
    steady_state_vector = Optim.minimizer(result)

    # return -
    return steady_state_vector
end
# -------------------------------------------------------------------------------------------------------------------------------- #

# - PRIVATE FUNCTIONS ------------------------------------------------------------------------------------------------------------ #
function evaluate_steady_state_balances(x,data_dictionary::Dict{String,Any})

    # pull the system array from the data_dictionary -
    G = build_system_gain_matrix(x,data_dictionary)

    # calculate the system rate vector -
    v = calculate_system_control_vector(x,data_dictionary)

    # what is the residual -
    residual = (x - G*v)

    # calculate the error -
    number_of_states = length(residual)
    transcript_range_array = data_dictionary["transcript_range_array"]
    W = Matrix{Float64}(I, number_of_states, number_of_states)
    for transcript_index in transcript_range_array
        W[transcript_index,transcript_index] = 10000
    end

    # compute the error, weighted by W matrix -
    error = transpose(residual)*W*residual

    # return -
    return error
end

function build_biophysical_dictionary()

    # ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)                   units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                             # mum
	mass_of_single_cell = 2.8e-13                   # g
	number_of_rnapII = 4600            	            # copies/cells
	number_of_ribosome = 50000         	            # copies/cells
	mRNA_half_life_TF = 0.016                       # hrs
	protein_half_life = 2.0                         # hrs
	infrastructure_half_life = 300					# hrs
	doubling_time_cell = 0.66                       # hrs
	max_translation_rate = 16.5                     # aa/sec
	max_transcription_rate = 60.0                   # nt/sec
	transcription_initiation_time_contstant = 400   # sec
	average_transcript_length = 1000   	            # nt
	average_protein_length = 330       	            # aa
	fraction_nucleus = 0.0             	            # dimensionless
	av_number = 6.02e23                             # number/mol
	avg_gene_number = 2                             # number of copies of a gene
	polysome_number = 4					            # number of ribsomoses per transcript
    number_of_basal_TF = 5700                       # copies per cell
	# ------------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # umol/gdw
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # umol/gdw

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                           # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                        # hr^-1
	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(0.5)			# hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                          # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*1e6             # umol/gdw

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 91*(1/av_number)*(1/mass_of_single_cell)*1e6               # umol/gdw
	saturation_translation = 160000*(1/av_number)*(1/mass_of_single_cell)*1e6                # umol/gdw

    # what is the sigma70 concentration?
    basal_TF_concentration = number_of_basal_TF*(1/av_number)*(1/mass_of_single_cell)*1e6    # umol/gdw
	# -------------------------------------------------------------------------------------------------#

	# Alias the txtl parameters -
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

    # return -
    return txtl_parameter_dictionary
end

function build_transcription_gain_matrix(x,data_dictionary::Dict{String,Any})

    # I need to sort this out ... but for now, let's just use hardcoded values -
    biophysical_dictionary = build_biophysical_dictionary()
    kcat_transcription_initiation = biophysical_dictionary["kcat_transcription_initiation"]
    avg_gene_concentration = biophysical_dictionary["avg_gene_concentration"]
    saturation_transcription = biophysical_dictionary["saturation_constant_transcription"]
    rnapII_concentration = biophysical_dictionary["rnapII_concentration"]
    basal_TF_concentration = biophysical_dictionary["basal_TF_concentration"]
    mugmax = biophysical_dictionary["maximum_specific_growth_rate"]
    kcat_transcription = biophysical_dictionary["kcat_transcription"]
    kdTX = biophysical_dictionary["degradation_constant_mRNA"]
    avg_nt_length = biophysical_dictionary["average_transcript_length"]

    # get gene models -
    control_structure_dictionary = data_dictionary["control_structure_dictionary"]

    # get the species array -
    species_array = data_dictionary["gene_species_array"]
    number_of_genes = length(species_array)
    tmp_array = Float64[]
    for species_id in species_array

        # control model -
        control_model = control_structure_dictionary[species_id]

        # compute length factor -
        length_factor = (avg_nt_length/control_model.nt_seq_length)

        # calculate the gain -
        kcat = (kcat_transcription*length_factor*kcat_transcription_initiation)/((kcat_transcription*length_factor+mugmax)*(mugmax+kdTX))
        gain_TX = kcat*(rnapII_concentration)*((avg_gene_concentration)/(saturation_transcription+avg_gene_concentration))

        # cache the gain -
        push!(tmp_array,gain_TX)
    end

    # compute return the gain array -
    return (Matrix{Float64}(I, number_of_genes, number_of_genes)).*tmp_array
end

function build_translation_gain_matrix(x,data_dictionary::Dict{String,Any})

    # I need to sort this out ... but for now, let's just use hardcoded values -
    biophysical_dictionary = build_biophysical_dictionary()
    kcat_translation = biophysical_dictionary["kcat_translation"]
    saturation_translation = biophysical_dictionary["saturation_constant_translation"]
    ribosome_concentration = biophysical_dictionary["ribosome_concentration"]
    mugmax = biophysical_dictionary["maximum_specific_growth_rate"]
    kcat_translation = biophysical_dictionary["kcat_translation"]
    kdTL = biophysical_dictionary["degradation_constant_protein"]
    avg_aa_seq_len = biophysical_dictionary["average_protein_length"]

    # get the range of mRNA -
    mRNA_range = data_dictionary["transcript_range_array"]
    mRNA_concentration_vector = x[collect(mRNA_range)]

    # get gene models -
    control_structure_dictionary = data_dictionary["control_structure_dictionary"]

    # initialize tmp array -
    tmp_array = Float64[]
    number_of_mRNA = length(mRNA_range)
    species_array = data_dictionary["gene_species_array"]
    for species_index = 1:number_of_mRNA

        # control model -
        control_model = control_structure_dictionary[species_array[species_index]]

        # compute the length factor -
        length_factor = (avg_aa_seq_len/control_model.aa_seq_length)

        # compute the mRNA saturation term -
        mhat = (mRNA_concentration_vector[species_index])/(saturation_translation+mRNA_concentration_vector[species_index])

        # compute the gain -
        gain_TL = 1.0*((kcat_translation*mhat*length_factor*ribosome_concentration)/(mugmax+kdTL))

        # cache -
        push!(tmp_array,gain_TL)
    end

    # compute return the gain array -
    return (Matrix{Float64}(I, number_of_mRNA, number_of_mRNA)).*tmp_array
end

function build_system_gain_matrix(x,data_dictionary::Dict{String,Any})

    # calculate the TX gain matrix -
    transcription_gain_matrix = build_transcription_gain_matrix(x,data_dictionary)
    translation_gain_matrix = build_translation_gain_matrix(x,data_dictionary)

    # what is the size of the TL/TX array's?
    (number_of_genes,n_TX_col) = size(transcription_gain_matrix)
    (number_of_proteins,n_TL_col) = size(translation_gain_matrix)

    # create a ZERO block to go after the TX and TL arrays -
    TX_zero_block = zeros(number_of_genes,n_TL_col)
    TL_zero_block = zeros(number_of_proteins,n_TX_col)

    # create the system gain matrix -
    system_gain_matrix = [
        transcription_gain_matrix TX_zero_block ;
        TL_zero_block translation_gain_matrix ;
    ]

    # return -
    return system_gain_matrix
end
# -------------------------------------------------------------------------------------------------------------------------------- #
