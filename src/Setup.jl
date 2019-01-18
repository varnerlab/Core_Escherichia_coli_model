include("Include.jl")


function main()

    # internal functions -
    isemptydir(dir::AbstractString) = isempty(readdir(dir))

    # load the data dictionary -
    path_to_cobra_mat_file::String = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
    model_name::String = "modelReg"
    dd = generate_default_data_dictionary(path_to_cobra_mat_file,model_name)

    # get the cobra dictionary -
    cobra_dictionary = dd["cobra_dictionary"]

    # Download a bunch of stuff from KEGG -
    kegg_organism_code::String = "eco"
    path_to_input_gene_file::String = "$(pwd())/cobra/config/data/gene_list.dat"

    # Step 1: Write the gene symbols to disk -
    path_to_output_gene_file::String = "$(pwd())/cobra/config/data"
    if (filesize("$(path_to_output_gene_file)/gene_list.dat")==0)
        write_gene_symbols_to_disk(cobra_dictionary,path_to_output_gene_file)
    end

    # Step 2: Write the EC numbers to disk -
    path_to_output_ec_file::String = "$(pwd())/cobra/config/data"
    if (filesize("$(path_to_output_ec_file)/ec_numbers.dat")==0)
        write_ec_numbers_to_disk(path_to_input_gene_file,path_to_output_ec_file,kegg_organism_code)
    end

    # Step 3: Write gene sequence to disk -
    path_to_output_gene_seq_files::String = "$(pwd())/cobra/config/data/gene_seq"
    if (isemptydir("$(path_to_output_gene_seq_files)") == true)
        write_gene_sequences_to_disk(path_to_input_gene_file,path_to_output_gene_seq_files,kegg_organism_code)
    end

    # Step 4: Write the protein sequence to disk -
    path_to_output_protein_seq_files::String = "$(pwd())/cobra/config/data/protein_seq"
    if (isemptydir("$(path_to_output_protein_seq_files)") == true)
        write_protein_sequences_to_disk(path_to_input_gene_file,path_to_output_protein_seq_files,kegg_organism_code)
    end

    # Step 5: Write the rules function -
    path_to_output_rules_function = "$(pwd())/flux/src"
    if (filesize("$(path_to_output_rules_function)/Rules.jl") == 0)
        generate_rules_function(cobra_dictionary, path_to_output_rules_function)
    end
end

main()
