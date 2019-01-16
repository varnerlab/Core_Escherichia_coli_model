function load_default_data_dictionary_test(organism_id)

    # setup the path -
    path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
    model_name = "modelReg"

    # load the data data_dictionary -
    dd = generate_default_data_dictionary(path_to_cobra_mat_file,model_name,organism_id)

    # return dd -
    return dd
end
