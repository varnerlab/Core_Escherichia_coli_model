function test_rules_generator(organism_id::Symbol)

    # setup the path -
    path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
    model_name = "modelReg"

    # load the data data_dictionary -
    dd = generate_default_data_dictionary(path_to_cobra_mat_file,model_name,organism_id)

    # get the cobra dictionary -
    cdb = dd["cobra_dictionary"]

    # generate rules -
    generate_rules_function(cdb,"./test/tmp/Rules.jl")

end
