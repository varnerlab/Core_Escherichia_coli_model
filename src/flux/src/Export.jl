function write_reactions_to_file(data_dictionary::Dict{String,Any}, path_to_reaction_file::String)

    reaction_list = data_dictionary["list_of_chemical_reaction_strings"]
    writedlm(path_to_reaction_file,reaction_list)

end
