function find_species_in_compartment(list_of_species,compartment_string)

    # initialize -
    index_array = Int64[]

    # main -
    counter = 1
    for species_symbol in list_of_species

        if (occursin(compartment_string,species_symbol) == true)
            push!(index_array,counter)
        end

        # update the counter -
        counter = counter + 1
    end

    # return -
    return index_array
end

function reconstruct_reaction_string_list(cobra_dictionary)

    # initialize -
    reaction_string_buffer = String[]

    # get the stoichiometric array -
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    list_of_reaction_tags = cobra_dictionary["rxns"]
    list_of_species = cobra_dictionary["mets"]
    list_of_reversible_reactions = cobra_dictionary["rev"]

    # what is the size?
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)
    for reaction_index = 1:number_of_reactions

        # initialize empty buffer -
        buffer = ""

        # get the reaction tag -
        reaction_tag_string = list_of_reaction_tags[reaction_index]

        # add the tag to the buffer -
        buffer *= "$(reaction_index),$(reaction_tag_string),"

        # find the reactants -
        idx_reactants = findall(stoichiometric_matrix[:,reaction_index].<0.0)
        if (isempty(idx_reactants) == true)
            buffer *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_reactants)
            counter = 1
            for index in idx_reactants

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != -1.0)
                    # add data to the buffer -
                    buffer *= "$(abs(stcoeff))*$(metabolite_string)"
                else
                    # add data to the buffer -
                    buffer *= "$(metabolite_string)"
                end



                # do we have more?
                if (counter < number_of_species)
                    buffer *= "+"
                else
                    buffer *= ","
                end

                counter = counter + 1
            end
        end

        # find the products -
        idx_products = findall(stoichiometric_matrix[:,reaction_index].>0.0)
        if (isempty(idx_products) == true)
            buffer *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_products)
            counter = 1
            for index in idx_products

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != 1.0)
                    # add data to the buffer -
                    buffer *= "$(stcoeff)*$(metabolite_string)"
                else
                    # add data to the buffer -
                    buffer *= "$(metabolite_string)"
                end

                # do we have more?
                if (counter < number_of_species)
                    buffer *= "+"
                else
                    buffer *= ","
                end

                counter = counter + 1
            end
        end

        # is this reaction reversible?
        rev_value = list_of_reversible_reactions[reaction_index]
        if (rev_value == 1.0)
            buffer *= "-inf,inf"
        else
            buffer *= "0,inf"
        end

        # add buffer to list of strings -
        push!(reaction_string_buffer,buffer)
    end

    # return reaction_string_buffer
    return reaction_string_buffer
end
