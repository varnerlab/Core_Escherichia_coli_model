using JSON

# script to parse Brenda files -
function parse_brenda_file(path_to_brenda_file)

    # iniialize -
    brenda_records_array = Dict{String,Any}[]

    # process the file -
    open(path_to_brenda_file) do f

        # do stuff with the open file
        buffer = readstring(f)
        list_of_records = split(buffer,'!')


        # how many records do we have?
        number_of_records = length(list_of_records)
        for record_index = 1:number_of_records

            # initalize empty dictionary -
            record_dictionary = Dict{String,Any}()

            # convert the record to a string -
            record = string(list_of_records[record_index])

            # split into fields -
            field_array = split(record,'#')
            number_of_fields = length(field_array)
            for field_index = 1:number_of_fields

                # split into key*value pairs -
                key_value_array = split(field_array[field_index],'*')

                # split again if we have two elements -
                if (length(key_value_array) == 2)

                    # split into key,value pair
                    key = string(key_value_array[1])
                    value = string(key_value_array[2])

                    # cache the values -
                    record_dictionary[key] = value

                end
            end

            # add this record to the array -
            push!(brenda_records_array,record_dictionary)
        end
    end

    # return the array -
    return brenda_records_array
end

function main(path_to_model_files,path_to_json_file)

    # initialize -
    brenda_dictionary = Dict{String,Any}()
    brenda_records_array = Dict{String,Any}[]

    # brenda files end w/.dat
    for file in filter(x -> endswith(x, "dat"), readdir(path_to_model_files))

        # path to specific file -
        path_to_brenda_file = "$(path_to_model_files)/$(file)"

        # parse the brenda file,
        brenda_records_array = parse_brenda_file(path_to_brenda_file)

        # grab the first record, and get the ecnumber -
        ec_number = brenda_records_array[1]["ecNumber"]
        ec_number_key = "ec:$(ec_number)"

        # cache this records collection -
        brenda_dictionary[ec_number_key] = brenda_records_array
    end

    # write dictionary to JSON -
    file = open(path_to_json_file,"w")
    JSON.print(file,brenda_dictionary,6)
    close(file)
end

# build the dictionary -
main("./brenda_data","./ECNDatabase-CCM-Palsson-SciReports-2017.json")
