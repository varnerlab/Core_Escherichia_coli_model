module VLCobraModule

    # include -
    top_level_path = pwd()
    include("$(top_level_path)/cobra/src/Include.jl")

    # exports -
    export download_reactions_from_kegg
    export write_ec_numbers_to_disk
    export write_protein_sequences_to_disk
    export write_gene_sequences_to_disk
    export write_gene_symbols_to_disk
    export generate_rules_function
end
