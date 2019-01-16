mutable struct Gene

    # information about a Gene -
    symbol::AbstractString
    organism_code::AbstractString


    function Gene()
        this = new()
    end
end

mutable struct Promoter

    # information about a Gene -
    symbol::AbstractString
    organism_code::AbstractString


    function Promoter()
        this = new()
    end
end
