# - LOAD SYSTEM PACKAGES ------------------------------------------- #
using Optim
using LinearAlgebra
using JSON
using DelimitedFiles
# ------------------------------------------------------------------ #

# - INCLUDE CODE IN THIS MODULE ------------------------------------ #
top_level_path = pwd()
include("$(top_level_path)/gene_expression/src/Control.jl")
include("$(top_level_path)/gene_expression/src/Types.jl")
include("$(top_level_path)/gene_expression/src/Data.jl")
include("$(top_level_path)/gene_expression/src/Model.jl")
# ------------------------------------------------------------------ #
