# - LOAD SYSTEM PACKAGES --------------------------------- #
using VLFluxModule
# -------------------------------------------------------- #

# - LOAD CODE IN THIS MODULE ----------------------------- #
top_level_path = pwd()
include("$(top_level_path)/test/src/DataDictionaryTests.jl")
include("$(top_level_path)/test/src/FluxTests.jl")
# -------------------------------------------------------- #
