# - LOAD SYSTEM PACKAGES --------------------------------- #
using VLFluxModule
using VLCobraModule
# -------------------------------------------------------- #

# - LOAD CODE IN THIS MODULE ----------------------------- #
top_level_path = pwd()
include("$(top_level_path)/test/src/DataDictionaryTests.jl")
include("$(top_level_path)/test/src/FluxTests.jl")
include("$(top_level_path)/test/src/RulesGeneratorTests.jl")
# -------------------------------------------------------- #
