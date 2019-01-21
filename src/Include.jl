# - LOAD SYSTEM PACKAGES --------------------------------- #
using Printf
using ProgressMeter
using Statistics
# -------------------------------------------------------- #

# - LOAD MY MODULES AND CODE ----------------------------- #
top_level_path = pwd()

# load the VLCobra module -
path_to_cobra_module = "$(top_level_path)/cobra"
push!(LOAD_PATH, path_to_cobra_module)
using VLCobraModule

# Load the VLGeneExpressionModule -
path_to_gex_module = "$(top_level_path)/gene_expression"
push!(LOAD_PATH, path_to_gex_module)
using VLGeneExpressionModule

# load the VLFlux module -
path_to_flux_module = "$(top_level_path)/flux"
push!(LOAD_PATH, path_to_flux_module)
using VLFluxModule

# load the VLTest module -
path_to_test_module = "$(top_level_path)/test"
push!(LOAD_PATH, path_to_test_module)
using VLTestModule
# -------------------------------------------------------- #
