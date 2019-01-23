# ---------------------------------------------------------#
# check - do we have all the required packages installed?
using Pkg

# which packages are installed?
installed_packages_dict = Pkg.installed()

# First things first, do we have JSON installed?
if (haskey(installed_packages_dict, "JSON") == false)
    Pkg.add("JSON")
end

# ok, load the Required.json file -
using JSON
tmp_json_blob = JSON.parsefile("../Required.json")
required_packages_array = tmp_json_blob["required_packages_array"]
for required_packages_dict in required_packages_array

    # package name -
    package_name = required_packages_dict["name"]

    # do we have the package installed?
    if (haskey(installed_packages_dict, package_name) == false)
        Pkg.add(package_name)
    end
end
# ---------------------------------------------------------#

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
