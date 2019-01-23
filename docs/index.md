[Julia](https://julialang.org) implementation of the [Core _Escherichia coli_ model of Palsson and coworkers](https://www.ncbi.nlm.nih.gov/pubmed/26443778).

### Requirements and Installation
To interact with this model code, you need to have [Julia](https://julialang.org) version 1.0 or greater installed on your local machine.
Alternatively, you can use [JuliaBox](https://juliabox.com) in the cloud. Download the code (either as a zip) or by cloning from GitHub:

  git clone https://github.com/varnerlab/Core_Escherichia_coli_model.git

Once the model code is own your machine, to execute model calculations, start the [Julia](https://julialang.org) REPL in the `src` subdirectory of the model and
execute the command:

  julia>include("Solve.jl")
