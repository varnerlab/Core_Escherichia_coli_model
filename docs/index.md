### What is this?
This repository holds a [Julia](https://julialang.org) implementation, developed
by the [Varnerlab](http://www.varnerlab.org) in the [Robert Frederick Smith School of Chemical and Biomolecular Engineering, Cornell University](https://www.cheme.cornell.edu/cbe),
of the [Core _Escherichia coli_ model of Palsson and coworkers](https://www.ncbi.nlm.nih.gov/pubmed/26443778).
We use this model as a teaching tool in [CHEME 7770 Advanced Biomolecular Engineering](https://varnerlab.github.io/CHEME-7770-Cornell-S19/) and to test new ways of thinking about metabolic modeling and data integration.

### Requirements and Installation
To use this model, you need a couple of things. First, [Julia](https://julialang.org) version 1.0 or greater must be installed on your local machine.
Alternatively, you can use [JuliaBox](https://juliabox.com) in the cloud. The code solves the Linear Programming (LP) problem associated with [Flux Balance Analysis (FBA)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/) using the
[GNU Linear Programming Kit (GLPK)](https://www.gnu.org/software/glpk/) and performs [MOMA calculations](https://www.pnas.org/content/99/23/15112) using the [Gurobi](http://www.gurobi.com) solver.
Both the GLPK and Gurobi libraries must be installed, and the paths configured, or the model calculations will not run properly.
Next, you'll need the model code. You can download the code either as a zip file, or by cloning directly from GitHub:

    git clone https://github.com/varnerlab/Core_Escherichia_coli_model.git

Once the model code is own your machine (or on your [JuliaBox](https://juliabox.com) instance in the cloud), to test the installation, start the [Julia](https://julialang.org) REPL in the ``src`` subdirectory and
execute the command:

    julia> include("Test.jl")

#### What does Test.jl do?
``Test.jl`` downloads (and installs) required packages that are not already installed on your machine. It then maximizes the aerobic growth of _E. coli_ on glucose with some example external measurement distrubutions.
