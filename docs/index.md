### What is this?
This repository holds a [Julia](https://julialang.org) implementation of the [Core _Escherichia coli_ model of Palsson and coworkers](https://www.ncbi.nlm.nih.gov/pubmed/26443778) developed
by the [Varnerlab](http://www.varnerlab.org) in the [Robert Frederick Smith School of Chemical and Biomolecular Engineering, Cornell University](https://www.cheme.cornell.edu/cbe).
We use this model as a teaching tool in [CHEME 7770 Advanced Biomolecular Engineering](https://varnerlab.github.io/CHEME-7770-Cornell-S19/).

### Requirements and Installation
To use this model, you need a couple of things. First, [Julia](https://julialang.org) version 1.0 or greater must be installed on your local machine.
Alternatively, you can use [JuliaBox](https://juliabox.com) in the cloud. Next, you'll need the model code. You can download the code either as a zip file, or by cloning directly from GitHub:

    git clone https://github.com/varnerlab/Core_Escherichia_coli_model.git

Once the model code is own your machine (or on your [JuliaBox](https://juliabox.com) instance in the cloud), to execute model calculations, start the [Julia](https://julialang.org) REPL in the ``src`` subdirectory and
execute a sample command:

    julia> include("Test.jl")

#### What does Test.jl do?
``Test.jl`` downloads (and installs) required packages that are not already installed on your machine. It then maximizes the aerobic growth of _E. coli_ on glucose with some example external measurement distrubutions.
