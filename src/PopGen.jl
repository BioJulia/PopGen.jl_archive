module PopGen
using DataFrames, PlotlyJS, Statistics

export PopObj,
    show,
    csv,
    genepop,
    indnames,
    loci,
    locations,
    locations!,
    popid,
    popid!,
    missing,
    plotmissing






    ##############################################################################
    ##
    ## Load files
    ##
    ##############################################################################

    include("/home/pdimens/PopGen/src/PopObj.jl")
    include("/home/pdimens/PopGen/src/Read.jl")
    include("/home/pdimens/PopGen/src/Manipulate.jl")
    include("/home/pdimens/PopGen/src/Plotting.jl")


end # module PopGen
