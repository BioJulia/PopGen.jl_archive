function loadgenepop(infile::String; ploidy::Int64 = 2, numpop::Int64)
    gpop = split(open(readlines,infile)[2:end], "POP")
    if length(gpop)-1 != numpop
        error("incorrect number of populations detected \n expected : $(length(gpop)-1) \n detected : $numpop \n see docstring to verify that your infile follows standard Genepop format ")
    end
    d = Dict()
    locinames = gpop[1]
    popid = []
    indnames = []
    for i in 2:length(gpop)
        append!(popid, fill(i-1,length(gpop[i])))
        for j in 1:length(gpop[i])
            push!(indnames, split( strip(gpop[i][j]), r"\,|\t")[1] )
            d[ last(indnames) ] = split( strip(gpop[i][j]), r"\s|\t" )[2:end] |> Array{String,1}
       end
    end
    ## print some basic information ##
    println("\n", "Input File : ", joinpath(@__DIR__, infile) )
    println( "Number of Individuals : ", length(indnames) )
    println( "Number of Loci : ", length(locinames) )
    println( "Number of Populations : ", maximum(popid) )
    println( "\t", "Pop | #Inds " )
    println( "\t", "----------- " )
    popcounts = hcat(unique(popid),[sum(popid .== i) for i in unique(popid)])
    for eachpop in 1:length(popcounts)÷2
        println("\t", popcounts[eachpop], "   |   ", popcounts[eachpop,2])
    end
    println()
    PopObj(indnames,
          popid,
          locinames,
          ploidy,
          d,
          fill( 0, length(indnames) ),
          fill( 0, length(indnames) )
          ) ;
end

# load in CSV files
function loadcsv(infile::String; delim::Union{Char,String,Regex}, ploidy::Int64 = 2)
    gpop = open(readlines,infile)
    locinames = split(gpop[1], delim)
    popid = []
    indnames = []
    d = Dict()
    for i in gpop[2:end]
        tmp = split(i, delim) |> Array{String,1}
        d[tmp[1]] = tmp[3:end]
        push!(indnames, tmp[1])
        push!(popid, parse(Int64,tmp[2]))
    end
    ## print some basic information ##
    println("\n", "Input File : ", joinpath(@__DIR__, infile) )
    println( "Number of Individuals : ", length(indnames) )
    println( "Number of Loci : ", length(locinames) )
    println( "Number of Populations : ", maximum(popid) )
    println( "\t", "Pop | #Inds " )
    println( "\t", "----------- " )
    popcounts = hcat(unique(popid),[sum(popid .== i) for i in unique(popid)])
    for eachpop in 1:length(popcounts)÷2
        println("\t", popcounts[eachpop], "   |   ", popcounts[eachpop,2])
    end
    println()
    PopObj(indnames,
          popid,
          locinames,
          ploidy,
          d,
          fill( 0, length(indnames) ),
          fill( 0, length(indnames) )
          ) ;
end
