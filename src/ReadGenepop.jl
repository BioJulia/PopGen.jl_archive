
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genos::Dict
    xloc::Array{Union{Int64,Float64},1}
    yloc::Array{Union{Int64,Float64},1}
end

"""
    loadgenepop(infile::String; ploidy::Int64 = 2, popsep::Any = "POP", numpop::Int64)
Load a Genepop format file into memory as a PopObj object.
- `infile` : path to Genepop file
- `ploidy` : ploidy of the organism
- `popsep` : word that separates populations in `infile` (e.g. "POP")
- `numpop` : number of populations in `infile` (used for checking parser)
File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after comment as one-per-line without commas
- A line with a particular keyword (e.g. "POP") must delimit populations
- File is tab or space delimted

usage: `waspsNY = loadgenepop("wasp_hive.gen", ploidy = 2, popsep = "POP", numpop = 2);`

Genepop file example:  \n
---------------------
Wasp populations in New York \n
Locus1 \n
Locus2 \n
Locus3 \n
POP \n
Oneida_01,  250230 564568 110100  \n
Oneida_02,  252238 568558 100120  \n
Oneida_03,  254230 564558 090100  \n
POP \n
Newcomb_01,  254230 564558 080100 \n
Newcomb_02,  000230 564558 090080 \n
Newcomb_03,  254230 000000 090100 \n
Newcomb_04,  254230 564000 090120 \n
---------------------
"""
function loadgenepop(infile::String; ploidy::Int64 = 2, popsep::Any = "POP", numpop::Int64)
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

"""
    loadcsv(infile::String; delim::Union{Char,String,Regex}, ploidy::Int64 = 2)
Load a CSV-type file into memory as a PopObj object
- `infile` : path to Genepop file
- `delim` values can be space (" "), comma (","), tab ("\\t"), etc.
- `ploidy` : ploidy of the organism
File formatting:
- Loci names must be first row
- Individuals names must be first value in row
- Population ID's must be second value in row

example: `lizardsCA = loadcsv("CA_lizards.csv", delim = ",", ploidy = 2);`

Formatting example:  \n
---------------------  \n
Locus1,Locus2,Locus3   \n
sierra_01,1,001001,002002,001001   \n
sierra_02,1,001001,001001,001002   \n
snbarb_03,2,001001,001001,001002 \n
snbarb_02,2,001001,001001,001001 \n
snbarb_03,2,001002,001001,001001 \n
---------------------
"""
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
