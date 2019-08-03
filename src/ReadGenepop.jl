mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genos::Dict
    xloc::Array{Union{Int64,Float64},1}
    yloc::Array{Union{Int64,Float64},1}
end

function readgenepop(infile::String, ploidy::Int64 = 2)
    gpop = split(open(readlines,infile)[2:end],"POP");
    d = Dict()
    locinames = gpop[1]
    popid = []
    indnames = []
    for i in 2:length(gpop)
        append!(popid, fill(i-1,length(gpop[i])))
        for j in 1:length(gpop[i])
            push!(indnames, split( gpop[i][j], r"\,|\t")[1] )
            d[ indnames[j] ] = split( strip(gpop[i][j]), r"\s|\t" )[2:end] |> Array{String,1}
       end
    end
    return PopObj(indnames,
                  popid,
                  locinames,
                  ploidy,
                  d,
                  fill( 0, length(indnames) ),
                  fill( 0, length(indnames) )
                  )
end
