## create struct for individuals

mutable struct GenoInd
    name::String
    genotypes::Array{Array{String,1},1}
    popid::Union{Int64,String}
    locx::Union{Int64,Float64}
    locy::Union{Int64,Float64}
end

# create struct for entire object. Similar-ish to a genind from R
mutable struct GenoLite
    indiv::Dict
    loci::Array{String,1}
    ploidy::Int64
    #other::Type{Any}
end

function readgenepop(infile::String, ploidy::Int64 = 2)
    seq = open(readlines,infile);
    d = Dict()
    locnames = Array{String,1}()
    popid = 0
    i = 2
    while uppercase(seq[i]) != "POP"
        push!(locnames,seq[i])
        i += 1
    end
    while i <= length(seq)
        if uppercase(seq[i]) == "POP"
            popid += 1
            i += 1
        else
            tmparray = []
            genoarray = split(seq[i],"\t")[2:end-1]
            genodigits = length(genoarray[1]) ÷ ploidy
            #strip all commas from infile earlier to avoid this
            indname = split(seq[i],",")[1]
            for locus in 1:length(genoarray)
                push!(tmparray,[String(genoarray[locus][1:genodigits]),String(genoarray[locus][genodigits+1:end])])
            end
            #push!(d,GenoInd(indname,tmparray, popid, 0, 0))
            d[indname] = GenoInd(indname,tmparray, popid, 0, 0)
            i += 1
            if i > length(seq)
                return GenoLite(d,locnames,ploidy)
            end
        end
    end
end

# usage: variablename = readgenepop("infile", 2)

# basic function to see indiv names
function names(data::GenoLite)
    keys(data.indiv) |> collect |> sort
end

# basic function to view indiv name along with their population ID
function popid(data::GenoLite)
    for i in keys(data.indiv)
        # maybe makes more sense to output a tuple? Array{Array{Any,1}}?
        println(data.indiv[i].name, "  ", data.indiv[i].popid)
    end
end

# must create another function to edit popID
# a find [key] replace-like function?
