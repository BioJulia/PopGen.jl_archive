"""
    indnames(x::PopObj)
View individual/sample names in a `PopObj`

Equivalent to `PopObj.inds`
"""
indnames(x::PopObj) = x.ind


"""
    loci(x::PopObj)
View loci names in a `PopObj`

Equivalent to `PopObj.loci`
"""
loci(x::PopObj) = x.loci


"""
    locations(x::PopObj)
View location data (`.longitude` and `.latitude`) in a `PopObj`

Use `locations!` to add spatial data to a `PopObj`
"""
locations(x::PopObj) = hcat(x.longitude, x.latitude)


"""
    locations!(x::PopObj; xloc::Array, yloc::Array)
Add location data (longitude `xloc`, latitude `yloc`) to `PopObj`. Takes decimal
degrees or decimal minutes format. **Must** use `-` symbol instead of cardinal directions.
Location data must be in order of `ind`. Replaces existing `PopObj` location data.
- Decimal Degrees : `-11.431`
- Decimal Minutes : `"-11 43.11"` (must use space and double-quotes)

If conversion is not necessary, can directly assign `PopObj.xloc` and `PopObj.yloc`
"""
function locations!(x::PopObj; xloc::Array, yloc::Array)
    # test for decimal degrees vs decimal minutes
    if occursin(" ", string(xloc[1])) == false && occursin(" ", string(yloc[1])) == false
        x.longitude = xloc ;
        x.latitude = yloc ;
        println("  xloc    yloc")
        hcat(x.xloc, x.yloc)
    else
        # make sure decimal minutes are Strings
        if typeof(xloc) != Array{String,1}
            xloc = string.(xloc)
        end
        if typeof(yloc) != Array{String,1}
            yloc = string.(yloc)
        end
        # convert xloc to decimal degrees
        xlocConverted = []
        for value in xloc
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0   # if negative, subtract
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(xlocConverted, decideg)
            else                           # if positive, add
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(xlocConverted, decideg)
            end
        end
        # convert yloc to decimal degrees
        ylocConverted = []
        for value in yloc
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(ylocConverted, decideg)
            else
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(ylocConverted, decideg)
            end
        end
        x.longitude = xlocConverted
        x.latitude = ylocConverted
        println("  xloc    yloc")
        hcat(x.xloc, x.yloc)
    end
end


"""
    popid(x::PopObj; listall::Bool = false)
View unique population ID's in a `PopObj`.

`listall = true`, displays `ind` and their `popid` instead (default = `false`).
"""
function popid(x::PopObj; listall::Bool = false)
    if listall == true
        hcat(x.ind, x.popid)
    else
        println( "   ", "#Inds | Pop " )
        println( "   ", "--------------" )
        popcounts = hcat([sum(x.popid .== i) for i in unique(x.popid)],unique(x.popid))
        for eachpop in 1:length(popcounts)÷2
            println("\t", popcounts[eachpop], "\t", " |", "\t", popcounts[eachpop,2])
        end
    end
end

"""
    popid!(x::PopObj; rename::Dict)
Rename the population ID's of `PopObj.popid`.

Uses a `Dict` of `[popid] => replacement` to rename

Example:

potatopops = Dict(1 => "Idaho", 2 => "Russet")

popid!(potatoes, rename = potatopops)
"""
function popid!(x::PopObj; rename::Dict)
    for eachkey in keys(rename)
        replace!(x.popid, eachkey => rename[eachkey])
    end
    popid(x, listall = true)
end


"""
    genotypes(x::PopObj; loci::Array{String,1})
Get the per-individual genotypes of an array of `loci` in a `PopObj`.
- Locus names are case insensitive, but must be in quotes
- Query a single locus as a 1-dimensional Array

Examples:

genotypes(eggplant, loci = ["contig_001", "contig_101", "contig_5150"])

genotypes(eggplant, loci = ["contig_099"])
"""
function genotypes(x::PopObj; loci::Array{String,1})
    positions = []
    for each in loci
        lowercase(each) ∉ lowercase.(x.loci) && error("$each not found in PopObj")
        geno_position = findall(i->i==lowercase(each), lowercase.(x.loci))
        length(geno_position) > 1 && error("More than one instance of $each in PopObj.genotypes. Please check data")
        push!(positions, geno_position[1])
    end
    returnarray = []
    for indiv in sort(collect(keys(x.genotypes)))
        tmp = []
        for posit in positions
            if length(tmp) == 0
                tmp = x.genotypes[indiv][posit]
            else
                tmp = hcat(tmp, x.genotypes[indiv][posit])

            end
        end
        if length(returnarray) == 0
            returnarray = tmp
        else
            returnarray = cat(returnarray,tmp, dims = 1)
        end
    end
    hcat(sort(collect(keys(x.genotypes))), returnarray)
end


"""
    Base.missing(x::PopObj; plot::Bool = false)
Identify and count missing loci in each individual of a `PopObj`. Will print out
counts per individual but return a `Dict` of which loci are missing for each
individual. Use `plot = true` to only print counts with a complementary plot.

Example:

`aardvark = genepop("aardvark.gen", numpop = 5)`  # load file to PopObj

`misscounts = missing(aardvark) ;`
"""
function Base.missing(x::PopObj; plot::Bool = false)
    d = Dict()
    nmissing = []
    for each in sort(collect(keys(x.genotypes)))
        missloci = x.loci[findall(i->i==(0,0), x.genotypes[each])]
        push!(nmissing, length(missloci))
        if length(missloci) == 0
            println(each, "  ", 0)
        else
            d[string(each)] =  missloci
            println(each, "  ", length(missloci))
        end
    end
    if plot == true
        summ = sortslices(hcat(nmissing,x.ind), dims = 1)
        missingplot = bar(x=summ[1:end,1],
                          y = summ[1:end,2],
                          name = "Missing Data",
                          orientation = :h,
                          )
        layout = Layout(title = "Missing loci per individual",
                        xaxis_title = "# missing loci",
                        left_margin = 110
                        )
        return plot(missingplot, layout)
    else
        return d
    end
end

