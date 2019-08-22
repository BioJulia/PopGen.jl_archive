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
        println( "   ", " #Inds | Pop " )
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
function genotypes(x::PopObj; loci::Union{String, Array})
    positions = []
    if typeof(loci) == String
        lowercase(loci) ∉ lowercase.(x.loci) && error("$loci not found in PopObj")
        geno_position = findall(i->i==lowercase(loci), lowercase.(x.loci))
        length(geno_position) > 1 && error("More than one instance of $loci in PopObj.genotypes. Please check data")
        push!(positions, geno_position[1])
    else
    for each in loci
            lowercase(each) ∉ lowercase.(x.loci) && error("$each not found in PopObj")
            geno_position = findall(i->i==lowercase(each), lowercase.(x.loci))
            length(geno_position) > 1 && error("More than one instance of $each in PopObj.genotypes. Please check data")
            push!(positions, geno_position[1])
        end
    end
    returnarray = []
    for indiv in x.ind
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
    hcat(x.ind, returnarray)
end


"""
    Base.missing(x::PopObj; plot::Bool = false)
Identify and count missing loci in each individual of a `PopObj`. Returns a tuple
of `DataFrames`: loci per individual, number per loci.

Example:

`aardvark = genepop("aardvark.gen", numpop = 5)`  # load file to PopObj

`missing_ind,missing_loc = missing(aardvark)`
"""
function Base.missing(x::PopObj)
    # by individuals
    nmissing = []
    misslociarray = []
    for each in x.ind
        missloci = x.loci[findall(i->i==(0,0), x.genotypes[each])]
        push!(nmissing, length(missloci))
        if length(missloci) == 0
            push!(misslociarray, [])
        else
            push!(misslociarray, missloci)
        end
    end
    ind_df = DataFrame(ind = x.ind, population = string.(x.popid),  nmissing = nmissing, loci = misslociarray)

    # count missing loci
    d = Dict()
    for b in ind_df[!, :4]
        for c in b
            if c ∉ keys(d)
                d[c] = 1
            else
                d[c] += 1
            end
        end
    end

    # add loci without any missing
    countarray = []
    for locus in x.loci
        if locus ∉ keys(d)
            d[locus] = 0
        end
        push!(countarray, d[locus])
    end

    #convert to DF and sort
    loci_df = hcat(x.loci, countarray) |> DataFrame ; sort!(loci_df, 2)
    return (ind_df, loci_df)
end


"""
    plot_missing(x::PopObj; color = false)
Return an interactive plot of the number of missing loci in individuals of a
`PopObj`, along with the number of missing individuals per locus. To set a
custom color palette, use `color = [color1, color2, etc.]`
"""
function plot_missing(x::PopObj; color = false)
    by_ind,by_loci = missing(x);
    ys = Array[subdf[!, :nmissing] for subdf in groupby(by_ind[!, 1:3], :population)]
    texts = Array[subdf[!, :ind] for subdf in groupby(by_ind[!, 1:3], :population)]
    popnum = length(by_ind[!, :population] |> unique)
    if color == false
        colors = ["hsl($i, 50%, 50%)" for i in range(0, stop=300, length=popnum)]
    else
        colors = color
    end
    returnplots = []

    byind = [box(y=y,
            marker_color=mc,
            name=name,
            text = text,
            jitter = 1,
            pointpos = 0,
            marker_size = 8.5,
            boxpoints = "all",
            marker=attr(line=attr(width=0.75)),
            )
            for (y, mc, name, text) in zip(ys, colors, unique(by_ind[!, :population]), texts)]

    bylocus = bar(;x = by_loci[!, :1],
                  y = by_loci[!, :2],
                  marker=attr(color="rgb(146,134,184)"),
                  name = "# missing data",
                  )

    layout_ind = Layout(title = "Number of missing loci per population",
                    hovermode = "closest",
                    yaxis = attr(title = "# missing loci", zeroline = false),
                    xaxis = attr(title = "Population", zeroline = false)
                    )
    layout_loci = Layout(title = "Missing data per locus",
                    hovermode = "closest",
                    bargap = 0.05,
                    yaxis = attr(title = "# missing", zeroline = false),
                    xaxis = attr(title = "loci", zeroline = false, showticklabels = false),
                    )
    ind_plot =  plot(byind, layout_ind)
    loci_plot = plot(bylocus, layout_loci)

    return [ind_plot loci_plot]
end


function remove_ind!(x::PopObj, inds::Union{String, Array})
    # get individuals indices
    if typeof(inds) == String
        inds ∉ x.ind && error("individual \"$inds\" not found")
        idx = findfirst(i -> i == inds, x.ind)
    else
        idx = []
        for each in inds
            each ∉ x.ind && error("individual \"$each\" not found")
            end
            push!(idx, findfirst(i -> i == each, x.ind) )
        end
    end
    deleteat!(x.ind, idx)  # delete name(s)
    deleteat!(x.popid, idx)    # delete popid(s)
    if length(x.longitude) != 0 && length(x.latitude) != 0
        deleteat!(x.longitude, idx)    # delete xloc(s)
        deleteat!(x.longitude, idx)    # delete yloc(s)
    end
    for each in inds
        delete!(x.genotypes, each)  # delete genotypes
    end
    return x
end

function remove_loci!(x::PopObj, loci::Union{String, Array})
    # get loci indices
    if typeof(loci) == String
        loci ∉ x.loci && error("locus \"$loci\" not found")
        idx = findfirst(i -> i == loci, x.loci)
    else
        idx = []
        for locus in loci
            locus ∉ x.loci && error("locus \"$locus\" not found")
            push!(idx, findfirst(i -> i == locus, x.loci) )
        end
    end
    deleteat!(x.loci, idx)  # delete loci names from list
    # remove genotypes of loci from all individuals
    for each in x.ind
        deleteat!(x.genotypes[each], idx)
    end
    return x
end
