"""
    convert_coord(coordinate::String)
Takes non-decimal-degree format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of coordinate strings to convert them.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
### Example
```
julia> convert_coord("-41 31.52")
-41.5253f0

julia> convert_coord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```
"""
function convert_coord(coordinate::String)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = parse.(Float32, split(coord_strip, r"\s|:"))
    split_coord[2] /=60.0
    if length(split_coord) == 3
        split_coord[3] /=3600.0
    end
    conv = sum(abs.(split_coord))
    # N + E are positive | S + W are negative
    if split_coord[1] < 0 || occursin(r"[SW]", uppercase(coordinate))
        # negative
        return round(conv * -1, digits = 4)
    else
        # positive
        return round(conv, digits = 4)
    end
end

"""
    multitest_missing(pvals::Vector{T}, correction::String) where T <: Union{Missing, <:AbstractFloat}
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. Missing values are first removed from the array, the appropriate
correction made, then missing values are re-added to the array at their original
positions. See MultipleTesting.jl docs for full more detailed information.
#### example
`multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")`

### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
"""
@inline function multitest_missing(pvals::Vector{T}, correction::String) where T <: Union{Missing, <:AbstractFloat}
    # get indices of where original missing are
    miss_idx = findall(i -> i === missing, pvals)

    # make seperate array for non-missing P vals
    p_no_miss = skipmissing(pvals) |> collect

    # make a dict of all possible tests and their respective functions
    d = Dict(
        "bonferroni" => Bonferroni(),
        "holm" => Holm(),
        "hochberg" => Hochberg(),
        "bh" => BenjaminiHochberg(),
        "by" => BenjaminiYekutieli(),
        "bl" => BenjaminiLiu(),
        "hommel" => Hommel(),
        "sidak" => Sidak(),
        "forwardstop" => ForwardStop(),
        "fs" => ForwardStop(),
        "bc" => BarberCandes(),
    )

    correct = adjust(p_no_miss, d[lowercase(correction)]) |> Vector{Union{Missing, Float64}}

    # re-add missing to original positions
    @inbounds for i in miss_idx
        @inbounds insert!(correct, i, missing)
    end
    return correct
end

"""
    nonmissing(vec::T) where T<:AbstractArray
Convenience function to count the number of non-`missing` values
in a vector.
"""
function nonmissing(vec::T) where T<:AbstractArray
    count(!ismissing, vec)
end


#TODO add to docs/API/Utils
"""
    permute_loci!(data::PopData)
Edits `PopData` in place with loci permuted across populations within
the `.loci` dataframe.
"""
function permute_loci!(data::PopData)
    loci_df = groupby(data.loci, :locus)
    for locus in loci_df
        shuffle!(locus.population)
    end
end

#TODO add to docs/API/Utils
"""
    permute_samples!(data::PopData; meta::Bool = false)
Edits `PopData` in place with samples permuted across populations within
the `.loci` dataframe. Since performance is important for many permutation,
the default is to only edit the `.loci` table in place; use `meta = true`
if you also require the `.meta` dataframe edited in place.
"""
function permute_samples!(data::PopData; meta::Bool = false)
    pops = shuffle(data.meta.population)

    if meta == true
        meta_pops = deepcopy(pops)
        meta_df = groupby(data.meta, :name)
        for name in meta_df
            name.population .= pop!(meta_pops)
        end
    end

    loci_df = groupby(data.loci, :name)
    for name in loci_df
        name.population .= pop!(pops)
    end
end

"""
    reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    !iszero(num) ? 1.0/float(num) : 0
end