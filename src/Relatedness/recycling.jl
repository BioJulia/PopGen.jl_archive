

#### Rework for efficiency ####

function QuellerGoodnight(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2

    n1 = sum((a == c, a == d, b == c, b == d)) - 2.0 * (alleles[loc][a] + alleles[loc][b])
    n2 = sum((a == c, a == d, b == c, b == d)) - 2.0 * (alleles[loc][c] + alleles[loc][d])

    d1 = 2.0 * (1.0 + (a==b) - alleles[loc][a] - alleles[loc][b])
    d2 = 2.0 * (1.0 + (c==d) - alleles[loc][c] - alleles[loc][d])
    return (n1, d1, n2, d2)
end

function Ritland(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2

    A = ((alleles[loc] |> length) - 1)

    R = 0.0
    for allele in keys(alleles[loc])
        # Individual locus relatedness value (eq 7 in paper)
        R += ((((a == allele) + (b == allele)) * ((c == allele) + (d == allele))) / (4.0 * alleles[loc][allele]))
    end
    R = (2 / A) * (R - 1.0)
    # return numerator,denominator
    return ((R * A), A, 0.0, 0.0)
end

function LynchRitland(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2
    A = ((alleles[loc] |> length) - 1)

    n1 = alleles[loc][a] * ((b == c) + (b == d)) + alleles[loc][b] * ((a == c) + (a == d)) - 4.0 * alleles[loc][a] * alleles[loc][b]
    n2 = alleles[loc][c] * ((d == a) + (d == b)) + alleles[loc][d] * ((c == a) + (c == b)) - 4.0 * alleles[loc][c] * alleles[loc][d]

    d1 = 2.0 * (1.0 + (a == b)) * (alleles[loc][a] + alleles[loc][b]) - 8.0 * alleles[loc][a] * alleles[loc][b]
    d2 = 2.0 * (1.0 + (c == d)) * (alleles[loc][c] + alleles[loc][d]) - 8.0 * alleles[loc][c] * alleles[loc][d]

    RL = (n1 / d1) + (n2 / d2)
    # return numerator, denominator
    return (RL, A, 0.0, 0.0)
end

"""
    relatedness_moment(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 10 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function relatedness_moment(data::PopData, ind1::String, ind2::String; alleles::T, method::Vector{Function}) where T <: NamedTuple
    #NEED TO CHECK TO CONFIRM EQUATIONS
    #Extract the pair of interest's genotypes
    gen1 = get_genotypes(data, ind1)
    gen2 = get_genotypes(data, ind2)

    d = Dict{Symbol,Vector{Float64}}()

    for (loc,geno1,geno2) in zip(skipmissings(Symbol.(loci(data)), gen1, gen2)...)
        for mthd in method
            nu, denm, nu2, denm2 = mthd(loc, geno1, geno2, alleles)
            get!(d, Symbol(mthd), zeros(4))[1] += nu
            get!(d, Symbol(mthd), zeros(4))[2] += denm
            get!(d, Symbol(mthd), zeros(4))[3] += nu2
            get!(d, Symbol(mthd), zeros(4))[4] += denm2
        end
    end
    if haskey(d, :QuellerGoodnight)
        qg = d[:QuellerGoodnight]
        d[:QuellerGoodnight][1] = (qg[1]/qg[2]) + (qg[3]/qg[4])
        d[:QuellerGoodnight][2] = 2.0
    end
    return NamedTuple{Tuple(keys(d))}(getindex.(values(d), 1) ./ getindex.(values(d), 2))
end