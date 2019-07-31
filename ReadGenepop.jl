function readgenepop(infile::String)
    seq = open(readlines,infile);
    d = Dict()
    locnames = Array{String,1}()
    popid = 0
    i = 2
    while seq[i] != "POP"
        push!(locnames,seq[i])
        i += 1
    end
    while i <= length(seq)
        if seq[i] == "POP"
            popid += 1
            i += 1
        else
            mergearray = []
            push!(mergearray, locnames)
            push!(mergearray, split(seq[i],"\t")[2:end-1])
            genodigits = length(mergearray[2][1]) ÷ 2
            #strip all commas from infile earlier to avoid this
            indname = split(seq[i],",")[1]
            d[indname] = Dict()
            d[indname]["pop"] = [popid]
            for locus in 1:length(mergearray[1])
                d[indname][mergearray[1][locus]] = [String(mergearray[2][locus][1:genodigits]),String(mergearray[2][locus][genodigits+1:end])]
            end
            i += 1
        end
    end
    return d
end
