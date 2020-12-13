function perm_samp(data::PopData)
    popcounts = countmap(data.meta.population)
    pops, counts = keys(popcounts), values(popcounts)
    gdf = groupby(data.loci, :name)
    perm_idx = shuffle(Xoroshiro128Star(), keys(gdf))
    
    #return gdf

end