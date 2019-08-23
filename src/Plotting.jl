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


"""
    plot_locations(x::PopObj; region::String = "world", projection::String = "mercator")

Returns a simple low resolution interactive scatterplot of the individuals in a
`PopObj`. Default `region` and `projection` are "world" and "mercator", respectively.

Example:

`plot_locations(manatees, region = "usa", projection = "albers usa")`

[regions]

"usa","europe", "asia", "africa", "north america", "south america"

[projections]

"equirectangular", "mercator", "orthographic", "natural earth",
"kavrayskiy7", "miller", "robinson", "eckert4",
"azimuthal equal area", "azimuthal equidistant",
"conic equal area", "conic conformal" , "conic equidistant",
"gnomonic", "stereographic", "mollweide", "hammer",
"transverse mercator", "albers usa", "winkel tripel",
"aitoff", "sinusoidal"
"""
function plot_locations(x::PopObj; region::String = "world", projection::String = "mercator")
    if projection ∉ ["equirectangular", "mercator", "orthographic", "natural earth",
                    "kavrayskiy7", "miller", "robinson", "eckert4",
                    "azimuthal equal area", "azimuthal equidistant",
                    "conic equal area", "conic conformal" , "conic equidistant",
                    "gnomonic", "stereographic", "mollweide", "hammer",
                    "transverse mercator", "albers usa", "winkel tripel",
                    "aitoff", "sinusoidal"]
        error("Projection not recognized. Please see the help doc for full list of projection options")
    end
    locs = scattergeo(lat=x.latitude,
                       lon=x.longitude,
                       hoverinfo="text",
                       text=["$i: $j" for (i,j) in zip(x.popid, x.ind)],
                       marker_line_color="black", marker_line_width=2
                       )
    if region == "world"
        geo = attr(scope = region,
                   projection_type = projection,
                   showcoastlines = false,
                   showcountries = true,
                   countrywidth = 0.75,
                   countrycolor = "rgb(255,255,255)",
                   subunitcolor = "rgb(255,255,255)",
                   showland = true,
                   landcolor = "rgb(217, 217, 217)",
                   )
   elseif region ∈ ["usa","europe", "asia", "africa", "north america", "south america"]
       geo = attr(scope = region,
                  showlakes = true,
                  lakecolor = "#a5a5a5",
                  showrivers = true,
                  rivercolor = "#a5a5a5",
                  showsubunits = true,
                  projection_type = projection,
                  showcoastlines = false,
                  showcountries = true,
                  countrywidth = 0.75,
                  countrycolor = "rgb(255,255,255)",
                  subunitcolor = "rgb(255,255,255)",
                  showland = true,
                  landcolor = "rgb(217, 217, 217)",
                  )
    else
        error("Leave \"region =\" empty or use one of: usa, europe, asia, africa, north america, south america")
    end
   layout = Layout(;title="Sample locations", showlegend=false, geo=geo)
   plot(locs, layout)
end
