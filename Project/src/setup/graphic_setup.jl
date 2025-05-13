#!/usr/bin/julia

# Set default options for plotting

using Plots; pgfplotsx()
using ColorSchemes
using LaTeXStrings

# Get color scheme
MyColors = ColorSchemes.tab10

Plots.default(
    size = (440, 320),
    bgcolor = :white,
    palette = :seaborn_colorblind,
    
    #gridlinewidth = 0.5,
    #gridstyle = :dot,
    #gridcolor = :gray,
    
    linewidth = 1.0,
    
    foreground_color = :black,
    foreground_color_axis = :black,
    
    tick_direction = :in,
    minorticks = true,
    
    legend_background_color = nothing,
    legend=:topright,
    legend_foreground_color= nothing,
    legend_font_halign=:center,

    titlefontsize=12,
    legendfontsize = 9,
    guidefontsize = 10,
    tickfontsize = 10,

    fontfamily = "Computer Modern",
    legendfontfamily = "Computer Modern",
    titlefontfamily = "Computer Modern",
    guidefontfamily = "Computer Modern",
    tickfontfamily = "Computer Modern",

    markerstrokewidth=0,
    #margin = 1Plots.mm,
    framestyle = :box,
)

