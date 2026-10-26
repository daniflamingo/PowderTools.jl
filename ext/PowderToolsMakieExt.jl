module PowderToolsMakieExt

using PowderTools
using PowderTools: Cut, PowderMap, PowderData
using Makie

import PowderTools: plot_cut, plot_cut!, plot_powder

_xlabel(c::Cut) = c.dir === :E ? "Energy" : "|Q|"

function plot_cut!(ax, c::Cut; label=nothing, color=nothing, errorbars::Bool=true, kwargs...)
    good = c.counts .> 0
    x, y, e = c.x[good], c.I[good], c.err[good]
    color = something(color, Makie.wong_colors()[1])
    scatterlines!(ax, x, y; label, color, kwargs...)
    if errorbars && any(isfinite, e)
        ok = isfinite.(e)
        errorbars!(ax, x[ok], y[ok], e[ok]; color)
    end
    return ax
end

function plot_cut(cuts::Cut...; labels=nothing, xlabel=nothing, ylabel="Intensity", title="")
    isempty(cuts) && throw(ArgumentError("give at least one Cut"))
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=something(xlabel, _xlabel(first(cuts))), ylabel, title)
    for (i, c) in enumerate(cuts)
        label = isnothing(labels) ? nothing : labels[i]
        plot_cut!(ax, c; label, color=Makie.wong_colors()[mod1(i, 7)])
    end
    isnothing(labels) || axislegend(ax)
    return fig
end

function plot_powder(m::PowderMap; colorrange=nothing, colormap=:viridis, title="")
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="|Q|", ylabel="Energy", title)
    hm = heatmap!(ax, m.radii, m.energies, m.data;
                  colormap, colorrange=something(colorrange, Makie.automatic))
    Colorbar(fig[1, 2], hm)
    return fig
end

function plot_powder(d::PowderData; colorrange=nothing, colormap=:viridis, markersize=6, title="")
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="|Q|", ylabel="Energy", title)
    sc = scatter!(ax, d.Q, d.E; color=d.I, colormap, markersize,
                  colorrange=something(colorrange, Makie.automatic))
    Colorbar(fig[1, 2], sc)
    return fig
end

end # module
