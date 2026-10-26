# Implemented in ext/PowderToolsMakieExt.jl; available after `using GLMakie`
# or `using CairoMakie`.

"""
    plot_cut(cuts...; labels=nothing, xlabel=nothing, ylabel="Intensity", title="")

Plot one or more [`Cut`](@ref)s (e.g. a data cut and a model cut) with error
bars where available. Returns the `Figure`. Requires a Makie backend.
"""
function plot_cut end

"""
    plot_cut!(ax, cut; label=nothing, color=nothing, errorbars=true, kwargs...)

Add a [`Cut`](@ref) to an existing Makie axis. Requires a Makie backend.
"""
function plot_cut! end

"""
    plot_powder(map::PowderMap; colorrange=nothing, colormap=:viridis, title="")
    plot_powder(data::PowderData; colorrange=nothing, colormap=:viridis,
                markersize=6, title="")

Heatmap of a powder map, or scatter plot of scattered powder data, with `|Q|`
on the x-axis and energy on the y-axis. Returns the `Figure`. Requires a
Makie backend.
"""
function plot_powder end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f in (plot_cut, plot_cut!, plot_powder)
            print(io, "\nPlotting requires a Makie backend: run `using GLMakie` or `using CairoMakie` first.")
        end
    end
end
