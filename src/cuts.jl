# One-dimensional cuts through scattered (q, E, I) data.
#
# To make a cut: choose the cut direction (the x-axis of the 1D plot), keep the
# points whose *other* coordinate lies in a window around a centre value, then
# bin along the cut direction and combine the intensities in each bin.

"""
    Cut

Result of [`make_cut`](@ref). Fields:

- `x`: bin centres along the cut direction.
- `I`: combined intensity per bin (`NaN` for empty bins).
- `err`: standard error per bin from the scatter of the points
  (`NaN` with fewer than 2 points).
- `counts`: number of data points per bin.
- `dir`: `:E` (intensity vs energy) or `:q` (intensity vs `|Q|`).
"""
struct Cut
    x::Vector{Float64}
    I::Vector{Float64}
    err::Vector{Float64}
    counts::Vector{Int}
    dir::Symbol
end

"""
    make_cut(q, E, I; dir, center, width, bins=:auto, reduce=:mean) -> Cut
    make_cut(map::PowderMap; dir, center, width, bins=:auto, reduce=:mean) -> Cut

One-dimensional cut through scattered data given as equal-length arrays of
coordinates `q`, `E` and intensities `I` (any shape; `NaN` intensities are
ignored).

- `dir = :E`: intensity vs energy at fixed `q` (constant-Q cut);
  `dir = :q`: intensity vs `q` at fixed energy (constant-E cut).
- `center`, `width`: the fixed coordinate is integrated over
  `|fixed - center| ≤ width/2`.
- `bins`: `:auto` (about √n bins), an integer number of equal bins spanning
  the selected data, or a range of bin *edges* such as `0:0.05:2`.
- `reduce`: `:mean` (robust to uneven point density) or `:sum`.

The second form cuts a model map the same way, so model and data cuts can be
compared directly.
"""
function make_cut(q::AbstractArray, E::AbstractArray, I::AbstractArray;
                  dir::Symbol, center::Real, width::Real, bins=:auto, reduce::Symbol=:mean)
    length(q) == length(E) == length(I) || throw(ArgumentError("q, E, I must have equal length"))
    dir in (:E, :q) || throw(ArgumentError("dir must be :E or :q"))
    reduce in (:mean, :sum) || throw(ArgumentError("reduce must be :mean or :sum"))
    width > 0 || throw(ArgumentError("width must be positive"))

    x, fixed = dir === :E ? (vec(E), vec(q)) : (vec(q), vec(E))
    Iv = vec(I)
    keep = (abs.(fixed .- center) .≤ width / 2) .& .!isnan.(Iv)
    any(keep) || throw(ArgumentError("no data points within $center ± $(width/2); " *
                                     "widen `width` or check `center`"))
    xin, Iin = x[keep], Iv[keep]

    edges = _bin_edges(xin, bins)
    nb = length(edges) - 1
    buckets = [Float64[] for _ in 1:nb]
    for (xi, v) in zip(xin, Iin)
        (xi < edges[1] || xi > edges[end]) && continue
        b = min(searchsortedlast(edges, xi), nb)      # right edge belongs to last bin
        push!(buckets[b], v)
    end

    Iout, err = fill(NaN, nb), fill(NaN, nb)
    counts = length.(buckets)
    for (b, vals) in enumerate(buckets)
        n = length(vals)
        n == 0 && continue
        Iout[b] = reduce === :mean ? mean(vals) : sum(vals)
        if n ≥ 2
            se = std(vals) / sqrt(n)
            err[b] = reduce === :mean ? se : se * n
        end
    end
    centres = (edges[1:end-1] .+ edges[2:end]) ./ 2
    return Cut(centres, Iout, err, counts, dir)
end

function make_cut(m::PowderMap; kwargs...)
    q = [Q for Q in m.radii, _ in m.energies]
    E = [ω for _ in m.radii, ω in m.energies]
    return make_cut(q, E, m.data; kwargs...)
end

function _bin_edges(x, bins)
    bins isa AbstractRange && return collect(float.(bins))
    lo, hi = extrema(x)
    lo == hi && (lo -= 0.5; hi += 0.5)
    n = if bins === :auto
        clamp(round(Int, sqrt(length(x))), 5, 200)
    elseif bins isa Integer
        bins > 0 ? Int(bins) : throw(ArgumentError("bins must be positive"))
    else
        throw(ArgumentError("bins must be :auto, an Integer, or a range of edges"))
    end
    return collect(range(lo, hi; length=n + 1))
end
