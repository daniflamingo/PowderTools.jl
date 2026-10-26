"""
    PowderMap(radii, energies, data)

Powder intensity map on a regular grid: `data[iQ, iE]` is the intensity at
`|Q| = radii[iQ]` and energy `energies[iE]`. Returned by
[`powder_one_magnon`](@ref) and [`powder_two_magnon`](@ref).

`PowderMap(res)` converts the result of Sunny's `powder_average` (whose
`data` is indexed `[iE, iQ]`). Maps on the same grid can be added and scaled,
e.g. `one + 0.3 * two`.
"""
struct PowderMap
    radii::Vector{Float64}
    energies::Vector{Float64}
    data::Matrix{Float64}
    function PowderMap(radii, energies, data)
        size(data) == (length(radii), length(energies)) ||
            throw(DimensionMismatch("data must be (length(radii), length(energies)) = " *
                                    "$((length(radii), length(energies))), got $(size(data))"))
        new(collect(Float64, radii), collect(Float64, energies), Matrix{Float64}(data))
    end
end

PowderMap(res::Sunny.PowderIntensities) = PowderMap(res.radii, res.energies, permutedims(real.(res.data)))

function _same_grid(a::PowderMap, b::PowderMap)
    a.radii == b.radii && a.energies == b.energies ||
        throw(ArgumentError("PowderMaps are on different (|Q|, E) grids"))
end

Base.:+(a::PowderMap, b::PowderMap) = (_same_grid(a, b); PowderMap(a.radii, a.energies, a.data + b.data))
Base.:-(a::PowderMap, b::PowderMap) = (_same_grid(a, b); PowderMap(a.radii, a.energies, a.data - b.data))
Base.:*(s::Real, a::PowderMap) = PowderMap(a.radii, a.energies, s * a.data)
Base.:*(a::PowderMap, s::Real) = s * a

function Base.show(io::IO, m::PowderMap)
    print(io, "PowderMap(", length(m.radii), " |Q| × ", length(m.energies), " E; |Q| ∈ ",
          extrema(m.radii), ", E ∈ ", extrema(m.energies), ")")
end

# Bilinear interpolation of the map at (Q, E); 0 outside the grid.
function _interp(m::PowderMap, Q::Real, E::Real)
    (; radii, energies, data) = m
    (radii[1] ≤ Q ≤ radii[end] && energies[1] ≤ E ≤ energies[end]) || return 0.0
    i = clamp(searchsortedlast(radii, Q), 1, max(length(radii) - 1, 1))
    j = clamp(searchsortedlast(energies, E), 1, max(length(energies) - 1, 1))
    (length(radii) == 1 || length(energies) == 1) && return data[i, j]
    tq = (Q - radii[i]) / (radii[i+1] - radii[i])
    te = (E - energies[j]) / (energies[j+1] - energies[j])
    return (1 - tq) * (1 - te) * data[i, j] + tq * (1 - te) * data[i+1, j] +
           (1 - tq) * te * data[i, j+1] + tq * te * data[i+1, j+1]
end
