# Directions on the sphere and preferred-orientation (texture) weights.

"""
    fibonacci_sphere(n; rotate=0.0) -> Vector{Vec3}

`n` quasi-uniform, equal-weight unit vectors covering the full sphere
(spherical Fibonacci lattice). `rotate` shifts the azimuth of all points.
"""
function fibonacci_sphere(n::Integer; rotate::Real=0.0)
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1"))
    golden = π * (3 - sqrt(5))
    return map(0:n-1) do i
        z = 1 - 2(i + 0.5) / n
        r = sqrt(max(0.0, 1 - z^2))
        θ = golden * i + rotate
        Vec3(r * cos(θ), r * sin(θ), z)
    end
end

"""
    random_rotation(rng=Random.default_rng()) -> Matrix{Float64}

Haar-random proper rotation (3×3, determinant +1). Used to decorrelate the
Fibonacci points between `|Q|` shells.
"""
function random_rotation(rng::AbstractRNG=Random.default_rng())
    F = qr(randn(rng, 3, 3))
    R = Matrix(F.Q) * Diagonal(sign.(diag(F.R)))
    det(R) < 0 && (R[:, 1] .*= -1)
    return R
end

"""
    march_dollase(r; axis=(0,0,1)) -> d -> W(d)

March–Dollase preferred-orientation weight as a function of a unit
wavevector direction `d` (global Cartesian):

    W(α) = (r² cos²α + sin²α / r)^(-3/2)

where `α` is the angle between `d` and the texture `axis` (e.g. the plate
normal; for the crystal `c` axis use `axis = cryst.latvecs[:, 3]`). `r = 1` is
isotropic, `r < 1` enhances `Q ∥ axis` (platelets), `r > 1` enhances
`Q ⊥ axis` (needles). `W` averages to 1 over the sphere for any `r`.

Pass the result as `texture` to [`powder_one_magnon`](@ref) or
[`powder_two_magnon`](@ref).
"""
function march_dollase(r::Real; axis=(0, 0, 1))
    r > 0 || throw(ArgumentError("March–Dollase r must be positive (got $r)"))
    t = normalize(Vec3(axis...))
    return d -> begin
        c2 = dot(Vec3(d...), t)^2 / dot(d, d)
        (r^2 * c2 + (1 - c2) / r)^(-3/2)
    end
end

# Texture weight of direction d (`nothing` = isotropic).
_weight(texture, d) = isnothing(texture) ? 1.0 : float(texture(d))
