# Two-magnon continuum in linear spin-wave theory.
#
# All wavevectors are in RLU of the original crystal, as in Sunny. Internal
# magnon momenta k are sampled uniformly over one Brillouin zone of the
# magnetic cell used by `swt`, so the results do not depend on how the system
# was reshaped.

const _SV3C = Sunny.SVector{3, ComplexF64}

# Dipole modes only (:dipole, and :dipole_large_s, which Sunny ≥ 0.7.3 stores
# as :dipole_uncorrected): one boson per site, so band s ↔ site s.
function _check_dipole(swt)
    swt.sys.mode == :SUN && throw(ArgumentError(
        "two-magnon functions require a system in a dipole mode, not :SUN"))
end

# Ordered-moment unit vectors n̂_s and positions r_s (global) of the magnetic cell.
function _frames(swt)
    sys = swt.sys
    L = Sunny.nbands(swt)
    ns = [normalize(sys.dipoles[1, 1, 1, s]) for s in 1:L]
    rs = [global_position(sys, (1, 1, 1, s)) for s in 1:L]
    return ns, rs
end

_random_k(rng, C) = C * Vec3(rand(rng), rand(rng), rand(rng))

# Neutron polarization projector δ - Q̂Q̂ for a global wavevector (identity at Q = 0).
function _polarization(Qg)
    n = norm(Qg)
    n < 1e-12 && return Sunny.Mat3(I)
    Qh = Vec3(Qg) / n
    return Sunny.Mat3(I) - Qh * Qh'
end

# Longitudinal two-magnon amplitudes A[a, b] (3-vectors, global frame) for
# creating magnons (k, band a) and (q - k, band b). The dynamics use the
# reduced wavevectors; the intra-cell phases use the full momentum transfer
# `Qg` (global Cartesian), which matters in the extended zone.
function _amplitudes!(A, ws1, ws2, swt, ns, rs, q, k, Qg)
    L = Sunny.nbands(swt)
    E1 = _modes!(ws1, swt, k)
    E2 = _modes!(ws2, swt, q - k)
    Tk, Tp = ws1.T, ws2.T
    for a in 1:L, b in 1:L
        acc = zero(_SV3C)
        for s in 1:L
            amp = cis(-dot(Qg, rs[s])) * (Tk[s, a] * Tp[s+L, b] + Tp[s, b] * Tk[s+L, a])
            acc += ns[s] * amp
        end
        A[a, b] = acc
    end
    return E1, E2
end

"""
    two_magnon_jdos(swt, q, energies; fwhm, nsamples=200_000, rng=Random.Xoshiro(1))

Two-magnon joint density of states at wavevector `q` (RLU),

    D(q, ω) = Σ_ab ⟨δ(ω - E_a(k) - E_b(q - k))⟩_k

with the matrix element set to 1, from a Monte Carlo average over the internal
momentum `k` in the magnetic Brillouin zone. The δ-function is a Gaussian of
width `fwhm`, so `∫ D dω = L²` (number of band pairs) when `energies` spans
the continuum.
"""
function two_magnon_jdos(swt::SpinWaveTheory, q, energies; fwhm::Real,
                         nsamples::Integer=200_000, rng::AbstractRNG=Random.Xoshiro(1))
    g = _gaussian(fwhm)
    C = _magnetic_to_orig(swt)
    ws1, ws2 = _Modes(swt), _Modes(swt)
    qv = Vec3(q...)
    hist = zeros(length(energies))
    for _ in 1:nsamples
        k = _random_k(rng, C)
        E1 = _modes!(ws1, swt, k)
        E2 = _modes!(ws2, swt, qv - k)
        for Ea in E1, Eb in E2, (i, ω) in enumerate(energies)
            hist[i] += g(ω - Ea - Eb)
        end
    end
    return hist ./ nsamples
end

"""
    two_magnon_edges(swt, q; ngrid=30) -> (lo, hi)

Lower and upper edges of the two-magnon continuum at `q` (RLU): the minimum
and maximum of `E_a(k) + E_b(q - k)` over all band pairs and an
`ngrid × ngrid × ngrid` grid of `k` in the magnetic Brillouin zone (which
includes `k = 0`).
"""
function two_magnon_edges(swt::SpinWaveTheory, q; ngrid::Integer=30)
    C = _magnetic_to_orig(swt)
    ws1, ws2 = _Modes(swt), _Modes(swt)
    qv = Vec3(q...)
    lo, hi = Inf, -Inf
    for i in 0:ngrid-1, j in 0:ngrid-1, l in 0:ngrid-1
        k = C * Vec3(i, j, l) / ngrid
        E1 = _modes!(ws1, swt, k)
        E2 = _modes!(ws2, swt, qv - k)
        lo = min(lo, minimum(E1) + minimum(E2))
        hi = max(hi, maximum(E1) + maximum(E2))
    end
    return lo, hi
end

"""
    two_magnon_intensity(swt, q, energies; fwhm, nsamples=200_000,
                         rng=Random.Xoshiro(1))

Longitudinal two-magnon scattering intensity at wavevector `q` (RLU), in
arbitrary units: Monte Carlo over the internal momentum `k` with the
amplitude for creating two magnons from the longitudinal spin fluctuations
`δS∥ ∝ a†a`, projected with the neutron polarization factor `δ - Q̂Q̂`.
A factor ½ corrects for counting each magnon pair twice. Requires a dipole
mode (`:dipole` or `:dipole_large_s`/`:dipole_uncorrected`), not `:SUN`.
Collinear ferromagnets give zero (no quantum fluctuations in the ground
state).
"""
function two_magnon_intensity(swt::SpinWaveTheory, q, energies; fwhm::Real,
                              nsamples::Integer=200_000, rng::AbstractRNG=Random.Xoshiro(1))
    _check_dipole(swt)
    g = _gaussian(fwhm)
    L = Sunny.nbands(swt)
    C = _magnetic_to_orig(swt)
    ns, rs = _frames(swt)
    ws1, ws2 = _Modes(swt), _Modes(swt)
    A = Matrix{_SV3C}(undef, L, L)
    qv = Vec3(q...)
    Qg = Sunny.orig_crystal(swt.sys).recipvecs * qv
    P = _polarization(Qg)
    hist = zeros(length(energies))
    for _ in 1:nsamples
        k = _random_k(rng, C)
        E1, E2 = _amplitudes!(A, ws1, ws2, swt, ns, rs, qv, k, Qg)
        for a in 1:L, b in 1:L
            w = real(dot(A[a, b], P * A[a, b]))
            Ω = E1[a] + E2[b]
            for (i, ω) in enumerate(energies)
                hist[i] += w * g(ω - Ω)
            end
        end
    end
    return hist .* (0.5 / nsamples)
end

"""
    powder_two_magnon(swt, radii, energies; fwhm, n_sphere=200, n_k=200,
                      texture=nothing, formfactor=nothing,
                      rng=Random.Xoshiro(1), randomize_shells=true) -> PowderMap

Powder average of [`two_magnon_intensity`](@ref): for each `|Q|` shell,
`n_sphere` Fibonacci directions (optionally texture-weighted, see
[`march_dollase`](@ref)) times `n_k` Monte Carlo samples of the internal
momentum. The full sphere is used (no symmetry folding), which matters when
the magnetic structure breaks inversion. The magnon energies and Bogoliubov
vectors use the reduced wavevector, while the intra-cell phases,
polarization factor and form factor use the full `Q`. Arbitrary units.

`formfactor` and `randomize_shells` work as in [`powder_one_magnon`](@ref).
"""
function powder_two_magnon(swt::SpinWaveTheory, radii, energies; fwhm::Real,
                           n_sphere::Integer=200, n_k::Integer=200, texture=nothing,
                           formfactor=nothing, rng::AbstractRNG=Random.Xoshiro(1),
                           randomize_shells::Bool=true)
    _check_dipole(swt)
    g = _gaussian(fwhm)
    L = Sunny.nbands(swt)
    C = _magnetic_to_orig(swt)
    recip = Sunny.orig_crystal(swt.sys).recipvecs
    ns, rs = _frames(swt)
    ws1, ws2 = _Modes(swt), _Modes(swt)
    A = Matrix{_SV3C}(undef, L, L)
    base = fibonacci_sphere(n_sphere)
    out = zeros(length(radii), length(energies))
    for (iQ, Q) in enumerate(radii)
        R = randomize_shells ? random_rotation(rng) : Matrix{Float64}(I, 3, 3)
        acc = zeros(length(energies))
        Wsum = 0.0
        for d0 in base
            d = Vec3(R * d0)
            W = _weight(texture, d)
            Wsum += W
            Qg = Q * d
            P = _polarization(Qg)
            q = recip \ Qg
            for _ in 1:n_k
                k = _random_k(rng, C)
                E1, E2 = _amplitudes!(A, ws1, ws2, swt, ns, rs, q, k, Qg)
                for a in 1:L, b in 1:L
                    w = W * real(dot(A[a, b], P * A[a, b]))
                    Ω = E1[a] + E2[b]
                    for (i, ω) in enumerate(energies)
                        acc[i] += w * g(ω - Ω)
                    end
                end
            end
        end
        out[iQ, :] .= _formfactor(formfactor, Q)^2 .* acc .* (0.5 / (n_k * Wsum))
    end
    return PowderMap(radii, energies, out)
end
