# Berry phases, Wilson loops, Berry curvature and Chern numbers of magnon bands.
#
# Magnons are bosonic Bogoliubov quasiparticles: their eigenvectors in the
# 2L-dimensional Nambu space are para-unitary, T† Ĩ T = Ĩ with
# Ĩ = diag(1, …, 1, -1, …, -1). The Berry connection must use this metric,
# A_n(q) = i ⟨u_n(q)| Ĩ |∂_q u_n(q)⟩ (Shindou et al., PRB 87, 174427 (2013));
# the plain inner product gives wrong results.
#
# Loops that close by a reciprocal lattice vector G reuse the eigenvector of the
# starting point, which is valid when the dynamical matrix is periodic,
# H(q + G) = H(q), in Sunny's gauge (bond phases exp(2πi q⋅n) only). This holds
# for reciprocal lattice vectors of the magnetic cell and, for couplings that
# only connect sites within planes, also for in-plane vectors of the original
# cell. It is checked numerically and a warning is given otherwise.
#
# All wavevectors are 3-vectors in RLU of the original crystal.

_bdg_metric(L::Integer) = Diagonal([ones(L); -ones(L)])

"""
    band_vectors(swt, q) -> (energies, U)

The `L` positive magnon energies at `q` (RLU), in Sunny's band order
(descending), and a `2L × L` matrix whose columns are the corresponding
para-unitarily normalized Bogoliubov vectors.
"""
function band_vectors(swt::SpinWaveTheory, q)
    L = Sunny.nbands(swt)
    E, T = excitations(swt, Vec3(q...))
    return E[1:L], T[:, 1:L]
end

# Is H(q + G) = H(q) for Sunny's dynamical matrix? Tested at a generic q.
function _is_period(swt, G; rtol=1e-8)
    norm(G) < 1e-12 && return true
    L = Sunny.nbands(swt)
    H1, H2 = zeros(ComplexF64, 2L, 2L), zeros(ComplexF64, 2L, 2L)
    q0 = Vec3(0.1234, 0.2718, 0.0577)
    Sunny.dynamical_matrix!(H1, swt, Sunny.to_reshaped_rlu(swt.sys, q0))
    Sunny.dynamical_matrix!(H2, swt, Sunny.to_reshaped_rlu(swt.sys, q0 + Vec3(G...)))
    return norm(H1 - H2) ≤ rtol * max(norm(H1), 1.0)
end

_is_integer(v; atol=1e-8) = all(x -> abs(x - round(x)) < atol, v)

# Loop points without a duplicated end point. If the last point equals the
# first up to a reciprocal lattice vector, it is dropped and the loop is closed
# with the starting eigenvector (checking periodicity); otherwise the loop is
# closed by a straight link from the last point back to the first.
function _loop_points(swt, qloop)
    qs = [Vec3(q...) for q in qloop]
    length(qs) ≥ 3 || throw(ArgumentError("a loop needs at least 3 points"))
    Δ = qs[end] - qs[1]
    if _is_integer(Δ) || _is_integer(Sunny.to_reshaped_rlu(swt.sys, Δ))
        _is_period(swt, Δ) || @warn "The loop closes by Δq = $(Tuple(Δ)), but the spin-wave " *
            "Hamiltonian is not periodic under this shift in Sunny's gauge. The Berry phase " *
            "will not be quantized; use a reciprocal lattice vector of the magnetic cell."
        return qs[1:end-1]
    end
    steps = norm.(diff(qs))
    norm(Δ) > 10 * median(steps) && @warn "The last loop point is far from the first " *
        "(|Δq| = $(round(norm(Δ); sigdigits=3)) vs typical step $(round(median(steps); sigdigits=3))); " *
        "the loop is closed by one long straight link. Check that `qloop` is a closed loop."
    return qs
end

function _check_gap(E, bands, gap_tol, where)
    others = setdiff(axes(E, 1), bands)
    isempty(others) && return Inf
    gap = minimum(abs(E[i, j] - E[o, j]) for i in bands, o in others, j in axes(E, 2))
    gap < gap_tol && @warn "Band(s) $(collect(bands)) come within $(round(gap; sigdigits=3)) of " *
        "another band $where; the result is not well defined. Use a loop that avoids the " *
        "crossing, or include all touching bands in `berry_phase_multiband`."
    return gap
end

"""
    berry_phase(swt, qloop, n; wrap=false, gap_tol=1e-6) -> γ

Discrete (Abelian) Berry phase in radians of band `n` around the loop sampled
by `qloop` (a vector of RLU 3-vectors), using the bosonic metric. The result
is gauge invariant.

The loop is closed automatically from the last point back to the first. If
the last point equals the first up to a reciprocal lattice vector (a Zak
phase across the Brillouin zone), it is identified with the first point.
Sample densely enough that each link is short, and keep band `n` gapped from
the others along the loop (a warning is given otherwise; use
[`berry_phase_multiband`](@ref) for touching bands).

With `wrap = false` the signed sum of link phases is returned (stable near
±π); `wrap = true` reduces it to `(-π, π]`. For loops that wind around the
Brillouin zone the Zak phase depends on the unit-cell origin; here it is the
value in Sunny's gauge, without intra-cell positions.
"""
function berry_phase(swt::SpinWaveTheory, qloop::AbstractVector, n::Integer;
                     wrap::Bool=false, gap_tol::Real=1e-6)
    L = Sunny.nbands(swt)
    1 ≤ n ≤ L || throw(ArgumentError("band index must be in 1:$L"))
    M = _bdg_metric(L)
    qs = _loop_points(swt, qloop)
    modes = [band_vectors(swt, q) for q in qs]
    _check_gap(reduce(hcat, first.(modes)), n:n, gap_tol, "along the loop")
    vs = [U[:, n] for (_, U) in modes]
    push!(vs, vs[1])
    γ = -sum(angle(vs[k]' * M * vs[k+1]) for k in 1:length(vs)-1)
    return wrap ? rem2pi(γ, RoundNearest) : γ
end

"""
    berry_phase_multiband(swt, qloop, bands; gap_tol=1e-6) -> (total, eigphases)

Non-Abelian Wilson loop of the subspace spanned by `bands` (e.g. `[1, 2]`)
around `qloop`, for bands that are degenerate or cross along the loop.
Returns the total Berry phase `-Im log det W` and the Wilson-loop eigenphases,
both in `(-π, π]`. Loop closing works as in [`berry_phase`](@ref); the
subspace must stay gapped from the remaining bands.
"""
function berry_phase_multiband(swt::SpinWaveTheory, qloop::AbstractVector,
                               bands::AbstractVector{<:Integer}; gap_tol::Real=1e-6)
    L = Sunny.nbands(swt)
    all(b -> 1 ≤ b ≤ L, bands) || throw(ArgumentError("band indices must be in 1:$L"))
    M = _bdg_metric(L)
    qs = _loop_points(swt, qloop)
    modes = [band_vectors(swt, q) for q in qs]
    _check_gap(reduce(hcat, first.(modes)), bands, gap_tol, "along the loop")
    Vs = [U[:, bands] for (_, U) in modes]
    push!(Vs, Vs[1])
    W = Matrix{ComplexF64}(I, length(bands), length(bands))
    for k in 1:length(Vs)-1
        W = W * (Vs[k]' * M * Vs[k+1])
    end
    eigphases = sort(rem2pi.(-angle.(eigvals(W)), RoundNearest))
    total = rem2pi(-imag(log(complex(det(W)))), RoundNearest)
    return total, eigphases
end

"""
    loop_around(center, u, v; n=60, r=0.05) -> Vector

Closed parallelogram loop of RLU points around `center`, spanned by the RLU
directions `u` and `v` with half-extent `r`, with `n` points per side (the
end point is not repeated). Useful to enclose a band-touching point, around
which a Dirac magnon gives a Berry phase of ±π.
"""
function loop_around(center, u, v; n::Integer=60, r::Real=0.05)
    c, U, V = Vec3(center...), r * Vec3(u...), r * Vec3(v...)
    corners = [c - U - V, c + U - V, c + U + V, c - U + V]
    pts = Vec3[]
    for s in 1:4
        a, b = corners[s], corners[mod1(s + 1, 4)]
        for t in range(0, 1; length=n + 1)[1:end-1]
            push!(pts, a + t * (b - a))
        end
    end
    return pts
end

# Band-n vectors and energies on the (N+1)×(N+1) grid origin + (i g1 + j g2)/N.
function _curvature(swt, n, N, origin, g1, g2)
    L = Sunny.nbands(swt)
    1 ≤ n ≤ L || throw(ArgumentError("band index must be in 1:$L"))
    o, G1, G2 = Vec3(origin...), Vec3(g1...), Vec3(g2...)
    for (name, G) in (("g1", G1), ("g2", G2))
        _is_period(swt, G) || @warn "The spin-wave Hamiltonian is not periodic under " *
            "$name = $(Tuple(G)) in Sunny's gauge, so the region is not a closed torus and the " *
            "total flux is not quantized. Use reciprocal lattice vectors of the magnetic cell."
    end
    V = Matrix{Vector{ComplexF64}}(undef, N + 1, N + 1)
    E = zeros(L, (N + 1)^2)
    for i in 0:N, j in 0:N
        Ek, U = band_vectors(swt, o + (i / N) * G1 + (j / N) * G2)
        V[i+1, j+1] = U[:, n]
        E[:, i * (N + 1) + j + 1] .= Ek
    end
    M = _bdg_metric(L)
    link(a, b) = (z = a' * M * b; z / abs(z))
    F = zeros(N, N)
    for i in 1:N, j in 1:N
        F[i, j] = angle(link(V[i, j], V[i+1, j]) * link(V[i+1, j], V[i+1, j+1]) *
                        link(V[i+1, j+1], V[i, j+1]) * link(V[i, j+1], V[i, j]))
    end
    centres = ((1:N) .- 0.5) ./ N
    return F, centres, centres, E
end

"""
    berry_curvature(swt, n; N=80, origin=(0,0,0), g1=(1,0,0), g2=(0,1,0))
        -> (F, s1, s2)

Lattice Berry curvature of band `n` over the parallelogram
`origin + s g1 + t g2`, `s, t ∈ [0, 1]` (RLU), on an `N × N` grid of plaquettes
(Fukui–Hatsugai–Suzuki). `F[i, j]` is the Berry flux through plaquette
`(i, j)`, in `(-π, π]`, and `s1`, `s2` are the fractional plaquette centres
along `g1` and `g2`. Over a full Brillouin zone, `sum(F) = 2π C`.

The defaults cover the in-plane `a*`–`b*` cell of the original crystal.
A warning is given if `g1` or `g2` is not a period of the spin-wave
Hamiltonian.
"""
function berry_curvature(swt::SpinWaveTheory, n::Integer; N::Integer=80,
                         origin=(0, 0, 0), g1=(1, 0, 0), g2=(0, 1, 0))
    F, s1, s2, _ = _curvature(swt, n, N, origin, g1, g2)
    return F, s1, s2
end

"""
    chern_number(swt, n; N=80, origin=(0,0,0), g1=(1,0,0), g2=(0,1,0),
                 gap_tol=1e-6) -> Float64

Chern number `C = Σ F / 2π` of band `n` over the region of
[`berry_curvature`](@ref). For a band that is gapped everywhere and a region
that is a full period, `C` is an integer for any reasonable `N`. Warnings are
given when the band touches a neighbour on the grid (then `C` is not defined;
gap it, e.g. with a DM interaction), when the region is not a period, or when
`C` is not close to an integer (increase `N`).
"""
function chern_number(swt::SpinWaveTheory, n::Integer; N::Integer=80, origin=(0, 0, 0),
                      g1=(1, 0, 0), g2=(0, 1, 0), gap_tol::Real=1e-6)
    F, _, _, E = _curvature(swt, n, N, origin, g1, g2)
    gap = _check_gap(E, n:n, gap_tol, "on the grid")
    C = sum(F) / 2π
    if gap ≥ gap_tol && abs(C - round(C)) > 1e-3
        @warn "Chern number C = $(round(C; digits=4)) is not close to an integer; increase N " *
              "or check that the region covers a full Brillouin zone."
    end
    return C
end
