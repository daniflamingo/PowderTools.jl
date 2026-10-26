using Test
using LinearAlgebra
using Random
using Sunny
using PowderTools

# Honeycomb ferromagnet (moments along -z in a field along +z, g = 2),
# optionally with a Kane–Mele-type DM interaction on second neighbours.
function honeycomb_fm(; J=-1.0, B=0.1, D=0.0, s=1/2, dims=(1, 1, 1))
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[1/3, 2/3, 0]], 191)
    sys = System(cryst, [1 => Moment(; s, g=2)], :dipole; dims, seed=0)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    iszero(D) || set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [1, 0, 0]))
    set_field!(sys, [0, 0, B])
    polarize_spins!(sys, [0, 0, -1])
    return sys
end

# Easy-axis honeycomb Néel antiferromagnet (gapped, with quantum fluctuations).
function honeycomb_afm(; J=1.0, A=0.2)
    cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[1/3, 2/3, 0]], 191)
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; seed=0)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_onsite_coupling!(sys, S -> -A * S[3]^2, 1)
    set_dipole!(sys, [0, 0, 1], (1, 1, 1, 1))
    set_dipole!(sys, [0, 0, -1], (1, 1, 1, 2))
    return sys
end

swt_of(sys) = SpinWaveTheory(sys; measure=ssf_perp(sys))

@testset "PowderTools" begin

    @testset "sampling and texture" begin
        pts = fibonacci_sphere(500)
        @test all(p -> norm(p) ≈ 1, pts)
        @test norm(sum(pts)) / 500 < 1e-2
        R = random_rotation(Xoshiro(3))
        @test R' * R ≈ I
        @test det(R) ≈ 1
        dense = fibonacci_sphere(4000)
        for r in (0.6, 1.7)
            W = march_dollase(r; axis=(1, 1, 0))
            @test sum(W, dense) / length(dense) ≈ 1 rtol=1e-2
        end
        @test march_dollase(1.0)([0.3, 0.4, 0.5]) ≈ 1
        @test_throws ArgumentError march_dollase(-1)
    end

    @testset "form factors" begin
        @test form_factor_ce3(0) ≈ 1 atol=1e-3
        @test issorted([form_factor_ce3(Q) for Q in 0:0.5:5]; rev=true)
        @test PowderTools._formfactor(nothing, 2.0) == 1
        @test PowderTools._formfactor(Q -> 0.5, 2.0) == 0.5
        @test PowderTools._formfactor(FormFactor("Fe2"), 0.0) ≈ 1 atol=1e-2
    end

    @testset "cuts" begin
        q = [x for x in 0:0.01:1, _ in 0:0.1:2]
        E = [y for _ in 0:0.01:1, y in 0:0.1:2]
        I = copy(E)
        I[1, 1] = NaN                                   # ignored
        c = make_cut(q, E, I; dir=:E, center=0.5, width=0.105, bins=-0.05:0.1:2.05)
        @test c.x ≈ 0:0.1:2
        @test c.I ≈ 0:0.1:2
        @test all(==(11), c.counts)
        @test all(e -> abs(e) < 1e-12, c.err)          # identical values in each bin
        cs = make_cut(q, E, I; dir=:E, center=0.5, width=0.105, bins=-0.05:0.1:2.05, reduce=:sum)
        @test cs.I ≈ 11 .* (0:0.1:2)
        cq = make_cut(q, E, q; dir=:q, center=1.0, width=0.05, bins=10)
        @test all(abs.(cq.I .- cq.x) .< 0.01)           # each bin mean within 0.01 of its centre
        @test length(make_cut(q, E, I; dir=:E, center=0.5, width=0.2).x) ≥ 5
        @test_throws ArgumentError make_cut(q, E, I; dir=:E, center=5.0, width=0.1)
        @test_throws ArgumentError make_cut(q, E, I; dir=:x, center=0.5, width=0.1)

        m = PowderMap(0:0.01:1, 0:0.1:2, E)
        cm = make_cut(m; dir=:E, center=0.5, width=0.105, bins=-0.05:0.1:2.05)
        @test cm.I ≈ 0:0.1:2
    end

    @testset "PowderMap and interpolation" begin
        radii, energies = 0:0.5:2, 0:1:3
        m = PowderMap(radii, energies, [2r + 3e for r in radii, e in energies])
        data = PowderData([0.3, 1.7, 5.0], [0.5, 2.2, 1.0], [0.0, 0.0, 0.0])
        @test model_at(m, data) ≈ [2 * 0.3 + 3 * 0.5, 2 * 1.7 + 3 * 2.2, 0.0]
        @test (m + m).data ≈ 2 .* m.data
        @test (0.5 * m).data ≈ 0.5 .* m.data
        @test_throws DimensionMismatch PowderMap(radii, energies, zeros(2, 2))
        @test_throws ArgumentError m + PowderMap(radii, 0:1:2, zeros(5, 3))
    end

    @testset "one-magnon powder" begin
        swt = swt_of(honeycomb_fm())
        radii, energies = range(0.5, 3.0, 6), 0:0.05:4
        m = powder_one_magnon(swt, radii, energies; fwhm=0.2, n_sphere=300)
        @test size(m.data) == (6, length(energies))
        @test all(≥(-1e-12), m.data)

        # Agrees with Sunny's powder_average (different random rotations)
        res = powder_average(Sunny.orig_crystal(swt.sys), radii, 300) do qs
            intensities(swt, qs; energies, kernel=gaussian(fwhm=0.2))
        end
        @test sum(m.data) ≈ sum(PowderMap(res).data) rtol=0.05

        # Isotropic texture and form factor scaling are exact
        mt = powder_one_magnon(swt, radii, energies; fwhm=0.2, n_sphere=300,
                               texture=march_dollase(1.0))
        @test mt.data ≈ m.data
        mf = powder_one_magnon(swt, radii, energies; fwhm=0.2, n_sphere=300, formfactor=Q -> 0.5)
        @test mf.data ≈ 0.25 .* m.data
    end

    @testset "two-magnon continuum" begin
        sys = honeycomb_fm()
        swt = swt_of(sys)
        EΓ = dispersion(swt, [[0, 0, 0]])[:, 1]

        # Edges at Γ: both extremes are reached at k = 0
        lo, hi = two_magnon_edges(swt, [0, 0, 0]; ngrid=6)
        @test lo ≈ 2 * minimum(EΓ)
        @test hi ≈ 2 * maximum(EΓ)

        # Same result for a magnetic supercell: wavevectors are in original RLU
        swt2 = swt_of(honeycomb_fm(dims=(2, 1, 1)))
        @test Sunny.nbands(swt2) == 4
        lo2, hi2 = two_magnon_edges(swt2, [0, 0, 0]; ngrid=6)
        @test lo2 ≈ lo
        @test hi2 ≈ hi

        # JDOS sum rule: ∫ D dω = L²
        ωs = 0:0.005:7.5
        D = two_magnon_jdos(swt, [1/3, 1/3, 0], ωs; fwhm=0.05, nsamples=400)
        @test sum(D) * step(ωs) ≈ 4 rtol=1e-3
        D2 = two_magnon_jdos(swt2, [1/3, 1/3, 0], ωs; fwhm=0.05, nsamples=200)
        @test sum(D2) * step(ωs) ≈ 16 rtol=1e-3

        # A collinear ferromagnet has no two-magnon (longitudinal) intensity
        I2 = two_magnon_intensity(swt, [0.2, 0.1, 0], ωs; fwhm=0.05, nsamples=100)
        @test maximum(abs, I2) < 1e-10

        # An antiferromagnet does, and it is non-negative
        swta = swt_of(honeycomb_afm())
        Ia = two_magnon_intensity(swta, [0.2, 0.1, 0], 0:0.02:12; fwhm=0.1, nsamples=200)
        @test all(≥(-1e-12), Ia)
        @test sum(Ia) > 0

        radii, energies = [0.8, 1.6], 0:0.1:12
        P = powder_two_magnon(swta, radii, energies; fwhm=0.3, n_sphere=20, n_k=20)
        @test size(P.data) == (2, length(energies))
        @test all(≥(-1e-12), P.data) && sum(P.data) > 0
        Pt = powder_two_magnon(swta, radii, energies; fwhm=0.3, n_sphere=20, n_k=20,
                               texture=march_dollase(1.0))
        @test Pt.data ≈ P.data
        Pf = powder_two_magnon(swta, radii, energies; fwhm=0.3, n_sphere=20, n_k=20,
                               formfactor=form_factor_ce3)
        @test Pf.data[1, :] ≈ form_factor_ce3(0.8)^2 .* P.data[1, :]

        sun = System(Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]]),
                     [1 => Moment(s=1, g=2)], :SUN)
        set_field!(sun, [0, 0, 1])
        polarize_spins!(sun, [0, 0, -1])
        @test_throws ArgumentError two_magnon_intensity(swt_of(sun), [0.1, 0, 0], 0:1; fwhm=0.1, nsamples=1)
    end

    @testset "fitting helpers" begin
        model = [0.1, 0.5, 0.2, 0.9, 0.4]
        data = PowderData(1:5, 1:5, 3 .* model .+ 0.1)
        @test shape_cost(model, data) ≈ 0 atol=1e-12
        @test shape_cost(-model, data) ≈ 2 atol=1e-12
        s, b = scale_background(model, data)
        @test s ≈ 3 && b ≈ 0.1
        s, b = scale_background(-model, data; clamp_positive=true)
        @test s == 0 && b ≈ sum(data.I) / 5

        d = PowderData([0.5 1.0; 1.5 2.0], [1 2; 3 4], [1.0 NaN; 3.0 4.0]; Erange=(0, 3.5))
        @test length(d) == 2                         # NaN and E = 4 dropped

        # Unstable parameters get the worst-case cost; other errors propagate
        unstable = p -> error("Not an energy-minimum; wavevector q = [0, 0, 0] unstable.")
        bad = @test_logs (:warn, r"unstable") fit_powder(unstable, [1.0], data; iterations=2)
        @test bad.cost == 2.0 && isnothing(bad.model)
        @test_throws ErrorException fit_powder(p -> error("bug"), [1.0], data; iterations=2)
    end

    @testset "synthetic powder fit" begin
        build(p) = honeycomb_fm(J=-p[1], B=p[2])
        model = powder_model(build; radii=range(0.2, 3.0, 40), energies=range(0, 4.5, 90),
                             fwhm=0.2, n_sphere=60, restarts=0)
        p_true = [1.0, 0.3]
        rng = Xoshiro(7)
        Q, E = 0.3 .+ 2.5 .* rand(rng, 600), 0.1 .+ 4.0 .* rand(rng, 600)
        clean = model_at(model(p_true), PowderData(Q, E, zeros(600)))
        I = 5 .* clean .+ 0.2 .+ 0.01 * maximum(clean) .* randn(rng, 600)
        data = PowderData(Q, E, I)

        fit = fit_powder(model, [0.8, 0.45], data; iterations=150)
        @test fit.params ≈ p_true rtol=0.05
        @test fit.cost < 0.05
        @test fit.scale ≈ 5 rtol=0.1
    end

    @testset "topology" begin
        swt = swt_of(honeycomb_fm())
        E, U = band_vectors(swt, [0.1, 0.2, 0])
        L = Sunny.nbands(swt)
        M = Diagonal([ones(L); -ones(L)])
        @test U' * M * U ≈ I

        # Dirac magnons: Berry phase ±π around K
        K = [1/3, 1/3, 0]
        loop = loop_around(K, [1, 0, 0], [0, 1, 0]; n=40, r=0.05)
        γ1 = berry_phase(swt, loop, 1; wrap=true)
        γ2 = berry_phase(swt, loop, 2; wrap=true)
        @test abs(γ1) ≈ π atol=1e-3
        @test abs(γ2) ≈ π atol=1e-3
        total, eig = berry_phase_multiband(swt, loop, [1, 2])
        @test abs(rem2pi(total, RoundNearest)) < 1e-6
        @test length(eig) == 2

        # Explicitly repeated end point, also shifted by a reciprocal lattice vector
        @test berry_phase(swt, [loop; [loop[1]]], 1; wrap=true) ≈ γ1
        @test berry_phase(swt, [loop; [loop[1] + [1, 0, 0]]], 1; wrap=true) ≈ γ1

        # Without DM the bands touch at K: the Chern number is not defined
        @test_logs (:warn, r"come within") match_mode=:any chern_number(swt, 1; N=12)

        # Kane–Mele-type DM gaps the Dirac points: Chern numbers ±1
        swtD = swt_of(honeycomb_fm(D=0.1))
        C1 = chern_number(swtD, 1; N=24)
        C2 = chern_number(swtD, 2; N=24)
        @test abs(C1 - round(C1)) < 1e-6
        @test abs(round(C1)) == 1
        @test round(C1) == -round(C2)
        F, s1, s2 = berry_curvature(swtD, 1; N=24)
        @test sum(F) / 2π ≈ C1
        @test length(s1) == length(s2) == 24

        # Half a Brillouin zone is not a period
        @test_logs (:warn, r"not periodic") match_mode=:any chern_number(swtD, 1; N=12, g1=(1/2, 0, 0))
    end
end
