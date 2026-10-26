# Dirac magnons and Chern bands on the honeycomb ferromagnet.
#
#   julia --project=. examples/topology.jl

using Sunny, PowderTools

cryst = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[1/3, 2/3, 0]], 191)
function honeycomb(D)
    sys = System(cryst, [1 => Moment(s=1/2, g=2)], :dipole; seed=0)
    set_exchange!(sys, -1.0, Bond(1, 2, [0, 0, 0]))
    iszero(D) || set_exchange!(sys, dmvec([0, 0, D]), Bond(1, 1, [1, 0, 0]))
    set_field!(sys, [0, 0, 0.1])
    polarize_spins!(sys, [0, 0, -1])
    return SpinWaveTheory(sys; measure=nothing)
end

# Without DM, the two bands touch at K with a Berry phase of ±π
swt = honeycomb(0.0)
loop = loop_around([1/3, 1/3, 0], [1, 0, 0], [0, 1, 0]; n=60, r=0.05)
for n in 1:2
    γ = berry_phase(swt, loop, n; wrap=true)
    println("D = 0:   band $n, Berry phase around K = $(round(γ / π; digits=4)) π")
end

# A Kane–Mele-type DM interaction gaps the Dirac points; the bands become
# topological with Chern numbers ±1
swt = honeycomb(0.1)
for n in 1:2
    println("D = 0.1: band $n, Chern number = $(round(chern_number(swt, n; N=60); digits=4))")
end

const HAVE_GLMAKIE = Base.find_package("GLMakie") !== nothing
HAVE_GLMAKIE && @eval using GLMakie

if HAVE_GLMAKIE
    F, s1, s2 = berry_curvature(swt, 1; N=120)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="s₁ (along a*)", ylabel="s₂ (along b*)", aspect=DataAspect(),
              title="Berry flux per plaquette, band 1")
    hm = heatmap!(ax, s1, s2, F; colormap=:balance)
    Colorbar(fig[1, 2], hm)
    display(fig)
end
