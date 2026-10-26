# One- and two-magnon powder spectra of an easy-axis honeycomb antiferromagnet,
# with preferred orientation, and cuts through the result.
#
#   julia --project=. examples/powder_spectra.jl
# (add GLMakie to your environment to see the plots)

using Sunny, PowderTools, LinearAlgebra

cryst = Crystal(lattice_vectors(3.5, 3.5, 12, 90, 90, 120), [[1/3, 2/3, 0]], 191)
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole; seed=0)
set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
set_onsite_coupling!(sys, S -> -0.2 * S[3]^2, 1)
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))

radii = range(0.1, 3.0, 60)        # |Q| in 1/Å
energies = range(0, 8, 160)        # meV

# Platelets with the c axis as plate normal (r < 1 favours Q ∥ c)
texture = march_dollase(0.7; axis=cryst.latvecs[:, 3])

one = powder_one_magnon(swt, radii, energies; fwhm=0.15, texture)
two = powder_two_magnon(swt, radii, energies; fwhm=0.15, texture, n_sphere=60, n_k=60)
println(one)
println(two)

# Continuum edges along a path
for q in ([0, 0, 0], [1/3, 1/3, 0], [1/2, 0, 0])
    lo, hi = two_magnon_edges(swt, q; ngrid=20)
    println("q = $q: two-magnon continuum from $(round(lo; digits=3)) to $(round(hi; digits=3)) meV")
end

# Constant-|Q| cut through the combined spectrum
total = one + 0.5 * two
cut = make_cut(total; dir=:E, center=1.5, width=0.1, bins=0:0.05:8)

const HAVE_GLMAKIE = Base.find_package("GLMakie") !== nothing
HAVE_GLMAKIE && @eval using GLMakie

if HAVE_GLMAKIE
    display(plot_powder(one; title="one-magnon"))
    display(plot_powder(two; title="two-magnon"))
    display(plot_cut(cut, make_cut(one; dir=:E, center=1.5, width=0.1, bins=0:0.05:8);
                     labels=["one + two", "one only"]))
end
