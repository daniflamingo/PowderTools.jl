"""
    PowderTools

Powder inelastic neutron scattering with [Sunny.jl](https://github.com/SunnySuite/Sunny.jl):
texture-weighted one- and two-magnon powder spectra, cuts through scattered
`(|Q|, E, I)` data, fitting of spin-wave models to powder data, and magnon
band topology (Berry phases, Berry curvature, Chern numbers).

Wavevectors are in reciprocal lattice units (RLU) of the crystal the system
was built from, as in Sunny. `|Q|` is in inverse length, energies in the
system's energy units.

Plotting helpers become available after loading a Makie backend
(`using GLMakie` or `using CairoMakie`).
"""
module PowderTools

using Sunny
using LinearAlgebra
using Random
using Statistics
import Optim

const Vec3 = Sunny.Vec3

include("common.jl")
include("sampling.jl")
include("formfactors.jl")
include("powdermap.jl")
include("onemagnon.jl")
include("twomagnon.jl")
include("cuts.jl")
include("fitting.jl")
include("topology.jl")
include("plotting.jl")

# Sampling and texture
export fibonacci_sphere, random_rotation, march_dollase

# Form factors
export form_factor_ce3

# Powder spectra
export PowderMap, powder_one_magnon
export two_magnon_jdos, two_magnon_edges, two_magnon_intensity, powder_two_magnon

# Cuts
export Cut, make_cut

# Fitting
export PowderData, model_at, shape_cost, scale_background, powder_model, fit_powder

# Topology
export band_vectors, berry_phase, berry_phase_multiband, loop_around
export berry_curvature, chern_number

# Plotting (implemented in the Makie extension)
export plot_cut, plot_cut!, plot_powder

end # module
