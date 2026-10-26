# PowderTools.jl

Tools for powder inelastic neutron scattering with
[Sunny.jl](https://github.com/SunnySuite/Sunny.jl), plus magnon band topology.

Main features:

- **Powder spectra**: one-magnon powder averages with preferred orientation
  (March–Dollase texture), and the two-magnon continuum: joint density of
  states, continuum edges, and the longitudinal two-magnon intensity for
  single crystals and powders.
- **Cuts**: constant-|Q| and constant-E cuts through scattered `(|Q|, E, I)`
  data or model maps, with error bars.
- **Fitting**: fit spin-wave models to powder data by comparing the shape of
  the spectrum (scale and background are fitted automatically).
- **Topology**: Berry phases, Wilson loops, Berry curvature and Chern numbers
  of magnon bands, with the correct bosonic (para-unitary) metric.
- **Plotting** (with a Makie backend): powder maps, scattered data and cuts.

## Installation

The package is not registered. From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/daniflamingo/PowderTools.jl")
```

or, from a local copy, `Pkg.develop(path="path/to/PowderTools")`.

Requires Julia ≥ 1.10 and Sunny 0.7 or 0.8.

## Quick example

```julia
using Sunny, PowderTools

# Easy-axis honeycomb antiferromagnet
cryst = Crystal(lattice_vectors(3.5, 3.5, 12, 90, 90, 120), [[1/3, 2/3, 0]], 191)
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
set_onsite_coupling!(sys, S -> -0.2 * S[3]^2, 1)
randomize_spins!(sys); minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))

radii, energies = range(0.1, 3.0, 60), range(0, 8, 160)
texture = march_dollase(0.7; axis=cryst.latvecs[:, 3])

one = powder_one_magnon(swt, radii, energies; fwhm=0.15, texture)
two = powder_two_magnon(swt, radii, energies; fwhm=0.15, texture)
cut = make_cut(one + two; dir=:E, center=1.5, width=0.1)

# Fit a model p -> System to measured (|Q|, E, I) arrays
data = PowderData(Q, E, I; Erange=(0.3, 4.0))
model = powder_model(build_system; radii, energies, fwhm=0.15)
fit = fit_powder(model, p0, data)
```

More in [`examples/`](examples): powder spectra and cuts
(`powder_spectra.jl`), a fit to synthetic data (`fit_synthetic.jl`), and
Dirac and Chern magnons on the honeycomb lattice (`topology.jl`).

## Overview

| Function | Purpose |
|:--|:--|
| `powder_one_magnon(swt, radii, energies; fwhm, texture, formfactor)` | Texture-weighted one-magnon powder spectrum |
| `powder_two_magnon(swt, radii, energies; fwhm, texture, formfactor)` | Powder-averaged longitudinal two-magnon spectrum |
| `two_magnon_intensity(swt, q, energies; fwhm)` | Two-magnon spectrum at a single wavevector |
| `two_magnon_jdos(swt, q, energies; fwhm)` | Two-magnon joint density of states |
| `two_magnon_edges(swt, q)` | Lower and upper edges of the two-magnon continuum |
| `march_dollase(r; axis)`, `fibonacci_sphere(n)` | Preferred orientation weight, sphere sampling |
| `PowderMap` | Intensity on a (\|Q\|, E) grid; maps can be added and scaled |
| `make_cut(q, E, I; dir, center, width, bins)` | 1D cut through scattered data or a `PowderMap` |
| `PowderData(Q, E, I; Qrange, Erange)` | Scattered data with a fit window |
| `powder_model(build_system; radii, energies, fwhm)` | Turn `p -> System` into a powder model |
| `fit_powder(model, p0, data)` | Fit by shape correlation (Nelder–Mead) |
| `model_at`, `shape_cost`, `scale_background` | Building blocks of the fit |
| `band_vectors(swt, q)` | Magnon energies and Bogoliubov vectors |
| `berry_phase`, `berry_phase_multiband`, `loop_around` | Berry phase / Wilson loop around a closed loop |
| `berry_curvature`, `chern_number` | Berry curvature and Chern number (Fukui–Hatsugai–Suzuki) |
| `plot_powder`, `plot_cut`, `plot_cut!` | Plotting (Makie) |

Every function has a docstring: type `?powder_two_magnon` in the REPL.
Wavevectors are in reciprocal lattice units of the crystal the system was
built from, as in Sunny.

## Running the tests

```julia
using Pkg
Pkg.test("PowderTools")
```

## Acknowledgements

Built on [Sunny.jl](https://github.com/SunnySuite/Sunny.jl). If you use this
package, please also cite Sunny:
D. Dahlbom *et al.*, *Sunny.jl: A Julia Package for Spin Dynamics*,
J. Open Source Softw. **10**, 8138 (2025),
[doi:10.21105/joss.08138](https://doi.org/10.21105/joss.08138).
The Berry connection follows R. Shindou *et al.*, Phys. Rev. B **87**, 174427
(2013), and the Chern number T. Fukui, Y. Hatsugai and H. Suzuki,
J. Phys. Soc. Jpn. **74**, 1674 (2005).

## License

MIT. See [LICENSE](LICENSE).
