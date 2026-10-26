"""
    powder_one_magnon(swt, radii, energies; fwhm, texture=nothing,
                      formfactor=nothing, n_sphere=200,
                      rng=Random.Xoshiro(1), randomize_shells=true) -> PowderMap

Texture-weighted one-magnon powder spectrum

    I(|Q|, ω) = f(|Q|)² Σ_d W(d) S(|Q| d, ω) / Σ_d W(d)

using Sunny's `intensities` (dispersion, polarization factor and g-tensors as
set by the `measure` of `swt`) on `n_sphere` Fibonacci directions `d` per
shell. `radii` are `|Q|` values in inverse length.

- `fwhm`: Gaussian energy resolution.
- `texture`: orientation weight `W(d)`, e.g. [`march_dollase`](@ref);
  `nothing` is an isotropic powder.
- `formfactor`: `nothing`, a function of `|Q|` (e.g. [`form_factor_ce3`](@ref))
  or a Sunny `FormFactor`; applied as `f(|Q|)²` per shell. Leave it `nothing`
  if the `measure` already includes form factors.
- `randomize_shells`: apply a fresh random rotation to the Fibonacci points
  of each shell, which removes `|Q|`-correlated sampling streaks.

Without texture this is equivalent to Sunny's `powder_average`.
"""
function powder_one_magnon(swt::SpinWaveTheory, radii, energies; fwhm::Real,
                           texture=nothing, formfactor=nothing, n_sphere::Integer=200,
                           rng::AbstractRNG=Random.Xoshiro(1), randomize_shells::Bool=true)
    recip = Sunny.orig_crystal(swt.sys).recipvecs
    base = fibonacci_sphere(n_sphere)
    kernel = gaussian(; fwhm)
    energies = collect(Float64, energies)
    out = zeros(length(radii), length(energies))
    for (iQ, Q) in enumerate(radii)
        R = randomize_shells ? random_rotation(rng) : Matrix{Float64}(I, 3, 3)
        dirs = [Vec3(R * d) for d in base]
        W = [_weight(texture, d) for d in dirs]
        qs = [recip \ (Q * d) for d in dirs]
        res = intensities(swt, qs; energies, kernel)        # data[iE, iq]
        out[iQ, :] .= _formfactor(formfactor, Q)^2 .* (real.(res.data) * W) ./ sum(W)
    end
    return PowderMap(radii, energies, out)
end
