# Magnetic form factors as functions of |Q| (inverse Å).

# Dipole approximation ⟨j₀⟩ for Ce³⁺:
#   f(s) = A exp(-a s²) + B exp(-b s²) + C exp(-c s²) + D exp(-d s²) + E,  s = |Q|/4π
# Coefficients as tabulated for the "J" form factors of the `periodictable`
# package (Lisher & Forsyth, Acta Cryst. A27, 545 (1971)).
const _CE3_J0 = (A=0.051183, a=6.115375, B=0.277738, b=7.952485,
                 C=0.654079, c=2.287000, D=0.0, d=0.0, E=0.016355)

"""
    form_factor_ce3(Q)

Magnetic form factor of Ce³⁺ in the dipole approximation, `⟨j₀(|Q|)⟩`, for
`|Q|` in inverse Å. (Sunny ≤ 0.7 has no Ce³⁺ entry; Sunny 0.8 has one with
different coefficients.) Pass as `formfactor = form_factor_ce3`.
"""
function form_factor_ce3(Q::Real)
    (; A, a, B, b, C, c, D, d, E) = _CE3_J0
    s2 = (Q / 4π)^2
    return A * exp(-a * s2) + B * exp(-b * s2) + C * exp(-c * s2) + D * exp(-d * s2) + E
end

# The `formfactor` argument accepted throughout: `nothing` (no form factor), a
# function of |Q|, or a Sunny `FormFactor` (e.g. `FormFactor("Fe2")`).
_formfactor(::Nothing, Q) = 1.0
_formfactor(f::Sunny.FormFactor, Q) = Sunny.compute_form_factor(f, Float64(Q)^2)
_formfactor(f, Q) = float(f(Q))
