# Internal helpers shared by several files.

# Gaussian broadening with full width at half maximum `fwhm`, normalized to
# unit area. Matches Sunny's `gaussian(; fwhm)`.
function _gaussian(fwhm)
    fwhm > 0 || throw(ArgumentError("fwhm must be positive (got $fwhm)"))
    σ = fwhm / (2 * sqrt(2 * log(2)))
    return x -> exp(-x^2 / (2σ^2)) / (σ * sqrt(2π))
end

# Matrix mapping RLU of the magnetic cell used by `swt` (Sunny's "reshaped"
# cell) to RLU of the original crystal: q_orig = C * q_magnetic.
_magnetic_to_orig(swt) = Sunny.orig_crystal(swt.sys).recipvecs \ swt.sys.crystal.recipvecs

# Workspace for repeated Bogoliubov diagonalizations without allocation.
struct _Modes
    T::Matrix{ComplexF64}
    H::Matrix{ComplexF64}
end
_Modes(swt) = (L = Sunny.nbands(swt); _Modes(zeros(ComplexF64, 2L, 2L), zeros(ComplexF64, 2L, 2L)))

# Positive-branch energies at `q` (original RLU). Afterwards `ws.T` holds the
# para-unitary Bogoliubov matrix, [a; a†] = T [b; b†].
function _modes!(ws::_Modes, swt, q)
    E = excitations!(ws.T, ws.H, swt, Vec3(q...))
    return view(E, 1:Sunny.nbands(swt))
end

# Sunny signals an unstable classical state with an `ErrorException`; the
# wording changed in 0.7.3. Other errors must not be mistaken for it.
_is_instability(e) = e isa ErrorException &&
    (occursin("Instability at wavevector", e.msg) || occursin("Not an energy-minimum", e.msg))
