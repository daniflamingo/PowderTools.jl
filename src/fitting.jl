# Fitting spin-wave models to powder data.
#
# The model intensities are in arbitrary units, so the default cost compares
# only the *shape* of the model and the data: 1 - (weighted Pearson
# correlation). An overall scale and a flat background are nuisance parameters
# that the correlation removes automatically; `scale_background` recovers them
# afterwards for plotting and χ².

"""
    PowderData(Q, E, I; weights=nothing, Qrange=(-Inf, Inf), Erange=(-Inf, Inf))

Scattered powder data: equal-size arrays of `|Q|`, energy and intensity
(any shape, e.g. the matrices of a time-of-flight data set). Points outside
`Qrange`/`Erange`, with non-finite values, or with zero weight are dropped.
`weights` defaults to 1; for Poisson-like statistics use e.g. `1 ./ (I .+ ε)`.
"""
struct PowderData
    Q::Vector{Float64}
    E::Vector{Float64}
    I::Vector{Float64}
    w::Vector{Float64}
end

function PowderData(Q::AbstractArray, E::AbstractArray, I::AbstractArray; weights=nothing,
                    Qrange=(-Inf, Inf), Erange=(-Inf, Inf))
    size(Q) == size(E) == size(I) || throw(DimensionMismatch("Q, E, I must have the same size"))
    w = isnothing(weights) ? ones(length(I)) : vec(collect(Float64, weights))
    length(w) == length(I) || throw(DimensionMismatch("weights must match the data size"))
    Q, E, I = vec(collect(Float64, Q)), vec(collect(Float64, E)), vec(collect(Float64, I))
    keep = isfinite.(Q) .& isfinite.(E) .& isfinite.(I) .& (w .> 0) .&
           (Qrange[1] .≤ Q .≤ Qrange[2]) .& (Erange[1] .≤ E .≤ Erange[2])
    any(keep) || throw(ArgumentError("no data points left after applying the windows"))
    return PowderData(Q[keep], E[keep], I[keep], w[keep])
end

Base.length(d::PowderData) = length(d.I)
Base.show(io::IO, d::PowderData) =
    print(io, "PowderData(", length(d), " points; |Q| ∈ ", extrema(d.Q), ", E ∈ ", extrema(d.E), ")")

"""
    model_at(map::PowderMap, data::PowderData) -> Vector

Model intensities at the data points, by bilinear interpolation of `map`
(zero outside its grid).
"""
model_at(m::PowderMap, data::PowderData) = [_interp(m, Q, E) for (Q, E) in zip(data.Q, data.E)]

function _weighted_cor(x, y, w)
    ws = sum(w)
    mx, my = sum(w .* x) / ws, sum(w .* y) / ws
    cxy = sum(w .* (x .- mx) .* (y .- my))
    cxx = sum(w .* (x .- mx) .^ 2)
    cyy = sum(w .* (y .- my) .^ 2)
    (cxx ≤ 0 || cyy ≤ 0) && return 0.0      # flat model or data
    return cxy / sqrt(cxx * cyy)
end

"""
    shape_cost(model, data::PowderData) -> Float64

`1 - ρ`, with `ρ` the weighted Pearson correlation between model and measured
intensities: 0 for a perfect match up to scale and background, 1 for no
correlation, 2 for perfect anti-correlation. `model` is a `PowderMap` or a
vector of intensities at the data points.
"""
shape_cost(model::AbstractVector, data::PowderData) = 1 - _weighted_cor(model, data.I, data.w)
shape_cost(m::PowderMap, data::PowderData) = shape_cost(model_at(m, data), data)

"""
    scale_background(model, data::PowderData; clamp_positive=false) -> (scale, background)

Weighted least-squares fit `data ≈ scale * model + background`. With
`clamp_positive = true`, a negative scale (model anti-correlated with the
data) is replaced by the background-only fit, `scale = 0`.
"""
function scale_background(model::AbstractVector, data::PowderData; clamp_positive::Bool=false)
    X = hcat(model, ones(length(model)))
    W = Diagonal(data.w)
    s, b = (X' * W * X) \ (X' * W * data.I)
    if clamp_positive && s < 0
        s, b = 0.0, sum(data.w .* data.I) / sum(data.w)
    end
    return s, b
end
scale_background(m::PowderMap, data::PowderData; kwargs...) = scale_background(model_at(m, data), data; kwargs...)

# Lowest-energy classical state over `restarts` random starts. Each restart
# uses a fresh system from `build_system(p)`.
function _ground_state(build_system, p; restarts::Integer, minimize)
    restarts == 0 && return build_system(p)
    best, Ebest = nothing, Inf
    for _ in 1:restarts
        sys = build_system(p)
        randomize_spins!(sys)
        minimize_energy!(sys; minimize...)
        E = energy_per_site(sys)
        E < Ebest && ((best, Ebest) = (sys, E))
    end
    return best
end

"""
    powder_model(build_system; radii, energies, fwhm, n_sphere=100, restarts=3,
                 measure=ssf_perp, texture=nothing, formfactor=nothing,
                 seed=0, minimize=(;)) -> (p -> PowderMap)

Turn a function `build_system(p)`, which returns a Sunny `System` with
couplings set from the parameter vector `p`, into a powder model for
[`fit_powder`](@ref). For each `p` the model

1. finds the classical ground state from `restarts` random starts
   (`restarts = 0` uses the system as returned, e.g. if `build_system` already
   sets the ordered state; `minimize` is passed to `minimize_energy!`),
2. builds `SpinWaveTheory(sys; measure=measure(sys))`, and
3. powder-averages on the `radii × energies` grid with `n_sphere` directions
   per shell: Sunny's `powder_average` for an isotropic sample, or
   [`powder_one_magnon`](@ref) when a `texture` is given. A `formfactor`
   (see [`powder_one_magnon`](@ref)) is applied as `f(|Q|)²`.

`seed` fixes the sphere sampling, so the model is a deterministic function of
`p`, which the optimizer needs.
"""
function powder_model(build_system; radii, energies, fwhm::Real, n_sphere::Integer=100,
                      restarts::Integer=3, measure=ssf_perp, texture=nothing,
                      formfactor=nothing, seed::Integer=0, minimize=(;))
    radii, energies = collect(Float64, radii), collect(Float64, energies)
    return function (p)
        sys = _ground_state(build_system, p; restarts, minimize)
        swt = SpinWaveTheory(sys; measure=measure(sys))
        if isnothing(texture)
            kernel = gaussian(; fwhm)
            res = powder_average(Sunny.orig_crystal(sys), radii, n_sphere; seed) do qs
                intensities(swt, qs; energies, kernel)
            end
            m = PowderMap(res)
            f2 = [_formfactor(formfactor, Q)^2 for Q in radii]
            return PowderMap(radii, energies, f2 .* m.data)
        else
            return powder_one_magnon(swt, radii, energies; fwhm, texture, formfactor,
                                     n_sphere, rng=Random.Xoshiro(seed))
        end
    end
end

# Cost assigned to parameters whose classical state is unstable in LSWT.
const COST_MAX = 2.0

"""
    fit_powder(model, p0, data::PowderData; cost=shape_cost,
               method=Optim.NelderMead(), iterations=200, show_trace=false)

Fit model parameters to powder data. `model(p)` returns a [`PowderMap`](@ref)
(typically from [`powder_model`](@ref)); `cost(map, data)` defaults to
[`shape_cost`](@ref). Parameters for which linear spin-wave theory is unstable
get the worst-case cost (2); any other error is raised. The default
derivative-free Nelder–Mead suits the noisy, possibly discontinuous cost
(ground states can switch between parameter sets).

Returns a named tuple `(; params, cost, scale, background, chi2, model, result)`
with the best-fit parameters, their cost, the scale and background of
[`scale_background`](@ref), the weighted χ² of `scale * model + background`,
the best-fit `PowderMap` (`nothing`, with `NaN` scale, background and χ², if
the model is unstable there), and the full `Optim` result.
"""
function fit_powder(model, p0, data::PowderData; cost=shape_cost,
                    method=Optim.NelderMead(), iterations::Integer=200, show_trace::Bool=false)
    function objective(p)
        try
            return float(cost(model(p), data))
        catch e
            _is_instability(e) || rethrow()
            return COST_MAX
        end
    end
    result = Optim.optimize(objective, collect(Float64, p0), method,
                            Optim.Options(; iterations, show_trace))
    params = Optim.minimizer(result)
    best = try
        model(params)
    catch e
        _is_instability(e) || rethrow()
        nothing
    end
    if isnothing(best)
        @warn "Linear spin-wave theory is unstable at the best parameters found; " *
              "try another starting point."
        scale = background = chi2 = NaN
    else
        vals = model_at(best, data)
        scale, background = scale_background(vals, data; clamp_positive=true)
        chi2 = sum(data.w .* (scale .* vals .+ background .- data.I) .^ 2)
    end
    return (; params, cost=Optim.minimum(result), scale, background, chi2, model=best, result)
end
