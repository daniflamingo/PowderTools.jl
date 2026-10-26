# Fit a spin-wave model to (synthetic) powder data.
#
#   julia --project=. examples/fit_synthetic.jl
#
# With real data, replace the "synthetic data" block by your |Q|, E and I
# arrays, e.g. read from CSV files:
#   using DelimitedFiles
#   Q, E, I = readdlm("Q.csv", ','), readdlm("E.csv", ','), readdlm("I.csv", ',')

using Sunny, PowderTools, Random

# Model: honeycomb ferromagnet in a field, parameters p = (|J|, B)
cryst = Crystal(lattice_vectors(3.5, 3.5, 12, 90, 90, 120), [[1/3, 2/3, 0]], 191)
function build_system(p)
    J, B = p
    sys = System(cryst, [1 => Moment(s=1/2, g=2)], :dipole; seed=0)
    set_exchange!(sys, -J, Bond(1, 2, [0, 0, 0]))
    set_field!(sys, [0, 0, B])
    return sys
end

model = powder_model(build_system; radii=range(0.1, 2.5, 50), energies=range(0, 4.5, 90),
                     fwhm=0.15, n_sphere=80, restarts=2)

# Synthetic data: scattered (|Q|, E) points, arbitrary scale, background, noise
p_true = [1.0, 0.3]
rng = Xoshiro(1)
Q, E = 0.2 .+ 2.2 .* rand(rng, 2000), 0.2 .+ 3.8 .* rand(rng, 2000)
clean = model_at(model(p_true), PowderData(Q, E, zeros(2000)))
I = 40 .* clean .+ 0.5 .+ 0.05 * maximum(40 .* clean) .* randn(rng, 2000)

# Fit window excludes the elastic line
data = PowderData(Q, E, I; Erange=(0.3, 4.0))
fit = fit_powder(model, [0.7, 0.5], data; iterations=200, show_trace=false)
println("true parameters:   ", p_true)
println("fitted parameters: ", round.(fit.params; digits=4))
println("correlation = ", round(1 - fit.cost; digits=4),
        ", scale = ", round(fit.scale; digits=3), ", background = ", round(fit.background; digits=3))

const HAVE_GLMAKIE = Base.find_package("GLMakie") !== nothing
HAVE_GLMAKIE && @eval using GLMakie

if HAVE_GLMAKIE
    display(plot_powder(data; title="data"))
    display(plot_powder(fit.scale * fit.model; title="best fit"))
    cut_data = make_cut(data.Q, data.E, data.I; dir=:E, center=1.2, width=0.2, bins=0.3:0.1:4)
    cut_fit = make_cut(fit.scale * fit.model; dir=:E, center=1.2, width=0.2, bins=0.3:0.1:4)
    cut_fit = Cut(cut_fit.x, cut_fit.I .+ fit.background, cut_fit.err, cut_fit.counts, cut_fit.dir)
    display(plot_cut(cut_data, cut_fit; labels=["data", "fit"]))
end
