cd(@__DIR__)
using Pkg
Pkg.activate("../")

using Plots

include("nufft_tests.jl")
include("data_samples.jl")

# %%
# Setup
ns = 2 .^ (4:16)
tol = 1e-12
f(x) = cos(x^2)

# %%
# Run codes
errs = zeros(Float64, length(ns), 3)
for i in 1:length(ns)
    generate_data(f, ns[i], ns[i], tol, 1)
    tol, n, data, m = load_data()
    nu_xs = data[:, 1]
    nu_ys = data[:, 2]

    truth, _ = direct1dt1()
    truth ./= maximum(abs.(truth))

    gauss, _ = nufft1d_gaussian()
    gauss ./= maximum(abs.(gauss))
    top = sum(abs.(truth .- gauss) .^ 2)
    bot = sum(abs.(truth) .^ 2)
    e1 = sqrt(top / bot)

    interp = nufft1d_interpolation(nu_xs, nu_ys, m)
    interp ./= maximum(abs.(interp))
    top = sum(abs.(truth .- interp) .^ 2)
    bot = sum(abs.(truth) .^ 2)
    e2 = sqrt(top / bot)

    finu = nufft1d_es(nu_xs, nu_ys, m, tol)
    finu ./= maximum(abs.(finu))
    top = sum(abs.(truth .- finu) .^ 2)
    bot = sum(abs.(truth) .^ 2)
    e3 = sqrt(top / bot)

    errs[i, 1] = e1
    errs[i, 2] = e2
    errs[i, 3] = e3
end

# %%
P = plot(ns, errs[1:end, 1], label="Gaussian", xlabel="Number of points (N)", ylabel="Time (seconds)")
# plot!(ns, errs[1:end, 2], label="Interpolation")
plot!(ns, errs[1:end, 3], label="Exponential-Semicircle")
# savefig(P, "../latex/presentation/images/n_vs_time.png")
display(P)