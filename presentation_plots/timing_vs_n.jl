cd(@__DIR__)
using Pkg
Pkg.activate("../")

using Plots

include("nufft_tests.jl")
include("data_samples.jl")

# %%
# Setup
ns = 2 .^ (4:16)
tol = 1e-10
f(x) = cos(x^2)

# %%
# Run codes
times = zeros(Float64, length(ns), 4)
for i in 1:length(ns)
    generate_data(f, ns[i], ns[i], tol, 1)
    tol, n, data, m = load_data()
    nu_xs = data[:, 1]
    nu_ys = data[:, 2]
    time1 = time()
    _, t1 = direct1dt1()
    time2 = time()
    _, t2 = nufft1d_gaussian()
    time3 = time()
    nufft1d_interpolation(nu_xs, nu_ys, m)
    time4 = time()
    nufft1d_es(nu_xs, nu_ys, m, tol)
    time5 = time()

    times[i, 1] = t1
    times[i, 2] = t2
    times[i, 3] = time4-time3
    times[i, 4] = time5-time4
end

# %%
P = plot(ns, times[1:end, 1], yaxis=:log, xaxis=:log, label="Direct", xlabel="Number of points (N)", ylabel="Time (seconds)")
plot!(ns, times[1:end, 2], yaxis=:log, xaxis=:log, label="Gaussian")
plot!(ns, times[1:end, 3], yaxis=:log, xaxis=:log, label="Interpolation")
plot!(ns, times[1:end, 4], yaxis=:log, xaxis=:log, label="Exponential-Semicircle")
savefig(P, "../latex/presentation/images/n_vs_time.png")
display(P)