cd(@__DIR__)
using Pkg
Pkg.activate("../")
include("../src/nufft.jl")
using CairoMakie

# %%
# Look at window
ts = range(0, 2pi, length=1000)
tau = .01
center = 2
tol = 1e-5
ys = nufft.gaussian_window.(ts, center, tau)
zs = float.(ys .< tol)
lines(ts, ys, label="Window function")
lines!(ts, zs, label="Less than $tol")
axislegend()
current_figure()

# %%
# Smoothing for some nonuniform data
ts = sort(rand(20) .* 2pi)
ys = sin.(ts)

# Data
lines(ts, ys)
scatter!(ts, ys, label="Nonuniform points")
for tau in [.00001, .001, .01, .1]
    full_ts = range(0, 2pi, length=1000)
    full_ys = zeros(size(full_ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.gaussian_window.(full_ts, t, tau)
    end
    lines!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="Gaussian tau=$tau")
end
axislegend()
current_figure()

# %%
# Smoothing for dense nonuniform data
ts = sort(rand(2000) .* 2pi)
ys = sin.(ts)

# Data
fig = Figure()
Axis(fig[1,1])
for tau in [.01, .1, 1, 2, 3]
    # full_ts = range(0, 2pi, length=1000)
    full_ys = zeros(size(ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.gaussian_window.(ts, t, tau)
    end
    lines!(ts, abs.(ys .- full_ys ./ maximum(abs.(full_ys))), label="Difference tau=$tau")
end
axislegend()
current_figure()