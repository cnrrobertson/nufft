cd(@__DIR__)
using Pkg
Pkg.activate("../../")
include("../../src/nufft.jl")
using Plots
using FFTW
using Interpolations
using Random

# %%
# Look at gaussian window function
ts = range(0, 2pi, length=1000)
tau = .01
center = 2
tol = 1e-5
ys = nufft.gaussian_window.(ts, center, tau)
zs = float.(ys .< tol)
P = plot(ts, ys, label="Gaussian window function")
display(P)

# %%
# Smoothing for some nonuniform data
rng = MersenneTwister(1)
ts = sort(rand(rng, 64) .* 2pi)
ys = sin.(ts)

P = scatter(ts, ys, label=nothing)
savefig(P, "../../latex/writeup/images/nu_points.png")
display(P)
# %%

# Data
P = plot(ts, ys)
scatter!(ts, ys, label="Nonuniform points")
for tau in [.00001, .001, .01, .1]
    full_ts = range(0, 2pi, length=1000)
    full_ys = zeros(size(full_ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.gaussian_window.(full_ts, t, tau)
    end
    plot!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="Gaussian tau=$tau")
end
display(P)

# %%
# Smoothing for dense nonuniform data
ts = sort(rand(1028) .* 2pi)
ys = sin.(ts)

# Data
P = plot(ts, sin.(ts), label="True")
for tau in [.01, .1, 1, 2, 3]
    full_ts = range(0, 2pi, length=length(ts)*2)
    full_ys = zeros(size(full_ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.gaussian_window.(full_ts, t, tau)
    end
    # lines!(full_ts, abs.(sin.(full_ts) .- full_ys ./ maximum(abs.(full_ys))), label="Difference tau=$tau")
    plot!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="tau=$tau")
end
display(P)

# %%
# FFT on smoothed dense nonuniform data
M = 1028
y(t) = sin(t) + cos(50t) - 10*sin(100t)
ts = sort(rand(M) .* 2pi)
ys = y.(ts)

# Data
R = 2
tau = .001
full_ts = range(0, 2pi, length=R*M)
full_ys = zeros(size(full_ts)...)
for (t, y) in zip(ts, ys)
    full_ys .+= (y .* nufft.gaussian_window.(full_ts, t, tau))
end
# full_ys ./= maximum(abs.(full_ys))

uys = y.(range(0, 2pi, length=M))
yhats = fftshift(fft(uys))
ks = 1:length(yhats)
# ys = real.(ifft(yhats))

full_yhats = fftshift(fft(full_ys))
full_ks = 1:length(full_yhats)
nu_ks = -(R*M ÷ 2):((R*M ÷ 2) - 1)
u_ks = -(M ÷ 2):((M ÷ 2) - 1)
full_yhats .*= sqrt(pi/tau) .* exp.((nu_ks .^ 2) .* tau)
full_yhats = fftshift(full_yhats)
# full_ys = real.(ifft(full_yhats))

p1 = plot(u_ks, real.(yhats))
plot!(u_ks, imag.(yhats))
p2 = plot(nu_ks, real.(full_yhats))
plot!(nu_ks, imag.(full_yhats))
P = plot(p1, p2, layout=(2, 1))
display(P)
