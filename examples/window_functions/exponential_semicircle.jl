cd(@__DIR__)
using Pkg
Pkg.activate("../../")
include("../../src/nufft.jl")
using Plots
using FFTW
using Interpolations
using Random

# %%
# Look at exponential semicircle function
ts = range(0, 2pi, length=1000)
beta = 10
center = 2
tol = 1e-5
ys = nufft.exponential_semicircle.(ts, center, beta)
zs = float.(ys .< tol)
P = plot(ts, ys, label="Exponential semicircle window function")
display(P)

# %%
# Smoothing for some nonuniform data
rng = MersenneTwister(1)
ts = sort(rand(rng, 64) .* 2pi)
ys = sin.(ts)

P = scatter(ts, ys, label=nothing)
display(P)
# %%
P = scatter(ts, ys, label="Original")
beta = 10
full_ts = range(0, 2pi, length=2*length(ts))
full_ys = zeros(size(full_ts)...)
for (t, y) in zip(ts, ys)
    full_ys .+= y .* nufft.exponential_semicircle.(full_ts, t, beta)
end
# plot!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="ES beta=$beta")
plot!(full_ts, full_ys ./ 6, label="Exponential of semicircle")
# savefig(P, "../../latex/writeup/images/es.png")
display(P)

# %%
# Data
# P = plot(ts, ys)
Pscatter!(ts, ys, label="Nonuniform points")
for beta in [1, 2, 3, 4, 10]
    full_ts = range(0, 2pi, length=1000)
    full_ys = zeros(size(full_ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.exponential_semicircle.(full_ts, t, beta)
    end
    plot!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="ES beta=$beta")
end
display(P)

# %%
# Smoothing for dense nonuniform data
ts = sort(rand(1028) .* 2pi)
ys = sin.(ts)

# Data
P = plot(ts, sin.(ts), label="True")
for beta in [.01, .1, 1, 2, 3]
    full_ts = range(0, 2pi, length=length(ts)*2)
    full_ys = zeros(size(full_ts)...)
    for (t, y) in zip(ts, ys)
        full_ys .+= y .* nufft.exponential_semicircle.(full_ts, t, beta)
    end
    # lines!(full_ts, abs.(sin.(full_ts) .- full_ys ./ maximum(abs.(full_ys))), label="Difference beta=$beta")
    plot!(full_ts, full_ys ./ maximum(abs.(full_ys)), label="beta=$beta")
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
beta = .001
full_ts = range(0, 2pi, length=R*M)
full_ys = zeros(size(full_ts)...)
for (t, y) in zip(ts, ys)
    full_ys .+= (y .* nufft.exponential_semicircle.(full_ts, t, beta))
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
full_yhats .*= sqrt(pi/beta) .* exp.((nu_ks .^ 2) .* beta)
full_yhats = fftshift(full_yhats)
# full_ys = real.(ifft(full_yhats))

p1 = plot(u_ks, real.(yhats))
plot!(u_ks, imag.(yhats))
p2 = plot(nu_ks, real.(full_yhats))
plot!(nu_ks, imag.(full_yhats))
P = plot(p1, p2, layout=(2, 1))
display(P)
