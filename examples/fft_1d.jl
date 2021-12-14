cd(@__DIR__)
using Pkg
Pkg.activate("../")
include("../src/nufft.jl")

using Plots

# %%
sr = 256
ts = 1.0/(sr-1)
t = 0:ts:1

freq = 1.0
x = 3 .* sin.((2pi*freq) .* t)

freq = 4.0
x .+= sin.((2pi*freq) .* t)

freq = 7.0
x .+= 0.5 .* sin.((2pi*freq) .* t)

X = nufft.fft(x)
sizes = reverse(sortperm(abs.(X)))

# %%
# Plot result
N = length(X)
n = 0:(N-1)
T = N/(sr-1)
freq = n ./ T
P = plot(freq, abs.(X), line=:stem)
display(P)

# %%
# Plot data
P = plot(t, x, label="True")
approx = real.(nufft.ifft(X))
plot!(t, approx, label="ifft(fft(true))")
display(P)

# %%
# Plot difference
approx = real.(nufft.ifft(X))
P = plot(t, abs.(x .- approx), label="abs difference")
display(P)