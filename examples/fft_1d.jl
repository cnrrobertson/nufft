cd(@__DIR__)
using Pkg
Pkg.activate("../")
include("../src/nufft.jl")

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

# For plotting
using CairoMakie

# Plot result
N = length(X)
n = 0:(N-1)
T = N/(sr-1)
freq = n ./ T
stem(freq, abs.(X))

# Plot data
lines(t, x, label="True")
approx = real.(nufft.ifft(X))
lines!(t, approx, label="ifft(fft(true))")
axislegend()
current_figure()

# Plot difference
approx = real.(nufft.ifft(X))
lines(t, abs.(x .- approx), label="abs difference")
axislegend()
current_figure()