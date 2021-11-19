# cd(@__DIR__)
# using Pkg
# Pkg.activate("../")
# include("../src/nufft.jl")

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
@test sizes[1] == 2 || sizes[2] == 2
@test sizes[3] == 5 || sizes[4] == 5
@test sizes[5] == 8 || sizes[6] == 8

# # For plotting
# using CairoMakie

# # Plot data
# lines(t, x)

# # Plot result
# N = length(X)
# n = 0:(N-1)
# T = N/(sr-1)
# freq = n ./ T
# stem(freq, abs.(X))