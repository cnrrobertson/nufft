cd(@__DIR__)
using Pkg
Pkg.activate("../")

using Plots
using Random

# Define function
f(x) = sin(x)

# Generate nonuniform points with random seed
M = 64
seed = 1
rng = MersenneTwister(seed)
nu_ts = sort(rand(rng, M) .* 2pi)
nu_ys = f.(nu_ts)

# Generate uniform points with same length
u_ts = range(0, 2pi, length=M)
u_ys = f.(u_ts)

# %%
# Plot points
P = scatter(nu_ts, nu_ys, label="Nonuniform")
scatter!(u_ts, u_ys, label="Uniform")
display(P)