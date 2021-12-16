cd(@__DIR__)
using Pkg
Pkg.activate("../../../")

using Plots

# %%
#####################################################################
# Generating input
#####################################################################
input_file = "input_1dt1.dat"

# Type-1 f(x[i]): nonuniform input, F(k[j]): uniform output
tolerance = 1e-12
n = 256
m = 90
t = rand(Float64, (n,1)) # random numbers in [0,1)
# scale xs between [-\pi, \pi]
freq = 1.0
c = 3.0
x = (2pi .* t) .- pi
sort!(x, dims=1)
fx = c .* sin.(freq .* x)
open(input_file, "w") do io
    println(io, tolerance)
    println(io, n)
    for i in 1:n
        println(io, x[i], ' ', fx[i])
    end

    # frequency
    println(io, m)
end

#####################################################################
# Running fortran script
#####################################################################
output = read(`./nufft1dt1`, String)

#####################################################################
# Parse output
#####################################################################
output = split(strip(output), "\n")
ms = parse.(Int, output[1])
ier = output[end-2]
err = output[end-1]
output = output[2:end-3]
fk0 = zeros(ComplexF64, ms)
fk1 = zeros(ComplexF64, ms)
for i in 1:ms
    out = parse.(Float64, split(strip(output[i])))
    fk0[i] = out[1] + 1im*out[2]
end
for i in (ms+2):2ms
    out = parse.(Float64, split(strip(output[i])))
    fk1[i-ms-1] = out[1] + 1im*out[2]
end

# %%
# Plot results
P = plot(1:ms, abs.(fk0))
plot!(1:ms, abs.(fk1))
display(P)