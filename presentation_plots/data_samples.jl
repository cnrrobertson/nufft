cd(@__DIR__)
using Pkg
Pkg.activate("../")

using Random

# n = 1028
# m = 256
# tolerance = 1e-12
# seed = 1

# Function choices
# freq = 1.0
# c = 3.0
# f(x) = c * sin(freq * x)
# f(x) = cos(x^2)

function generate_data(f, n, m, tolerance, seed)
    input_file = "input_1dt1.dat"

    rng = MersenneTwister(seed)
    t = rand(rng, Float64, (n,1)) # random numbers in [0,1)
    # scale xs between [-\pi, \pi]
    x = (2pi .* t) .- pi
    sort!(x, dims=1)
    fx = f.(x)
    open(input_file, "w") do io
        println(io, tolerance)
        println(io, n)
        for i in 1:n
            println(io, x[i], ' ', fx[i])
        end

        # frequency
        println(io, m)
    end
    # println("Data saved at $input_file")
end

function load_data(input_file="input_1dt1.dat")
    input = read(input_file, String)
    input = split(strip(input), "\n")
    tol = parse(Float64, input[1])
    n = parse(Int, input[2])
    data = zeros(Float64, n, 2)
    for i in 1:n
        nums = split(strip(input[i+2]), ' ')
        data[i, 1] = parse(Float64, nums[1])
        data[i, 2] = parse(Float64, nums[2])
    end
    m = parse(Int, input[end])
    return tol, n, data, m
end