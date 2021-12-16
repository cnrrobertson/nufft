cd(@__DIR__)
using Pkg
Pkg.activate("..")

using FFTW
using Interpolations
using FINUFFT

function direct1dt1()
    output = read(`./direct1dt1`, String)
    output = split(strip(output), "\n")
    ms = parse.(Int, output[1])
    fk0 = zeros(ComplexF64, ms)
    for i in 2:ms+1
        out = parse.(Float64, split(strip(output[i])))
        fk0[i-1] = out[1] + 1im*out[2]
    end
    ttime = parse.(Float64, output[end])
    return fk0, ttime
end

function nufft1d_gaussian()
    output = read(`./nufft1dt1`, String)
    output = split(strip(output), "\n")
    ms = parse.(Int, output[1])
    fk0 = zeros(ComplexF64, ms)
    for i in 2:ms+1
        out = parse.(Float64, split(strip(output[i])))
        fk0[i-1] = out[1] + 1im*out[2]
    end
    ttime = parse.(Float64, output[end])
    return fk0, ttime
end

function nufft1d_interpolation(nu_xs, nu_ys, m)
    ys_fit = LinearInterpolation(nu_xs, nu_ys, extrapolation_bc=Line())
    x_min = minimum(nu_xs)
    x_max = maximum(nu_xs)
    u_xs = range(x_min, x_max, m)
    u_ys = ys_fit.(u_xs)

    return fftshift(fft(u_ys))
end

function nufft1d_es(nu_xs, nu_ys, m, tol)
    u_fk = nufft1d1(nu_xs, nu_ys .+ 1im*0, -1, tol, m)
    return u_fk
end