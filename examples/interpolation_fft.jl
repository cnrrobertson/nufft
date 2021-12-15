cd(@__DIR__)
using Pkg
Pkg.activate("../")

using FFTW
using Plots
using Interpolations
using Random
using FINUFFT
using LaTeXStrings

# %%
# Nonuniformly and uniformly sampled sin function data
M = 64
rng = MersenneTwister(1)
nu_ts = sort(rand(M) .* 2pi)
u_ts = range(0, 2pi, length=M)
nu_ys = sin.(nu_ts)
u_ys = sin.(u_ts)

# Linear interpolation of nonuniform points
ys_fit = LinearInterpolation(nu_ts, nu_ys, extrapolation_bc=Line())

# Comparison of data and interpolation
P = scatter(nu_ts, nu_ys, label="Original")
# P = scatter(u_ts, ys_fit.(u_ts), label="Original")
plot!(u_ts, ys_fit.(u_ts), label="Interpolation")
savefig(P, "../latex/writeup/images/interp.png")
display(P)

# %%
# Compare fft of uniform data with fft of resampled points
true_fk = fft(u_ys)
approx_fk = fft(ys_fit.(u_ys))
ks = -(M ÷ 2):((M ÷ 2) - 1)

p1 = plot(ks, fftshift(abs.(true_fk)), label="FFT of uniform points")
p2 = plot(ks, fftshift(abs.(approx_fk)), label="FFT of interpolated points")
P = plot(p1, p2, layout=(2, 1))
display(P)

# %%
# Compare with FINUFFT
true_fk = fft(u_ys)
approx_fk = fft(ys_fit.(u_ys))
approx_fk2 = nufft1d1(nu_ts, nu_ys .+ 1im*0, -1, 1e-10, 2*M)
ks = -(M ÷ 2):((M ÷ 2) - 1)
ks2 = -(M):(M - 1)

p1 = plot(ks, fftshift(abs.(true_fk)), label="FFT of uniform points")
p2 = plot(ks, fftshift(abs.(approx_fk)), label="FFT of interpolated points")
p3 = plot(ks2[M÷2:end-(M÷2)], abs.(approx_fk2)[M÷2:end-(M÷2)], label="FFT of smoothed points", xlabel=L"k", ylabel=L"|\hat{f}\;|")
P = plot(p1, p2, p3, layout=(3, 1))
display(P)
savefig(P, "../latex/writeup/images/fft_comparison.png")

# %%
# Speed plot
Ms = 2 .^ (6:18)
itimes = []
estimes = []
for M in Ms
    time1 = time()
    rng = MersenneTwister(1)
    nu_ts = sort(rand(M) .* 2pi)
    u_ts = range(0, 2pi, length=M)
    nu_ys = sin.(nu_ts)
    u_ys = sin.(u_ts)

    # Linear interpolation of nonuniform points
    ys_fit = LinearInterpolation(nu_ts, nu_ys, extrapolation_bc=Line())
    approx_fk = fft(ys_fit.(u_ys))
    time2 = time()
    push!(itimes, time2-time1)

    # FINUFFT of nonuniform points
    time1 = time()
    approx_fk = nufft1d1(nu_ts, nu_ys .+ 1im*0, -1, 1e-10, M)
    time2 = time()
    push!(estimes, time2-time1)
end

P = plot(Ms, itimes, label="Interpolation time")
plot!(Ms, estimes, label="Kernel time", xlabel="Number of samples", ylabel="Time (seconds)")
savefig(P, "../latex/writeup/images/speed.png")
display(P)