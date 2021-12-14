cd(@__DIR__)
using Pkg
Pkg.activate("../")

using FFTW
using Plots
using Interpolations

# %%
# Nonuniformly and uniformly sampled sin function data
M = 1024
nu_ts = sort(rand(M) .* 2pi)
u_ts = range(0, 2pi, length=M)
nu_ys = sin.(nu_ts)
u_ys = sin.(u_ts)

# Linear interpolation of nonuniform points
ys_fit = LinearInterpolation(nu_ts, nu_ys, extrapolation_bc=Line())

# Comparison of data and interpolation
P = plot(u_ts, u_ys, label="Original")
plot!(u_ts, ys_fit.(u_ts), label="Interpolation")
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