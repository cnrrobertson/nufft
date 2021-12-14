cd(@__DIR__)
using Pkg
Pkg.activate("../../")

using FINUFFT
using FFTW
using Plots
using GaussQuadrature

# %%
# Nonuniformly sampled points
M = 128
nu_ts = sort(rand(M) .* 2pi)
u_ts = range(0, 2pi, M)
nu_ys = sin.(nu_ts)
u_ys = sin.(u_ts)
P = scatter(nu_ts, nu_ys)
display(P)

# %%
# Type 1 nufft - from nonuniform samples to uniform frequencies (2M)
tol = 1e-12
nu_fk = nufft1d1(nu_ts, nu_ys .+ 1im*0, 1, tol, 2M)
nu_ks = -M:(M-1)
p1 = plot(nu_ks, abs.(nu_fk), label="NUFFT type 1")

# Compare with fft of uniform points
u_fk = fft(u_ys)
u_ks = -(M ÷ 2):((M ÷ 2) - 1)
p2 = plot(u_ks, fftshift(abs.(u_fk)), label="FFT of uniform points", xlabel="k")
P = plot(p1, p2, layout=(2, 1))
display(P)

# %%
# Try to do inufft using type 2 nufft
nu_fs = nufft1d2(nu_ts, -1, tol, nu_fk)
P = plot(nu_ts, real.(nu_fs))
display(P)

# %%
# Use legendre nonuniform points
f(x) = sin(x)
M = 128
kmax = 2M
nu_ts,nu_ws = legendre(M)
nu_ts .*= pi # rescale
u_ts = range(0, 2pi, M)
nu_ys = f.(nu_ts)
u_ys = f.(u_ts)

P = scatter(nu_ts, nu_ys, label="sin at legendre points")
display(P)

# %%
tol = 1e-10
# nu_fk = nufft1d1(nu_ts .*dk , f.(nu_ts) .* nu_ws .+ 1im*0, 1, tol, kmax)
nu_fk = nufft1d1(nu_ts, nu_ys .+ 1im*0, 1, tol, kmax)
ks = -(kmax ÷ 2):((kmax ÷ 2) - 1)
p1 = plot(ks, abs.(nu_fk), label="FFT of sin at legendre points")
display(P)

u_fk = fft(u_ys)
u_ks = -(M ÷ 2):((M ÷ 2) - 1)
p2 = plot(u_ks, fftshift(abs.(u_fk)), label="FFT of uniform points", xlabel="k")
P = plot(p1, p2, layout=(2, 1))
display(P)