cd(@__DIR__)
using Pkg
Pkg.activate("../")

using Plots
using Interpolations
using LaTeXStrings
using FINUFFT

include("../src/nufft.jl")
include("data_samples.jl")
include("nufft_tests.jl")

# %%
f(x) = cos(x^2)
n = 256
m = n
tol = 1e-12
seed = 2
generate_data(f, n, m, tol, seed)
_, _, data, _ = load_data()
nu_xs = data[1:end, 1]
nu_ys = data[1:end, 2]
u_xs = range(minimum(nu_xs), maximum(nu_xs), length=n)
u_ys = f.(u_xs)

# Fit and smooth
ys_fit = LinearInterpolation(nu_xs, nu_ys, extrapolation_bc=Line())
interp = ys_fit.(u_xs)

es = zeros(size(u_xs)...)
beta = .01
for (x, y) in zip(nu_xs, nu_ys)
    es .+= y .* nufft.exponential_semicircle.(u_xs, x, beta)
end
es ./= maximum(abs.(es))

gaussian = zeros(size(u_xs)...)
tau = .002
for (x, y) in zip(nu_xs, nu_ys)
    gaussian .+= y .* nufft.gaussian_window.(u_xs, x, tau)
end
gaussian ./= maximum(abs.(gaussian))

P = scatter(nu_xs, nu_ys, label="Nonuniform points", xlabel=L"x", ylabel=L"f(x)")
plot!(u_xs, interp, linestyle=:dash, label="Interpolation")
# plot!(u_xs, es, label="Exponential-Semicircle")
plot!(u_xs, gaussian, linestyle=:dash, label="Gaussian")
savefig(P, "../latex/presentation/images/conv_vs_interp.png")
display(P)

# %%
# Compare outputs
fk_i = nufft1d_interpolation(nu_xs, nu_ys, m)
fk_i ./= maximum(abs.(fk_i))
fk_es = nufft1d_es(nu_xs, nu_ys, m, tol)
fk_es ./= maximum(abs.(fk_es))
fk_gaussian, _ = nufft1d_gaussian()
fk_gaussian ./= maximum(abs.(fk_gaussian))
fk_dir, _ = direct1dt1()
fk_dir ./= maximum(abs.(fk_dir))

ks = -(m ÷ 2):((m ÷ 2) - 1)
p1 = plot(ks, abs.(fk_dir), label="Directe", xlabel=L"k", ylabel=L"\hat{f}(x)")
p2 = plot(ks, abs.(fk_i), label="Interpolation", xlabel=L"k", ylabel=L"\hat{f}(x)")
p3 = plot(ks, abs.(fk_es), label="Exponential-Semicircle", xlabel=L"k", ylabel=L"\hat{f}(x)")
p4 = plot(ks, abs.(fk_gaussian), label="Gaussian", xlabel=L"k", ylabel=L"\hat{f}(x)")
savefig(P, "../latex/presentation/images/spectrum.png")
P = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600))
display(P)