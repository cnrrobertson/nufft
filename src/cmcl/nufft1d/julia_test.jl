cd(@__DIR__)
using Pkg
Pkg.activate("../../../")

nj = 32
xj = sort(rand(nj)) .* 2pi
cj = sin.(xj) .+ 1im*0
iflag = -1
ms = 16
fk0 = zeros(ComplexF64, ms)

ccall((:dirft1d1_, "./nufft1d.so"), Vector{ComplexF64}, (Int64, Vector{Float64}, Vector{ComplexF64}, Int64, Int64, Vector{ComplexF64}), nj, xj, cj, iflag, ms, fk0)
