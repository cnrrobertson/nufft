# Type-2 F(k[j]): uniform input, f(x[i]: nonuniform output 
tolerance = 1e-12
n = 256
t = rand(Float64, (n,1)) # random numbers in [0,1)
# scale xs between [-\pi, \pi]
x = Array{Float64}(undef,n,1)
for i in 1:n
  x[i] = (t[i]*2.0 - 1.0)*pi
end
sort!(x, dims=1)

m = 90
F = rand(Float64, (m,1)) # random numbers in [0,1)

println(tolerance)
println(n)
for i in 1:n
  println(x[i])
end
println(m)
for i in 1:m
  println(F[i])
end
