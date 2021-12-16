# Type-1 f(x[i]): nonuniform input, F(k[j]): uniform output
tolerance = 1e-12
n = 256
t = rand(Float64, (n,1)) # random numbers in [0,1)
# scale xs between [-\pi, \pi]
x = Array{Float64}(undef,n,1)
fx = Array{Float64}(undef,n,1)
freq = 1.0
c = 3.0
for i in 1:n
  x[i] = (t[i]*2.0 - 1.0)*pi
end
sort!(x, dims=1)
for i in 1:n
  fx[i] = c .* sin.(freq .* x[i])
end

println(tolerance)
println(n)
for i in 1:n
  println(x[i], ' ', fx[i])
end

# frequency
m = 90
println(m)
#k = Int(-m/2):1:Int(floor((m-1)/2))
