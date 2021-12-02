sr = 256
ts = 1.0/(sr-1)
t = 0:ts:1

freq = 1.0
x = 3 .* sin.((2pi*freq) .* t)

freq = 4.0
x .+= sin.((2pi*freq) .* t)

freq = 7.0
x .+= 0.5 .* sin.((2pi*freq) .* t)

# Test fft
X = nufft.fft(x)
sizes = reverse(sortperm(abs.(X)))
@test sizes[1] == 2 || sizes[2] == 2
@test sizes[3] == 5 || sizes[4] == 5
@test sizes[5] == 8 || sizes[6] == 8

# Test ifft
X = nufft.fft(x)
new_x = real.(nufft.ifft(X))
@test maximum(abs.(x .- new_x)) < 1e-12