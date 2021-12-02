function fft(x)
    """
    Recursive implementation from:
    https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.03-Fast-Fourier-Transform.html
    """
    N = length(x)
    if N == 1
        return x
    else
        X_even = fft(x[1:2:N])
        X_odd = fft(x[2:2:N])
        factor = exp.((-2im*pi) .* ((0:(N-1)) ./ N))

        N2 = N ÷ 2
        X = [X_even .+ (factor[1:N2] .* X_odd); X_even .+ (factor[N2+1:N] .* X_odd)]
        return X
    end
end

function ifft(x)
    N = length(x)
    if N == 1
        return x
    else
        X_even = ifft(x[1:2:N])
        X_odd = ifft(x[2:2:N])
        factor = exp.((2im*pi) .* ((0:(N-1)) ./ N))

        N2 = N ÷ 2
        X = [X_even .+ (factor[1:N2] .* X_odd); X_even .+ (factor[N2+1:N] .* X_odd)]
        return X ./ 2
    end
end