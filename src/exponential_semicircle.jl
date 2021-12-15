function exponential_semicircle(x, center, beta)
    if abs(x-center) <= 1
        return exp(beta*(sqrt(1-(x-center)^2) - 1))
    else
        return 0
    end
end