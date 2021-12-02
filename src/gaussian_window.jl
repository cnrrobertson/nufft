function gaussian_window(x, center, tau)
    y = 0
    for l in -1:1
        denom = 4*tau
        num = ((x - center) - 2pi*l)^2
        y += exp(-num/denom)
    end
    return y
end