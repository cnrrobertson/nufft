function kaiser_bessel(x, center, beta, I0=1)
    if abs(x-center) <= 1
        return I0*beta*sqrt(1 - (x-center)^2)
    else
        return 0
    end
end

function kaiser_bessel_ft(x, center, beta)
    if abs(x-center) <= 1
        sq_b_x = sqrt(beta^2 - (x-center)^2)
        return 2*(sinh(sq_b_x) / sq_b_x)
    else
        return 0
    end
end
