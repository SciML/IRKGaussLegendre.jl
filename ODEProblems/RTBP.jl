function RTBPEnergy(u, mu)
    """
         RTBP-Nbody problem
    """

    @inbounds begin
        x = u[1]   # x
        y = u[2]   # y
        px = u[3]  # px
        py = u[4]  # py

        umu = 1 - mu

        r1 = ((x + mu)^2 + y^2)^(1 / 2)
        r2 = ((x - umu)^2 + y^2)^(1 / 2)

        Energy = (px * px + py * py) / 2 + px * y - py * x - umu / r1 - mu / r2 -
            mu * umu / 2

        return (Energy)
    end
end

f = (du, u, p, t) -> begin
    μ = p[1]
    uμ = 1 - μ
    @inbounds begin
        # 1 = y₁
        # 2 = y₂
        # 3 = y₁'
        # 4 = y₂'
        D₁ = ((u[1] + μ)^2 + u[2]^2)^(3 / 2)
        D₂ = ((u[1] - uμ)^2 + u[2]^2)^(3 / 2)
        du[1] = u[3] + u[2]
        du[2] = u[4] - u[1]
        du[3] = u[4] - uμ * (u[1] + μ) / D₁ - μ * (u[1] - uμ) / D₂
        du[4] = -u[3] - uμ * u[2] / D₁ - μ * u[2] / D₂
    end
end
