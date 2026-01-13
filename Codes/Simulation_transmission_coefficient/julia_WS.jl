using GeneralizedTransferMatrixMethod
using LinearAlgebra
using MAT
using Unzip
using Plots
import PhysicalConstants.CODATA2022: c_0, ħ, k_B, e, m_e, ε_0
eps_0 = ε_0


function eps_diagonal(no=1,ne=1)
    eps = Diagonal([no^2,ne^2,no^2])
    return eps
end

function eps_Weyl(eps_d=1,eps_a=0)
    eps = eps_d * I(3) + [[0,0,-1im * eps_a] [0,0,0] [1im * eps_a,0,0]]
    return eps
end

function reflection(wvl, eps_mat, d, α, ψ)
    layers = []
    push!(layers, eps_mat(d=d,ψ=-deg2rad(ψ)))

    @permittivity "Superstrate" λ -> eps_diagonal()
    @permittivity "Substrate" λ -> eps_diagonal(1,1)
    structure = LayeredStructure(
                superstrate = Superstrate(),
                layers = layers,
                substrate = Substrate()
            )

    R = calculate_reflection.(wvl, deg2rad(α), structure)
    R_pp, R_ss, R_ps, R_sp = unzip(R)
    return R_pp, R_ss, R_ps, R_sp
end

function reflection_analytic(wvl,eps_mat, d, α, ψ)
    kz0 = cosd(α)
    kr = sind(α)
    kx = sind(α) * cosd(ψ)
    eps = eps_mat(d=d).ϵ(wvl)
    eps_d = eps[1,1]
    eps_xz = eps[1,3]
    eps_v = eps_d - abs(eps_xz)^2 / eps_d
    # if ψ==90
    #     # yz plane
    #     kz1 = sqrt(eps_d + eps_v - 2*kx^2 + sqrt(eps_d - eps_v) * sqrt(eps_d - eps_v + 4*kx^2))
    #     R_pp = abs((kz0*eps_v - kz1 + kx * eps_xz/eps_d)/(kz0*eps_v + kz1 - kx * eps_xz/eps_d))^2
    #     kz1 = sqrt(eps_d + eps_v - 2*kx^2 - sqrt(eps_d - eps_v) * sqrt(eps_d - eps_v + 4*kx^2))
    #     R_ss = abs((kz1 - kz0)/(kz1 + kz0))^2
    # else
    # xz plane
    kz1 = sqrt(eps_v - kx^2)
    R_pp = abs((kz0*eps_v - kz1 + kx * eps_xz/eps_d)/(kz0*eps_v + kz1 - kx * eps_xz/eps_d))^2
    kz1 = sqrt(eps_d - kx^2)
    R_ss = abs((kz1 - kz0)/(kz1 + kz0))^2
    #end
    return R_pp, R_ss
end


function compute_reflection(alpha, phi)
    d = 10e-6

    @permittivity "AnisotropicMat" λ -> eps_Weyl(9+0.3im,9)

    R_pp, R_ss, R_ps, R_sp = reflection([big(1e-6)], AnisotropicMat, d, alpha, phi)
    R_p = R_sp + R_pp
    R_s = R_ps + R_ss
    return (R_p[1] + R_s[1])/2
end

#println((compute_reflection(10, 0)))