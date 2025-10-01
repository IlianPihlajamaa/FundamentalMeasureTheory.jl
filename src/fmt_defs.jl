

function find_analytical_C_PY(k, η)
    A = -(1 - η)^-4 * (1 + 2η)^2
    B = (1 - η)^-4 * 6η * (1 + η / 2)^2
    D = -(1 - η)^-4 * 1 / 2 * η * (1 + 2η)^2
    Cₖ = @. 4π / k^6 *
            (
        24 * D - 2 * B * k^2 - (24 * D - 2 * (B + 6 * D) * k^2 + (A + B + D) * k^4) * cos(k)
        +
        k * (-24 * D + (A + 2 * B + 4 * D) * k^2) * sin(k)
    )
    return Cₖ
end


function compute_pair_direct_correlation_function(Φ, n, ω, homogeneity_substitutions)
    χ2 = Matrix{Any}(undef, length(n), length(n))
    for i in eachindex(n)
        for j in eachindex(n)
            χ2[i, j] = substitute(expand_derivatives(Differential(n[j])(Differential(n[i])(Φ))), homogeneity_substitutions)
        end
    end

    c2_matrix = Matrix{Any}(undef, length(n), length(n))
    for i in eachindex(n)
        for j in eachindex(n)
            c2_matrix[i, j] = -χ2[i, j] * ω[i](k1x, k1y, k1z) * ω[j](-k1x, -k1y, -k1z)
        end
    end
    c2 = sum(c2_matrix)

    c2_final = substitute(c2, Dict(k1x => 0, Pi => pi, k1y => 0, k1z => k1l, R => 0.5))
    f_expr = build_function(c2_final.re, k1l, ρ)
    direct_cor = eval(f_expr)

    # check symmetry
    for i in eachindex(n)
        for j in eachindex(n)
            t1 = substitute(χ2[i, j], Dict(R => 0.5, ρ => 0.5, Pi => pi))
            t2 = substitute(χ2[j, i], Dict(R => 0.5, ρ => 0.5, Pi => pi))
            if !(t1 ≈ t2)
                println("Not Symmetric")
            end
        end
    end

    return direct_cor
end


function compute_triplet_direct_correlation_function(Φ, n, ω, homogeneity_substitutions)
    # println("Constructing triplet direct correlation function")

    begin
        χ1 = Array{Any,1}(undef, length(n))
        for i in eachindex(n)
            χ1[i] = expand_derivatives(Differential(n[i])(Φ))
        end
        χ2 = Array{Any,2}(undef, length(n), length(n))
        for i in eachindex(n)
            for j in eachindex(n)
                χ2[i, j] = expand_derivatives(Differential(n[j])(χ1[i]))
            end
        end
        χ3 = Array{Any,3}(undef, length(n), length(n), length(n))
        for i in eachindex(n)
            for j in eachindex(n)
                for k in eachindex(n)
                    χ3[i, j, k] = substitute(expand_derivatives(Differential(n[k])(χ2[i, j])), homogeneity_substitutions)
                end
            end
        end
    end

    c3_matrix = Array{Any,3}(undef, length(n), length(n), length(n))
    for i in eachindex(n)
        for j in eachindex(n)
            for k in eachindex(n)
                c3_matrix[i, j, k] = -χ3[i, j, k] * ω[i](0, 0, k1l) * ω[j](k2l * sin(θ), 0, k2l * cos(θ)) * ω[k](-0 - k2l * sin(θ), -0 - 0, -k1l - k2l * cos(θ))
            end
        end
    end
    c3_final = sum(c3_matrix)
    f_expr3 = build_function(c3_final.re, k1l, k2l, θ, ρ)
    direct_cor3 = eval(f_expr3)

    return direct_cor3
end


function compute_quadruplet_direct_correlation_function(Φ, n, ω, homogeneity_substitutions)
    println("Constructing quadruplet direct correlation function")
    χ4 = Array{Any,4}(undef, length(n), length(n), length(n), length(n))


    @time begin
        χ1 = Array{Any,1}(undef, length(n))
        for i in eachindex(n)
            χ1[i] = expand_derivatives(Differential(n[i])(Φ))
        end
        χ2 = Array{Any,2}(undef, length(n), length(n))
        for i in eachindex(n)
            for j in eachindex(n)
                χ2[i, j] = expand_derivatives(Differential(n[j])(χ1[i]))
            end
        end
        χ3 = Array{Any,3}(undef, length(n), length(n), length(n))
        for i in eachindex(n)
            for j in eachindex(n)
                for k in eachindex(n)
                    χ3[i, j, k] = expand_derivatives(Differential(n[k])(χ2[i, j]))
                end
            end
        end
        for i in eachindex(n)
            for j in eachindex(n)
                for k in eachindex(n)
                    for l in eachindex(n)
                        χ4_inhom = expand_derivatives(Differential(n[l])(χ3[i, j, k]))
                        χ4[i, j, k, l] = substitute(χ4_inhom, homogeneity_substitutions)
                    end
                end
            end
        end
    end

    c4_final = Num(0)
    @time for i in eachindex(n)
        w1 = ω[i](0, 0, k1l)
        for j in eachindex(n)
            w2 = ω[j](k2l * sin(θ12), 0, k2l * cos(θ12))
            for k in eachindex(n)
                w3 = ω[k](k3l * sin(θ13) * cos(ϕ23), k3l * sin(θ13) * sin(ϕ23), k3l * cos(θ13))
                for l in eachindex(n)
                    w4 = ω[l](-(0 + k2l * sin(θ12) + k3l * sin(θ13) * cos(ϕ23)), -(0 + 0 + k3l * sin(θ13) * sin(ϕ23)), -(k1l + k2l * cos(θ12) + k3l * cos(θ13)))
                    c4_final += -χ4[i, j, k, l] * w1 * w2 * w3 * w4
                end
            end
        end
    end
    f_expr4 = build_function(c4_final.re, k1l, k2l, k3l, θ12, θ13, ϕ23, ρ)
    direct_cor4 = eval(f_expr4)


    return direct_cor4
end

function get_S2(direct_cor2, k, ρ)
    c2 = direct_cor2(k, ρ)
    S2 = 1 + ρ * c2 / (1 - ρ * c2)
    return S2
end

function get_S3(direct_cor2, direct_cor3, k1, k2, θ, ρ)
    k3 = sqrt(k1^2 + k2^2 + 2k1 * k2 * cos(θ))

    S2_k1 = get_S2(direct_cor2, k1, ρ)
    S2_k2 = get_S2(direct_cor2, k2, ρ)
    S2_k3 = get_S2(direct_cor2, k3, ρ)

    S3 = S2_k1 * S2_k2 * S2_k3 * (1 + ρ^2 * direct_cor3(k1, k2, θ, ρ))
    return S3
end

function get_c3(direct_cor3, k1vec, k2vec, ρ)
    k1 = norm(k1vec)
    k2 = norm(k2vec)
    θ = acos(dot(k1vec, k2vec) / (k1 * k2))
    return direct_cor3(k1, k2, θ, ρ)
end

function get_S3(direct_cor2, direct_cor3, k1vec, k2vec, ρ)
    k1 = norm(k1vec)
    k2 = norm(k2vec)
    θ = acos(dot(k1vec, k2vec) / (k1 * k2))
    return get_S3(direct_cor2, direct_cor3, k1, k2, θ, ρ)
end

function get_S4(direct_cor2, direct_cor3, direct_cor4, k1, k2, k3, θ12, θ13, ϕ23, ρ)
    K1 = [0, 0, k1]
    K2 = [k2 * sin(θ12), 0, k2 * cos(θ12)]
    K3 = [k3 * sin(θ13) * cos(ϕ23), k3 * sin(θ13) * sin(ϕ23), k3 * cos(θ13)]
    K4 = -K1 - K2 - K3
    k4 = norm(K4)

    S2_k1 = get_S2(direct_cor2, k1, ρ)
    S2_k2 = get_S2(direct_cor2, k2, ρ)
    S2_k3 = get_S2(direct_cor2, k3, ρ)
    S2_k4 = get_S2(direct_cor2, k4, ρ)

    S3_2_4 = get_S3(direct_cor2, direct_cor3, K2, K4, ρ)
    S3_2_3 = get_S3(direct_cor2, direct_cor3, K2, K3, ρ)
    S3_2_1 = get_S3(direct_cor2, direct_cor3, K2, K1, ρ)

    term1 = S2_k1 * S2_k2 * S2_k3 * S2_k4 * (-2 + ρ^3 * direct_cor4(k1, k2, k3, θ12, θ13, ϕ23, ρ))
    term2 = S2_k1 * S2_k3 * S3_2_4 * (1 + ρ^2 * get_c3(direct_cor3, K1, K3, ρ))
    term3 = S2_k1 * S2_k4 * S3_2_3 * (1 + ρ^2 * get_c3(direct_cor3, K1, K4, ρ))
    term4 = S2_k3 * S2_k4 * S3_2_1 * (1 + ρ^2 * get_c3(direct_cor3, K3, K4, ρ))
    S4 = term1 + term2 + term3 + term4
    return S4
end


j0(x) = sin(x) / x
j1(x) = sin(x) / x^2 - cos(x) / x
j2(x) = 3 * sin(x) / x^3 - 3 * cos(x) / x^2 - sin(x) / x

K(kx, ky, kz) = sqrt(kx^2 + ky^2 + kz^2)
ω0(kx, ky, kz) = j0(K(kx, ky, kz) * R)
ω1(kx, ky, kz) = R * j0(K(kx, ky, kz) * R)
ω2(kx, ky, kz) = 4 * Pi * R^2 * j0(K(kx, ky, kz) * R)
ω3(kx, ky, kz) = 4Pi * R^3 / (R * K(kx, ky, kz)) * j1(K(kx, ky, kz) * R)

ω2x(kx, ky, kz) = -im * kx * ω3(kx, ky, kz)
ω2y(kx, ky, kz) = -im * ky * ω3(kx, ky, kz)
ω2z(kx, ky, kz) = -im * kz * ω3(kx, ky, kz)

ω1x(kx, ky, kz) = ω2x(kx, ky, kz) / (4 * Pi * R)
ω1y(kx, ky, kz) = ω2y(kx, ky, kz) / (4 * Pi * R)
ω1z(kx, ky, kz) = ω2z(kx, ky, kz) / (4 * Pi * R)

ω0_KR(kx, ky, kz) = cos(K(kx, ky, kz) * R) + K(kx, ky, kz) * R * sin(K(kx, ky, kz) * R) / 2
ω1_KR(kx, ky, kz) = (sin(K(kx, ky, kz) * R) + K(kx, ky, kz) * R * cos(K(kx, ky, kz) * R)) / (2 * K(kx, ky, kz))
ω2_KR(kx, ky, kz) = 4 * Pi * R^2 * j0(K(kx, ky, kz) * R) + 0 * im
ω3_KR(kx, ky, kz) = 4Pi * R^3 / (R * K(kx, ky, kz)) * j1(K(kx, ky, kz) * R)

#tarazona's functions 
traceless = true
ω2xx(kx, ky, kz) = -4 * Pi * R^2 * (kx * kx / K(kx, ky, kz)^2 - traceless / 3) * j2(K(kx, ky, kz) * R)
ω2xy(kx, ky, kz) = -4 * Pi * R^2 * kx * ky / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2xz(kx, ky, kz) = -4 * Pi * R^2 * kx * kz / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2yx(kx, ky, kz) = -4 * Pi * R^2 * ky * kx / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2yy(kx, ky, kz) = -4 * Pi * R^2 * (ky * ky / K(kx, ky, kz)^2 - traceless / 3) * j2(K(kx, ky, kz) * R)
ω2yz(kx, ky, kz) = -4 * Pi * R^2 * ky * kz / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2zx(kx, ky, kz) = -4 * Pi * R^2 * kz * kx / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2zy(kx, ky, kz) = -4 * Pi * R^2 * kz * ky / K(kx, ky, kz)^2 * j2(K(kx, ky, kz) * R)
ω2zz(kx, ky, kz) = -4 * Pi * R^2 * (kz * kz / K(kx, ky, kz)^2 - traceless / 3) * j2(K(kx, ky, kz) * R)
@variables n0 n1 n2 n3 n1x n2x n1y n2y n1z n2z
@variables n2xx n2xy n2xz n2yx n2yy n2yz n2zx n2zy n2zz
@variables k1x k2x k1y k2y k1z k2z k3x k3y k3z
@variables k1l k2l k3l θ θ12 θ13 ϕ23
@variables ρ

Pi::Float64 = Float64(π)
R::Float64 = 0.5

n = [n0, n1, n2, n3, n1x, n1y, n1z, n2x, n2y, n2z]
ω = [ω0, ω1, ω2, ω3, ω1x, ω1y, ω1z, ω2x, ω2y, ω2z]

n_KR = [n0, n1, n2, n3]
ω_KR = [ω0_KR, ω1_KR, ω2_KR, ω3_KR]

n_T = [n0, n1, n2, n3, n1x, n1y, n1z, n2x, n2y, n2z, n2xx, n2xy, n2xz, n2yx, n2yy, n2yz, n2zx, n2zy, n2zz]
ω_T = [ω0, ω1, ω2, ω3, ω1x, ω1y, ω1z, ω2x, ω2y, ω2z, ω2xx, ω2xy, ω2xz, ω2yx, ω2yy, ω2yz, ω2zx, ω2zy, ω2zz]

n1dotn2 = n1x * n2x + n1y * n2y + n1z * n2z
n2dotn2 = n2x^2 + n2y^2 + n2z^2
n_tensor = [n2xx n2xy n2xz; n2yx n2yy n2yz; n2zx n2zy n2zz]
n_vector = [n2x, n2y, n2z]

vTv = (reshape(n_vector, (1, 3))*n_tensor*n_vector)[1]
n_squared = n_tensor * n_tensor
n_cubed = n_tensor * n_tensor * n_tensor


Φ_RF = -n0 * log(1 - n3) +
       (n1 * n2 - n1dotn2) / (1 - n3) +
       (n2^3 - 3 * n2 * n2dotn2) / (24 * Pi * (1 - n3)^2)

Φ_WB = -n0 * log(1 - n3) +
       (n1 * n2 - n1dotn2) / (1 - n3) +
       (n2^3 - 3 * n2 * n2dotn2) / (36 * Pi * n3^2 * (1 - n3)^2) * (n3 + (1 - n3)^2 * log(1 - n3))

ϕ2WB2(n3) = (6n3 - 3n3^2 + 6 * (1 - n3) * log(1 - n3)) / n3^3
ϕ3WB2(n3) = (6n3 - 9n3^2 + 6n3^3 + 6 * (1 - n3)^2 * log(1 - n3)) / (4n3^3)

Φ_WB2 = -n0 * log(1 - n3) +
        (1 + n3^2 * ϕ2WB2(n3) / 9) * (n1 * n2 - n1dotn2) / (1 - n3) +
        (1 - 4 * n3 * ϕ3WB2(n3) / 9) * (n2^3 - 3 * n2 * n2dotn2) / (24 * Pi * (1 - n3)^2)

Φ_M = -n0 * log(1 - n3) +
      (n1 * n2 - n1dotn2) / (1 - n3) +
      (n2^3 - 3 * n2 * n2dotn2) / (108 * Pi * n3^2 * (1 - n3)^2) * (8 * (1 - n3)^2 * log(1 - n3) + 8n3 - 15n3^2 / 2 + 2n3^3)

Φ_KR = -n0 * log(1 - n3) + (n1 * n2) / (1 - n3) + (n2^3) / (24 * Pi * (1 - n3)^2)

φ3_tensor = (n2^3 - 3 * n2 * n2dotn2 + 9 / 2 * (vTv - tr(n_cubed))) / (24 * Pi * (1 - n3)^2) # roth
# φ3_tensor =  3*(-n2 * n2dotn2 + vTv * n2*tr(n_squared) - tr(n_cubed)) / (16*Pi*(1-n3)^2) # tarazona wrong: should use non-traceless tensor

Φ_T = -n0 * log(1 - n3) + (n1 * n2 - n1dotn2) / (1 - n3) + φ3_tensor

Φ_WB2t = -n0 * log(1 - n3) + (1 + n3^2 * ϕ2WB2(n3) / 9) * (n1 * n2 - n1dotn2) / (1 - n3) +
         (1 - 4 * n3 * ϕ3WB2(n3) / 9) * φ3_tensor

#Lutsko’s Phi3
Φ3_L_func(A, B) = ((8A + 2B) / 9 * n2^3 - 2 * A * n2 * n2dotn2 + 3 * A * vTv - (A + B) * n2 * tr(n_squared) + (2B - A) * tr(n_cubed)) / (24 * Pi * (1 - n3)^2)

# Tarazona is A=3/2, B=-3/2: fill in random numbers to check
@assert abs(substitute(φ3_tensor, Dict(n3 => 0.1523, n2 => 0.53, n2x => 0.1, n2y => 0.2, n2z => 0.3, n2xx => 0.01, n2xy => 0.02, n2xz => 0.03, n2yx => 0.02, n2yy => 0.04, n2yz => 0.05, n2zx => 0.03, n2zy => 0.05, n2zz => 0.06)) -
            substitute(Φ3_L_func(3 / 2, -3 / 2), Dict(n2 => 0.53, n2x => 0.1, n2y => 0.2, n2z => 0.3, n2xx => 0.01, n2xy => 0.02, n2xz => 0.03, n2yx => 0.02, n2yy => 0.04, n2yz => 0.05, n2zx => 0.03, n2zy => 0.05, n2zz => 0.06, n3 => 0.1523,))) < 1e-17

# Lutsko is A=1, B=0
Φ_L = -n0 * log(1 - n3) + (n1 * n2 - n1dotn2) / (1 - n3) + Φ3_L_func(1, 0)
# Gül is A = 1.3, B = -1
Φ_G = -n0 * log(1 - n3) + (n1 * n2 - n1dotn2) / (1 - n3) + Φ3_L_func(1.3, -1)

homogeneity_substitutions = Dict(
    Pi => pi,
    n0 => ρ,
    n1 => ρ * R,
    n2 => 4 * Pi * ρ * R^2,
    n3 => 4 * Pi * R^3 * ρ / 3,
    n1x => 0, n1y => 0, n1z => 0,
    n2x => 0, n2y => 0, n2z => 0,
    n2xx => 0, n2xy => 0, n2xz => 0, n2yx => 0, n2yy => 0, n2yz => 0, n2zx => 0, n2zy => 0, n2zz => 0
)

function c_from_s(Sk, rho)
    ρCk = @. (Sk - 1) / Sk
    return ρCk / rho
end

function δc_from_s(Sk, δSk, rho)
    ρCk = @. 1 - 1 / Sk
    δρCk = @. δSk / (Sk^2)
    return δρCk / rho
end