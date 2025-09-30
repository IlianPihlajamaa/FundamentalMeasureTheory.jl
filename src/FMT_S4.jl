println("Starting the script")
using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie
using LaTeXStrings
using HDF5
using Dierckx

mytheme = Theme( Axis = (xtickalign = 1,
xminortickalign = 1,
xticksmirrored = true,
ytickalign = 1,
yminortickalign = 1,
yticksmirrored = true, 
xgridvisible = false, 
ygridvisible = false,))
set_theme!(merge(mytheme, theme_latexfonts()))
CairoMakie.theme(:palette).color[] = ColorSchemes.:ColorSchemes.:seaborn_muted6



function find_analytical_C_k(k, η)
    A = -(1 - η)^-4 *(1 + 2η)^2
    B = (1 - η)^-4*  6η*(1 + η/2)^2
    D = -(1 - η)^-4 * 1/2 * η*(1 + 2η)^2
    Cₖ = @. 4π/k^6 * 
    (
        24*D - 2*B * k^2 - (24*D - 2 * (B + 6*D) * k^2 + (A + B + D) * k^4) * cos(k)
     + k * (-24*D + (A + 2*B + 4*D) * k^2) * sin(k)
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
            c2_matrix[i,j] = - χ2[i, j] * ω[i](k1x, k1y, k1z) * ω[j](-k1x, -k1y, -k1z)
        end
    end
    c2 = sum(c2_matrix)

    c2_final = substitute(c2, Dict(k1x=>0, Pi=>pi, k1y=>0, k1z=>k1l, R=>0.5))
    f_expr = build_function(c2_final.re, k1l, ρ)
    direct_cor = eval(f_expr)

    # check symmetry
    for i in eachindex(n)
        for j in eachindex(n)
            t1 = substitute(χ2[i,j], Dict(R=>0.5, ρ=>0.5, Pi=>pi))
            t2 = substitute(χ2[j,i], Dict(R=>0.5, ρ=>0.5, Pi=>pi))
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
                    χ3[i, j, k] = substitute(expand_derivatives(Differential(n[k])(χ2[i,j])), homogeneity_substitutions)
                end
            end
        end
    end

    c3_matrix = Array{Any,3}(undef, length(n), length(n), length(n))
    for i in eachindex(n)
        for j in eachindex(n)
            for k in eachindex(n)
                c3_matrix[i,j,k] = - χ3[i, j, k] * ω[i](0, 0, k1l) * ω[j](k2l*sin(θ), 0, k2l*cos(θ)) * ω[k](-0-k2l*sin(θ), -0-0, -k1l-k2l*cos(θ))
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
                    χ3[i, j, k] = expand_derivatives(Differential(n[k])(χ2[i,j]))
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
            w2 = ω[j](k2l*sin(θ12), 0, k2l*cos(θ12))
            for k in eachindex(n)
                w3 = ω[k](k3l*sin(θ13)*cos(ϕ23), k3l*sin(θ13)*sin(ϕ23), k3l*cos(θ13))
                for l in eachindex(n)
                    w4 = ω[l](-(0+k2l*sin(θ12)+k3l*sin(θ13)*cos(ϕ23)), -(0+0+k3l*sin(θ13)*sin(ϕ23)), -(k1l+k2l*cos(θ12)+k3l*cos(θ13)))
                    c4_final += - χ4[i, j, k, l] * w1 * w2 * w3 * w4
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
    k3 = sqrt(k1^2 + k2^2 + 2k1*k2*cos(θ))

    S2_k1 = get_S2(direct_cor2, k1, ρ)
    S2_k2 = get_S2(direct_cor2, k2, ρ)
    S2_k3 = get_S2(direct_cor2, k3, ρ)
    
    S3 = S2_k1*S2_k2*S2_k3 * (1 + ρ^2 * direct_cor3(k1, k2, θ, ρ))
    return S3
end

function get_c3(direct_cor3, k1vec, k2vec, ρ)
    k1 = norm(k1vec)
    k2 = norm(k2vec)
    θ = acos(dot(k1vec, k2vec)/(k1*k2))
    return direct_cor3(k1, k2, θ, ρ)
end

function get_S3(direct_cor2, direct_cor3, k1vec, k2vec, ρ)
    k1 = norm(k1vec)
    k2 = norm(k2vec)
    θ = acos(dot(k1vec, k2vec)/(k1*k2))
    return get_S3(direct_cor2, direct_cor3, k1, k2, θ, ρ)
end

function get_S4(direct_cor2, direct_cor3, direct_cor4, k1, k2, k3, θ12, θ13, ϕ23, ρ)
    K1 = [0, 0, k1]
    K2 = [k2*sin(θ12), 0, k2*cos(θ12)]
    K3 = [k3*sin(θ13)*cos(ϕ23), k3*sin(θ13)*sin(ϕ23), k3*cos(θ13)]
    K4 = -K1 - K2 - K3
    k4 = norm(K4)

    S2_k1 = get_S2(direct_cor2, k1, ρ)
    S2_k2 = get_S2(direct_cor2, k2, ρ)
    S2_k3 = get_S2(direct_cor2, k3, ρ)
    S2_k4 = get_S2(direct_cor2, k4, ρ)

    S3_2_4 = get_S3(direct_cor2, direct_cor3, K2, K4, ρ)
    S3_2_3 = get_S3(direct_cor2, direct_cor3, K2, K3, ρ)
    S3_2_1 = get_S3(direct_cor2, direct_cor3, K2, K1, ρ)

    term1 = S2_k1*S2_k2*S2_k3*S2_k4*(-2 + ρ^3 * direct_cor4(k1, k2, k3, θ12, θ13, ϕ23, ρ))
    term2 = S2_k1*S2_k3*S3_2_4*(1 + ρ^2 * get_c3(direct_cor3, K1, K3, ρ))
    term3 = S2_k1*S2_k4*S3_2_3*(1 + ρ^2 * get_c3(direct_cor3, K1, K4, ρ))
    term4 = S2_k3*S2_k4*S3_2_1*(1 + ρ^2 * get_c3(direct_cor3, K3, K4, ρ))
    S4 = term1 + term2 + term3 + term4
    return S4
end


j0(x) = sin(x)/x
j1(x) = sin(x)/x^2 - cos(x)/x
j2(x) = 3*sin(x)/x^3 - 3*cos(x)/x^2 - sin(x)/x

K(kx, ky, kz) = sqrt(kx^2 + ky^2 + kz^2)
ω0(kx, ky, kz) = j0(K(kx, ky, kz)*R)
ω1(kx, ky, kz) = R*j0(K(kx, ky, kz)*R)
ω2(kx, ky, kz) = 4*Pi*R^2*j0(K(kx, ky, kz)*R)
ω3(kx, ky, kz) = 4Pi*R^3/(R*K(kx, ky, kz))*j1(K(kx, ky, kz)*R)

ω2x(kx, ky, kz) = -im * kx * ω3(kx, ky, kz)
ω2y(kx, ky, kz) = -im * ky * ω3(kx, ky, kz)
ω2z(kx, ky, kz) = -im * kz * ω3(kx, ky, kz)

ω1x(kx, ky, kz) = ω2x(kx, ky, kz) / (4*Pi*R)
ω1y(kx, ky, kz) = ω2y(kx, ky, kz) / (4*Pi*R)
ω1z(kx, ky, kz) = ω2z(kx, ky, kz) / (4*Pi*R)

ω0_KR(kx, ky, kz) = cos(K(kx, ky, kz)*R) + K(kx, ky, kz)*R*sin(K(kx, ky, kz)*R)/2
ω1_KR(kx, ky, kz) = (sin(K(kx, ky, kz)*R) + K(kx, ky, kz)*R*cos(K(kx, ky, kz)*R))/(2*K(kx, ky, kz))
ω2_KR(kx, ky, kz) = 4*Pi*R^2*j0(K(kx, ky, kz)*R) + 0*im
ω3_KR(kx, ky, kz) = 4Pi*R^3/(R*K(kx, ky, kz))*j1(K(kx, ky, kz)*R)

#tarazona's functions 
traceless = true
ω2xx(kx, ky, kz) = -4*Pi*R^2 * (kx*kx/K(kx, ky, kz)^2 - traceless/3)  * j2(K(kx, ky, kz)*R)
ω2xy(kx, ky, kz) = -4*Pi*R^2 *  kx*ky/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R) 
ω2xz(kx, ky, kz) = -4*Pi*R^2 *  kx*kz/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R) 
ω2yx(kx, ky, kz) = -4*Pi*R^2 *  ky*kx/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R) 
ω2yy(kx, ky, kz) = -4*Pi*R^2 * (ky*ky/K(kx, ky, kz)^2 - traceless/3)  * j2(K(kx, ky, kz)*R)
ω2yz(kx, ky, kz) = -4*Pi*R^2 *  ky*kz/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R)  
ω2zx(kx, ky, kz) = -4*Pi*R^2 *  kz*kx/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R) 
ω2zy(kx, ky, kz) = -4*Pi*R^2 *  kz*ky/K(kx, ky, kz)^2                 * j2(K(kx, ky, kz)*R) 
ω2zz(kx, ky, kz) = -4*Pi*R^2 * (kz*kz/K(kx, ky, kz)^2 - traceless/3)  * j2(K(kx, ky, kz)*R)
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

n1dotn2 = n1x*n2x + n1y*n2y + n1z*n2z
n2dotn2 = n2x^2 + n2y^2 + n2z^2
n_tensor = [n2xx n2xy n2xz; n2yx n2yy n2yz; n2zx n2zy n2zz]
n_vector = [n2x, n2y, n2z]

vTv = (reshape(n_vector,(1, 3)) * n_tensor * n_vector)[1]
n_squared = n_tensor * n_tensor
n_cubed = n_tensor * n_tensor * n_tensor


Φ_RF = -n0*log(1-n3) +
        (n1*n2-n1dotn2)/(1-n3) +
        (n2^3 - 3*n2 * n2dotn2)/(24*Pi*(1-n3)^2)

Φ_WB = -n0*log(1-n3) +
        (n1*n2-n1dotn2)/(1-n3) +
        (n2^3 - 3*n2 * n2dotn2)/(36*Pi*n3^2*(1-n3)^2) * (n3 + (1-n3)^2*log(1-n3))

ϕ2WB2(n3) = (6n3 - 3n3^2 + 6*(1-n3)*log(1-n3)) / n3^3
ϕ3WB2(n3) = (6n3 - 9n3^2 + 6n3^3 + 6*(1-n3)^2*log(1-n3)) / (4n3^3)

Φ_WB2 = -n0*log(1-n3) +
        (1+n3^2*ϕ2WB2(n3)/9)*(n1*n2-n1dotn2)/(1-n3) +
        (1-4*n3*ϕ3WB2(n3)/9) * (n2^3 - 3*n2 * n2dotn2)/(24*Pi*(1-n3)^2)  

Φ_M = -n0*log(1-n3) +
        (n1*n2-n1dotn2)/(1-n3) +
        (n2^3 - 3*n2 * n2dotn2)/(108*Pi*n3^2*(1-n3)^2)*(8*(1-n3)^2*log(1-n3)+8n3-15n3^2/2+2n3^3)

Φ_KR = -n0*log(1-n3) + (n1*n2)/(1-n3) + (n2^3)/(24*Pi*(1-n3)^2)

φ3_tensor = (n2^3 - 3*n2 * n2dotn2 + 9/2*(vTv - tr(n_cubed)))/(24*Pi*(1-n3)^2) # roth
# φ3_tensor =  3*(-n2 * n2dotn2 + vTv * n2*tr(n_squared) - tr(n_cubed)) / (16*Pi*(1-n3)^2) # tarazona wrong: should use non-traceless tensor

Φ_T = -n0*log(1-n3) + (n1*n2-n1dotn2)/(1-n3) + φ3_tensor

Φ_WB2t = - n0*log(1-n3) +  (1+n3^2*ϕ2WB2(n3)/9)*(n1*n2-n1dotn2)/(1-n3) +
        (1-4*n3*ϕ3WB2(n3)/9) * φ3_tensor

#Lutsko’s Phi3
Φ3_L_func(A, B) = ((8A + 2B)/9*n2^3 - 2*A*n2*n2dotn2 +3*A*vTv - (A+B)*n2*tr(n_squared) + (2B-A)*tr(n_cubed))/(24*Pi*(1-n3)^2)

# Tarazona is A=3/2, B=-3/2: fill in random numbers to check
@assert abs(substitute(φ3_tensor, Dict(n3 =>0.1523, n2=>0.53, n2x=>0.1, n2y=>0.2, n2z=>0.3, n2xx=>0.01, n2xy=>0.02, n2xz=>0.03, n2yx=>0.02, n2yy=>0.04, n2yz=>0.05, n2zx=>0.03, n2zy=>0.05, n2zz=>0.06)) - 
substitute(Φ3_L_func(3/2, -3/2), Dict(n2=>0.53, n2x=>0.1, n2y=>0.2, n2z=>0.3, n2xx=>0.01, n2xy=>0.02, n2xz=>0.03, n2yx=>0.02, n2yy=>0.04, n2yz=>0.05, n2zx=>0.03, n2zy=>0.05, n2zz=>0.06, n3 =>0.1523, ))) < 1e-17

# Lutsko is A=1, B=0
Φ_L = - n0*log(1-n3) +  (n1*n2-n1dotn2)/(1-n3) + Φ3_L_func(1, 0)
# Gül is A = 1.3, B = -1
Φ_G = - n0*log(1-n3) +  (n1*n2-n1dotn2)/(1-n3) + Φ3_L_func(1.3, -1)

homogeneity_substitutions = Dict(
    Pi=>pi,
    n0=>ρ, 
    n1=>ρ*R, 
    n2=>4*Pi*ρ*R^2, 
    n3=>4*Pi*R^3*ρ/3, 
    n1x=>0, n1y=>0, n1z=>0,
    n2x=>0, n2y=>0, n2z=>0,
    n2xx=>0, n2xy=>0, n2xz=>0, n2yx=>0, n2yy=>0, n2yz=>0, n2zx=>0, n2zy=>0, n2zz=>0
    )

direct_cor2_RF = compute_pair_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor2_WB = compute_pair_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor2_WB2 = compute_pair_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor2_M = compute_pair_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor2_KR = compute_pair_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor2_T = compute_pair_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor2_WB2t = compute_pair_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor2_L = compute_pair_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor2_G = compute_pair_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@show direct_cor2_RF(0.1, 0.94)*0.94
@show direct_cor2_WB(0.1, 0.94)*0.94
@show direct_cor2_WB2(0.1, 0.94)*0.94
@show direct_cor2_M(0.1, 0.94)*0.94
@show direct_cor2_KR(0.1, 0.94)*0.94
@show direct_cor2_T(0.1, 0.94)*0.94
@show direct_cor2_WB2t(0.1, 0.94)*0.94
@show direct_cor2_L(0.1, 0.94)*0.94
@show direct_cor2_G(0.1, 0.94)*0.94
@assert abs(direct_cor2_RF(0.1, 0.94) - find_analytical_C_k(0.1, 0.94/6*π)) < 1e-5
@assert abs(direct_cor2_KR(0.1, 0.94) - find_analytical_C_k(0.1, 0.94/6*π)) < 1e-5
# @assert abs(direct_cor2_T(0.1, 0.94) - find_analytical_C_k(0.1, 0.94/6*π)) < 1e-5


direct_cor3_RF = compute_triplet_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor3_WB = compute_triplet_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor3_WB2 = compute_triplet_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor3_M = compute_triplet_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor3_KR = compute_triplet_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor3_T = compute_triplet_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor3_WB2t = compute_triplet_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor3_L = compute_triplet_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor3_G = compute_triplet_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@show direct_cor3_RF(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_WB(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_WB2(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_M(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_KR(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_T(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_WB2t(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_L(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
@show direct_cor3_G(0.1, 0.1, acos(-0.175), 0.94)*0.94^2

direct_cor4_RF = compute_quadruplet_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor4_WB = compute_quadruplet_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor4_WB2 = compute_quadruplet_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor4_M = compute_quadruplet_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor4_KR = compute_quadruplet_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor4_T = compute_quadruplet_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor4_WB2t = compute_quadruplet_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor4_L = compute_quadruplet_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor4_G = compute_quadruplet_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@show direct_cor4_RF(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_WB(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_WB2(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_M(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_KR(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_T(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_WB2t(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_L(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3
@show direct_cor4_G(2.2, 2.2, 2.2, π/4, acos(0.025), 0.5, 0.94)*0.94^3

@assert direct_cor2_RF(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB(0.1, 0.94) * 0.94 ≈ -51.88510371086945
@assert direct_cor2_WB2(0.1, 0.94) * 0.94 ≈ -51.885375608479336
@assert direct_cor2_M(0.1, 0.94) * 0.94 ≈ -51.61924863988694
@assert direct_cor2_KR(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_T(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB2t(0.1, 0.94) * 0.94 ≈ -51.885375608479336

@assert direct_cor3_RF(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -229.81290257878078
@assert direct_cor3_WB(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -193.7786362758227
@assert direct_cor3_WB2(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -193.78046452495985
@assert direct_cor3_M(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -190.14206661188794
@assert direct_cor3_KR(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -229.81290257878078
@assert direct_cor3_T(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -229.81290257878078
@assert direct_cor3_WB2t(0.1, 0.1, acos(-0.175), 0.94) * 0.94 ^ 2 ≈ -193.78046452495985

@assert direct_cor4_RF(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -253.09097694589278
@assert direct_cor4_WB(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -204.86292780419802
@assert direct_cor4_WB2(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -207.35891643000647
@assert direct_cor4_M(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -196.93451533988554
@assert direct_cor4_KR(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -253.0909769458927
# @assert direct_cor4_T(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -254.58630508349358
# @assert direct_cor4_WB2t(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94 ^ 3 ≈ -208.41390338083065





function c_from_s(Sk, rho)
    ρCk = @. (Sk-1)/Sk
    return ρCk/rho
end




begin 
    fig = Figure(size = (800, 1200))
    ax_RF       = Axis(fig[1, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}", title=L"S^{(4)}(k, k, k, \theta_{12}, \theta_{13}, \phi_{23})")
    ax2_RF      = Axis(fig[1, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}", title=L"S^{(4)}(k, k, 1.0, \theta_{12}, \theta_{13}, \phi_{23})")
    ax_WB       = Axis(fig[2, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_WB      = Axis(fig[2, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax_WB2      = Axis(fig[3, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_WB2     = Axis(fig[3, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax_M        = Axis(fig[4, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_M       = Axis(fig[4, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax_T        = Axis(fig[5, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_T       = Axis(fig[5, 3], ylabel = L"\cos
    \theta_{13}", xlabel = L"\phi_{23}")
    ax_WB2t    = Axis(fig[6, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_WB2t   = Axis(fig[6, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax_conv     = Axis(fig[7, 1], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")
    ax2_conv    = Axis(fig[7, 3], ylabel = L"\cos \theta_{13}", xlabel = L"\phi_{23}")


    ρval = 0.94
    kval = 7.2
    θ12 = 2pi/3  #acos(-0.5)
    cosθ13_arr = range(-0.99, 0.99, length = 50)
    ϕ23_arr = range(0.01, π-0.01, length = 40)
    S4_RF = @time [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB = @time  [get_S4(direct_cor2_WB, direct_cor3_WB, direct_cor4_WB, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB2 =  @time [get_S4(direct_cor2_WB2, direct_cor3_WB2, direct_cor4_WB2, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_M =  @time [get_S4(direct_cor2_M, direct_cor3_M, direct_cor4_M, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_T =  @time [get_S4(direct_cor2_T, direct_cor3_T, direct_cor4_T, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB2t = @time [get_S4(direct_cor2_WB2t, direct_cor3_WB2t, direct_cor4_WB2t, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    
    foo(x...) = 0.0 
    S4_conv = @time [get_S4(direct_cor2_WB2, foo, foo, kval, kval, kval, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]

    heatmap!(ax_RF,  ϕ23_arr, cosθ13_arr, S4_RF', colorrange=(-maximum(abs.(S4_RF)),maximum(abs.(S4_RF))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_WB,  ϕ23_arr, cosθ13_arr, S4_WB', colorrange=(-maximum(abs.(S4_WB)),maximum(abs.(S4_WB))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_WB2,  ϕ23_arr, cosθ13_arr, S4_WB2', colorrange=(-maximum(abs.(S4_WB2)),maximum(abs.(S4_WB2))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_M,  ϕ23_arr, cosθ13_arr, S4_M', colorrange=(-maximum(abs.(S4_M)),maximum(abs.(S4_M))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_T,  ϕ23_arr, cosθ13_arr, S4_T',  colorrange=(-maximum(abs.(S4_T)),maximum(abs.(S4_T))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_WB2t,  ϕ23_arr, cosθ13_arr, S4_WB2t', colorrange=(-maximum(abs.(S4_WB2t)),maximum(abs.(S4_WB2t))), colormap=ColorSchemes.:bwr)
    heatmap!(ax_conv,  ϕ23_arr, cosθ13_arr, S4_conv',  colorrange=(-maximum(abs.(S4_conv)),maximum(abs.(S4_conv))), colormap=ColorSchemes.:bwr)

    Colorbar(fig[1, 2], limits=(-maximum(abs.(S4_RF)),maximum(abs.(S4_RF))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[2, 2], limits=(-maximum(abs.(S4_WB)),maximum(abs.(S4_WB))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[3, 2], limits=(-maximum(abs.(S4_WB2)),maximum(abs.(S4_WB2))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[4, 2], limits=(-maximum(abs.(S4_M)),maximum(abs.(S4_M))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[5, 2], limits=(-maximum(abs.(S4_T)),maximum(abs.(S4_T))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[6, 2], limits=(-maximum(abs.(S4_WB2t)),maximum(abs.(S4_WB2t))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[7, 2], limits=(-maximum(abs.(S4_conv)),maximum(abs.(S4_conv))), colormap=ColorSchemes.:bwr)



    S4_RF =  @time  [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB =  @time [get_S4(direct_cor2_WB, direct_cor3_WB, direct_cor4_WB, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB2 = @time  [get_S4(direct_cor2_WB2, direct_cor3_WB2, direct_cor4_WB2, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_M = @time  [get_S4(direct_cor2_M, direct_cor3_M, direct_cor4_M, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_T = @time  [get_S4(direct_cor2_T, direct_cor3_T, direct_cor4_T, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_WB2t = @time [get_S4(direct_cor2_WB2t, direct_cor3_WB2t, direct_cor4_WB2t, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    S4_conv = @time [get_S4(direct_cor2_WB2, foo, foo, kval, kval, 2.0, θ12, acos(cosθ13), ϕ23, ρval) for cosθ13 in cosθ13_arr, ϕ23 in ϕ23_arr]
    heatmap!(ax2_RF,  ϕ23_arr, cosθ13_arr, S4_RF',colorrange=(-maximum(abs.(S4_RF)),maximum(abs.(S4_RF))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_WB,  ϕ23_arr, cosθ13_arr, S4_WB',  colorrange=(-maximum(abs.(S4_WB)),maximum(abs.(S4_WB))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_WB2,  ϕ23_arr, cosθ13_arr, S4_WB2', colorrange=(-maximum(abs.(S4_WB2)),maximum(abs.(S4_WB2))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_M,  ϕ23_arr, cosθ13_arr, S4_M', colorrange=(-maximum(abs.(S4_M)),maximum(abs.(S4_M))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_T,  ϕ23_arr, cosθ13_arr, S4_T', colorrange=(-maximum(abs.(S4_T)),maximum(abs.(S4_T))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_WB2t,  ϕ23_arr, cosθ13_arr, S4_WB2t', colorrange=(-maximum(abs.(S4_WB2t)),maximum(abs.(S4_WB2t))), colormap=ColorSchemes.:bwr)
    heatmap!(ax2_conv,  ϕ23_arr, cosθ13_arr, S4_conv', colorrange=(-maximum(abs.(S4_conv)),maximum(abs.(S4_conv))), colormap=ColorSchemes.:bwr)

    Colorbar(fig[1,4], limits=(-maximum(abs.(S4_RF)),maximum(abs.(S4_RF))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[2, 4], limits=(-maximum(abs.(S4_WB)),maximum(abs.(S4_WB))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[3, 4], limits=(-maximum(abs.(S4_WB2)),maximum(abs.(S4_WB2))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[4, 4], limits=(-maximum(abs.(S4_M)),maximum(abs.(S4_M))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[5, 4], limits=(-maximum(abs.(S4_T)),maximum(abs.(S4_T))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[6, 4], limits=(-maximum(abs.(S4_WB2t)),maximum(abs.(S4_WB2t))), colormap=ColorSchemes.:bwr)
    Colorbar(fig[7, 4], limits=(-maximum(abs.(S4_conv)),maximum(abs.(S4_conv))), colormap=ColorSchemes.:bwr)


    display(fig)
end


begin
    ρval = 0.94
    k1 = k2 = 7.2
    plot_idx = 1
    for cosθ12 = ([cos(π/6), cos(1π/3), cos(π/4), cos(π/2), cos(2π/3), cos(5pi/6), cos(π)])

        for cosθ13 = ([cos(π/6), cos(1π/3), cos(π/4), cos(π/2), cos(2π/3), cos(5pi/6), cos(π)])
            ϕ23_arr = range(0.01, π-0.01, length = 400)

            fig = Figure(size = (800, 400))


            for (j, k3) in enumerate([2.0, 7.0])

                f = h5open("Processed_Data/S4/Mean/S4_rho_$(ρval)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(cosθ12).h5", "r")

                costheta13_data = f["costheta13"][]
                ϕ23_data = f["phi23"][]
                S4_data = f["S4"][]
                S4_err = f["S4_stderr"][]
        
                close(f)


                if cosθ12 ≈ 1.0
                    θ12 = acos(1.0 - 0.025/2)
                elseif cosθ12 ≈ -1.0
                    θ12 = acos(-1.0 + 0.025/2)
                else
                    θ12 = acos(cosθ12)
                end

                if cosθ13 ≈ 1.0
                    θ13 = acos(1.0 - 0.025/2)
                elseif cosθ13 ≈ -1.0
                    θ13 = acos(-1.0 + 0.025/2)
                else
                    θ13 = acos(cosθ13)
                end

                @show k1, k2, k3, θ12, θ13
                ax =  Axis(fig[1, j], ylabel = L"S_4", xlabel = L"\phi_{23}", title=L"S^{(4)}(k_1=7, k_2=7, k_3=%$(k3), \theta_{12}=%$(round(θ12/pi, digits=2))\pi, \theta_{13}=%$(round(acos(cosθ13)/pi, digits=2))\pi, \phi_{23})")
                S4_RF = @time [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_WB = @time [get_S4(direct_cor2_WB, direct_cor3_WB, direct_cor4_WB, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_WB2 = @time [get_S4(direct_cor2_WB2, direct_cor3_WB2, direct_cor4_WB2, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_M = @time [get_S4(direct_cor2_M, direct_cor3_M, direct_cor4_M, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_T = @time [get_S4(direct_cor2_T, direct_cor3_T, direct_cor4_T, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_WB2t = @time [get_S4(direct_cor2_WB2t, direct_cor3_WB2t, direct_cor4_WB2t, k1, k2, k3, θ12, θ13, 
                ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_L = @time [get_S4(direct_cor2_L, direct_cor3_L, direct_cor4_L, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_G = @time [get_S4(direct_cor2_G, direct_cor3_G, direct_cor4_G, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_conv2 = @time [get_S4(direct_cor2_WB2, foo, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
                S4_conv3 = @time [get_S4(direct_cor2_RF, direct_cor3_RF, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]



                lines!(ax, ϕ23_arr, S4_RF, linewidth = 3, label = "Rosenfeld")
                lines!(ax, ϕ23_arr, S4_WB, linewidth = 3, label = "White-Bear")
                lines!(ax, ϕ23_arr, S4_WB2, linewidth = 3, label = "White-Bear II")
                lines!(ax, ϕ23_arr, S4_M,  linewidth = 3, label = "Malijevský", linestyle = :dash)
                lines!(ax, ϕ23_arr, S4_T, linewidth = 3, label = "Tarazona", linestyle = :dash)
                lines!(ax, ϕ23_arr, S4_WB2t, linewidth = 3, label = "White-Bear II (tensor)", linestyle = :dash)
                lines!(ax, ϕ23_arr, S4_conv2, linewidth = 3, label = "Convolution 2-body (RF)", linestyle = :dot)
                lines!(ax, ϕ23_arr, S4_L, linewidth = 2, label = "Lutsko", linestyle= :dashdot)
                lines!(ax, ϕ23_arr, S4_G, linewidth = 2, label = "Gül", linestyle= :dashdot)
                
                # lines!(ax, ϕ23_arr, S4_conv3, linewidth = 2, label = "Convolution 3-body (RF)", linestyle= :dash, color= :blue)
                


                spl = Spline2D(costheta13_data, ϕ23_data, S4_data)
                S4_data = [spl(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]
                spl_err = Spline2D(costheta13_data, ϕ23_data, S4_err)
                S4_err = [spl_err(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]

                scatter!(ax, ϕ23_data, S4_data, color = :black, label = "MC Data", markersize=10)
                errorbars!(ax, ϕ23_data, S4_data, S4_err, color = :black, whiskerwidth=8)
                
                if j == 2
                    global plot_idx
                    ax3 = Axis(fig[1, 3], title="Plot $(plot_idx)")
                    hidedecorations!(ax3)
                    hidespines!(ax3)
                    Legend(fig[1,3], ax, framevisible=false)
                    
                    plot_idx += 1
                end
                
            end
            display(fig)
            save("Plots/S4_2_7_costheta12_$(cosθ12)_costheta13_$(cosθ13).png", fig, px_per_unit = 5)
        end
    end
end


begin
    theta12_mat = [
        π/3 π/3
        π/3 π/3
        2π/3 2π/3
    ]
    theta13_mat = [
        π/3 π/3
        2π/3 2π/3
        2π/3 2π/3
    ]
    k3_mat = [
        2.0 7.0
        2.0 7.0
        2.0 7.0
    ]

    ρval = 0.94
    k1 = k2 = 7.2

    fig = Figure(size = (500, 750))
    ax11 = Axis(fig[1, 1], ylabel = L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
    ax12 = Axis(fig[1, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
    ax21 = Axis(fig[2, 1], ylabel = L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
    ax22 = Axis(fig[2, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
    ax31 = Axis(fig[3, 1], ylabel = L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xlabel= L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0,π,nothing,nothing))
    ax32 = Axis(fig[3, 2], xlabel= L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0,π,nothing,nothing))
    axs_mat = [
        ax11 ax12
        ax21 ax22
        ax31 ax32
    ]
    for row = 1:3
        for col = 1:2
            ax = axs_mat[row, col]
            θ12 = theta12_mat[row, col]
            θ13 = theta13_mat[row, col]
            k3 = k3_mat[row, col]
            
            ϕ23_arr = range(0.01, π-0.01, length = 400)

            S4_RF = [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB = [get_S4(direct_cor2_WB, direct_cor3_WB, direct_cor4_WB, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB2 = [get_S4(direct_cor2_WB2, direct_cor3_WB2, direct_cor4_WB2, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_M = [get_S4(direct_cor2_M, direct_cor3_M, direct_cor4_M, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_T = [get_S4(direct_cor2_T, direct_cor3_T, direct_cor4_T, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB2t = [get_S4(direct_cor2_WB2t,direct_cor3_WB2t,direct_cor4_WB2t,k1,k2,k3 ,θ12 ,θ13 ,ϕ23 ,ρval) for ϕ23 in ϕ23_arr]
            S4_L = [get_S4(direct_cor2_L, direct_cor3_L, direct_cor4_L, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_G = [get_S4(direct_cor2_G, direct_cor3_G, direct_cor4_G, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv = [get_S4(direct_cor2_RF,foo,foo,k1,k2,k3 ,θ12 ,θ13 ,ϕ23 ,ρval) for ϕ23 in ϕ23_arr]

            lines!(ax, ϕ23_arr, S4_RF, linewidth = 3, label = "R")
            lines!(ax, ϕ23_arr, S4_WB, linewidth = 3, label = "WBI")
            lines!(ax, ϕ23_arr, S4_WB2, linewidth = 3, label = "WBII")
            lines!(ax, ϕ23_arr, S4_M,  linewidth = 3, label = "M", linestyle = :dash)
            lines!(ax, ϕ23_arr, S4_T, linewidth = 3, label = "T", linestyle = :dash)
            lines!(ax, ϕ23_arr, S4_WB2t, linewidth = 3, label = "WBIIt", linestyle = :dash)
            lines!(ax, ϕ23_arr, S4_L, linewidth = 2, label = "L", linestyle= :dashdot)
            lines!(ax, ϕ23_arr, S4_G, linewidth = 2, label = "G", linestyle= :dashdot)
            lines!(ax, ϕ23_arr, S4_conv, linewidth = 3, label = "conv", linestyle = :dot)

            f = h5open("Processed_Data/S4/Mean/S4_rho_$(ρval)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(cos(θ12)).h5", "r")
            costheta13_data = f["costheta13"][]
            ϕ23_data = f["phi23"][]
            S4_data = f["S4"][]
            S4_err = f["S4_stderr"][]
            close(f)


            cosθ13 = cos(θ13)
            spl = Spline2D(costheta13_data, ϕ23_data, S4_data)
            S4_data = [spl(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]
            spl_err = Spline2D(costheta13_data, ϕ23_data, S4_err)
            S4_err = [spl_err(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]

            scatter!(ax, ϕ23_data, S4_data, color = :black, label = "MC", markersize=10)
            errorbars!(ax, ϕ23_data, S4_data, S4_err, color = :black, whiskerwidth=8)

            int12 = round(Int, 3θ12/pi)
            int12str =  int12 == 1 ? "" : int12 == 2 ? "2" : string(int12)
            int13 = round(Int, 3θ13/pi)
            int13str = int13 == 1 ? "" : int13 == 2 ? "2" : string(int13)

            if col == 1 && row == 2
                text!(ax, 1.0, 1.0, text=L"k_3 = %$(k3)", align=(:right, :top), offset=(-10,-27), space=:relative)
                text!(ax, 1.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:right, :top), offset=(-10,-44), space=:relative)
                text!(ax, 1.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:right, :top), offset=(-10,-61), space=:relative)
            else
                text!(ax, 0.0, 1.0, text=L"k_3 = %$(k3)", align=(:left, :top), offset=(10,-27), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:left, :top), offset=(10,-44), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:left, :top), offset=(10,-61), space=:relative)
            end
        end
    end
    for (ax, label) in zip([ax11, ax12, ax21, ax22, ax31, ax32], ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"])
        text!(ax, 0.0, 1.0, text=label, align=(:left, :top), offset=(10,-10), space=:relative)
    end

    axislegend(ax12, framevisible=false, position=(0.7, 1.0), rowgap=-3)

    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 10)


    display(fig)
    save("Plots/S4_all.pdf", fig, pt_per_unit = 5)
    save("Plots/S4_all.png", fig)

end



begin
    theta12_mat = [
        π/3 π/3
        2π/3 2π/3
    ]
    theta13_mat = [
        π/3 π/3
        2π/3 2π/3
    ]
    k3_mat = [
        2.0 7.0
        2.0 7.0
    ]

    ρval = 0.94
    k1 = k2 = 7.2

    fig = Figure(size = (500, 500))
    ax11 = Axis(fig[1, 1], ylabel = L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
    ax12 = Axis(fig[1, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0,π,nothing,nothing))
        ax31 = Axis(fig[2, 1], ylabel = L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xlabel= L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0,π,nothing,nothing))
    ax32 = Axis(fig[2, 2], xlabel= L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0,π,nothing,nothing))
    axs_mat = [
        ax11 ax12
        ax31 ax32
    ]
    for row = 1:2
        for col = 1:2
            ax = axs_mat[row, col]
            θ12 = theta12_mat[row, col]
            θ13 = theta13_mat[row, col]
            k3 = k3_mat[row, col]
            
            ϕ23_arr = range(0.01, π-0.01, length = 400)

            S4_RF = [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv3 = [get_S4(direct_cor2_RF, direct_cor3_RF, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv2 = [get_S4(direct_cor2_RF, foo, foo, k1,k2,k3 ,θ12 ,θ13 ,ϕ23 ,ρval) for ϕ23 in ϕ23_arr]

            lines!(ax, ϕ23_arr, S4_RF, linewidth = 3, label = "R")
            lines!(ax, ϕ23_arr, S4_conv3, linewidth = 3, label = "3-conv", color= :green)
            lines!(ax, ϕ23_arr, S4_conv2, linewidth = 3, label = "2-conv", color= :red)


            int12 = round(Int, 3θ12/pi)
            int12str =  int12 == 1 ? "" : int12 == 2 ? "2" : string(int12)
            int13 = round(Int, 3θ13/pi)
            int13str = int13 == 1 ? "" : int13 == 2 ? "2" : string(int13)

            # if col == 1 && row == 2
                # text!(ax, 1.0, 1.0, text=L"k_3 = %$(k3)", align=(:right, :top), offset=(-10,-10), space=:relative)
                # text!(ax, 1.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:right, :top), offset=(-10,-27), space=:relative)
                # text!(ax, 1.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:right, :top), offset=(-10,-44), space=:relative)
            # else
                text!(ax, 0.0, 1.0, text=L"k_3 = %$(k3)", align=(:left, :top), offset=(10,-27), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:left, :top), offset=(10,-44), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:left, :top), offset=(10,-61), space=:relative)
            # end
        end
    end
    for (ax, label) in zip([ax11, ax12,  ax31, ax32], ["(a)", "(b)", "(c)", "(d)"])
        text!(ax, 0.0, 1.0, text=label, align=(:left, :top), offset=(10,-10), space=:relative)
    end

    axislegend(ax12, framevisible=false, position=(0.8, 1.0), rowgap=-3)

    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 10)


    display(fig)
    save("Plots/S4_convs.pdf", fig, pt_per_unit = 5)
    save("Plots/S4_convs.png", fig)
end