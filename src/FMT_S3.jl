println("Starting the script")


using Pkg; Pkg.activate(".")

using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie, LaTeXStrings, Dierckx

mytheme = Theme( Axis = (xtickalign = 1,
xminortickalign = 1,
xticksmirrored = true,
ytickalign = 1,
yminortickalign = 1,
yticksmirrored = true, 
xgridvisible = false, 
ygridvisible = false,))
set_theme!(merge(mytheme, theme_latexfonts()))
Makie.theme(:palette).color[] = ColorSchemes.:ColorSchemes.:seaborn_muted6

colors = [ColorSchemes.:Dark2_8[(i-1)/7] for i in 1:8]
linestyles = [ :solid, (:dash, :dense), :dash, (:dash, :loose), :dashdot, (:dot, :dense), :dot, (:dot, :loose)]



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

# @show direct_cor2_RF(0.1, 0.94)*0.94
# @show direct_cor2_WB(0.1, 0.94)*0.94
# @show direct_cor2_WB2(0.1, 0.94)*0.94
# @show direct_cor2_M(0.1, 0.94)*0.94
# @show direct_cor2_KR(0.1, 0.94)*0.94
# @show direct_cor2_T(0.1, 0.94)*0.94
# @show direct_cor2_WB2t(0.1, 0.94)*0.94
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

# @show direct_cor3_RF(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_WB(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_WB2(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_M(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_KR(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_T(0.1, 0.1, acos(-0.175), 0.94)*0.94^2
# @show direct_cor3_WB2t(0.1, 0.1, acos(-0.175), 0.94)*0.94^2

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




Nsweeps = 20



for rho in [0.8, 0.94]
    for sweep = 1:Nsweeps
        if !isfile("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
            println("File not found for rho $(rho) and sweep $(sweep)")
            continue
        end
        data = readdlm("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
        k1 = data[:, 1]
        k2 = data[:, 2]

        cos_theta = data[:, 3]
        cos_theta[cos_theta .≈ 1.0] .= 0.99
        cos_theta[cos_theta .≈ -1.0] .= -0.99

        
        S3_mean = data[:, 4]
        S3_std_err = data[:, 5]
        t = eachindex(S3_mean)

        t_interpl = 1:0.01:length(S3_mean)
        k1_interpl = Spline1D(t, k1).(t_interpl)
        k2_interpl = Spline1D(t, k2).(t_interpl)
        cos_theta_interpl = Spline1D(t, cos_theta).(t_interpl)
        

        S3_RF = @. get_S3(direct_cor2_RF, direct_cor3_RF, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB = @. get_S3(direct_cor2_WB, direct_cor3_WB, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB2 = @. get_S3(direct_cor2_WB2, direct_cor3_WB2, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_M = @. get_S3(direct_cor2_M, direct_cor3_M, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_T = @. get_S3(direct_cor2_T, direct_cor3_T, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB2t = @. get_S3(direct_cor2_WB2t, direct_cor3_WB2t, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        # S3_KR = @. get_S3(direct_cor2_KR, direct_cor3_KR, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_L = @. get_S3(direct_cor2_L, direct_cor3_L, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_G = @. get_S3(direct_cor2_G, direct_cor3_G, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
 

        # find which type of sweep it is by seeing if costheta is constant
        f = Figure(size=(400, 400))

        if cos_theta[1] == cos_theta[end]

            #check if k2 is constant
            cosθ = round(cos_theta[1], digits=5)
            θdivπ = round(acos(cosθ)/π, digits=3)
            if k2[1] == k2[end]
                println("Sweep $(sweep) is a k1 sweep")
                ax = Axis(f[1, 1], xlabel=L"k_1", ylabel=L"S_3(k_1,k_2=%$(k2[1]),\,\theta=%$(θdivπ)\pi)")
            else
                println("Sweeping k1=k2=k")
                ax = Axis(f[1, 1], xlabel=L"k", ylabel=L"S_3(k_1=k_2=k,\,\theta=%$(θdivπ)\pi)")
            end
            lines!(ax, k1_interpl, S3_RF, linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, k1_interpl, S3_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, k1_interpl, S3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, k1_interpl, S3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, k1_interpl, S3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, k1_interpl, S3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, k1_interpl, S3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, k1_interpl, S3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])
            scatter!(ax, k1, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, k1, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
        else
            #sweeping costheta at fixed k1 and k2
            println("Sweep $(sweep) is a costheta sweep")
            k1val = k1[1]
            k2val = k2[1]
            @assert k1val == k2val
            ax = Axis(f[1, 1], xlabel=L"\cos(\theta)", ylabel=L"S_.3(k_1=k_2=%$(k1val),\,\cos\,\theta)")

            lines!(ax, cos_theta_interpl, S3_RF,linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, cos_theta_interpl, S3_WB, linewidth = 3, label = "WBI" , color=colors[2], linestyle=linestyles[2])
            lines!(ax, cos_theta_interpl, S3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, cos_theta_interpl, S3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, cos_theta_interpl, S3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, cos_theta_interpl, S3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, cos_theta_interpl, S3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, cos_theta_interpl, S3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])
            scatter!(ax, cos_theta, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, cos_theta, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
        end
        
        axislegend(ax, position = :rt,framevisible = false)
        # f[1,2] = Legend(f, ax, framevisible=false)

        display(f)
        save("Plots/S3_sweep_$(sweep)_rho_$(rho).pdf", f)

    end

end

function transpose_matrix(matrix)
    return [matrix[j, i] for i in 1:size(matrix, 2), j in 1:size(matrix, 1)]
end


begin
fig = Figure(size=(800, 600))
ax11 = Axis(fig[1, 1], xlabel=L"k_1", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(5,15,nothing,nothing))
ax12 = Axis(fig[1, 2], xlabel=L"k_1", )
ax13 = Axis(fig[1, 3], xlabel=L"k_1", limits=(4,15,nothing,nothing))
# ax21 = Axis(fig[2, 1], xlabel=L"k", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(4,13,nothing,nothing))
# ax22 = Axis(fig[2, 2], xlabel=L"k", limits=(4,11,nothing,nothing) )
# ax23 = Axis(fig[2, 3], xlabel=L"k", limits=(4,14,nothing,nothing))
ax31 = Axis(fig[2, 1], xlabel=L"\cos(\theta)", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(nothing, nothing, nothing, 0.001))
ax32 = Axis(fig[2, 2], xlabel=L"\cos(\theta)", )
ax33 = Axis(fig[2, 3], xlabel=L"\cos(\theta)", )

sweeps = [  2 4 5;
            # 15 17 18;
            8 11 13;    

            ]
axs = [ax11 ax12 ax13;
    #    ax21 ax22 ax23;
       ax31 ax32 ax33]
for i in 1:2
    for j in 1:3
        ax = axs[i, j]
        sweep = sweeps[i, j]
        rho = 0.94
        if !isfile("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
            println("File not found for rho $(rho) and sweep $(sweep)")
            continue
        end
        data = readdlm("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
        k1 = data[:, 1]
        k2 = data[:, 2]

        cos_theta = data[:, 3]
        cos_theta[cos_theta .≈ 1.0] .= 0.99
        cos_theta[cos_theta .≈ -1.0] .= -0.99

        
        S3_mean = data[:, 4]
        S3_std_err = data[:, 5]
        t = eachindex(S3_mean)
        t_interpl = 1:0.1:length(S3_mean)
        k1_interpl = Spline1D(t, k1).(t_interpl)
        k2_interpl = Spline1D(t, k2).(t_interpl)
        cos_theta_interpl = Spline1D(t, cos_theta).(t_interpl)

        S3_RF = @. get_S3(direct_cor2_RF, direct_cor3_RF, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB = @. get_S3(direct_cor2_WB, direct_cor3_WB, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB2 = @. get_S3(direct_cor2_WB2, direct_cor3_WB2, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_M = @. get_S3(direct_cor2_M, direct_cor3_M, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_T = @. get_S3(direct_cor2_T, direct_cor3_T, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_WB2t = @. get_S3(direct_cor2_WB2t, direct_cor3_WB2t, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        # S3_KR = @. get_S3(direct_cor2_KR, direct_cor3_KR, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_L = @. get_S3(direct_cor2_L, direct_cor3_L, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_G = @. get_S3(direct_cor2_G, direct_cor3_G, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        foo(x...) = 0.0 
        S3_conv = @. get_S3(direct_cor2_RF, foo, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        j==1&&i==3 && @show S3_conv
        if i == 2
            #sweeping costheta at fixed k1 and k2

            k1val = k1[1]
            k2val = k2[1]
            @assert k1val == k2val
            lines!(ax, cos_theta_interpl, S3_RF,linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, cos_theta_interpl, S3_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, cos_theta_interpl, S3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, cos_theta_interpl, S3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, cos_theta_interpl, S3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, cos_theta_interpl, S3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, cos_theta_interpl, S3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, cos_theta_interpl, S3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])
            lines!(ax, cos_theta_interpl, S3_conv, linewidth = 3, label = "conv", color=:black, linestyle=:dash)
            scatter!(ax, cos_theta, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, cos_theta, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
            text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-50,-10), space=:relative)
        else
            # sweeping k1 at fixed k2 or k1 and k2
            lines!(ax, k1_interpl, S3_RF,linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, k1_interpl, S3_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, k1_interpl, S3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, k1_interpl, S3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, k1_interpl, S3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, k1_interpl, S3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, k1_interpl, S3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, k1_interpl, S3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])
            lines!(ax, k1_interpl, S3_conv, linewidth = 3, label = "conv", color=:black, linestyle=:dash)

            scatter!(ax, k1, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, k1, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
            text!(ax, 1, 1, text=L"\theta=%$(round(acos(cos_theta[1])/π, digits=3))\pi", align=(:right, :top), offset=(-50,-10), space=:relative)
        end
        if i==1
            text!(ax, 1, 1, text=L"k_2=%$(k2[1])", align=(:right, :top), offset=(-50,-30), space=:relative)
        end
            

    end
end
for (ax,label) in zip(transpose_matrix(axs)[:], ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"])
    text!(ax, 1, 1, text=label, align=(:right, :top), offset=(-10,-10), space=:relative)
end
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)
axislegend(ax13, position = :rc, framevisible = false, patchsize=(30,20), rowgap=-2)


display(fig)
save("Plots/S3_all_sweeps.pdf", fig)
save("Plots/S3_all_sweeps.png", fig)

end


begin
    fig = Figure(size=(800, 250), figure_padding=(5, 5, 5, 5))
    # ax11 = Axis(fig[1, 1], xlabel=L"k_1", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(5,15,nothing,nothing))
    # ax12 = Axis(fig[1, 2], xlabel=L"k_1", )
    # ax13 = Axis(fig[1, 3], xlabel=L"k_1", limits=(4,15,nothing,nothing))
    ax21 = Axis(fig[1, 1], xlabel=L"k", ylabel=L"\rho^2c_3(k_1,\, k_2,\,\theta)", limits=(0,5,nothing,nothing))
    ax22 = Axis(fig[1, 2], xlabel=L"\cos(\theta)",  limits=(-1, 1,nothing, nothing))
    ax23 = Axis(fig[1, 3], xlabel=L"\cos(\theta)", limits=(-1, 1,nothing, nothing))
    ax24 = Axis(fig[1, 4], xlabel=L"\cos(\theta)", limits=(-1, 1,nothing, nothing))
    # ax31 = Axis(fig[2, 1], xlabel=L"\cos(\theta)", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(nothing, nothing, nothing, 0.001))
    # ax32 = Axis(fig[2, 2], xlabel=L"\cos(\theta)", )
    # ax33 = Axis(fig[2, 3], xlabel=L"\cos(\theta)", )
    
    sweeps = [ 15 8 11 13;   
                ]
    axs = [
           ax21 ax22 ax23 ax24;
           ]
    for i in 1:1
        for j in 1:4
            ax = axs[i, j]
            sweep = sweeps[i, j]
            rho = 0.94
            if !isfile("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
                println("File not found for rho $(rho) and sweep $(sweep)")
                continue
            end
            data = readdlm("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt")
            k1 = data[:, 1]
            k2 = data[:, 2]
    
            cos_theta = data[:, 3]
            cos_theta[cos_theta .≈ 1.0] .= 0.99
            cos_theta[cos_theta .≈ -1.0] .= -0.99
    
            
            S3_mean = data[:, 4]
            S3_std_err = data[:, 5]


            S2data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")

            S2_k = S2data[:, 1]
            S2_mean = S2data[:, 2]
            mask = isfinite.(S2_mean)
            Sk_spl = Spline1D(S2_k[mask], S2_mean[mask], k=1)


            SSS = @. Sk_spl(k1) * Sk_spl(k2) * Sk_spl(sqrt(k1^2 + k2^2 + 2*k1*k2*cos_theta))

            ρ²C3 = @. S3_mean / SSS - 1 
            ρ²C3_err = @. S3_std_err / SSS 


            t = eachindex(S3_mean)
            t_interpl = 1:0.1:length(S3_mean)
            k1_interpl = Spline1D(t, k1).(t_interpl)
            k2_interpl = Spline1D(t, k2).(t_interpl)
            cos_theta_interpl = Spline1D(t, cos_theta).(t_interpl)
    
            C3_RF = @. direct_cor3_RF(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_WB = @. direct_cor3_WB(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_WB2 = @. direct_cor3_WB2(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_M = @. direct_cor3_M(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_T = @. direct_cor3_T(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_WB2t = @. direct_cor3_WB2t(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_L = @. direct_cor3_L(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2
            C3_G = @. direct_cor3_G(k1_interpl, k2_interpl, acos(cos_theta_interpl), rho) * rho^2



            foo(x...) = 0.0 
            j==1&&i==3 && @show S3_conv
            if j == 3 || j == 2 || j == 4
                #sweeping costheta at fixed k1 and k2
    
                k1val = k1[1]
                k2val = k2[1]
                @assert k1val == k2val
                lines!(ax, cos_theta_interpl, C3_RF,linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, cos_theta_interpl, C3_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, cos_theta_interpl, C3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, cos_theta_interpl, C3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, cos_theta_interpl, C3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, cos_theta_interpl, C3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, cos_theta_interpl, C3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, cos_theta_interpl, C3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])

                scatter!(ax, cos_theta, ρ²C3, label="MC Data", color=:black)
                errorbars!(ax, cos_theta, ρ²C3, ρ²C3_err, color=:black, whiskerwidth=4)
                # text!(ax, 1, 0, text=L"k_1=k_2=%$(k1val)", align=(:right, :bottom), offset=(-50,10), space=:relative)
            else
                # sweeping k1 at fixed k2 or k1 and k2
                lines!(ax, k1_interpl, C3_RF,linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, k1_interpl, C3_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, k1_interpl, C3_WB2, linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, k1_interpl, C3_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, k1_interpl, C3_T, linewidth = 3, label = "T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, k1_interpl, C3_WB2t, linewidth = 3, label = "WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, k1_interpl, C3_L, linewidth = 3, label = "L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, k1_interpl, C3_G, linewidth = 3, label = "G", color=colors[8], linestyle=linestyles[8])
                scatter!(ax, k1[1:2:end], ρ²C3[1:2:end], label="MC Data", color=:black)
                errorbars!(ax, k1[1:2:end], ρ²C3[1:2:end], ρ²C3_err[1:2:end], color=:black, whiskerwidth=4)
            end
            if j==1 
                text!(ax, 0, 1, text=L"k_1=k_2=k", align=(:left, :top), offset=(10,-27), space=:relative)
                text!(ax, 0, 1, text=L"\theta=%$(round(acos(cos_theta[1])/π, digits=3))\pi", align=(:left, :top), offset=(10,-10), space=:relative)
                text!(ax, 1, 0, text="(a)", align=(:right, :bottom), offset=(-10,10), space=:relative)

            elseif j==2
                text!(ax, 0, 1, text=L"k_1=k_2=%$(k1val)", align=(:left, :top), offset=(10,-10), space=:relative)
                text!(ax, 1, 0, text="(b)", align=(:right, :bottom), offset=(-10,10), space=:relative)
            elseif j==3 
                text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-40,-10), space=:relative)
                text!(ax, 1, 1, text="(c)", align=(:right, :top), offset=(-10,-10), space=:relative)
            elseif j==4
                text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-40,-10), space=:relative)
                text!(ax, 1, 1, text="(d)", align=(:right, :top), offset=(-10,-10), space=:relative)
            end
        end
    end






    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    
    Legend(fig[1, 5], axs[1, 2], framevisible=false, tellwidth=true, patchsize=(30,20), rowgap=-2)    
    
    display(fig)
    save("Plots/C3_all_sweeps.pdf", fig)
    save("Plots/C3_all_sweeps.png", fig)
    
    end