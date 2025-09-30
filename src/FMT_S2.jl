println("Starting the script")
using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie, LaTeXStrings
using OrnsteinZernike

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

function find_analytical_C_PY(k, η)
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

colors = [ColorSchemes.:Dark2_8[(i-1)/7] for i in 1:8]
linestyles = [ :solid, (:dash, :dense), :dash, (:dash, :loose), :dashdot, (:dot, :dense), :dot, (:dot, :loose)]


function get_S2(direct_cor2, k, ρ)
    c2 = direct_cor2(k, ρ)
    S2 = 1 + ρ * c2 / (1 - ρ * c2)
    return S2
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
# φ3_tensor =  3*(-n2 * n2dotn2 + vTv * n2*tr(n_squared) - tr(n_cubed)) / (16*Pi*(1-n3)^2) # tarazona: dont use: should use traceless tensor

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
@assert abs(direct_cor2_RF(0.1, 0.94) - find_analytical_C_PY(0.1, 0.94/6*π)) < 1e-5
@assert abs(direct_cor2_KR(0.1, 0.94) - find_analytical_C_PY(0.1, 0.94/6*π)) < 1e-5
# @assert abs(direct_cor2_T(0.1, 0.94) - find_analytical_C_k(0.1, 0.94/6*π)) < 1e-5




@assert direct_cor2_RF(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB(0.1, 0.94) * 0.94 ≈ -51.88510371086945
@assert direct_cor2_WB2(0.1, 0.94) * 0.94 ≈ -51.885375608479336
@assert direct_cor2_M(0.1, 0.94) * 0.94 ≈ -51.61924863988694
@assert direct_cor2_KR(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_T(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB2t(0.1, 0.94) * 0.94 ≈ -51.885375608479336





function c_from_s(Sk, rho)
    ρCk = @. (Sk-1)/Sk
    return ρCk/rho
end
ρval = 0.94

dims = 3; kBT = 1.0
potential = HardSpheres(1.0)
system = SimpleLiquid(dims, ρval, kBT, potential)
closure1 = Verlet()
method = NgIteration(M=10000, dr=0.01)
sol1 = @time solve(system, closure1, method);
closure2 = ModifiedHypernettedChain(ρval/6*π)
sol2 = @time solve(system, closure2, method);

begin 
    fig = Figure(size = (500, 300))
    ax11 = Axis(fig[1, 1], xlabel = L"k", ylabel = L"\rho c^{(2)}(k)", limits=(0,10,nothing,nothing))
    ax12 = Axis(fig[1, 2], xlabel = L"k", ylabel = L"\rho c^{(2)}(k)", limits=(0,2,-53,-45))


    for ax in [ax11, ax12]
        k = range(0.01, 15, length = 200)
        c2_RF = [direct_cor2_RF(ki, ρval)*ρval for ki in k]
        c2_WB = [direct_cor2_WB(ki, ρval)*ρval for ki in k]
        c2_WB2 = [direct_cor2_WB2(ki, ρval)*ρval for ki in k]
        c2_M = [direct_cor2_M(ki, ρval)*ρval for ki in k]
        c2_L = [direct_cor2_L(ki, ρval)*ρval for ki in k]
        c2_G = [direct_cor2_G(ki, ρval)*ρval for ki in k]

        lines!(ax, k, c2_RF,  linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
        lines!(ax, k, c2_WB, linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
        lines!(ax, k, c2_WB2,  linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
        lines!(ax, k, c2_M, linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
        lines!(ax, k, c2_L, linewidth = 3, label = "L", color=colors[5], linestyle=linestyles[5])
        lines!(ax, k, c2_G, linewidth = 3, label = "G", color=colors[6], linestyle=linestyles[6])
        lines!(ax, sol1.k, sol1.ck,  linewidth = 3, label = "V", color=colors[7], linestyle=linestyles[7])
        lines!(ax, sol2.k, sol2.ck,  linewidth = 3, label = "MHNC", color=colors[8], linestyle=linestyles[8])
    end
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(ρval).txt")
    k_data = data[:, 1]
    S2_data = data[:, 2]

    c2_data = c_from_s(S2_data, ρval)*ρval
    scatter!(ax11, k_data, c2_data, color = :black, label = "MC Data", markersize=6)
    scatter!(ax12, k_data, c2_data, color = :black, label = "MC Data", markersize=10)

    axislegend(ax11, position = :rb, framevisible=false,patchsize=(30,20), rowgap=-2)

    for (ax, label) in zip([ax11, ax12], ["(a)", "(b)"])
        text!(ax,
        0, 1,
        text= label,
        align=(:left, :top),
        offset=(10, -10),
        space=:relative
        )
    end


    display(fig)
    save("Plots/ck_rho_0.94.pdf", fig)
    save("Plots/ck_rho_0.94.png", fig)
end




fig = Figure(size = (500, 500))
ax1 = Axis(fig[1, 1:2], xlabel = L"k", ylabel = L"S(k)", limits=(0,20, nothing,nothing))
ax2 = Axis(fig[2, 1], xlabel = L"k", ylabel = L"S(k)", limits=(6.5,7.5, 2,3.3))
ax3 = Axis(fig[2, 2], xlabel = L"k", ylabel = L"S(k)", limits=(12,14, 1.2,1.45), ylabelvisible=false)
k = range(0.01, 20, length = 1000)
S2_RF = [get_S2(direct_cor2_RF, ki, ρval) for ki in k]
S2_WB = [get_S2(direct_cor2_WB, ki, ρval) for ki in k]
S2_WB2 = [get_S2(direct_cor2_WB2, ki, ρval) for ki in k]
S2_M = [get_S2(direct_cor2_M, ki, ρval) for ki in k]
S2_L = [get_S2(direct_cor2_L, ki, ρval) for ki in k]
S2_G = [get_S2(direct_cor2_G, ki, ρval) for ki in k]
data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(ρval).txt")
k_data = data[:, 1]
S2_data = data[:, 2]

for ax in [ax1, ax2, ax3]
    lines!(ax, k, S2_RF,  linewidth = 3, label = "R", color=colors[1], linestyle=linestyles[1])
    lines!(ax, k, S2_WB,  linewidth = 3, label = "WBI", color=colors[2], linestyle=linestyles[2])
    lines!(ax, k, S2_WB2,  linewidth = 3, label = "WBII", color=colors[3], linestyle=linestyles[3])
    lines!(ax, k, S2_M,  linewidth = 3, label = "M", color=colors[4], linestyle=linestyles[4])
    lines!(ax, k, S2_L,  linewidth = 3, label = "L", color=colors[5], linestyle=linestyles[5])
    lines!(ax, k, S2_G,  linewidth = 3, label = "G", color=colors[6], linestyle=linestyles[6])
    lines!(ax, sol1.k, sol1.Sk,  linewidth = 3, label = "V", color=colors[7], linestyle=linestyles[7])
    lines!(ax, sol2.k, sol2.Sk,  linewidth = 3, label = "MHNC", color=colors[8], linestyle=linestyles[8])
end
scatter!(ax1, k_data, S2_data, color = :black, label = "MC Data", markersize=6)
scatter!(ax2, k_data, S2_data, color = :black, label = "MC Data", markersize=10)
scatter!(ax3, k_data, S2_data, color = :black, label = "MC Data", markersize=10)
text!(ax2,
0.5, 0.0, 
text= "first peak",
align=(:center, :bottom),
offset=(0, 10),
space=:relative
)

text!(ax3,
0.5, 0.0,
text= "second peak",
align=(:center, :bottom),
offset=(0, 10),
space=:relative
)

for (ax, label) in zip([ax1, ax2, ax3], ["(a)", "(b)", "(c)"])
    text!(ax,
    0, 1,
    text= label,
    align=(:left, :top),
    offset=(10, -10),
    space=:relative
    )
end


# axislegend(ax, position = :rt)
axislegend(ax1, position = :rt, nbanks=2, framevisible=false, patchsize=(30,20), rowgap=-2)
display(fig)

save("Plots/Sk_rho_0.94.pdf", fig)
save("Plots/Sk_rho_0.94.png", fig)
