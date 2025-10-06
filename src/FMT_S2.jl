println("Starting the script")
using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie, LaTeXStrings
using OrnsteinZernike

include("plot_settings.jl")
include("fmt_defs.jl")


direct_cor2_RF = compute_pair_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor2_WB = compute_pair_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor2_WB2 = compute_pair_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor2_M = compute_pair_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor2_KR = compute_pair_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor2_T = compute_pair_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor2_WB2t = compute_pair_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor2_L = compute_pair_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor2_G = compute_pair_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@show direct_cor2_RF(0.1, 0.94) * 0.94
@show direct_cor2_WB(0.1, 0.94) * 0.94
@show direct_cor2_WB2(0.1, 0.94) * 0.94
@show direct_cor2_M(0.1, 0.94) * 0.94
@show direct_cor2_KR(0.1, 0.94) * 0.94
@show direct_cor2_T(0.1, 0.94) * 0.94
@show direct_cor2_WB2t(0.1, 0.94) * 0.94
@show direct_cor2_L(0.1, 0.94) * 0.94
@show direct_cor2_G(0.1, 0.94) * 0.94
@assert abs(direct_cor2_RF(0.1, 0.94) - find_analytical_C_PY(0.1, 0.94 / 6 * π)) < 1e-5
@assert abs(direct_cor2_KR(0.1, 0.94) - find_analytical_C_PY(0.1, 0.94 / 6 * π)) < 1e-5
# @assert abs(direct_cor2_T(0.1, 0.94) - find_analytical_C_k(0.1, 0.94/6*π)) < 1e-5

@assert direct_cor2_RF(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB(0.1, 0.94) * 0.94 ≈ -51.88510371086945
@assert direct_cor2_WB2(0.1, 0.94) * 0.94 ≈ -51.885375608479336
@assert direct_cor2_M(0.1, 0.94) * 0.94 ≈ -51.61924863988694
@assert direct_cor2_KR(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_T(0.1, 0.94) * 0.94 ≈ -58.16983018151516
@assert direct_cor2_WB2t(0.1, 0.94) * 0.94 ≈ -51.885375608479336

ρval = 0.94
dims = 3;
kBT = 1.0;
potential = HardSpheres(1.0)
system = SimpleLiquid(dims, ρval, kBT, potential)
closure1 = Verlet()
method = NgIteration(M=10000, dr=0.01)
sol1 = @time solve(system, closure1, method);
closure2 = ModifiedHypernettedChain(ρval / 6 * π)
sol2 = @time solve(system, closure2, method);

begin
    fig = Figure(size=(500, 300))
    ax11 = Axis(fig[1, 1], xlabel=L"k", ylabel=L"\rho c^{(2)}(k)", limits=(0, 10, nothing, nothing))
    ax12 = Axis(fig[1, 2], xlabel=L"k", ylabel=L"\rho c^{(2)}(k)", limits=(0, 2, -59, -45))


    for ax in [ax11, ax12]
        k = range(0.01, 15, length=200)
        c2_RF = [direct_cor2_RF(ki, ρval) * ρval for ki in k]
        c2_WB = [direct_cor2_WB(ki, ρval) * ρval for ki in k]
        c2_WB2 = [direct_cor2_WB2(ki, ρval) * ρval for ki in k]
        c2_M = [direct_cor2_M(ki, ρval) * ρval for ki in k]
        c2_L = [direct_cor2_L(ki, ρval) * ρval for ki in k]
        c2_G = [direct_cor2_G(ki, ρval) * ρval for ki in k]

        lines!(ax, k, c2_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
        lines!(ax, k, c2_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
        lines!(ax, k, c2_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
        lines!(ax, k, c2_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
        lines!(ax, k, c2_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
        lines!(ax, k, c2_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
        lines!(ax, sol1.k, sol1.ck, linewidth=3, label="V", color=:blue, linestyle=linestyles[5])
        lines!(ax, sol2.k, sol2.ck, linewidth=3, label="MHNC", color=:darkred, linestyle=linestyles[6])
    end
    data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(ρval).txt")
    k_data = data[:, 1]
    S2_data = data[:, 2]
    δS2_data = data[:, 3]

    c2_data = c_from_s(S2_data, ρval) * ρval
    δc2_data = δc_from_s(S2_data, δS2_data, ρval) * ρval
    mask = .!isnan.(δc2_data)
    errorbars!(ax11, k_data[mask], c2_data[mask], δc2_data[mask], color=:black, whiskerwidth=6)
    scatter!(ax11, k_data[mask], c2_data[mask], color=:black, label="MC Data", markersize=6)
    scatter!(ax12, k_data[mask], c2_data[mask], color=:black, label="MC Data", markersize=10)
    errorbars!(ax12, k_data[mask], c2_data[mask], δc2_data[mask], color=:black, whiskerwidth=10, label="MC Data")


    axislegend(ax11, position=:rb, framevisible=false, patchsize=(30, 20), rowgap=-2, merge=true)

    for (ax, label) in zip([ax11, ax12], ["(a)", "(b)"])
        text!(ax,
            0, 1,
            text=label,
            align=(:left, :top),
            offset=(10, -10),
            space=:relative
        )
    end

    colgap!(fig.layout, 5)
    display(fig)
    save("Plots/ck_rho_0.94.pdf", fig)
    save("Plots/ck_rho_0.94.png", fig)
end



begin
    fig = Figure(size=(500, 500))
    ax1 = Axis(fig[1, 1:2], xlabel=L"k", ylabel=L"S(k)", limits=(0, 20, nothing, nothing))
    ax2 = Axis(fig[2, 1], xlabel=L"k", ylabel=L"S(k)", limits=(6.5, 7.5, 2, 3.3))
    ax3 = Axis(fig[2, 2], xlabel=L"k", ylabel=L"S(k)", limits=(12, 14, 1.2, 1.45), ylabelvisible=false)
    k = range(0.01, 20, length=1000)
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
        lines!(ax, k, S2_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
        lines!(ax, k, S2_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
        lines!(ax, k, S2_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
        lines!(ax, k, S2_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
        lines!(ax, k, S2_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
        lines!(ax, k, S2_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
        lines!(ax, sol1.k, sol1.Sk, linewidth=3, label="V", color=:blue, linestyle=linestyles[5])
        lines!(ax, sol2.k, sol2.Sk, linewidth=3, label="MHNC", color=:darkred, linestyle=linestyles[6])
    end
    scatter!(ax1, k_data, S2_data, color=:black, label="MC Data", markersize=6)
    scatter!(ax2, k_data, S2_data, color=:black, label="MC Data", markersize=10)
    scatter!(ax3, k_data, S2_data, color=:black, label="MC Data", markersize=10)
    text!(ax2,
        0.5, 0.0,
        text="first peak",
        align=(:center, :bottom),
        offset=(0, 10),
        space=:relative
    )

    text!(ax3,
        0.5, 0.0,
        text="second peak",
        align=(:center, :bottom),
        offset=(0, 10),
        space=:relative
    )

    for (ax, label) in zip([ax1, ax2, ax3], ["(a)", "(b)", "(c)"])
        text!(ax,
            0, 1,
            text=label,
            align=(:left, :top),
            offset=(10, -10),
            space=:relative
        )
    end


    # axislegend(ax, position = :rt)
    axislegend(ax1, position=:rt, nbanks=2, framevisible=false, patchsize=(30, 20), rowgap=-2)
    display(fig)

    save("Plots/Sk_rho_0.94.pdf", fig)
    save("Plots/Sk_rho_0.94.png", fig)
end