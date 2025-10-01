println("Starting the script")


using Pkg;
Pkg.activate(".");

using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie, LaTeXStrings, Dierckx

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

direct_cor3_RF = compute_triplet_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor3_WB = compute_triplet_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor3_WB2 = compute_triplet_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor3_M = compute_triplet_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor3_KR = compute_triplet_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor3_T = compute_triplet_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor3_WB2t = compute_triplet_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor3_L = compute_triplet_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor3_G = compute_triplet_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@assert direct_cor3_RF(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -229.81290257878078
@assert direct_cor3_WB(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -193.7786362758227
@assert direct_cor3_WB2(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -193.78046452495985
@assert direct_cor3_M(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -190.14206661188794
@assert direct_cor3_KR(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -229.81290257878078
@assert direct_cor3_T(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -229.81290257878078
@assert direct_cor3_WB2t(0.1, 0.1, acos(-0.175), 0.94) * 0.94^2 ≈ -193.78046452495985

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
        cos_theta[cos_theta.≈1.0] .= 0.99
        cos_theta[cos_theta.≈-1.0] .= -0.99


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
        S3_L = @. get_S3(direct_cor2_L, direct_cor3_L, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)
        S3_G = @. get_S3(direct_cor2_G, direct_cor3_G, k1_interpl, k2_interpl, acos(cos_theta_interpl), rho)

        # find which type of sweep it is by seeing if costheta is constant
        f = Figure(size=(400, 400))

        if cos_theta[1] == cos_theta[end]
            #check if k2 is constant
            cosθ = round(cos_theta[1], digits=5)
            θdivπ = round(acos(cosθ) / π, digits=3)
            if k2[1] == k2[end]
                println("Sweep $(sweep) is a k1 sweep")
                ax = Axis(f[1, 1], xlabel=L"k_1", ylabel=L"S_3(k_1,k_2=%$(k2[1]),\,\theta=%$(θdivπ)\pi)")
            else
                println("Sweeping k1=k2=k")
                ax = Axis(f[1, 1], xlabel=L"k", ylabel=L"S_3(k_1=k_2=k,\,\theta=%$(θdivπ)\pi)")
            end
            lines!(ax, k1_interpl, S3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, k1_interpl, S3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, k1_interpl, S3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, k1_interpl, S3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, k1_interpl, S3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, k1_interpl, S3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, k1_interpl, S3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, k1_interpl, S3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
            scatter!(ax, k1, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, k1, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
        else
            #sweeping costheta at fixed k1 and k2
            println("Sweep $(sweep) is a costheta sweep")
            k1val = k1[1]
            k2val = k2[1]
            @assert k1val == k2val
            ax = Axis(f[1, 1], xlabel=L"\cos(\theta)", ylabel=L"S_.3(k_1=k_2=%$(k1val),\,\cos\,\theta)")

            lines!(ax, cos_theta_interpl, S3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, cos_theta_interpl, S3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, cos_theta_interpl, S3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, cos_theta_interpl, S3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, cos_theta_interpl, S3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, cos_theta_interpl, S3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, cos_theta_interpl, S3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, cos_theta_interpl, S3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
            scatter!(ax, cos_theta, S3_mean, label="MC Data", color=:black)
            errorbars!(ax, cos_theta, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
        end

        axislegend(ax, position=:rt, framevisible=false)

        display(f)
        save("Plots/S3_sweep_$(sweep)_rho_$(rho).pdf", f)
    end
end

function transpose_matrix(matrix)
    return [matrix[j, i] for i in 1:size(matrix, 2), j in 1:size(matrix, 1)]
end


begin
    fig = Figure(size=(800, 600))
    ax11 = Axis(fig[1, 1], xlabel=L"k_1", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(5, 15, 0, 14))
    ax12 = Axis(fig[1, 2], xlabel=L"k_1", limits=(0, 20, 0, 7))
    ax13 = Axis(fig[1, 3], xlabel=L"k_1", limits=(4, 14, 0, nothing))
    ax31 = Axis(fig[2, 1], xlabel=L"\cos(\theta)", ylabel=L"S_3(k_1,\, k_2,\,\theta)", limits=(nothing, nothing, nothing, 0.001))
    ax32 = Axis(fig[2, 2], xlabel=L"\cos(\theta)", limits=(-1, 1, 0, nothing))
    ax33 = Axis(fig[2, 3], xlabel=L"\cos(\theta)", limits=(-1, 1, 0, nothing))

    sweeps = [2 4 5;
        8 11 13]
    axs = [ax11 ax12 ax13;
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
            cos_theta[cos_theta.≈1.0] .= 0.99
            cos_theta[cos_theta.≈-1.0] .= -0.99

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
            j == 1 && i == 3 && @show S3_conv
            if i == 2
                #sweeping costheta at fixed k1 and k2
                k1val = k1[1]
                k2val = k2[1]
                @assert k1val == k2val
                lines!(ax, cos_theta_interpl, S3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, cos_theta_interpl, S3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, cos_theta_interpl, S3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, cos_theta_interpl, S3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, cos_theta_interpl, S3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, cos_theta_interpl, S3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, cos_theta_interpl, S3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, cos_theta_interpl, S3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
                lines!(ax, cos_theta_interpl, S3_conv, linewidth=3, label="conv", color=:black, linestyle=:dash)
                scatter!(ax, cos_theta, S3_mean, label="MC Data", color=:black)
                errorbars!(ax, cos_theta, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
                text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-50, -10), space=:relative)
            else
                # sweeping k1 at fixed k2 or k1 and k2
                lines!(ax, k1_interpl, S3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, k1_interpl, S3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, k1_interpl, S3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, k1_interpl, S3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, k1_interpl, S3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, k1_interpl, S3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, k1_interpl, S3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, k1_interpl, S3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
                lines!(ax, k1_interpl, S3_conv, linewidth=3, label="conv", color=:black, linestyle=:dash)

                scatter!(ax, k1, S3_mean, label="MC Data", color=:black)
                errorbars!(ax, k1, S3_mean, S3_std_err, color=:black, whiskerwidth=4)
                text!(ax, 1, 1, text=L"\theta=%$(round(acos(cos_theta[1])/π, digits=3))\pi", align=(:right, :top), offset=(-50, -10), space=:relative)
            end
            if i == 1
                text!(ax, 1, 1, text=L"k_2=%$(k2[1])", align=(:right, :top), offset=(-50, -30), space=:relative)
            end
        end
    end
    for (ax, label) in zip(transpose_matrix(axs)[:], ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"])
        text!(ax, 1, 1, text=label, align=(:right, :top), offset=(-10, -10), space=:relative)
    end
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    axislegend(ax13, position=:rc, framevisible=false, patchsize=(30, 20), rowgap=-5)


    display(fig)
    save("Plots/S3_all_sweeps.pdf", fig)
    save("Plots/S3_all_sweeps.png", fig)
end


begin
    fig = Figure(size=(800, 250), figure_padding=(5, 5, 5, 5))
    ax21 = Axis(fig[1, 1], xlabel=L"k", ylabel=L"\rho^2c_3(k_1,\, k_2,\,\theta)", limits=(0, 4, nothing, nothing))
    ax22 = Axis(fig[1, 2], xlabel=L"\cos(\theta)", limits=(-1, 1, -180, -70))
    ax23 = Axis(fig[1, 3], xlabel=L"\cos(\theta)", limits=(-1, 1, nothing, nothing))
    ax24 = Axis(fig[1, 4], xlabel=L"\cos(\theta)", limits=(-1, 1, -0.2, 0.6))

    sweeps = [15 8 11 13;
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
            cos_theta[cos_theta.≈1.0] .= 0.99
            cos_theta[cos_theta.≈-1.0] .= -0.99

            S3_mean = data[:, 4]
            S3_std_err = data[:, 5]

            S2data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")

            S2_k = S2data[:, 1]
            S2_mean = S2data[:, 2]
            mask = isfinite.(S2_mean)
            Sk_spl = Spline1D(S2_k[mask], S2_mean[mask], k=1)

            SSS = @. Sk_spl(k1) * Sk_spl(k2) * Sk_spl(sqrt(k1^2 + k2^2 + 2 * k1 * k2 * cos_theta))

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
            j == 1 && i == 3 && @show S3_conv
            if j == 3 || j == 2 || j == 4
                #sweeping costheta at fixed k1 and k2
                k1val = k1[1]
                k2val = k2[1]
                @assert k1val == k2val
                lines!(ax, cos_theta_interpl, C3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, cos_theta_interpl, C3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, cos_theta_interpl, C3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, cos_theta_interpl, C3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, cos_theta_interpl, C3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, cos_theta_interpl, C3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, cos_theta_interpl, C3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, cos_theta_interpl, C3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])

                scatter!(ax, cos_theta, ρ²C3, label="MC Data", color=:black)
                errorbars!(ax, cos_theta, ρ²C3, ρ²C3_err, color=:black, whiskerwidth=4)
            else
                # sweeping k1 at fixed k2 or k1 and k2
                lines!(ax, k1_interpl, C3_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
                lines!(ax, k1_interpl, C3_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
                lines!(ax, k1_interpl, C3_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
                lines!(ax, k1_interpl, C3_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
                lines!(ax, k1_interpl, C3_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
                lines!(ax, k1_interpl, C3_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
                lines!(ax, k1_interpl, C3_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
                lines!(ax, k1_interpl, C3_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
                scatter!(ax, k1[1:2:end], ρ²C3[1:2:end], label="MC Data", color=:black)
                errorbars!(ax, k1[1:2:end], ρ²C3[1:2:end], ρ²C3_err[1:2:end], color=:black, whiskerwidth=4)
            end
            if j == 1
                text!(ax, 0, 1, text=L"k_1=k_2=k", align=(:left, :top), offset=(10, -27), space=:relative)
                text!(ax, 0, 1, text=L"\theta=%$(round(acos(cos_theta[1])/π, digits=3))\pi", align=(:left, :top), offset=(10, -10), space=:relative)
                text!(ax, 1, 0, text="(a)", align=(:right, :bottom), offset=(-10, 10), space=:relative)

            elseif j == 2
                text!(ax, 0, 1, text=L"k_1=k_2=%$(k1val)", align=(:left, :top), offset=(10, -10), space=:relative)
                text!(ax, 1, 0, text="(b)", align=(:right, :bottom), offset=(-10, 10), space=:relative)
            elseif j == 3
                text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-40, -10), space=:relative)
                text!(ax, 1, 1, text="(c)", align=(:right, :top), offset=(-10, -10), space=:relative)
            elseif j == 4
                text!(ax, 1, 1, text=L"k_1=k_2=%$(k1val)", align=(:right, :top), offset=(-40, -10), space=:relative)
                text!(ax, 1, 1, text="(d)", align=(:right, :top), offset=(-10, -10), space=:relative)
            end
        end
    end

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)

    Legend(fig[1, 5], axs[1, 2], framevisible=false, width=85, patchsize=(30, 20))

    display(fig)
    save("Plots/C3_all_sweeps.pdf", fig)
    save("Plots/C3_all_sweeps.png", fig)

end