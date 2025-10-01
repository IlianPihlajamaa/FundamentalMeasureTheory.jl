println("Starting the script")
using Symbolics, LinearAlgebra
using ColorSchemes
using DelimitedFiles
using CairoMakie
using LaTeXStrings
using HDF5
using Dierckx

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

direct_cor4_RF = compute_quadruplet_direct_correlation_function(Φ_RF, n, ω, homogeneity_substitutions)
direct_cor4_WB = compute_quadruplet_direct_correlation_function(Φ_WB, n, ω, homogeneity_substitutions)
direct_cor4_WB2 = compute_quadruplet_direct_correlation_function(Φ_WB2, n, ω, homogeneity_substitutions)
direct_cor4_M = compute_quadruplet_direct_correlation_function(Φ_M, n, ω, homogeneity_substitutions)
direct_cor4_KR = compute_quadruplet_direct_correlation_function(Φ_KR, n_KR, ω_KR, homogeneity_substitutions)
direct_cor4_T = compute_quadruplet_direct_correlation_function(Φ_T, n_T, ω_T, homogeneity_substitutions)
direct_cor4_WB2t = compute_quadruplet_direct_correlation_function(Φ_WB2t, n_T, ω_T, homogeneity_substitutions)
direct_cor4_L = compute_quadruplet_direct_correlation_function(Φ_L, n_T, ω_T, homogeneity_substitutions)
direct_cor4_G = compute_quadruplet_direct_correlation_function(Φ_G, n_T, ω_T, homogeneity_substitutions)

@assert direct_cor4_RF(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94^3 ≈ -253.09097694589278
@assert direct_cor4_WB(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94^3 ≈ -204.86292780419802
@assert direct_cor4_WB2(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94^3 ≈ -207.35891643000647
@assert direct_cor4_M(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94^3 ≈ -196.93451533988554
@assert direct_cor4_KR(2.2, 2.2, 2.2, π / 4, acos(0.025), 0.5, 0.94) * 0.94^3 ≈ -253.0909769458927

foo = (x...)->0.0


begin
    ρval = 0.94
    k1 = k2 = 7.2
    plot_idx = 1
    for cosθ12 = ([cos(π / 6), cos(1π / 3), cos(π / 4), cos(π / 2), cos(2π / 3), cos(5pi / 6), cos(π)])
        for cosθ13 = ([cos(π / 6), cos(1π / 3), cos(π / 4), cos(π / 2), cos(2π / 3), cos(5pi / 6), cos(π)])
            ϕ23_arr = range(0.01, π - 0.01, length=400)
            fig = Figure(size=(800, 400))

            for (j, k3) in enumerate([2.0, 7.0])

                f = h5open("Processed_Data/S4/Mean/S4_rho_$(ρval)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(cosθ12).h5", "r")

                costheta13_data = f["costheta13"][]
                ϕ23_data = f["phi23"][]
                S4_data = f["S4"][]
                S4_err = f["S4_stderr"][]

                close(f)

                if cosθ12 ≈ 1.0
                    θ12 = acos(1.0 - 0.025 / 2)
                elseif cosθ12 ≈ -1.0
                    θ12 = acos(-1.0 + 0.025 / 2)
                else
                    θ12 = acos(cosθ12)
                end

                if cosθ13 ≈ 1.0
                    θ13 = acos(1.0 - 0.025 / 2)
                elseif cosθ13 ≈ -1.0
                    θ13 = acos(-1.0 + 0.025 / 2)
                else
                    θ13 = acos(cosθ13)
                end

                @show k1, k2, k3, θ12, θ13
                ax = Axis(fig[1, j], ylabel=L"S_4", xlabel=L"\phi_{23}", title=L"S^{(4)}(k_1=7, k_2=7, k_3=%$(k3), \theta_{12}=%$(round(θ12/pi, digits=2))\pi, \theta_{13}=%$(round(acos(cosθ13)/pi, digits=2))\pi, \phi_{23})", titlesize=12)
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

                lines!(ax, ϕ23_arr, S4_RF, linewidth=3, label="Rosenfeld", color=colors[1], linestyle=linestyles[1])
                lines!(ax, ϕ23_arr, S4_WB, linewidth=3, label="White-Bear", color=colors[2], linestyle=linestyles[2])
                lines!(ax, ϕ23_arr, S4_WB2, linewidth=3, label="White-Bear II", color=colors[3], linestyle=linestyles[3])
                lines!(ax, ϕ23_arr, S4_M, linewidth=3, label="Malijevský", color=colors[4], linestyle=linestyles[4])
                lines!(ax, ϕ23_arr, S4_T, linewidth=3, label="Tarazona", color=colors[5], linestyle=linestyles[5])
                lines!(ax, ϕ23_arr, S4_WB2t, linewidth=3, label="White-Bear II (tensor)", color=colors[6], linestyle=linestyles[6])
                lines!(ax, ϕ23_arr, S4_L, linewidth=3, label="Lutsko", color=colors[7], linestyle=linestyles[7])
                lines!(ax, ϕ23_arr, S4_G, linewidth=3, label="Gül", color=colors[8], linestyle=linestyles[8])
                lines!(ax, ϕ23_arr, S4_conv2, linewidth=3, label="Convolution 2-body (RF)", linestyle=:dash, color=:black)

                spl = Spline2D(costheta13_data, ϕ23_data, S4_data)
                S4_data = [spl(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]
                spl_err = Spline2D(costheta13_data, ϕ23_data, S4_err)
                S4_err = [spl_err(cosθ13, ϕ_23) for ϕ_23 in ϕ23_data]

                scatter!(ax, ϕ23_data, S4_data, color=:black, label="MC Data", markersize=10)
                errorbars!(ax, ϕ23_data, S4_data, S4_err, color=:black, whiskerwidth=8)

                if j == 2
                    global plot_idx
                    ax3 = Axis(fig[1, 3], title="Plot $(plot_idx)")
                    hidedecorations!(ax3)
                    hidespines!(ax3)
                    Legend(fig[1, 3], ax, framevisible=false, patchsize=(30, 20), rowgap=-2)

                    plot_idx += 1
                end

            end
            display(fig)
            save("Plots/S4_2_7_costheta12_$(cosθ12)_costheta13_$(cosθ13).png", fig, px_per_unit=5)
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

    fig = Figure(size=(500, 750))
    ax11 = Axis(fig[1, 1], ylabel=L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax12 = Axis(fig[1, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax21 = Axis(fig[2, 1], ylabel=L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax22 = Axis(fig[2, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax31 = Axis(fig[3, 1], ylabel=L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xlabel=L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0, π, nothing, nothing))
    ax32 = Axis(fig[3, 2], xlabel=L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0, π, nothing, nothing))
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

            ϕ23_arr = range(0.01, π - 0.01, length=400)

            S4_RF = [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB = [get_S4(direct_cor2_WB, direct_cor3_WB, direct_cor4_WB, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB2 = [get_S4(direct_cor2_WB2, direct_cor3_WB2, direct_cor4_WB2, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_M = [get_S4(direct_cor2_M, direct_cor3_M, direct_cor4_M, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_T = [get_S4(direct_cor2_T, direct_cor3_T, direct_cor4_T, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_WB2t = [get_S4(direct_cor2_WB2t, direct_cor3_WB2t, direct_cor4_WB2t, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_L = [get_S4(direct_cor2_L, direct_cor3_L, direct_cor4_L, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_G = [get_S4(direct_cor2_G, direct_cor3_G, direct_cor4_G, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv = [get_S4(direct_cor2_RF, foo, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]

            lines!(ax, ϕ23_arr, S4_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, ϕ23_arr, S4_WB, linewidth=3, label="WBI", color=colors[2], linestyle=linestyles[2])
            lines!(ax, ϕ23_arr, S4_WB2, linewidth=3, label="WBII", color=colors[3], linestyle=linestyles[3])
            lines!(ax, ϕ23_arr, S4_M, linewidth=3, label="M", color=colors[4], linestyle=linestyles[4])
            lines!(ax, ϕ23_arr, S4_T, linewidth=3, label="T", color=colors[5], linestyle=linestyles[5])
            lines!(ax, ϕ23_arr, S4_WB2t, linewidth=3, label="WBIIt", color=colors[6], linestyle=linestyles[6])
            lines!(ax, ϕ23_arr, S4_L, linewidth=3, label="L", color=colors[7], linestyle=linestyles[7])
            lines!(ax, ϕ23_arr, S4_G, linewidth=3, label="G", color=colors[8], linestyle=linestyles[8])
            lines!(ax, ϕ23_arr, S4_conv, linewidth=3, label="conv", linestyle=:dash, color=:black)

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

            scatter!(ax, ϕ23_data, S4_data, color=:black, label="MC", markersize=10)
            errorbars!(ax, ϕ23_data, S4_data, S4_err, color=:black, whiskerwidth=8)

            int12 = round(Int, 3θ12 / pi)
            int12str = int12 == 1 ? "" : int12 == 2 ? "2" : string(int12)
            int13 = round(Int, 3θ13 / pi)
            int13str = int13 == 1 ? "" : int13 == 2 ? "2" : string(int13)

            if col == 1 && row == 2
                text!(ax, 1.0, 1.0, text=L"k_3 = %$(k3)", align=(:right, :top), offset=(-10, -27), space=:relative)
                text!(ax, 1.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:right, :top), offset=(-10, -44), space=:relative)
                text!(ax, 1.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:right, :top), offset=(-10, -61), space=:relative)
            else
                text!(ax, 0.0, 1.0, text=L"k_3 = %$(k3)", align=(:left, :top), offset=(10, -27), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:left, :top), offset=(10, -44), space=:relative)
                text!(ax, 0.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:left, :top), offset=(10, -61), space=:relative)
            end
        end
    end
    for (ax, label) in zip([ax11, ax12, ax21, ax22, ax31, ax32], ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"])
        text!(ax, 0.0, 1.0, text=label, align=(:left, :top), offset=(10, -10), space=:relative)
    end

    axislegend(ax12, framevisible=false, position=(0.75, 1.15), patchsize=(30, 20), rowgap=-5)

    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 10)


    display(fig)
    save("Plots/S4_all.pdf", fig, pt_per_unit=5)
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

    fig = Figure(size=(500, 500))
    ax11 = Axis(fig[1, 1], ylabel=L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax12 = Axis(fig[1, 2], xticks=0:0.25π:π, xticklabelsvisible=false, limits=(0, π, nothing, nothing))
    ax31 = Axis(fig[2, 1], ylabel=L"S_4(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)", xlabel=L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0, π, nothing, nothing))
    ax32 = Axis(fig[2, 2], xlabel=L"\phi_{23}", xticks=(0:0.25π:π, ["0", "π/4", "π/2", "3π/4", "π"]), limits=(0, π, nothing, nothing))
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

            ϕ23_arr = range(0.01, π - 0.01, length=400)

            S4_RF = [get_S4(direct_cor2_RF, direct_cor3_RF, direct_cor4_RF, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv3 = [get_S4(direct_cor2_RF, direct_cor3_RF, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]
            S4_conv2 = [get_S4(direct_cor2_RF, foo, foo, k1, k2, k3, θ12, θ13, ϕ23, ρval) for ϕ23 in ϕ23_arr]

            lines!(ax, ϕ23_arr, S4_RF, linewidth=3, label="R", color=colors[1], linestyle=linestyles[1])
            lines!(ax, ϕ23_arr, S4_conv3, linewidth=3, label="3-conv", color=:black, linestyle=:dot)
            lines!(ax, ϕ23_arr, S4_conv2, linewidth=3, label="2-conv", color=:black, linestyle=:dash)

            int12 = round(Int, 3θ12 / pi)
            int12str = int12 == 1 ? "" : int12 == 2 ? "2" : string(int12)
            int13 = round(Int, 3θ13 / pi)
            int13str = int13 == 1 ? "" : int13 == 2 ? "2" : string(int13)

            text!(ax, 0.0, 1.0, text=L"k_3 = %$(k3)", align=(:left, :top), offset=(10, -27), space=:relative)
            text!(ax, 0.0, 1.0, text=L"θ_{12} = %$(int12str)π/3", align=(:left, :top), offset=(10, -44), space=:relative)
            text!(ax, 0.0, 1.0, text=L"θ_{13} = %$(int13str)π/3", align=(:left, :top), offset=(10, -61), space=:relative)
        end
    end
    for (ax, label) in zip([ax11, ax12, ax31, ax32], ["(a)", "(b)", "(c)", "(d)"])
        text!(ax, 0.0, 1.0, text=label, align=(:left, :top), offset=(10, -10), space=:relative)
    end

    axislegend(ax12, framevisible=false, position=(0.8, 1.0), patchsize=(30, 20), rowgap=-2)

    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 10)


    display(fig)
    save("Plots/S4_convs.pdf", fig, pt_per_unit=5)
    save("Plots/S4_convs.png", fig)
end