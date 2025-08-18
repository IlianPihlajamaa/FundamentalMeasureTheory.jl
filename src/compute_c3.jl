import Pkg
Pkg.activate(".")

using HDF5, Dierckx, DelimitedFiles, CairoMakie

rho_list = [0.8, 0.94]

for rho in rho_list
    Sk_data = readdlm("Processed_Data/Sk/mean/Sk_rho_$(rho).txt")
    k_Sk = Sk_data[:, 1]
    Sk = Sk_data[:, 2]
    mask = .!isnan.(Sk)
    k_Sk = k_Sk[mask]
    Sk = Sk[mask]
    @show extrema(k_Sk)
    S_interp = Spline1D(k_Sk, Sk, k=3)


    costheta_bin_centers, k_bin_centers, S3, S3_std = h5open("Processed_Data/S3/Mean/S3_rho_$rho.h5") do f
        k_bin_centers = read(f, "k_sample_arr")
        costheta_bin_centers = read(f, "costheta_bin_centers")
        S3 =  read(f, "S3")
        S3_std = read(f, "S3_std")
        (costheta_bin_centers, k_bin_centers, S3, S3_std)
    end
    ρ²C3 = zeros(length(costheta_bin_centers), length(k_bin_centers), length(k_bin_centers))
    ρ²C3_std = zeros(length(costheta_bin_centers), length(k_bin_centers), length(k_bin_centers))

    for (itheta, costheta) = enumerate(costheta_bin_centers)
        for (ik1, k1) = enumerate(k_bin_centers)
            for (ik2, k2) = enumerate(k_bin_centers)
                k3 = sqrt(k1^2 + k2^2 + 2k1*k2*costheta)
                S3_conv = S_interp(k3)*S_interp(k1)*S_interp(k2)
                ρ²C3[itheta, ik1, ik2] = S3[itheta, ik1, ik2] / S3_conv - 1
                ρ²C3_std[itheta, ik1, ik2] = S3_std[itheta, ik1, ik2]^2 / S3_conv
            end
        end
    end

    # plot ρ²C3

    for itheta = 1:length(costheta_bin_centers)
        fig = Figure(size = (1200, 600))
        ax = Axis(fig[1, 1], xlabel="k1", ylabel="k2", title="ρ²C3 at $(rho), costheta = $(costheta_bin_centers[itheta])", aspect=1)
        heatmap!(ax, k_bin_centers, k_bin_centers, ρ²C3[itheta, :, :])
        ax2 = Axis(fig[1, 2], xlabel="k", ylabel="ρ²C3(k,k,theta)", title="ρ²C3 at $(rho), costheta = $(costheta_bin_centers[itheta])", limits=(0,10,nothing, nothing))
        scatter!(ax2, k_bin_centers, [ρ²C3[itheta, ik, ik] for ik in 1:length(k_bin_centers)], label="ρ²C3(k,k,theta)")
        errorbars!(ax2, k_bin_centers, [ρ²C3[itheta, ik, ik] for ik in 1:length(k_bin_centers)], [ρ²C3_std[itheta, ik, ik] for ik in 1:length(k_bin_centers)], label="ρ²C3(k,k,theta) std", whiskerwidth = 10)
        ax3 = Axis(fig[1, 3], xlabel="k", ylabel="ρ²C3(k,k,theta)", title="ρ²C3 at $(rho), costheta = $(costheta_bin_centers[itheta])", limits=(0,10,-2,5))
        scatter!(ax3, k_bin_centers, [ρ²C3[itheta, ik, ik] for ik in 1:length(k_bin_centers)], label="ρ²C3(k,k,theta)")
        errorbars!(ax3, k_bin_centers, [ρ²C3[itheta, ik, ik] for ik in 1:length(k_bin_centers)], [ρ²C3_std[itheta, ik, ik] for ik in 1:length(k_bin_centers)], label="ρ²C3(k,k,theta) std", whiskerwidth = 10)
        display(fig)

    end

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel="costheta", ylabel="S3(k,k,costheta)", title="ρ²C3 at $(rho)")
    k = 3.0
    ik = argmin(abs.(k_bin_centers .- k))
    scatterlines!(ax, costheta_bin_centers, ρ²C3[:, ik, ik], label="k = $(k_bin_centers[ik])")
    display(fig)

end        

