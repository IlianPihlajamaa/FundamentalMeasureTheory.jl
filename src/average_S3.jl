import Pkg; Pkg.activate(".")


using DelimitedFiles
using CairoMakie, HDF5, Statistics


files = filter(x -> !(contains(x, "Equilibration")), readdir("Data"))
Nsweeps = 20


for rho in [0.8, 0.94]
    
    files_rho = filter(x -> contains(x, "rho_$(rho)"), files)
    for sweep in 1:Nsweeps
        S3_total = []
        k1 = []
        k2 = []
        cos_theta = []

        for file in files_rho
            seed = parse(Int, split(split(file, "seed_")[2], ".")[1])
            savefile = "Processed_Data/S3/S3_sweep$(sweep)_rho_$(rho)_seed_$(seed).h5"
            if isfile(savefile)
                file = h5open(savefile, "r")

                S3 = read(file["S3"])
                k1_new = read(file["k_sample_arr1"])
                k2_new = read(file["k_sample_arr2"])
                cos_theta_new = read(file["costheta_bin_centers"])
                
                push!(S3_total, S3)
                push!(k1, k1_new)
                push!(k2, k2_new)
                push!(cos_theta, cos_theta_new)

                # println("Loaded $(savefile)")
            else
                # println("File not found: $savefile")
            end
        end
        if length(S3_total) == 0
            println("No files found for rho $(rho) and sweep $(sweep)")
            continue
        end
        S3_total_mean = sum(S3_total) / length(S3_total)

        S3_std_err = (std(stack(S3_total), dims=2) / sqrt(length(S3_total)))[:]

        @assert allequal(k1)
        @assert allequal(k2)
        @assert allequal(cos_theta)
        k1_mean = k1[1]
        k2_mean = k2[1]
        cos_theta_mean = cos_theta[1]
        println("averaged S3 for ρ = $(rho) for $(length(S3_total)) seeds")
        
        open("Processed_Data/S3/Mean/S3_sweep$(sweep)_rho_$(rho).txt", "w") do io
            writedlm(io, [k1_mean k2_mean cos_theta_mean S3_total_mean S3_std_err])
        end
    end
end

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
        S3_mean = data[:, 4]
        S3_std_err = data[:, 5]

        # find which type of sweep it is by seeing if costheta is constant
        if cos_theta[1] == cos_theta[end]

            #check if k2 is constant
            cosθ = cos_theta[1]
            if k2[1] == k2[end]
                println("Sweep $(sweep) is a k1 sweep")
                f = Figure(size=(1000, 1000))
                ax = Axis(f[1, 1], xlabel="k1", ylabel="S3(k1,k2,costheta)", title="S3 for rho $(rho) and costheta=$(cosθ)")
                scatterlines!(ax, k1, S3_mean, label="S3 mean")
                errorbars!(ax, k1, S3_mean, S3_std_err, label="S3 std err")
                display(f)
            else
                println("Sweeping k1=k2=k")
                f = Figure(size=(1000, 1000))
                ax = Axis(f[1, 1], xlabel="k", ylabel="S3(k1=k2=k, costheta)", title="S3 for rho $(rho) and costheta=$(cosθ)")
                scatterlines!(ax, k1, S3_mean, label="S3 mean")
                errorbars!(ax, k1, S3_mean, S3_std_err, label="S3 std err")
                display(f)
            end
        else
            #sweeping costheta at fixed k1 and k2
            println("Sweep $(sweep) is a costheta sweep")
            f = Figure(size=(1000, 1000))
            ax = Axis(f[1, 1], xlabel="costheta", ylabel="S3(k1,k2,costheta)", title="S3 for rho $(rho) and k1=$(k1[1]) and k2=$(k2[1])")
            scatterlines!(ax, cos_theta, S3_mean, label="S3 mean")
            errorbars!(ax, cos_theta, S3_mean, S3_std_err, label="S3 std err")
            display(f)
        end

    end
end

