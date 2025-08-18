import Pkg
Pkg.activate(".")

using HDF5, Statistics, Dierckx

dir = "Processed_Data/S4"

files = readdir(dir, join=true)

rho_list = [0.8, 0.94]
k3_list = [2.0, 7.0]
costheta12_list = [cos(π/6), cos(π/4), cos(π/3),  cos(π/2), cos(2π/3), cos(5pi/6), cos(π)]
k1 = 7.2; k2 = 7.2
for rho in rho_list

    for k3 in k3_list, costheta12 in costheta12_list
        S4_list = []
        costheta13_array_list = []
        phi23_array_list = []
        for file in files
            if !occursin("rho_$(rho)_", file)
                continue
            end
            if !occursin("k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).h5", file)
                continue
            end
            if filesize(file) < 3000
                continue
            end


            h5open(file, "r") do file
                S4 = read(file, "S4")
                costheta13_array = read(file, "costheta13")
                phi23_array = read(file, "phi23")
                push!(S4_list, S4)
                push!(costheta13_array_list, costheta13_array)
                push!(phi23_array_list, phi23_array)
            end
        end
        println("Averaging S4 for rho = $(rho), k3 = $(k3), costheta12 = $(costheta12)")
        println("Number of files found: ", length(S4_list))

        S4 = mean(S4_list)
        S4_delta = [S4_list[i] - S4 for i in 1:length(S4_list)]
        S4_delta_squared = [(S4_delta[i]).^2 for i in 1:length(S4_list)]
        S4_std = sqrt.(sum(S4_delta_squared)/(length(S4_list) - 1))
        S4_std_of_mean = S4_std/sqrt(length(S4_list))
        costheta13_array = mean(costheta13_array_list)
        phi23_array = mean(phi23_array_list)
        @assert allequal(costheta13_array_list)
        @assert allequal(phi23_array_list)

        h5open("Processed_Data/S4/Mean/S4_rho_$(rho)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).h5", "w") do file
            write(file, "S4", S4)
            write(file, "S4_stderr", S4_std_of_mean)

            write(file, "costheta13", costheta13_array)
            write(file, "phi23", phi23_array)
        end

    end
end


# plot S4
using CairoMakie, DelimitedFiles

for rho in rho_list, k3 in k3_list, costheta12 in costheta12_list
    h5open("Processed_Data/S4/Mean/S4_rho_$(rho)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).h5", "r") do file
        S4 = read(file, "S4")
        costheta13_array = read(file, "costheta13")
        phi23_array = read(file, "phi23")

        fig = Figure()
        ax = Axis(fig[1, 1], ylabel="cos(theta13)", xlabel="phi23", title="S4 at $(rho), k3 = $(k3), theta12 = $(round(acos(costheta12)/π, digits=3))π")
        hm = heatmap!(ax, phi23_array, costheta13_array, S4', colormap=:bwr, colorrange=(-maximum(abs, S4), maximum(abs, S4)))
        Colorbar(fig[1, 2], hm, label="S4")
        display(fig)
    end
end

# for rho = 0.94, theta12 = 2π/3, k3 = 2.0, plot line for costheta13 = 0.0

begin 
    h5open("Processed_Data/S4/Mean/S4_rho_0.94_k1_7.2_k2_7.2_k3_7.0_costheta12_$(cos(2π/3)).h5", "r") do file
        S4 = read(file, "S4")
        costheta13_array = read(file, "costheta13")
        phi23_array = read(file, "phi23")

        spl = Spline2D(costheta13_array, phi23_array, S4)
        data = spl.(0.0, phi23_array)

        fig = Figure()
        ax = Axis(fig[1, 1], ylabel="S4", xlabel="phi23")
        itheta = findmin(abs.(costheta13_array .- 0.0))[2]
        lines!(ax, phi23_array, data, color=:black, label="costheta13 = $(round(costheta13_array[1], digits=3))")
        println("[")
        for i in 1:length(S4[itheta, :])
            println("$(phi23_array[i]) $(data[i])")
        end
        println("]")
        display(fig)
    end
end

# for rho = 0.94, theta12 = 2π/3, k3 = 7.0, plot line for costheta13 = 0.0
