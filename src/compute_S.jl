#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=phys.default.q
#SBATCH --error=Errors/slurm-%j.err
#SBATCH --output=Errors/slurm-%j.out
#=
julia -t ${SLURM_CPUS_PER_TASK} $0
exit
# =#

import Pkg; Pkg.activate(".")


import SimulationAnalysis
# using CairoMakie
using DelimitedFiles
# using HDF5
using LoopVectorization
using Random





files = filter(x -> !(contains(x, "Equilibration")), readdir("Data"))


function _find_density_modes!(Reρ, Imρ, r, kspace)
    Ndim, N, N_timesteps = size(r)
    Nk = kspace.Nk
    k_array = kspace.k_array
    if Ndim == 3
        for t = 1:N_timesteps
            @tturbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                kz = k_array[3, i_k]
                Reρkt = 0.0
                Imρkt = 0.0
                for particle = 1:N 
                    rx = r[1, particle,  t]
                    ry = r[2, particle,  t]
                    rz = r[3, particle,  t]
                    kr = kx*rx + ky*ry + kz*rz
                    sinkr, coskr = sincos(kr)
                    Reρkt += coskr
                    Imρkt += sinkr
                end
                Reρ[t, i_k] = Reρkt
                Imρ[t, i_k] = Imρkt
            end
        end
    elseif Ndim == 2
        for t = 1:N_timesteps
            @turbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                Reρkt = 0.0
                Imρkt = 0.0
                for particle = 1:N 
                    rx = r[1, particle, t]
                    ry = r[2, particle, t]
                    kr = kx*rx + ky*ry
                    sinkr, coskr = sincos(kr)
                    Reρkt += coskr
                    Imρkt += sinkr
                end
                Reρ[t, i_k] = Reρkt
                Imρ[t, i_k] = Imρkt
            end
        end
    else
        throw(ArgumentError("Only 2D and 3D simulations are supported"))
    end
end


function find_density_modes_i(s::SimulationAnalysis.SingleComponentSimulation, kspace::SimulationAnalysis.KSpace, i; verbose=true)
    Ndim, N, _ = size(s.r_array)
    Nk = kspace.Nk
    N_timesteps = length(i)
    Reρ = zeros(length(i), Nk)
    Imρ = zeros(length(i), Nk)
    if verbose
        println("Calculating density modes for $N particles at $N_timesteps time points for $Nk wave vectors")
        println("Memory usage: $(Base.format_bytes(2*Base.summarysize(Reρ)))")
        println("Based on 10 GFLOPS, this will take approximately $(round(Nk*s.N*N_timesteps*9/10^10, digits=1)) seconds.")
    end
    tstart = time()
    
    @views _find_density_modes!(Reρ, Imρ, s.r_array[:, :, i], kspace)
    if verbose
        tstop = time()
        println("Elapsed time: $(round(tstop-tstart,digits=3)) seconds")
        println("Achieved GFLOPS: $(round(Nk*s.N*N_timesteps*9/(tstop-tstart)/10^9, digits=3))")
    end
    return SimulationAnalysis.SingleComponentDensityModes(Reρ, Imρ)
end

# ifile = parse(Int, ARGS[1])

for filename = shuffle(files)#[files[ifile]]
    ρ = split(filename, '_')[3]
    seed = split(filename, '_')[5][1:end-3]

    filename = "Data/HS3_rho_$(ρ)_seed_$(seed).h5"
    if !isfile("Processed_Data/Sk/Sk_rho_$(ρ)_seed_$(seed).txt")

        # create empty Sk file
        open("Processed_Data/Sk/Sk_rho_$(ρ)_seed_$(seed).txt", "w") do f
            writedlm(f, [0.0 0.0])
        end

        s = SimulationAnalysis.read_monodisperse_hard_sphere_simulation(filename; original=false, velocities=false, forcestype=false, dtarr=false)
        k_sample_array = 0.1:0.1:79.9

    
        Sk = zeros(length(k_sample_array)-1)
        @time for i = reverse(1:length(k_sample_array)-1)
            k1 = k_sample_array[i]
            k2 = k_sample_array[i+1]
            kmean = (k1 + k2) / 2

            dk = 2pi/s.box_sizes[1]

            estimated_kpoints = 2pi*kmean^2 * (k2-k1)/dk^3
            max_kpoints = 1000
            if estimated_kpoints > max_kpoints
                #set delta_k such that estimated kpoints = max_kpoints
                deltak = dk^3 * max_kpoints / (2pi*kmean^2) 
                k1 = kmean - deltak/2
                k2 = kmean + deltak/2
                @time Sk[i] = SimulationAnalysis.find_structure_factor(s; kmin=k1, kmax=k2, kfactor=1)
            else
                @time Sk[i] = SimulationAnalysis.find_structure_factor(s; kmin=k1, kmax=k2, kfactor=1)
            end
            flush(stdout)
        end
            

        k_centers = (k_sample_array[1:end-1] .+ k_sample_array[2:end]) ./ 2


        open("Processed_Data/Sk/Sk_rho_$(ρ)_seed_$(seed).txt", "w") do f
            writedlm(f, [k_centers Sk])
        end

    else
        @show filename
        println("File exists. Skipping...")
    end
    # fig = Figure(size = (800, 600))
    # ax = Axis(fig[1, 1], xlabel = "k", ylabel = "S(k)")
    # lines!(ax, k_sample_array, Sk, color = :black)
    # save("Plots/Sk_rho_$(ρ)_seed_$(seed).png", fig)
    # display(fig)

end