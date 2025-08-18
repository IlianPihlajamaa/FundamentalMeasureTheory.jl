#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=phys.bigmem.q
#SBATCH --error=Errors/slurm-%j.err
#SBATCH --output=Errors/slurm-%j.out
#=
srun hostname
julia -t ${SLURM_CPUS_PER_TASK} $0
exit
# =#

import Pkg;
Pkg.activate(".");

import SimulationAnalysis
using DelimitedFiles
using HDF5
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




for filename = shuffle(files)
    ρ = split(filename, '_')[3]
    if ρ == "0.8"
        continue
    end
    seed = split(filename, '_')[5][1:end-3]

    file = "Data/HS3_rho_$(ρ)_seed_$(seed).h5"
    k1 = 7.2; k2 = 7.2; k3 = 7.2; costheta12 = cos(1π/3)
    if isfile("Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).txt")
        println("File exists. Skipping...")
        continue
    end
    open("Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).txt", "w") do f
        writedlm(f, [])
    end

    s = SimulationAnalysis.read_monodisperse_hard_sphere_simulation(file; original=false, velocities=false, forcestype=false, dtarr=false)
    kspace = SimulationAnalysis.construct_k_space(s, (0.0, 7.2*3+0.15); kfactor=1, negative=true, rectangular=true)
    indices = 1:length(s.t_array)
    # savefilename = "Processed Data/Sk/temp/Sk_rho_$(ρ)_seed_$(seed)_batch_$(i).txt"
    # if isfile(savefilename)
    #     print("File exists. Skipping...")
    #     continue
    # end
    ρkt = find_density_modes_i(s, kspace, indices; verbose=true)

    Ntheta = 24; Nphi = 40
     costheta13_array, phi23_array, S4, S4conv = SimulationAnalysis.find_S4_offiagonal(
         s, kspace, ρkt, Ntheta, Nphi; 
         q1=k1, dq1=0.05, q2=k2, dq2=0.05, 
         costheta12=costheta12, dcostheta12=0.025, 
         q3=k3, dq3=0.05, dcostheta13=0.025, 
         dphi23=0.05, maxsamples=10^9
         )

    savefile = "Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).txt"
    open(savefile, "w") do io
        writedlm(io, S4)
    end
    open("Processed_Data/S4/costheta13.txt", "w") do io
        writedlm(io, costheta13_array)
    end
    open("Processed_Data/S4/phi23.txt", "w") do io
        writedlm(io, phi23_array)
    end

    k1 = 7.2; k2 = 7.2; k3 = 2.0; costheta12 = cos(2π/3)
    costheta13_array, phi23_array, S4, S4conv = SimulationAnalysis.find_S4_offiagonal(
        s, kspace, ρkt, Ntheta, Nphi;
        q1=k1, dq1=0.05, q2=k2, dq2=0.05, 
        costheta12=costheta12, dcostheta12=0.025, 
        q3=k3, dq3=0.05, dcostheta13=0.025, 
        dphi23=0.05, maxsamples=10^9
        )
    
    savefile = "Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).txt"
    open(savefile, "w") do io
        writedlm(io, S4)
    end
end