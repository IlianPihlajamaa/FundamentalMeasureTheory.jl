#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=220G
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
using Base.Threads

files = filter(x -> !(contains(x, "Equilibration")), readdir("Data"))



function _find_density_modes!(Reρ, Imρ, r, kspace)
    Ndim, N, N_timesteps = size(r)
    Nk = kspace.Nk
    k_array = kspace.k_array
    t1 = time()
    count = Atomic{Int}(0)
    if Ndim == 3
        @threads for t = 1:N_timesteps

            @turbo for i_k = 1:Nk
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
            count[] += 1

            if threadid() == 1
                println("Done  = $(count[])/$(N_timesteps). Elapsed time: $(round(time()-t1, digits=3)) seconds")
                flush(stdout)
            end
        end
    elseif Ndim == 2
        @tturbo for t = 1:N_timesteps
             for i_k = 1:Nk
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
        println("Calculating density modes for $N particles at $N_timesteps time points for $Nk wave vectors with $(nthreads()) threads")
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
    flush(stdout)
    return SimulationAnalysis.SingleComponentDensityModes(Reρ, Imρ)
end


function compute_S3(ρkt, kspace, k_sample_arr, cosθ_sample_arr, L, N; nmax=10)
    N_timesteps, Nk = size(ρkt.Re)
    N_kbins = length(k_sample_arr)
    N_cosθbins = length(cosθ_sample_arr)
    cosθ_binwidth = cosθ_sample_arr[2] - cosθ_sample_arr[1]

    S3 = zeros(N_cosθbins-1, N_kbins-1, N_kbins-1)
    count_arr = zeros(Int, N_cosθbins-1, N_kbins-1, N_kbins-1)
    k_binwidth = k_sample_arr[2] - k_sample_arr[1]
    k_array = kspace.k_array
    k_lengths = kspace.k_lengths
    cart2lin = kspace.cartesian_to_linear
    dk = 2π/L * kspace.kfactor
    idone = 0
    for i_k1 = shuffle(1:Nk)
        idone += 1
        if idone % 100000 == 0
            println("i_k1 = $idone / $Nk, computed = $(sum(count_arr)), max = $(prod(size(count_arr))*nmax)")
            flush(stdout)
        end

        k1x =  k_array[1, i_k1]
        k1y =  k_array[2, i_k1]
        k1z =  k_array[3, i_k1]
        k1 = k_lengths[i_k1]
        if k1 == 0.0
            continue
        end
        if k1 > k_sample_arr[end]
            continue
        end
        k1_index = ceil(Int, k1 / k_binwidth)
        if all(count_arr[:, k1_index, :] .>= nmax)
            continue
        end
        for i_k2 = 1:Nk
            k2x =  k_array[1, i_k2]
            k2y =  k_array[2, i_k2]
            k2z =  k_array[3, i_k2]
            k2 = k_lengths[i_k2]
            if k2 == 0.0
                continue
            end
            if k2 > k_sample_arr[end]
                continue
            end
            k2_index = ceil(Int, k2 / k_binwidth)

            cosθ = (k1x*k2x + k1y*k2y + k1z*k2z) / (k1*k2)

            θ_index = ceil(Int, (cosθ + 1.0) / cosθ_binwidth)
            if θ_index == 0
                θ_index = 1
            end
            if count_arr[θ_index, k1_index, k2_index] >= nmax
                continue
            end

            i_k3 = cart2lin[round(Int, -(k1x+k2x)/dk), round(Int, -(k1y+k2y)/dk), round(Int, -(k1z+k2z)/dk)]
            k3 = k_lengths[i_k3]
            if k3 == 0.0
                continue
            end


            S3_new = 0.0
            for it = 1:N_timesteps
                Reρ1 = ρkt.Re[it, i_k1]
                Imρ1 = ρkt.Im[it, i_k1]
                Reρ2 = ρkt.Re[it, i_k2]
                Imρ2 = ρkt.Im[it, i_k2]
                Reρ3 = ρkt.Re[it, i_k3]
                Imρ3 = ρkt.Im[it, i_k3]

                S3_new += Reρ1*Reρ2*Reρ3 - Reρ1*Imρ2*Imρ3 - Imρ1*Imρ2*Reρ3 - Imρ1*Reρ2*Imρ3  # Re((Reρ1 + im*Imρ1)*(Reρ2 + im*Imρ2)*(Reρ3 + im*Imρ3))
            end
            S3_new /= N_timesteps
            S3[θ_index, k1_index, k2_index] += S3_new
            count_arr[θ_index, k1_index, k2_index] += 1
        end
    end

    S3 ./= count_arr*N

    return S3
end




for filename = shuffle(files)
    @show filename
    ρ = split(filename, '_')[3]

    seed = split(filename, '_')[5][1:end-3]

    if isfile("Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(7.2)_k2_$(7.2)_k3_$(7.0)_costheta12_$(cos(π)).h5")
        println("File exists. Skipping...")
        continue
    end
    h5open("Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(7.2)_k2_$(7.2)_k3_$(7.0)_costheta12_$(cos(π)).h5", "w") do f
        write(f, "seed", seed)
    end
    file = "Data/HS3_rho_$(ρ)_seed_$(seed).h5"

    s = SimulationAnalysis.read_monodisperse_hard_sphere_simulation(file; original=false, velocities=false, forcestype=false, dtarr=false)
    kspace = SimulationAnalysis.construct_k_space(s, (0.0, 7.2*3+0.15); kfactor=1, negative=true, rectangular=true)
    indices = 1:length(s.t_array)

    ρkt = find_density_modes_i(s, kspace, indices; verbose=true)

    
    k1 = 7.2; k2 = 7.2; 
    for k3 = [2.0, 7.0]
        for costheta12 = [cos(π/6), cos(1π/3), cos(π/4), cos(π/2), cos(2π/3), cos(5pi/6), cos(π)]
            Ntheta = 24; Nphi = 40
            costheta13_array, phi23_array, S4, S4conv = @time SimulationAnalysis.find_S4_offiagonal(
                s, kspace, ρkt, Ntheta, Nphi; 
                q1=k1, dq1=0.05, q2=k2, dq2=0.05, 
                costheta12=costheta12, dcostheta12=0.025, 
                q3=k3, dq3=0.05, dcostheta13=0.025, 
                dphi23=0.05, maxsamples=10^9
                )
            flush(stdout)
            savefile = "Processed_Data/S4/S4_rho_$(ρ)_seed_$(seed)_k1_$(k1)_k2_$(k2)_k3_$(k3)_costheta12_$(costheta12).h5"

            h5open(savefile, "w") do f
                write(f, "S4", S4)
                write(f, "costheta13", costheta13_array)
                write(f, "phi23", phi23_array)
            end
        end
    end


    # k_sample_arr = 0.0:0.2:10.0
    # cosθ_sample_arr = -1.0:0.05:1.0
    # S3 = @time compute_S3(ρkt, kspace, k_sample_arr, cosθ_sample_arr, s.box_sizes[1], s.N; nmax = 100000)
    # cosθ_bin_centers = (cosθ_sample_arr[1:end-1] + cosθ_sample_arr[2:end]) / 2

    # h5open(S3savefile, "w") do f
    #     for itheta = 1:length(cosθ_bin_centers)
    #         write(f, "S3_$(itheta)", S3[itheta, :, :])
    #     end
    #     write(f, "costheta_bin_centers", collect(cosθ_bin_centers))
    #     write(f, "k_sample_arr", collect(k_sample_arr))
    # end
    run(`sbatch compute_S34.jl`)
    exit()
end