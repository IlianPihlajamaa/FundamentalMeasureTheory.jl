#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
            @tturbo for i_k = 1:Nk
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

function compute_rhokt_single(r_array, kx, ky, kz, it)
    N = size(r_array, 2)
    Reρkt = 0.0
    Imρkt = 0.0
    @turbo for particle = 1:N 
        rx = r_array[1, particle, it]
        ry = r_array[2, particle, it]
        rz = r_array[3, particle, it]
        kr = kx*rx + ky*ry + kz*rz
        sinkr, coskr = sincos(kr)
        Reρkt += coskr
        Imρkt += sinkr
    end
    return Reρkt, Imρkt
end



function compute_S3(s, kspace, k1_indices, k2_indices, costhetamin, costhetamax, nmax)
    N_timesteps = length(s.t_array)
    S3 = 0.0
    count = 0

    k_array = kspace.k_array
    k_lengths = kspace.k_lengths
    r_array = s.r_array

    for i_k1 in k1_indices
        k1x =  k_array[1, i_k1]
        k1y =  k_array[2, i_k1]
        k1z =  k_array[3, i_k1]
        k1 = k_lengths[i_k1]
        if abs(k1) < 1e-6
            continue
        end
        if count >= nmax
            break
        end

        for i_k2 = k2_indices
            k2x =  k_array[1, i_k2]
            k2y =  k_array[2, i_k2]
            k2z =  k_array[3, i_k2]
            k2 = k_lengths[i_k2]
            if abs(k2) < 1e-6
                continue
            end
            cosθ = (k1x*k2x + k1y*k2y + k1z*k2z) / (k1*k2)

            if cosθ < costhetamin || cosθ > costhetamax
                continue
            end

            k3x = -k1x - k2x
            k3y = -k1y - k2y
            k3z = -k1z - k2z

            if k3x == 0.0 && k3y == 0.0 && k3z == 0.0
                continue
            end


            S3_new = 0.0
            for it = 1:N_timesteps
                # Reρ1 = ρkt.Re[it, i_k1]
                # Imρ1 = ρkt.Im[it, i_k1]
                # Reρ2 = ρkt.Re[it, i_k2]
                # Imρ2 = ρkt.Im[it, i_k2]
                # Reρ3 = ρkt.Re[it, i_k3]
                # Imρ3 = ρkt.Im[it, i_k3]
                Reρ1, Imρ1 = compute_rhokt_single(r_array, k1x, k1y, k1z, it)
                Reρ2, Imρ2 = compute_rhokt_single(r_array, k2x, k2y, k2z, it)
                Reρ3, Imρ3 = compute_rhokt_single(r_array, k3x, k3y, k3z, it)

                S3_new += Reρ1*Reρ2*Reρ3 - Reρ1*Imρ2*Imρ3 - Imρ1*Imρ2*Reρ3 - Imρ1*Reρ2*Imρ3  # Re((Reρ1 + im*Imρ1)*(Reρ2 + im*Imρ2)*(Reρ3 + im*Imρ3))
            end
            S3_new /= N_timesteps
            S3 += S3_new
            count += 1
            if count >= nmax
                break
            end
        end
    end
    return S3, count
end



function compute_S3_sampled(s, k_sample_arr1, k_sample_arr2, cosθ_sample_arr, dk, dcostheta; nmax=10)
    N = s.N

    @assert length(k_sample_arr1) == length(k_sample_arr2) == length(cosθ_sample_arr)
    N_samples = length(k_sample_arr1)


    S3 = zeros(N_samples)
    count_arr = zeros(Int, N_samples)
    kspace = @time SimulationAnalysis.construct_k_space(s, (0.0, maximum(k_sample_arr1) + maximum(k_sample_arr2)+ 2dk); kfactor=1, negative=true, rectangular=true)
    println("Memory use = $(Base.format_bytes(Base.summarysize(kspace)))")

    for (i, (k1, k2, costheta)) in enumerate(zip(k_sample_arr1, k_sample_arr2, cosθ_sample_arr))
        k1min = k1 - dk/2
        k1max = k1 + dk/2
    
        k1_indices = shuffle(findall(x -> (x > k1min && x <= k1max), kspace.k_lengths))
         
        println("k1 = $(round(k1, digits=4)), k2 = $(round(k2, digits=4)), costheta = $costheta, i = $i/$(length(k_sample_arr1))")
        k2min = k2 - dk/2
        k2max = k2 + dk/2

        k2_indices = shuffle(findall(x -> x > k2min && x <= k2max, kspace.k_lengths))

        costhetamin = costheta - dcostheta/2
        costhetamax = costheta + dcostheta/2
    
        S3_new, count = @time compute_S3(s, kspace, k1_indices, k2_indices, costhetamin, costhetamax, nmax)
        S3[i] = S3_new
        count_arr[i] = count
        flush(stdout)
    end


    S3 ./= count_arr*N

    return S3, count_arr
end

kmax = 20.0
for filename = shuffle(files)

    ρ = split(filename, '_')[3]
    seed = split(filename, '_')[5][1:end-3]

    if ρ == "0.8"
        continue
    end

    if isfile("Processed_Data/S3/S3_sweep$(1)_rho_$(ρ)_seed_$(seed).h5.h5")
        println("File exists. Skipping...")
        continue
    else
        h5open("Processed_Data/S3/S3_sweep$(1)_rho_$(ρ)_seed_$(seed).h5.h5", "w") do f
            write(f, "seed", seed)
        end
    end


    file = "Data/HS3_rho_$(ρ)_seed_$(seed).h5"
    @show file

    s = SimulationAnalysis.read_monodisperse_hard_sphere_simulation(file; original=false, velocities=false, forcestype=false, dtarr=false)
    println("Loaded simulation with $(s.N) particles and $(length(s.t_array)) time points")
    println("Memory use = $(Base.format_bytes(Base.summarysize(s)))")


    
    k_sample_arr1 = 0.0:0.1:5.0
    k_sample_arr2 = 5.0:0.2:10.0
    k_sample_arr3 = 10.0:0.4:20.0
    kcenters1 = k_sample_arr1[1:end-1] .+ diff(k_sample_arr1) / 2
    kcenters2 = k_sample_arr2[1:end-1] .+ diff(k_sample_arr2) / 2
    kcenters3 = k_sample_arr3[1:end-1] .+ diff(k_sample_arr3) / 2
    k_sample_arr = vcat(kcenters1, kcenters2, kcenters3)

    dk = 0.1
    dcostheta = 0.05

    # sweep k1 for constant k2
    k2 = 7.0

    sweepnumber = 1
    for costheta in [cos(π/6), cos(π/4), cos(π/3),  cos(π/2), cos(2π/3), cos(5pi/6), cos(π)]
        costheta_arr = [costheta for i in k_sample_arr]
        k1_arr = k_sample_arr
        k2_arr = [k2 for i in k_sample_arr] 
        
        S3, counts = compute_S3_sampled(s, k1_arr, k2_arr, costheta_arr, dk, dcostheta; nmax=100)

        h5open("Processed_Data/S3/S3_sweep$(sweepnumber)_rho_$(ρ)_seed_$(seed).h5", "w") do file
            write(file, "S3", S3)
            write(file, "counts", counts)
            write(file, "costheta_bin_centers", costheta_arr)
            write(file, "k_sample_arr1", k1_arr)
            write(file, "k_sample_arr2", k2_arr)
        end
        sweepnumber += 1
    end

    # sweep costheta for constant k1, k2
    costheta_arr = collect(-1:0.05:1)
    for k in [2.0, 3.0, 5.0, 7.0, 12.0, 14.0]
        k1_arr = [k for i in costheta_arr]
        k2_arr = [k for i in costheta_arr] 
        S3, counts = compute_S3_sampled(s, k1_arr, k2_arr, costheta_arr, dk, dcostheta; nmax=100)

        h5open("Processed_Data/S3/S3_sweep$(sweepnumber)_rho_$(ρ)_seed_$(seed).h5", "w") do file
            write(file, "S3", S3)
            write(file, "counts", counts)
            write(file, "costheta_bin_centers", costheta_arr)
            write(file, "k_sample_arr1", k1_arr)
            write(file, "k_sample_arr2", k2_arr)
        end
        sweepnumber += 1
    end

    # sweep k1, k2 for constant costheta

    for costheta in [cos(π/6), cos(π/4), cos(π/3),  cos(π/2), cos(2π/3), cos(5pi/6), cos(π)]
        costheta_arr = [costheta for i in k_sample_arr]
        k1_arr = k_sample_arr
        k2_arr = k_sample_arr
        
        S3, counts = compute_S3_sampled(s, k1_arr, k2_arr, costheta_arr, dk, dcostheta; nmax=100)

        h5open("Processed_Data/S3/S3_sweep$(sweepnumber)_rho_$(ρ)_seed_$(seed).h5", "w") do file
            write(file, "S3", S3)
            write(file, "counts", counts)
            write(file, "costheta_bin_centers", costheta_arr)
            write(file, "k_sample_arr1", k1_arr)
            write(file, "k_sample_arr2", k2_arr)
        end
        sweepnumber += 1
    end

end