#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=phys.default.q
#SBATCH --error=Errors/slurm-%j.err
#SBATCH --output=Errors/slurm-%j.out
#=
module load Julia/1.10.5
julia $0 $1
exit
# =#

cd("/vast.mnt/project/phys-glasscracking/Ilian/Hard_Spheres/")
import Pkg; Pkg.activate(".")

using Random
using StaticArrays, HDF5
using SimulationCode
dims = 3
kBT = 1.0 #parse(Float64, ARGS[1])
ρ = parse(Float64, ARGS[1])
interaction_potential = SimulationCode.HardSphere()

# System
Δt = 0.01                           # Time step
m = 1.0                             # mass
system = SimulationCode.Brownian(kBT, 0.0, 0.0, dims)


N = 30000                       # Number of particles
force_cutoff = 1.25             # Force cutoff
q_cutoff = 0.0                  # Size of the shell of the calculation of the Steinhardt order parameters
N_stepsMC = 10^6                # Number of MC steps to take
N_stepsMD = 0                # Number of time steps to take
swap_probability = 0.0          # Probability of choosing swap over displacement
max_MC_displacement = 0.1       # Maximal displacement in a displacement step in one direction
N_MD_equilibration_steps = 0   # Number of steps for short-time MD equilibration


random_seed = rand(1:10^9)      # Seed for the random number generator
box_size = (N/ρ)^(1/dims)          # Cubic box dimension
simulation_folder = "Data"      # Name of folder in which to store datafile
simulation_name = joinpath(simulation_folder, "HS$(dims)_rho_$(ρ)_seed_$(random_seed)")    # Name of the datafile
simulation_suffix = "_Equilibration.h5"
simulation_name_full = simulation_name*simulation_suffix     # Name of the datafile

# For neighbor lists
skin_distanceMC = 0.7           # Size of the verlet cells for swap MC
skin_distanceMD = 0.0          # Size of the verlet cells for MD

dump_info = SimulationCode.DumpInfo(
    true, #save
    simulation_name_full,
    [-1],#SimulationCode.create_when_to_save_array(N_stepsMD, 200), #when save
    true, #r
    false, #v
    false, #F
    false, #D
    false, #Epot      
)


########### Initialize structs and set random seed
cb = (arrays, parameters, output, neighborlist) -> nothing
parameters = SimulationCode.ParameterStruct(
    N_MD_equilibration_steps, random_seed, ρ, N, box_size, N_stepsMC, swap_probability, max_MC_displacement, 
    force_cutoff^2, q_cutoff, system, interaction_potential, dump_info, cb)
Random.seed!(random_seed)

arrays = SimulationCode.ArrayStruct(N, dims)
SimulationCode.generate_diameters!(arrays, parameters, interaction_potential)

output = SimulationCode.OutputStruct() 
println("Random seed = $random_seed") 

########### Long Equilibration ################

neighborlist = SimulationCode.initialize_neighbor_struct(skin_distanceMC, box_size, force_cutoff, N, arrays.D_array, dims; maxneighbours=100)
SimulationCode.find_random_initial_configuration!(arrays, parameters, output, neighborlist, steps=10^3)
println("\nStarting Long Equilibriation Procedure\n\n")
SimulationCode.perform_swap_monte_carlo!(arrays, parameters, output, neighborlist)

########### Production ################

simulation_suffix = ".h5"
simulation_name_full = simulation_name*simulation_suffix     # Name of the datafile
dump_info = SimulationCode.DumpInfo(
    true, #save
    simulation_name_full,
    0:1000:N_stepsMC,#SimulationCode.create_when_to_save_array(N_stepsMD, 200), #when save
    true, #r
    false, #v
    false, #F
    false, #D
    false, #Epot      
)
parameters = SimulationCode.ParameterStruct(
    N_MD_equilibration_steps, random_seed, ρ, N, box_size, N_stepsMC, swap_probability, max_MC_displacement, 
    force_cutoff^2, q_cutoff, system, interaction_potential, dump_info, cb)
output = SimulationCode.OutputStruct() 
SimulationCode.perform_swap_monte_carlo!(arrays, parameters, output, neighborlist)