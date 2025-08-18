# FundamentalMeasureTheory.jl

This repository contains the code used to generate the data and plots for the paper "Fundamental measure theory for predicting many-body correlation functions". The code computes the 2, 3, and 4-body direct correlation functions and structure factors for a hard-sphere fluid and compares them to various flavors of Fundamental Measure Theory (FMT).

## Repository Structure

- `src/`: Contains all the Julia scripts.
  - `simulationcode.jl`: Runs the Monte Carlo simulations.
  - `compute_*.jl`: Scripts for analyzing the raw simulation data.
  - `average_*.jl`: Scripts for averaging the results from multiple simulations.
  - `FMT_*.jl`: Scripts for comparing the averaged data with FMT and generating plots.
- `Processed_Data/`: Contains the processed and averaged data.
- `Plots/`: Contains the final plots.
- `Project.toml` and `Manifest.toml`: Julia package management files.

## Reproduction Steps

To reproduce the results, follow these steps:

### 1. Install Dependencies

The required Julia packages are listed in the `Project.toml` and `Manifest.toml` files. To install them, open a Julia REPL in the repository's root directory and run:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install the correct versions of all dependencies, including `SimulationCode.jl` from `https://github.com/IlianPihlajamaa/SimulationCode.jl` and `SimulationAnalysis.jl` from ``https://github.com/IlianPihlajamaa/SimulationAnalysis.jl``

### 2. Data Generation (Optional)

The raw data from the Monte Carlo simulations is not included in this repository due to its large size. You can generate it by running the `src/simulationcode.jl` script.

The script takes the density `ρ` as a command-line argument. For example:

```bash
julia src/simulationcode.jl 0.94
```

This will generate a data file in a `Data/` directory (which you may need to create). The paper uses 100 datasets for each density.

### 3. Data Analysis

The analysis is performed in three stages: computation, averaging, and analysis/plotting. The scripts should be run in this order.

#### a. Computation

These scripts process the raw data from the `Data/` directory and save the results in the `Processed_Data/` directory.

```bash
julia src/compute_S.jl
julia src/compute_S3.jl
julia src/compute_S4.jl
```

#### b. Averaging

These scripts average the results from the computation step over all the random seeds.

```bash
julia src/average_Sk.jl
julia src/average_S3.jl
julia src/average_S4.jl
julia src/compute_c3.jl # this uses the mean Sk data
```

#### c. Comparison with FMT and Plotting

These scripts compare the averaged data with the predictions from Fundamental Measure Theory and generate the final plots, which are saved in the `Plots/` directory.

```bash
julia src/FMT_S2.jl
julia src/FMT_S3.jl
julia src/FMT_S4.jl
```

After running these scripts, the `Plots/` directory will contain the figures from the paper and a some that arre not included.

## Citation

If you use this code in your research, please add a citation to the original paper. 