# Stochastic Optimization and Decomposition Methods Teaching Repository

This repository contains teaching materials and Julia implementations for decomposition methods in stochastic optimization, prepared for the course **Optimization Techniques** at **Universidad Pontificia Comillas**.

Prepared by **Diego Tejada**.

> **Disclaimer**  
> Part of code was created with the help of **GitHub Copilot**. The markdown notes and Julia code are original content created for this course, but they were generated with the assistance of AI tools to enhance clarity and completeness.

## What this repository is for

The goal of this repository is to help students understand, model, and solve:

1. The deterministic Fixed-Charge Transportation Problem (FCTP).
2. Its stochastic extension (SFCTP).
3. Advanced exact methods used in large-scale optimization:
   1. Benders Decomposition (single-cut and multicut).
   2. Column Generation.
   3. Branch-and-Price.

The repository combines:

1. Working Julia code with JuMP + HiGHS.
2. Didactic markdown notes with equations, diagrams, and interpretation.

## Repository contents (quick guide)

### Julia models and algorithms

- [FCTP.jl](FCTP.jl): Deterministic fixed-charge transportation model and scenario-by-scenario baseline solves.
- [SFCTP.jl](SFCTP.jl): Stochastic fixed-charge transportation model, stochastic metrics, and out-of-sample evaluation.
- [SFCTP-BD.jl](SFCTP-BD.jl): Benders decomposition (single-cut) for SFCTP.
- [SFCTP-BD-multicut.jl](SFCTP-BD-multicut.jl): Benders decomposition with multicut extension.
- [SFCTP-CG.jl](SFCTP-CG.jl): Column generation approach for SFCTP.
- [SFCTP-CG-B&P.jl](SFCTP-CG-B&P.jl): Branch-and-price implementation for SFCTP.

### Documentation files

- [intro_stochastic_optimization_fctp.md](intro_stochastic_optimization_fctp.md): Short introduction to stochastic optimization and deterministic vs stochastic FCTP, including standard stochastic metrics.
- [benders_decomposition_and_multicut_sfctp.md](benders_decomposition_and_multicut_sfctp.md): Practical guide to Benders decomposition and multicut acceleration for SFCTP.
- [column_generation_and_branch_and_price.md](column_generation_and_branch_and_price.md): Theory and implementation walkthrough for column generation and branch-and-price on SFCTP.

## The example used across the repository

The running example is a transportation network with:

1. Origins (supply nodes) with capacities.
2. Destinations (demand nodes).
3. Transportation arcs between origins and destinations.

Each arc has:

1. A fixed investment cost (pay once if the arc is built).
2. A variable per-unit transportation cost.

In the stochastic version, demand is uncertain and represented by scenarios with probabilities.

## Basic equations of the Fixed-Charge Transportation model

### Deterministic FCTP

Decision variables:

- $Y_{ij} \in \{0,1\}$: 1 if arc $(i,j)$ is built.
- $X_{ij} \ge 0$: flow sent on arc $(i,j)$.

Objective:

$$
\min \sum_{i,j} F_{ij} Y_{ij} + \sum_{i,j} C_{ij} X_{ij}
$$

Constraints:

$$
\sum_j X_{ij} \le A_i \quad \forall i
$$

$$
\sum_i X_{ij} \ge B_j \quad \forall j
$$

$$
X_{ij} \le U_{ij} Y_{ij} \quad \forall i,j
$$

where typically $U_{ij} = \min(A_i, B_j)$.

### Stochastic extension (SFCTP, two-stage)

First stage:

- Choose $Y_{ij}$ before demand realization.

Second stage (for each scenario $s$):

- Choose $X_{ij}^s$ and optional shortage $Z_j^s$.

Objective:

$$
\min \sum_{i,j} F_{ij} Y_{ij} + \sum_s p_s \left(\sum_{i,j} C_{ij} X_{ij}^s + \kappa \sum_j Z_j^s\right)
$$

Key linking constraint:

$$
X_{ij}^s \le U_{sij} Y_{ij} \quad \forall s,i,j
$$

This non-anticipative structure (same $Y$ for all scenarios) is the reason decomposition methods are effective in this repository.

## Suggested reading order

1. Start with [intro_stochastic_optimization_fctp.md](intro_stochastic_optimization_fctp.md).
2. Continue with [benders_decomposition_and_multicut_sfctp.md](benders_decomposition_and_multicut_sfctp.md).
3. Finish with [column_generation_and_branch_and_price.md](column_generation_and_branch_and_price.md).

Then run the Julia scripts to connect theory and implementation.

## Software stack

- Julia
- JuMP
- HiGHS
- DataFrames and Plots (for reporting and visualizations)

Check the [Project.toml](Project.toml) file for more information.

## How to use this repository

These steps are intended for students who want to download the repository, install the required Julia packages, and run the example files.

### 1. Clone the repository

Open a terminal and run:

```bash
git clone <repository-url>
cd stochastic-optimization-decomposition-methods-101
```

If you downloaded the repository as a ZIP file instead of cloning it, just extract it and open a terminal in the project folder.

### 2. Open Julia in the project folder

From the repository root, start Julia:

```bash
julia
```

### 3. Instantiate the project environment

Inside the Julia REPL, run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install the packages listed in [Project.toml](Project.toml) and create the [Manifest.toml](Manifest.toml) file.

### 4. Run any of the example scripts

Still inside Julia, you can run the scripts with:

```julia
include("FCTP.jl")
include("SFCTP.jl")
include("SFCTP-BD.jl")
include("SFCTP-BD-multicut.jl")
include("SFCTP-CG.jl")
include("SFCTP-CG-B&P.jl")
```

You can also run a file directly from the terminal:

```bash
julia --project=. FCTP.jl
julia --project=. SFCTP.jl
julia --project=. SFCTP-BD.jl
julia --project=. SFCTP-BD-multicut.jl
julia --project=. SFCTP-CG.jl
julia --project=. "SFCTP-CG-B&P.jl"
```

### 5. Suggested execution order for students

To follow the course material progressively, run the files in this order:

1. [FCTP.jl](FCTP.jl)
2. [SFCTP.jl](SFCTP.jl)
3. [SFCTP-BD.jl](SFCTP-BD.jl)
4. [SFCTP-BD-multicut.jl](SFCTP-BD-multicut.jl)
5. [SFCTP-CG.jl](SFCTP-CG.jl)
6. [SFCTP-CG-B&P.jl](SFCTP-CG-B&P.jl)

### 6. Notes for students

1. Some scripts generate plots, so running them in VS Code or a Julia environment with graphics support is recommended.
2. The first package installation may take a few minutes.
3. If package installation fails, run `Pkg.instantiate()` again.
4. If you want a clean Julia session between examples, restart Julia before running the next file.
