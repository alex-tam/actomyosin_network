# Perform actomyosin network simulations
# Alex Tam, 12/10/2020

# Load packages
using Revise # Modify code without restarting Julia
using Parameters # Tools for storing parameters in data structures
using LinearAlgebra # Matrix and vector operations
using StatsBase # Sampling without replacement
using Plots # Plotting library
using LaTeXStrings # Display LaTeX output
using ForwardDiff # Automatic differentiation tools
using Optim # Optimisation routines
using LineSearches # Line searches for optimisation
using Printf # C-style printing macros
using Distributions # Probability distributions
using DelimitedFiles # Read and write delimited data

# Include code from external files
include("model_parameters.jl")
include("actomyosin_network.jl")
include("State.jl")
include("initial_condition.jl")
include("Actin_Filament.jl")
include("periodic.jl")
include("Cross_Link.jl")
include("Myosin_Motor.jl")
include("intersection_search.jl")
include("spatial_statistics.jl")
include("network_force.jl")
include("draw_network.jl")
include("turnover.jl")
include("optimise_network.jl")
include("energy.jl")

function run_simulations()
    par_ref = 0.04;
    nTrials = 10;
    for i = 1:length(par_ref)
        par = par_ref[i]; # Select parameter
        bulk_stress = Vector{Float64}(); # Pre-allocate time-averaged bulk stress
        curvature = Vector{Float64}(); # Pre-allocate curvature
        index = Vector{Float64}(); # Pre-allocate integrated tension
        for j = 1:nTrials
            # Specify parameters
            parN = Numerical_Parameters(); # Initialise struct of numerical parameters
            parA = Actin_Properties(k_off = par); # Initialise struct of actin filament properties
            parM = Myosin_Properties(); # Initialise struct of myosin motor properties

            # Run simulations
            @time state, af, mm, xl, Bulk_Stress, Curvature, Index = actomyosin_network(parN, parA, parM, i, j);
            push!(bulk_stress, Bulk_Stress); # Store time-averaged net stress
            push!(curvature, mean(Curvature)/parN.dt); # Store time-averaged integrated curvature
            push!(index, mean(Index)/parN.dt); # Store time-averaged two-filament index
            writedlm("bulk_stress-$par.txt", bulk_stress); # Write bulk stress to file
            writedlm("curvature-$par.txt", curvature); # Write curvature to file
            writedlm("index-$par.txt", index); # Write two-filament index to file
        end
    end
end

run_simulations();
