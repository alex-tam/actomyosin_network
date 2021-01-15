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
    nTrials = 1;
    for i = 1:length(par_ref)
        par = par_ref[i]; # Select parameter
        tension = Vector{Float64}(); # Pre-allocate time-averaged tension
        curvature = Vector{Float64}(); # Pre-allocate curvature
        index = Vector{Float64}(); # Pre-allocate integrated tension
        for j = 1:nTrials
            # Specify parameters
            parN = Numerical_Parameters(nT = 101); # Initialise struct of numerical parameters
            parA = Actin_Properties(k_off = par); # Initialise struct of actin filament properties
            parM = Myosin_Properties(); # Initialise struct of myosin motor properties

            # Run simulations
            @time state, af, mm, xl, Tension, Curvature, Index = actomyosin_network(parN, parA, parM, i, j);
            push!(tension, Tension); # Store time-averaged net tension
            push!(curvature, mean(Curvature)/parN.dt); # Store time-averaged net tension
            push!(index, mean(Index)/parN.dt); # Store time-averaged net tension
            writedlm("tension-$par.txt", tension); # Write tension to file
            writedlm("curvature-$par.txt", curvature); # Write tension to file
            writedlm("index-$par.txt", index); # Write tension to file
        end
    end
end

run_simulations();
