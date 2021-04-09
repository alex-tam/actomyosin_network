# Perform actomyosin network simulations
# Alex Tam, 12/10/2020

# Load packages
using Revise # Modify code without restarting Julia
using Parameters # Tools for storing parameters in data structures
using LinearAlgebra # Matrix and vector operations
using StatsBase # Sampling without replacement
using Plots # Plotting library
using Plots.PlotMeasures # Enable re-sizing plot margins
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
include("thermal.jl")
include("network_statistics.jl")
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
        bulk_stress = Vector{Float64}(); # Pre-allocate time-averaged bulk stress
        curvature = Vector{Float64}(); # Pre-allocate curvature
        index = Vector{Float64}(); # Pre-allocate integrated tension
        filament_speed = Vector{Float64}(); # Pre-allocate actin filament node speed
        motor_speed = Vector{Float64}(); # Pre-allocate motor head speed
        angle_roc = Vector{Float64}(); # Pre-allocate rate of change of angles
        motor_pos = Vector{Float64}(); # Pre-allocate time-averaged motor position
        motor_angle = Vector{Float64}(); # Pre-allocate time-averaged motor angle
        for j = 1:nTrials
            # Specify parameters
            parN = Numerical_Parameters(); # Initialise struct of numerical parameters
            parA = Actin_Properties(k_off = par); # Initialise struct of actin filament properties
            parM = Myosin_Properties(); # Initialise struct of myosin motor properties

            # Run simulations
            @time state, af, mm, xl, Bulk_Stress, Curvature, Index, Filament_Speed, Motor_Speed, Angle_ROC, Motor_Pos, Motor_Angle, Bins, Hist = actomyosin_network(parN, parA, parM, i, j);
            push!(bulk_stress, Bulk_Stress); # Store time-averaged net stress
            push!(curvature, Curvature); # Store time-averaged integrated curvature
            push!(index, Index); # Store time-averaged two-filament index
            push!(filament_speed, Filament_Speed); # Store time-averaged filament node speed
            push!(motor_speed, Motor_Speed); # Store time-averaged motor head speed
            push!(angle_roc, Angle_ROC); # Store time-averaged rate of change of angle
            push!(motor_pos, Motor_Pos); # Store time-averaged motor position
            push!(motor_angle, Motor_Angle); # Store time-averaged motor angle
            writedlm("bulk_stress-$par.txt", bulk_stress); # Write time-averaged bulk stress per trial to file
            writedlm("curvature-$par.txt", curvature); # Write time-averaged curvature per trial to file
            writedlm("index-$par.txt", index); # Write time-averaged two-filament index per trial to file
            writedlm("filament_speed-$par.txt", filament_speed); # Write time-averaged actin filament node speed per trial to file
            writedlm("motor_speed-$par.txt", motor_speed); # Write time-averaged motor head speed per trial to file
            writedlm("angle_roc-$par.txt", angle_roc); # Write time-averaged angle rate of change to file
            writedlm("motor_pos-$par.txt", motor_pos); # Write time-averaged motor position to file
            writedlm("motor_angle-$par.txt", motor_angle); # Write time-averaged motor angle to file
            writedlm("pd_end_bins-par-$i-trial-$j.txt", Bins); # Write paried-distances bins to file
            writedlm("pd_end_counts-par-$i-trial-$j.txt", Hist); # Write paired-distances counts to file
        end
    end
end

run_simulations();
