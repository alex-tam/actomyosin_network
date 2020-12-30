# Data structures for model parameters
# Alex Tam, 13/10/2020

# Parameters
"Global numerical parameters"
@with_kw struct Numerical_Parameters
    nA::Int = 50 # [-] Number of actin filaments
    nM::Int = 10 # [-] Number of myosin motors
    nT::Int = 601 # [-] Number of time steps
    dt::Float64 = 0.1 # [s] Time step size
    lxx::Float64 = 2.5 # [μm] Reference domain width (x)
    lyy::Float64 = 2.5 # [μm] Reference domain width (y)
    xTol::Float64 = 1e-8 # [-] DOF tolerance for optimisation
    fTol::Float64 = 1e-8 # [-] Objective function tolerance for optimisation
    gTol::Float64 = 1e-8 # [-] Gradient tolerance for optimisation
end

"Actin filament properties"
@with_kw struct Actin_Properties
    nSeg::Int = 4; # [-] Number of segments
    LSeg::Float64 = 0.25 # [μm] Equilibrium segment length
    k::Float64 = 1000 # [pN/μm] Spring constant
    lambda_a::Float64 = 0.01 # [pN/(μm^2)*s] Actin-background drag coefficient
    kappa::Float64 = 0.073 # [pN*μm^2] Flexural rigidity
    lambda_xl::Float64 = 20.0 # [pN/μm*s] Cross-linker drag coefficient
    k_off::Float64 = 0.04 # [/filament/s] Turnover rate
    k_p::Float64 = 0.0 # [/filament/s] Polymerisation rate
end

"Myosin motor properties"
@with_kw struct Myosin_Properties
    k::Float64 = 1000 # [pN/μm] Spring constant
    Fs::Float64 = 5 # [pN] Stall force
    Vm::Float64 = 0.5 # [μm/s] Load-free velocity
    k_off::Float64 = 0.35 # [/motor/s] Reference off-rate in the absence of force
    F_ref::Float64 = 12.6 # [pN] Reference force for computing force-dependent off-rate
end