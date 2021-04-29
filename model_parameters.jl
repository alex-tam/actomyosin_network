# Data structures for model parameters
# Alex Tam, 13/10/2020

# Parameters
"Global numerical parameters"
@with_kw struct Numerical_Parameters
    nA::Int = 50 # [-] Number of actin filaments
    nM::Int = 10 # [-] Number of myosin motors
    nT::Int = 1201 # [-] Number of time steps
    dt::Float64 = 0.05 # [s] Time step size
    lxx::Float64 = 2.5 # [μm] Reference domain width (x)
    lyy::Float64 = 2.5 # [μm] Reference domain width (y)
    xTol::Float64 = 1e-8 # [-] DOF tolerance for optimisation
    fTol::Float64 = 1e-8 # [-] Objective function tolerance for optimisation
    gTol::Float64 = 1e-8 # [-] Gradient tolerance for optimisation
end

"Actin filament properties"
@with_kw struct Actin_Properties
    nSeg::Int = 5; # [-] Number of segments
    LSeg::Float64 = 0.2 # [μm] Equilibrium segment length
    k::Float64 = 1000 # [pN/μm] Spring constant
    lambda_a::Float64 = 0.05 # [pN/(μm^2)*s] Actin-background drag coefficient
    kappa::Float64 = 0.073 # [pN*μm^2] Flexural rigidity
    lambda_pf::Float64 = 30.0 # [pN/μm*s] Protein friction drag coefficient
    k_off::Float64 = 0.04 # [/filament/s] Turnover rate
    k_p::Float64 = 0.0 # [/filament/s] Polymerisation rate
    kb::Float64 = 1.380649e-5 # [μm*pN/K] Boltzmann constant
    T::Float64 = 298.15 # [K] Temperature
end

"Myosin motor properties"
@with_kw struct Myosin_Properties
    k::Float64 = 1000 # [pN/μm] Spring constant
    Fs::Float64 = 5 # [pN] Stall force
    Vm::Float64 = 0.5 # [μm/s] Load-free velocity
    k_off::Float64 = 0.35 # [/motor/s] Reference off-rate in the absence of force
    F_ref::Float64 = 12.6 # [pN] Reference force for computing force-dependent off-rate
end