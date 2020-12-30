# Simulate a 2D actomyosin network
# Alex Tam, 12/10/2020

# Simulate network
"Control function for actomyosin network simulations"
function actomyosin_network(parN, parA, parM, par, trial)
    # Specify domain width
    Lxx::Float64 = 2.5; Lxy::Float64 = 0; Lyx::Float64 = 0; Lyy::Float64 = 2.5; # Actual domain widths
    # Pre-allocate
    Force = [[0.0, 0.0, 0.0, 0.0] for idx in 1:parN.nT]; # [pN] Network force
    Curvature = [0.0 for idx in 1:parN.nT]; # Mean network curvature
    Dipole_Index = [0.0 for idx in 1:parN.nT]; # Mean dipole index
    # Generate initial conditions
    state = State{Float64}(Vector{Vector{Vector}}(), Vector{Vector}()); # Initialise empty State struct
    mm = Vector{Myosin_Motor}(); # Pre-allocate empty myosin motors
    af, state = actin_ic(state, parN, parA, Lxx, Lyy); # Initialise actin filaments
    xl = intersection_search(state, af, mm); # Initialise cross-links
    mm, xl, state = myosin_ic(state, mm, parN, xl, Lxx, Lxy, Lyx, Lyy); # Initialise myosin motors
    state_old = state; # Store initial state to compute force
    # Time-stepping
    animation = @animate for i = 1:parN.nT
        # Write DOF vector to file
        if i == parN.nT
            dof = build_dof(state);
            writedlm("dof-$par-trial-$trial-$i.csv", dof);
        end
        # Spatial statistics
        if parA.nSeg > 1
            Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy); # Compute curvature at current time step
        else
            Curvature[i] = 0;
        end
        savefig("actomyosin_curvature-$par-trial-$trial.png"); # Save histogram of filament mean curvature
        Dipole_Index[i] = dipole_index(mm, state, parN, Lxx, Lxy, Lyx, Lyy); # Compute mean dipole index at current time step
        savefig("actomyosin_dipole_index-$par-trial-$trial.png"); # Save histogram of filament dipole index
        pcf(af, state, Lxx, Lyy); # Compute paired distances between filament nodes
        savefig("actomyosin_pcf-$par-trial-$trial.png"); # Save plot of pair-correlation function
        # Compute force and draw network
        Force[i] = network_force(state, state_old, af, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy);
        draw_network(state, af, xl, mm, parN, parA, Force[i], Lxx, Lxy, Lyx, Lyy);
        if i == 1
            savefig("actomyosin_ic-par-$par-trial-$trial.png"); # Save image of initial condition
        end
        # Simulate one time step
        if i != parN.nT # Ensure correct looping sequence
            # Turnover and polymerisation
            af, mm, state = actin_turnover(state, af, mm, parN, parA, Lxx, Lyy);
            af = segment_translations(state, af); # Update filament translations
            xl = intersection_search(state, af, mm); # Update cross-links
            mm, xl, state = myosin_turnover(state, xl, mm, parN, parM, Lxx, Lxy, Lyx, Lyy); # Myosin binding/unbinding
            # Compute solution and store data
            state_old = state; # Store current state for energy functional
            new_dof = optimise_network(state_old, af, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy); # Compute network solution
            state = build_state(new_dof, af, mm); # Construct State from vector
        end
    end
    # Output
    gif(animation, "actomyosin-par-$par-trial-$trial.gif", fps = 10) # Save .gif of the network
    savefig("actomyosin_end-par-$par-trial-$trial.png");
    times, Tension, Tension_Int = draw_tension(parN, Force, parN.nT);
    @printf("Time-averaged net tension is %f pN/μm*s.\n", Tension_Int[end]/times[end])
    savefig("actomyosin_tension-par-$par-trial-$trial.png");
    Curvature_Int, Dipole_Int = draw_tension_spatial(parN, Force, parN.nT, Curvature, Dipole_Index)
    savefig("actomyosin_tension_spatial-par-$par-trial-$trial.png");
    writedlm("times.csv", times); writedlm("tension-par-$par-trial-$trial.csv", Tension);
    draw_force(parN, Force, parN.nT);
    savefig("actomyosin_force-par-$par-trial-$trial.png");
    return state, af, mm, xl, Tension_Int[end]/times[end], Curvature_Int[end]/times[end], Dipole_Int[end]/times[end]
end