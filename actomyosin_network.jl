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
    Index = [0.0 for idx in 1:parN.nT]; # Mean two-filament index
    Filament_Speed = [0.0 for idx in 1:parN.nT]; # [μm/s] Mean filament node speed
    Motor_Speed = [0.0 for idx in 1:parN.nT]; # [μm/s] Mean motor head speed
    Angle_ROC = [0.0 for idx in 1:parN.nT]; # [rad/s] Mean rate of change of angle
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
        # savefig("actomyosin_curvature-$par-trial-$trial.png"); # Save histogram of filament mean curvature
        Index[i] = two_filament_index(mm, state, parN, Lxx, Lxy, Lyx, Lyy); # Compute mean two-filament index at current time step
        # savefig("actomyosin_2f_index-$par-trial-$trial.png"); # Save histogram of two-filament index
        pcf(af, state, Lxx, Lyy); # Compute paired distances between filament nodes
        savefig("actomyosin_pcf-$par-trial-$trial.png"); # Save plot of pair-correlation function
        Filament_Speed[i] = filament_speed(parN, af, state, state_old, Lxx, Lxy, Lyx, Lyy);
        # savefig("actomyosin_filament_speed-$par-trial-$trial.png"); # Save histogram of filament node speeds
        Motor_Speed[i] = motor_speed(parN, mm, state, state_old, Lxx, Lxy, Lyx, Lyy);
        # savefig("actomyosin_motor_speed-$par-trial-$trial.png"); # Save histogram of motor head speeds
        Angle_ROC[i] = motor_angle_roc(parN, mm, state, state_old, Lxx, Lxy, Lyx, Lyy);
        # savefig("actomyosin_angle_roc-$par-trial-$trial.png"); # Save histogram of angle rate of change
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
    # Output .gif and image of the network
    gif(animation, "actomyosin-par-$par-trial-$trial.gif", fps = 10)
    savefig("actomyosin_end-par-$par-trial-$trial.png");
    # Draw force components
    draw_force(parN, Force, parN.nT);
    savefig("actomyosin_force-par-$par-trial-$trial.png");
    # Plot bulk stress versus time
    times, Bulk_Stress, Bulk_Stress_Int = draw_stress(parN, Force, parN.nT);
    @printf("Time-averaged net stress is %f pN/μm*s.\n", Bulk_Stress_Int[end]/times[end])
    savefig("actomyosin_stress-par-$par-trial-$trial.png");
    # Plot bulk stress with network statistics versus time
    Curvature_Int, Index_Int, Filament_Speed_Int, Motor_Speed_Int, Angle_ROC_Int = draw_stress_spatial(parN, Force, parN.nT, Curvature, Index, Filament_Speed, Motor_Speed, Angle_ROC)
    savefig("actomyosin_stress_spatial-par-$par-trial-$trial.png");
    # Save mean stress and spatial measures time series data to a file
    writedlm("times.csv", times);
    writedlm("stress-par-$par-trial-$trial.csv", Bulk_Stress);
    writedlm("curvature-par-$par-trial-$trial.csv", Curvature);
    writedlm("index-par-$par-trial-$trial.csv", Index);
    writedlm("filament_speed-par-$par-trial-$trial.csv", Filament_Speed);
    writedlm("motor_speed-par-$par-trial-$trial.csv", Motor_Speed);
    return state, af, mm, xl, Bulk_Stress_Int[end]/times[end], Curvature_Int[end]/times[end], Index_Int[end]/times[end], Filament_Speed_Int[end]/times[end], Motor_Speed_Int[end]/times[end], Angle_ROC_Int[end]/times[end]
end