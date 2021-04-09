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
    Motor_Pos = [0.0 for idx in 1:parN.nT]; # [--] Mean myosin motor head position
    Motor_Angle = [0.0 for idx in 1:parN.nT]; # [rad] Mean angle between filaments at motor binding site
    Bins = Vector{Float64}(); Hist = Vector{Int}(); # [μm, --] Pre-allocate paired-distance histogram
    # Generate initial conditions
    state = State{Float64}(Vector{Vector{Vector}}(), Vector{Vector}()); # Initialise empty State struct
    mm = Vector{Myosin_Motor}(); # Pre-allocate empty myosin motors
    af, state = actin_ic(state, parN, parA, Lxx, Lyy); # Initialise actin filaments
    xl = intersection_search(state, af, mm); # Initialise cross-links
    mm, xl, state = myosin_ic(state, mm, parN, xl, Lxx, Lxy, Lyx, Lyy); # Initialise myosin motors
    state_old = state; # Store initial state to compute force
    random = thermal(af, state); # Random variables for thermal motion
    # Time-stepping
    animation = @animate for i = 1:parN.nT
        # Write final DOF vector to file
        if i == parN.nT
            dof = build_dof(state);
            writedlm("dof-$par-trial-$trial-$i.csv", dof);
        end
        # Compute network characteristics
        Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_curvature-$par-trial-$trial.svg"); # Mean filament curvature
        Index[i] = two_filament_index(mm, state, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_2f_index-$par-trial-$trial.svg"); # Mean two-filament index
        Bins, Hist = pcf(af, state, Lxx, Lyy); savefig("actomyosin_pcf-$par-trial-$trial.svg"); # Paired distances between filament nodes
        Filament_Speed[i] = filament_speed(parN, af, state, state_old, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_filament_speed-$par-trial-$trial.svg"); # Filament node speeds
        Motor_Speed[i] = motor_speed(parN, mm, state, state_old, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_motor_speed-$par-trial-$trial.svg"); # Motor head speeds
        Angle_ROC[i] = motor_angle_roc(parN, mm, state, state_old, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_angle_roc-$par-trial-$trial.svg"); # Angle rate of change
        Motor_Pos[i] = motor_position(mm, state); # savefig("actomyosin_angle_roc-$par-trial-$trial.svg"); # Myosin motor head positions
        Motor_Angle[i] = motor_angle(mm, state, parN, Lxx, Lxy, Lyx, Lyy); # savefig("actomyosin_angle_roc-$par-trial-$trial.svg"); # Myosin motor angles
        Force[i] = network_force(state, state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy); # Force
        # Draw network
        draw_network(state, af, xl, mm, parN, parA, Force[i], Lxx, Lxy, Lyx, Lyy);
        if i == 1
            savefig("actomyosin_ic-par-$par-trial-$trial.svg"); # Save image of initial condition
        end
        # Simulate one time step
        if i != parN.nT # Ensure correct looping sequence
            # Turnover and polymerisation
            af, mm, state = actin_turnover(state, af, mm, parN, parA, Lxx, Lyy);
            af = segment_translations(state, af); # Update filament translations
            xl = intersection_search(state, af, mm); # Update cross-links
            mm, xl, state = myosin_turnover(state, xl, mm, parN, parM, Lxx, Lxy, Lyx, Lyy); # Myosin binding/unbinding
            random = thermal(af, state); # Random variables for thermal motion
            # Compute solution and store data
            state_old = state; # Store current state for energy functional
            new_dof = optimise_network(state_old, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy); # Compute network solution
            state = build_state(new_dof, af, mm); # Construct State from vector
        end
    end
    # Output .gif and image of the network
    gif(animation, "actomyosin-par-$par-trial-$trial.gif", fps = 10); savefig("actomyosin_end-par-$par-trial-$trial.svg");
    # Plot quantities versus time
    draw_stress(parN, Force, parN.nT, Lxx, Lyy); savefig("actomyosin_stress-par-$par-trial-$trial.svg"); # Stress components
    times, Bulk_Stress, Bulk_Stress_Int = draw_bulk_stress(parN, Force, parN.nT, Lxx, Lyy); savefig("actomyosin_bulk_stress-par-$par-trial-$trial.svg"); # Bulk stress
    @printf("Time-averaged bulk stress is %f pN/μm*s.\n", Bulk_Stress_Int[end]/times[end])
    # Compute time-integrated network statistics
    Curvature_Int, Index_Int, Filament_Speed_Int, Motor_Speed_Int, Angle_ROC_Int, Motor_Pos_Int, Motor_Angle_Int = integrated_statistics(parN, Curvature, Index, Filament_Speed, Motor_Speed, Angle_ROC, Motor_Pos, Motor_Angle)
    # Save mean stress and spatial measures time series data to files
    writedlm("times.csv", times);
    writedlm("stress-par-$par-trial-$trial.csv", Bulk_Stress);
    writedlm("curvature-par-$par-trial-$trial.csv", Curvature);
    writedlm("index-par-$par-trial-$trial.csv", Index);
    writedlm("filament_speed-par-$par-trial-$trial.csv", Filament_Speed);
    writedlm("motor_speed-par-$par-trial-$trial.csv", Motor_Speed);
    writedlm("motor_pos-par-$par-trial-$trial.csv", Motor_Pos);
    writedlm("motor_angle-par-$par-trial-$trial.csv", Motor_Angle);
    return state, af, mm, xl, Bulk_Stress_Int[end]/times[end], Curvature_Int[end]/times[end], Index_Int[end]/times[end], Filament_Speed_Int[end]/times[end], Motor_Speed_Int[end]/times[end], Angle_ROC_Int[end]/times[end], Motor_Pos_Int[end]/times[end], Motor_Angle_Int[end]/times[end], Bins, Hist
end