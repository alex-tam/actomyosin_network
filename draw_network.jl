# Draw actomyosin networks, force components, and stress
# Alex Tam, 12/10/2020

"Draw the network configuration"
function draw_network(s, af, xl, mm, parN, parA, Force, Lxx, Lxy, Lyx, Lyy)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    default(titlefont = (18, "times"), guidefont = (26, "times"), tickfont = (18, "times"))
    # Actin Filaments
    fx = zeros(Float64, parA.nSeg+1, length(af));
    fy = zeros(Float64, parA.nSeg+1, length(af)); # Pre-allocate un-translated filament node positions
    n1x = Vector{Float64}(); n1y = Vector{Float64}(); 
    n2x = Vector{Float64}(); n2y = Vector{Float64}(); # Pre-allocate translated node positions (non-plus end segments)
    p1x = Vector{Float64}(); p1y = Vector{Float64}(); 
    p2x = Vector{Float64}(); p2y = Vector{Float64}(); # Pre-allocate translated node positions (plus end segments)
    for f in af
        nx = Vector{}(); ny = Vector{}(); # Pre-allocate vectors of current filament nodes
        for seg in f.segments
            for j = 1:length(seg.t)
                if j == 1
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], Lxx, Lxy, Lyx, Lyy); # Get segment nodes
                    push!(nx, mx); push!(nx, px); # Store unprojected nodes of a segment (x)
                    push!(ny, my); push!(ny, py); # Store unprojected nodes of a segment (y)
                else
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], Lxx, Lxy, Lyx, Lyy); # Get segment nodes
                    if seg.index == length(f.segments)
                        push!(p1x, mx); push!(p2x, px); push!(p1y, my); push!(p2y, py); # Store nodes (plus-end)
                    else
                        push!(n1x, mx); push!(n2x, px); push!(n1y, my); push!(n2y, py); # Store nodes (non-plus-end)
                    end
                end
            end
        end
        nx = unique(nx); ny = unique(ny); # Remove duplicates
        fx[:, f.index] = nx; fy[:, f.index] = ny; # Store filament nodes in matrix
    end
    plot!(fx, fy, legend = false, color = "red", aspect_ratio=:equal, xlabel = L"$x$", ylabel = L"$y$", xlims = [0, parN.lxx], ylims = [0, parN.lyy]) # Plot un-translated segments
    scatter!(fx[1:end-1,:], fy[1:end-1,:], color = "red"); # Plot nodes as points
    scatter!(fx[end,:], fy[end,:], color = "black"); # Plot plus end nodes
    plot!(transpose(hcat(n1x, n2x)), transpose(hcat(n1y, n2y)), color = "green") # Plot translated segments (non-plus-end)
    plot!(transpose(hcat(p1x, p2x)), transpose(hcat(p1y, p2y)), color = "green") # Plot translated segments (plus end)
    scatter!([n1x, n2x], [n1y, n2y], color = "green"); # Plot translated nodes (non-plus-ends)
    scatter!([p1x], [p1y], color = "green"); # Plot translated nodes
    scatter!([p2x], [p2y], color = "black"); # Plot translated nodes (plus ends)
    # Cross-Links
    l1x = Vector{Float64}(); l1y = Vector{Float64}(); # Pre-allocate outside loop to plot all  cross-links at once
    l2x = Vector{Float64}(); l2y = Vector{Float64}();
    for l in xl
        # Obtain cross-link positions (dimensionless, un-translated)
        x1, y1, x2, y2 = get_xl_pos(l, s);
        # Convert to dimensional, translated positions
        lx1 = (x1 - l.t1[1])*Lxx + (y1 - l.t1[2])*Lyx;
        ly1 = (x1 - l.t1[1])*Lxy + (y1 - l.t1[2])*Lyy;
        lx2 = (x2 - l.t2[1])*Lxx + (y2 - l.t2[2])*Lyx;
        ly2 = (x2 - l.t2[1])*Lxy + (y2 - l.t2[2])*Lyy;
        push!(l1x, lx1); push!(l1y, ly1); push!(l2x, lx2); push!(l2y, ly2);
    end
    scatter!([l1x, l2x], [l1y, l2y], color = "pink") # Plot cross-link sites as points
    # plot!(transpose(hcat(l1x, l2x)), transpose(hcat(l1y, l2y)), color = "pink") # Draw cross-linkers
    # Myosin Motors
    m1x = Vector{Float64}(); tm1x = Vector{Float64}(); # Pre-allocate outside loop to plot all motors at once
    m1y = Vector{Float64}(); tm1y = Vector{Float64}();
    m2x = Vector{Float64}(); tm2x = Vector{Float64}();
    m2y = Vector{Float64}(); tm2y = Vector{Float64}();
    for m in mm
        x1, y1, x2, y2 = get_motor_pos(m, s, Lxx, Lxy, Lyx, Lyy); # Obtain motor positions (dimensionless, un-translated)
        # Convert to dimensional, translated positions
        mx1 = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx;
        my1 = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy;
        mx2 = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx;
        my2 = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy;
        # Compute translations for visualisation only
        myosin_translations = periodic([mx1/parN.lxx, my1/parN.lyy], [mx2/parN.lxx, my2/parN.lyy]);
        for j = 1:length(myosin_translations)
            if j == 1
                push!(m1x, mx1); push!(m1y, my1); push!(m2x, mx2); push!(m2y, my2);
            else
                t = myosin_translations[j];
                # Store dimensional, translated positions
                push!(tm1x, mx1 - (t[1]*Lxx + t[2]*Lyx)); 
                push!(tm1y, my1 - (t[1]*Lxy + t[2]*Lyy));
                push!(tm2x, mx2 - (t[1]*Lxx + t[2]*Lyx)); 
                push!(tm2y, my2 - (t[1]*Lxy + t[2]*Lyy));
            end
        end
    end
    scatter!([m1x, m2x], [m1y, m2y], color = "blue"); # Plot binding sites as points
    scatter!([tm1x, tm2x], [tm1y, tm2y], color = "orange"); # Plot binding sites as points
    plot!(transpose(hcat(m1x, m2x)), transpose(hcat(m1y, m2y)), color = "blue") # Draw motors
    plot!(transpose(hcat(tm1x, tm2x)), transpose(hcat(tm1y, tm2y)), color = "orange") # Draw projected motors
    # Principal stress directions
    Stress = [ Force[1]/Lyy Force[3]/Lxx ; Force[2]/Lyy Force[4]/Lxx ]; # Compute stress tensor
    evals = eigvals(Stress); evecs = eigvecs(Stress); # Compute eigenvalues and eigenvectors
    if evals[1] >= 0
        plot!([parN.lxx/2, parN.lxx/2 + evals[1]*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + evals[1]*evecs[1,2]], arrow = :arrow, color = "orange")
    else
        plot!([parN.lxx/2, parN.lxx/2 + abs(evals[1])*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[1])*evecs[1,2]], arrow = :arrow, color = "blue")
    end
    if evals[2] >= 0
        plot!([parN.lxx/2, parN.lxx/2 + evals[2]*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + evals[2]*evecs[2,2]], arrow = :arrow, color = "orange")
    else
        plot!([parN.lxx/2, parN.lxx/2 + abs(evals[2])*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[2])*evecs[2,2]], arrow = :arrow, color = "blue")
    end
end

"Plot force components"
function draw_force(parN, Force, time_step)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    default(titlefont = (18, "times"), guidefont = (26, "times"), tickfont = (18, "times"))
    # Pre-allocate
    cfxx = Vector{Float64}(); # Pre-allocate vector for force component (xx)
    cfxy = Vector{Float64}(); # Pre-allocate vector for force component (xy)
    cfyx = Vector{Float64}(); # Pre-allocate vector for force component (yx)
    cfyy = Vector{Float64}(); # Pre-allocate vector for force component (yy)
    cfxx_tot = Vector{Float64}(); # Pre-allocate vector for force component (xx)
    cfxy_tot = Vector{Float64}(); # Pre-allocate vector for force component (xy)
    cfyx_tot = Vector{Float64}(); # Pre-allocate vector for force component (yx)
    cfyy_tot = Vector{Float64}(); # Pre-allocate vector for force component (yy)
    times = (0:time_step-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:time_step
        push!(cfxx, Force[i][1]); # Construct vector for force component (xx)
        push!(cfxy, Force[i][2]); # Construct vector for force component (xy)
        push!(cfyx, Force[i][3]); # Construct vector for force component (yx)
        push!(cfyy, Force[i][4]); # Construct vector for force component (yy)
        push!(cfxx_tot, trapz(cfxx, times)); # Integrate force over time
        push!(cfxy_tot, trapz(cfxy, times)); # Integrate force over time
        push!(cfyx_tot, trapz(cfyx, times)); # Integrate force over time
        push!(cfyy_tot, trapz(cfyy, times)); # Integrate force over time
    end
    # Plot network force
    p1 = plot(times, cfxx, color = "orange", linewidth = 2, label=:false, xlabel = L"$t$", ylabel = "Force (pN)", legend=:bottomleft, size = (1000, 500));
    plot!(times, cfxy, color = "purple", linewidth = 2, label=:false);
    plot!(times, cfyx, color = "green", linewidth = 2, label=:false);
    plot!(times, cfyy, color = "blue", linewidth = 2, label=:false);
    plot!(times, movingaverage(cfxx, 100), color = "orange", linewidth = 2, linestyle = :dash, label = L"F_{xx}");
    plot!(times, movingaverage(cfxy, 100), color = "purple", linewidth = 2, linestyle = :dash, label = L"F_{xy}");
    plot!(times, movingaverage(cfyx, 100), color = "green", linewidth = 2, linestyle = :dash, label = L"F_{yx}");
    plot!(times, movingaverage(cfyy, 100), color = "blue", linewidth = 2, linestyle = :dash, label = L"F_{yy}");
    # Plot integrated force
    p2 = plot(times, cfxx_tot, color = "orange", linewidth=2, label = L"\int F_{xx} \; dt", xlabel = L"$t$", ylabel = "Integrated Force (pNs)", legend=:bottomleft, size = (1000, 500));
    plot!(times, cfxy_tot, color = "purple", linewidth=2, label = L"\int F_{xy} \;dt");
    plot!(times, cfyx_tot, color = "green", linewidth=2, label = L"\int F_{yx} \;dt");
    plot!(times, cfyy_tot, color = "blue", linewidth=2, label = L"\int F_{yy} \;dt");
    plot(p1, p2, layout = (1,2))
end

"Plot bulk stress and moving average"
function draw_bulk_stress(parN, Force, time_step, Lxx, Lyy)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    default(titlefont = (18, "times"), guidefont = (26, "times"), tickfont = (18, "times"))
    Bulk_Stress = Vector{Float64}(); Bulk_Stress_Int = Vector{Float64}(); # Pre-allocate
    times = (0:parN.nT-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:parN.nT
        Stress = [ Force[i][1]/Lyy Force[i][3]/Lxx ; Force[i][2]/Lyy Force[i][4]/Lxx ]; # Compute stress tensor
        push!(Bulk_Stress, sum(eigvals(Stress))/2);
        push!(Bulk_Stress_Int, trapz(Bulk_Stress, times)); # Integrate bulk stress over time
    end
    Bulk_Stress_Moving = movingaverage(Bulk_Stress[1:time_step], 100);
    plot(times[1:time_step], Bulk_Stress[1:time_step], linewidth = 2, label = "Bulk Stress", xlabel = L"$t$", ylabel = "Bulk Stress (pN/μm)", legend=:topleft);
    plot!(times[1:time_step], Bulk_Stress_Moving[1:time_step], linewidth = 2, linestyle = :dot, label = "Moving Average");
    return times, Bulk_Stress, Bulk_Stress_Int
end

"Plot bulk stress against spatial measures"
function draw_stress_spatial(parN, Force, time_step, Curvature, Index, Filament_Speed, Motor_Speed, Angle_ROC, Motor_Pos, Motor_Angle, Lxx, Lyy)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    default(titlefont = (18, "times"), guidefont = (26, "times"), tickfont = (18, "times"))
    Bulk_Stress = Vector{Float64}(); # Pre-allocate bulk stress
    Tension_Int = Vector{Float64}(); Curvature_Int = Vector{Float64}(); Index_Int = Vector{Float64}(); Filament_Speed_Int = Vector{Float64}(); 
    Motor_Speed_Int = Vector{Float64}(); Angle_ROC_Int = Vector{Float64}(); Motor_Pos_Int = Vector{Float64}(); Motor_Angle_Int = Vector{Float64}(); # Pre-allocate time-integrated variables
    times = (0:parN.nT-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:parN.nT
        Stress = [ Force[i][1]/Lyy Force[i][3]/Lxx ; Force[i][2]/Lyy Force[i][4]/Lxx ]; # Compute stress tensor
        push!(Bulk_Stress, sum(eigvals(Stress))/2);
        push!(Tension_Int, trapz(Bulk_Stress, times)); # Integrate force over time
        push!(Curvature_Int, trapz(Curvature, times)); # Integrate curvature over time
        push!(Index_Int, trapz(Index, times)); # Integrate two-filament index over time
        push!(Filament_Speed_Int, trapz(Filament_Speed, times)); # Integrate filament node speed over time
        push!(Motor_Speed_Int, trapz(Motor_Speed, times)); # Integrate motor head speed over time
        push!(Angle_ROC_Int, trapz(Angle_ROC, times)); # Integrate motor head speed over time
        push!(Motor_Pos_Int, trapz(Motor_Pos, times)); # Integrate motor head position over time
        push!(Motor_Angle_Int, trapz(Motor_Angle, times)); # Integrate motor angle over time
    end
    Tension_Moving = movingaverage(Bulk_Stress[1:time_step], 100);
    plot(times[1:time_step], Bulk_Stress[1:time_step], linewidth = 2, label = "Bulk Stress", xlabel = L"$t$", ylabel = "Bulk Stress (pN/μm)", legend=:topleft);
    plot!(times[1:time_step], Tension_Moving[1:time_step], linewidth = 2, linestyle = :dot, label = "Moving Average");
    # plot!(times[1:time_step], Curvature[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(Curvature), linewidth = 2, label = "Curvature");
    plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Index)), linewidth = 2, label = "Two-Filament Index");
    # plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Filament_Speed)), linewidth = 2, label = "Filament Speed");
    # plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Motor_Speed)), linewidth = 2, label = "Motor Speed");
    # plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Angle_ROC)), linewidth = 2, label = "Rate of Change of Angle");
    # plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Motor_Pos_Int)), linewidth = 2, label = "Motor Head Position");
    # plot!(times[1:time_step], Index[1:time_step]*maximum(abs.(Bulk_Stress))/maximum(abs.(Motor_Angle_Int)), linewidth = 2, label = "Motor Angle");
    return Curvature_Int, Index_Int, Filament_Speed_Int, Motor_Speed_Int, Angle_ROC_Int, Motor_Pos_Int, Motor_Angle_Int
end

"Integrate vector data using the trapezoidal rule"
function trapz(y, x)
    int = 0.0
    for i = 1:length(y)-1
        int += 0.5*(y[i] + y[i+1])*(x[i+1]-x[i]); # Composite trapezoidal rule
    end
    return int
end

"Calculate a moving average"
function movingaverage(x::Vector, numofele::Int)
    BackDelta = div(numofele,2)
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(x)
    y = similar(x)
    for n = 1:len
        lo = max(1, n-BackDelta)
        hi = min(len, n+ForwardDelta)
        y[n] = mean(x[lo:hi])
    end
    return y
end