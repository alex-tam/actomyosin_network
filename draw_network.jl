# Draw actomyosin networks, force components, and tension
# Alex Tam, 12/10/2020

"Draw the network configuration"
function draw_network(s, af, xl, mm, parN, parA, Force, Lxx, Lxy, Lyx, Lyy)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
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
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], parN, Lxx, Lxy, Lyx, Lyy); # Get segment nodes
                    push!(nx, mx); push!(nx, px); # Store unprojected nodes of a segment (x)
                    push!(ny, my); push!(ny, py); # Store unprojected nodes of a segment (y)
                else
                    mx, my, px, py = get_segment_nodes(f, seg, s, seg.t[j], parN, Lxx, Lxy, Lyx, Lyy); # Get segment nodes
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
        lx1 = (x1 - l.t1[1]*Lxx/parN.lxx - l.t1[2]*Lyx/parN.lxx)*Lxx + (y1 - l.t1[1]*Lxy/parN.lyy - l.t1[2]*Lyy/parN.lyy)*Lyx;
        ly1 = (y1 - l.t1[1]*Lxy/parN.lyy - l.t1[2]*Lyy/parN.lyy)*Lyy + (x1 - l.t1[1]*Lxx/parN.lxx - l.t1[2]*Lyx/parN.lxx)*Lxy;
        lx2 = (x2 - l.t2[1]*Lxx/parN.lxx - l.t2[2]*Lyx/parN.lxx)*Lxx + (y2 - l.t2[1]*Lxy/parN.lyy - l.t2[2]*Lyy/parN.lyy)*Lyx;
        ly2 = (y2 - l.t2[1]*Lxy/parN.lyy - l.t2[2]*Lyy/parN.lyy)*Lyy + (x2 - l.t2[1]*Lxx/parN.lxx - l.t2[2]*Lyx/parN.lxx)*Lxy;
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
        mx1 = (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxx + (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyx;
        my1 = (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyy + (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxy;
        mx2 = (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxx + (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyx;
        my2 = (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyy + (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxy;
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
    # # Principal stress directions
    # Stress = [ Force[1]/parN.lxx Force[2]/parN.lyy ; Force[3]/parN.lxx Force[4]/parN.lyy ]; # Compute stress tensor
    # evals = eigvals(Stress); evecs = eigvecs(Stress); # Compute eigenvalues and eigenvectors
    # if evals[1] >= 0
    #     plot!([parN.lxx/2, parN.lxx/2 + evals[1]*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + evals[1]*evecs[1,2]], arrow = :arrow, color = "orange")
    # else
    #     plot!([parN.lxx/2, parN.lxx/2 + abs(evals[1])*evecs[1,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[1])*evecs[1,2]], arrow = :arrow, color = "blue")
    # end
    # if evals[2] >= 0
    #     plot!([parN.lxx/2, parN.lxx/2 + evals[2]*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + evals[2]*evecs[2,2]], arrow = :arrow, color = "orange")
    # else
    #     plot!([parN.lxx/2, parN.lxx/2 + abs(evals[2])*evecs[2,1]], [parN.lyy/2 , parN.lyy/2 + abs(evals[2])*evecs[2,2]], arrow = :arrow, color = "blue")
    # end
end

"Plot force components"
function draw_force(parN, Force, time_step)
    cfxx = Vector{Float64}(); # Pre-allocate vector for force component (x)
    cfxy = Vector{Float64}(); # Pre-allocate vector for force component (y)
    cfyx = Vector{Float64}(); # Pre-allocate vector for force component (y)
    cfyy = Vector{Float64}(); # Pre-allocate vector for force component (
    cfxx_tot = Vector{Float64}(); # Pre-allocate vector for force component (x)
    cfxy_tot = Vector{Float64}(); # Pre-allocate vector for force component (y)
    cfyx_tot = Vector{Float64}(); # Pre-allocate vector for force component (x)
    cfyy_tot = Vector{Float64}(); # Pre-allocate vector for force component (x)
    times = (0:time_step-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:time_step
        push!(cfxx, Force[i][1]); # Construct vector for force component (x)
        push!(cfxy, Force[i][2]); # Construct vector for force component (y)
        push!(cfyx, Force[i][3]); # Construct vector for force component (y)
        push!(cfyy, Force[i][4]); # Construct vector for force component (y)
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

"Plot tension"
function draw_tension(parN, Force, time_step)
    Tension = Vector{Float64}(); Tension_Int = Vector{Float64}(); # Pre-allocate
    times = (0:parN.nT-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:parN.nT
        Stress = [ Force[i][1]/parN.lxx Force[i][2]/parN.lyy ; Force[i][3]/parN.lxx Force[i][4]/parN.lyy];
        push!(Tension, sum(eigvals(Stress))/2);
        push!(Tension_Int, trapz(Tension, times)); # Integrate force over time
    end
    Tension_Moving = movingaverage(Tension[1:time_step], 100);
    plot(times[1:time_step], Tension[1:time_step], linewidth = 2, label = "Bulk Stress", xlabel = L"$t$", ylabel = "Bulk Stress (pN/μm)", legend=:topleft);
    plot!(times[1:time_step], Tension_Moving[1:time_step], linewidth = 2, linestyle = :dot, label = "Moving Average");
    return times, Tension, Tension_Int
end

"Plot tension with spatial measures"
function draw_tension_spatial(parN, Force, time_step, Curvature, Dipole_Index)
    Tension = Vector{Float64}(); # Pre-allocate bulk stress
    Tension_Int = Vector{Float64}(); Curvature_Int = Vector{Float64}(); Dipole_Int = Vector{Float64}(); # Pre-allocate time-integrated variables
    times = (0:parN.nT-1).*parN.dt; # Vector of times at which we obtain measurements
    for i = 1:parN.nT
        Stress = [ Force[i][1]/parN.lxx Force[i][2]/parN.lyy ; Force[i][3]/parN.lxx Force[i][4]/parN.lyy];
        push!(Tension, sum(eigvals(Stress))/2);
        push!(Tension_Int, trapz(Tension, times)); # Integrate force over time
        push!(Curvature_Int, trapz(Curvature, times)); # Integrate force over time
        push!(Dipole_Int, trapz(Dipole_Index, times)); # Integrate force over time
    end
    Tension_Moving = movingaverage(Tension[1:time_step], 100);
    plot(times[1:time_step], Tension[1:time_step], linewidth = 2, label = "Bulk Stress", xlabel = L"$t$", ylabel = "Bulk Stress (pN/μm)", legend=:topleft);
    plot!(times[1:time_step], Tension_Moving[1:time_step], linewidth = 2, linestyle = :dot, label = "Moving Average");
    # plot!(times[1:time_step], Curvature[1:time_step]*maximum(abs.(Tension))/maximum(Curvature), linewidth = 2, label = "Curvature");
    plot!(times[1:time_step], Dipole_Index[1:time_step]*maximum(abs.(Tension))/maximum(abs.(Dipole_Index)), linewidth = 2, label = "Dipole Index");
    return Curvature_Int, Dipole_Int
end

"Integrate vector data using the trapezoidal rule"
function trapz(y, x)
    int = 0.0
    for i = 1:length(y)-1
        int += 0.5*(y[i] + y[i+1])*(x[i+1]-x[i]); # Composite trapezoidal rule
    end
    return int
end

"Moving average"
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