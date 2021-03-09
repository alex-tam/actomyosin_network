# Statistics for interpreting actomyosin simulation results
# Alex Tam, 23/11/20

"Actin filament curvature"
function curvature(af, s, Lxx, Lxy, Lyx, Lyy)
    filament_curvatures = Vector{Float64}(); # Pre-allocate mean curvature data
    # Loop over filaments
    for j = 1:length(af)
        f::Actin_Filament = af[j]; # Extract current filament
        seg_lengths = get_segment_lengths(f, s, Lxx, Lxy, Lyx, Lyy); # Dimensional segment lengths
        total_filament_curvature = 0.0; # Pre-allocate integrated curvature of current filament
        # Compute curvature at interior nodes
        for i = 1:length(s.an[f.index])-2 
            # Extract physical positions of plus node
            xp = s.an[f.index][i+2][1]*Lxx + s.an[f.index][i+2][2]*Lyx;
            yp = s.an[f.index][i+2][1]*Lxy + s.an[f.index][i+2][2]*Lyy;
            # Extract physical positions of centre node
            xc = s.an[f.index][i+1][1]*Lxx + s.an[f.index][i+1][2]*Lyx;
            yc = s.an[f.index][i+1][1]*Lxy + s.an[f.index][i+1][2]*Lyy;
            # Extract physical positions of minus node
            xm = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx;
            ym = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy;
            # Compute numerical second derivatives
            L01 = seg_lengths[i]; L12 = seg_lengths[i+1]; # Extract segment lengths
            avl = (L01 + L12)/2; # Average length
            ddfx = ((xp-xc)/L12 - (xc-xm)/L01)/avl; # Numerical second derivative (x)
            ddfy = ((yp-yc)/L12 - (yc-ym)/L01)/avl; # Numerical second derivative (y)
            # Add contribution to integrated curvature
            total_filament_curvature += sqrt(ddfx^2 + ddfy^2)*avl;
        end
        push!(filament_curvatures, total_filament_curvature);
    end
    # histogram(filament_curvatures)
    return mean(filament_curvatures)
end

"Two-filament index"
function two_filament_index(mm, s, Lxx, Lxy, Lyx, Lyy)
    index = Vector{Float64}(); # Pre-allocate index data
    # Loop over motors
    for i = 1:length(mm)
        m::Myosin_Motor = mm[i]; # Extract current motor
        f1::Actin_Filament = m.f1; f2::Actin_Filament = m.f2; # Extract filaments attached to current motor
        # Calculate filament lengths
        L1 = sum(get_segment_lengths(f1, s, Lxx, Lxy, Lyx, Lyy));
        L2 = sum(get_segment_lengths(m.f2, s, Lxx, Lxy, Lyx, Lyy));
        # Calculate motor positions
        x1, y1, x2, y2 = get_motor_pos(m, s, Lxx, Lxy, Lyx, Lyy);
        m1x = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx;
        m1y = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy;
        m2x = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx;
        m2y = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy;     
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s, Lxx, Lxy, Lyx, Lyy);
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s, m.t1, Lxx, Lxy, Lyx, Lyy);
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s, m.t2, Lxx, Lxy, Lyx, Lyy);
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1]; # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2]; # Vector between motor and plus end of filament 2
        cosine = (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2)); # Cosine of angle between vectors
        # Prevent spurious errors that occur in rare instances
        if abs(cosine) <= 1
            theta = acos(cosine); # Angle between vectors
        elseif cosine > 1
            theta = acos(1); # Angle between vectors
        else
            theta = acos(-1); # Angle between vectors
        end
        # Compute index
        if all(s.mp[i] .<=1)
            Id = ( 2*(L1*s.mp[i][1] + L2*s.mp[i][2])/(L1+L2) - 1)*(1-cos(theta/2)^2); # Heuristic index based on four filament branches
        else
            Id = 0; # Return zero if motor is detached
        end
        push!(index, Id);
    end
    # histogram(index)
    return mean(index)
end

"Actin filament speed"
function filament_speed(parN, af, s, s_old, Lxx, Lxy, Lyx, Lyy)
    node_speeds = Vector{Float64}(); # Pre-allocate vector of node speeds
    # Loop over filaments
    for f in af
        # Loop over nodes
        for i = 1:length(s.an[f.index])
            # Extract physical node position (new)
            xn = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx;
            yn = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy;
            # Extract physical node position (old)
            xo = s_old.an[f.index][i][1]*Lxx + s_old.an[f.index][i][2]*Lyx;
            yo = s_old.an[f.index][i][1]*Lxy + s_old.an[f.index][i][2]*Lyy;
            # Calculate speed
            speed = sqrt((xn-xo)^2 + (yn-yo)^2)/parN.dt;
            push!(node_speeds, speed)
        end
    end
    # histogram(node_speeds)
    return mean(node_speeds)
end

"Myosin motor speed"
function motor_speed(parN, mm, s, s_old, Lxx, Lxy, Lyx, Lyy)
    motor_speeds = Vector{Float64}(); # Pre-allocate vector of motor head speeds
    # Loop over motors
    for m in mm
        # Calculate (current) filament lengths
        L1 = sum(get_segment_lengths(m.f1, s, Lxx, Lxy, Lyx, Lyy));
        L2 = sum(get_segment_lengths(m.f2, s, Lxx, Lxy, Lyx, Lyy));
        # Calculate motor speeds
        s1 = L1*(s.mp[m.index][1]-s_old.mp[m.index][1])/parN.dt;
        s2 = L2*(s.mp[m.index][2]-s_old.mp[m.index][2])/parN.dt;
        push!(motor_speeds, s1); push!(motor_speeds, s2);
    end
    # histogram(motor_speeds)
    return mean(motor_speeds)
end

"Rate of change of angle between two filaments attached to a motor"
function motor_angle_roc(parN, mm, s, s_old, Lxx, Lxy, Lyx, Lyy)
    motor_angle_roc = Vector{Float64}();  # Pre-allocate vector of motor head speeds
    # Loop over motors
    for m in mm
        f1::Actin_Filament = m.f1; f2::Actin_Filament = m.f2; # Extract filaments attached to current motor
        # 1. Compute angles
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s, Lxx, Lxy, Lyx, Lyy);
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s, m.t1, Lxx, Lxy, Lyx, Lyy);
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s, m.t2, Lxx, Lxy, Lyx, Lyy);
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1]; # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2]; # Vector between motor and plus end of filament 2
        cosine = (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2)); # Cosine of angle between vectors
        # Prevent spurious errors that occur in rare instances
        if abs(cosine) <= 1
            theta = acos(cosine); # Angle between vectors
        elseif cosine > 1
            theta = acos(1); # Angle between vectors
        else
            theta = acos(-1); # Angle between vectors
        end
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s_old, Lxx, Lxy, Lyx, Lyy);
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s_old, m.t1, Lxx, Lxy, Lyx, Lyy);
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s_old, m.t2, Lxx, Lxy, Lyx, Lyy);
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1]; # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2]; # Vector between motor and plus end of filament 2
        cosine = (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2)); # Cosine of angle between vectors
        # Prevent spurious errors that occur in rare instances
        if abs(cosine) <= 1
            theta_old = acos(cosine); # Angle between vectors
        elseif cosine > 1
            theta_old = acos(1); # Angle between vectors
        else
            theta_old = acos(-1); # Angle between vectors
        end
        # 2. Calculate rate of change of angle
        rate = (theta - theta_old)/(parN.dt);
        push!(motor_angle_roc, rate);
    end
    # histogram(motor_angle_roc)
    return mean(motor_angle_roc)
end

"Myosin motor head positions"
function motor_position(mm, s)
    motor_pos = Vector{Float64}();  # Pre-allocate vector of motor head speeds
    # Loop over motors
    for m in mm
        push!(motor_pos, s.mp[m.index][1]);
        push!(motor_pos, s.mp[m.index][2]);
    end
    # histogram(motor_pos)
    return mean(motor_pos)
end

"Angle between two filaments attached to a motor"
function motor_angle(mm, s, parN, Lxx, Lxy, Lyx, Lyy)
    motor_angle = Vector{Float64}();  # Pre-allocate vector of motor head speeds
    # Loop over motors
    for m in mm
        f1::Actin_Filament = m.f1; f2::Actin_Filament = m.f2; # Extract filaments attached to current motor
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s, Lxx, Lxy, Lyx, Lyy);
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s, m.t1, Lxx, Lxy, Lyx, Lyy);
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s, m.t2, Lxx, Lxy, Lyx, Lyy);
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1]; # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2]; # Vector between motor and plus end of filament 2
        cosine = (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2)); # Cosine of angle between vectors
        # Prevent spurious errors that occur in rare instances
        if abs(cosine) <= 1
            theta = acos(cosine); # Angle between vectors
        elseif cosine > 1
            theta = acos(1); # Angle between vectors
        else
            theta = acos(-1); # Angle between vectors
        end
        push!(motor_angle, theta);
    end
    # histogram(motor_angle)
    return mean(motor_angle)
end

"Paired distances between nodes"
function pcf(af, s, Lxx, Lyy)
    distances = Vector{Float64}(); # Pre-allocate vector of all paired distances
    # Loop over all filament pairs
    for i = 1:length(af)
        for j = 1:length(af)
            if j > i # Prevent double counting
                # Loop over all nodes
                for k = 1:length(s.an[i])
                    n1x = Lxx*(s.an[i][k][1] - floor(s.an[i][k][1])); # In-domain node position on 1st filament (x)
                    n1y = Lyy*(s.an[i][k][2] - floor(s.an[i][k][2])); # In-domain node position on 1st filament (y)
                    for l = 1:length(s.an[j])
                        n2x = Lxx*(s.an[j][l][1] - floor(s.an[j][l][1])); # In-domain node position on 2nd filament (x)
                        n2y = Lyy*(s.an[j][l][2] - floor(s.an[j][l][2])); # In-domain node position on 2nd filament (y)
                        dist_candidates = Vector{Float64}(); # Pre-allocate possible shortest distances given periodic BCs
                        test = 0; #[-1, 0, 1]; # Translations to test for shortest distance
                        # Loop over translations
                        for t1x in test
                            for t1y in test
                                for t2x in test
                                    for t2y in test
                                        push!(dist_candidates, sqrt( (n1x - n2x - t1x*Lxx + t2x*Lxx)^2 + (n1y - n2y - t1y*Lyy + t2y*Lyy)^2 ));
                                    end
                                end
                            end
                        end
                        push!(distances, minimum(dist_candidates)); # Store shortest distance between nodes 
                    end
                end
            end
        end
    end
    # Generate normalised histogram
    nBins = 50;
    BinWidth = sqrt(Lxx^2 + Lyy^2)/nBins;
    Bins = Vector{Float64}(undef, nBins);
    Hist = Vector{Int}(undef, nBins); # Pre-allocate histogram
    for i = 1:nBins
        count::Int = 0;
        for d in distances
            if (d >= (i-1)*BinWidth) && (d < i*BinWidth)
                count += 1;
            end
        end
        Hist[i] = count;
        Bins[i] = i*BinWidth;
    end
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    default(guidefont = (18, "times"), tickfont = (12, "times"))
    bar(Bins, Hist/(sum(Hist)*BinWidth), label=:false, left_margin = 5mm, bottom_margin = 5mm, xlabel = L"d", ylabel = L"P(d)")
    x1::Vector{Float64} = 0.0:0.001:1;
    x2::Vector{Float64} = 1.0:(sqrt(2)-1)/1000:sqrt(2);
    p1 = 2*x1.*(x1.^2 .- 4*x1 .+ pi);
    p2 = 2*x2.*(4*sqrt.(x2.^2 .- 1) - (x2.^2 .+ 2 .- pi) - 4*atan.(sqrt.(x2.^2 .- 1)));
    plot!(2.5*x1, p1/2.5, linewidth=2, color = "black", label=:false);
    plot!(2.5*x2, p2/2.5, linewidth=2, color = "black", label=:false);
    default(titlefont = (18, "times"), guidefont = (26, "times"), tickfont = (18, "times"))
    return Bins, Hist
    # histogram(distances)
end