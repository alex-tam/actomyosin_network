# Spatial statistics for actomyosin networks
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
    histogram(filament_curvatures)
    return mean(filament_curvatures)
end

"Two-filament index"
function 2f_index(mm, s, parN, Lxx, Lxy, Lyx, Lyy)
    index = Vector{Float64}(); # Pre-allocate index data
    # Loop over motors
    for i = 1:length(mm)
        m::Myosin_Motor = mm[i]; # Extract current motor
        f1::Actin_Filament = m.f1; f2::Actin_Filament = m.f2; # Extract filaments attached to current motor
        # Calculate motor positions
        x1, y1, x2, y2 = get_motor_pos(m, s, Lxx, Lxy, Lyx, Lyy);
        m1x = (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxx + (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyx;
        m1y = (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyy + (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxy;
        m2x = (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxx + (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyx;
        m2y = (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyy + (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxy;       
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s, Lxx, Lxy, Lyx, Lyy);
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s, m.t1, parN, Lxx, Lxy, Lyx, Lyy);
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s, m.t2, parN, Lxx, Lxy, Lyx, Lyy);
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1]; # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2]; # Vector between motor and plus end of filament 2
        theta = acos( (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2) ) ); # Angle between vectors
        # Compute index
        if all(s.mp[i] .<=1)
            Id = (s.mp[i][1] + s.mp[i][2] - (1-s.mp[i][1]) - (1-s.mp[i][2]))*(1-cos(theta/2)); # Index based on four filament branches
        else
            Id = 0; # Return zero if motor is detached
        end
        push!(index, Id);
    end
    histogram(index)
    return mean(index)
end

"Pair correlation function"
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
                                        push!(dist_candidates, sqrt( (n1x-n2x-t1x*Lxx+t2x*Lxx)^2 + (n1y-n2y-t1y*Lyy+t2y*Lyy)^2 ));
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
    histogram(distances)
end