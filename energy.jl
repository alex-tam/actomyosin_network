# Actomyosin network energy functional
# Alex Tam, 12/10/2020

"Energy functional"
function energy_functional(x::Vector{T}, s_old::State{Float64}, af, xl, mm, random, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) where{T}
    s = build_state(x, af, mm); # Rebuild state from vector input
    energy = zero(T); # Pre-allocate energy
    for f in af
        energy += energy_actin_thermal(f, s, random, parN, parA, Lxx, Lxy, Lyx, Lyy);
        energy += energy_actin_spring(f, s, parA, Lxx, Lxy, Lyx, Lyy);
        energy += energy_actin_drag(f, s, s_old, parN, parA, Lxx, Lxy, Lyx, Lyy);
        energy += energy_actin_bending(f, s, parA, Lxx, Lxy, Lyx, Lyy);
    end
    for l in xl
        energy += energy_cross_link(l, s, parN, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for m in mm
        energy += energy_myosin_spring(m, s, parM, Lxx, Lxy, Lyx, Lyy);
        energy += energy_myosin_actin(m, s, s_old, parN, parM, Lxx, Lxy, Lyx, Lyy);
    end
    return energy
end

"Energy contribution of thermal forces on actin filaments"
function energy_actin_thermal(f::Actin_Filament, s::State{T}, random, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    seg_lengths = get_segment_lengths(f, s, Lxx, Lxy, Lyx, Lyy); # Dimensional
    for i = 1:length(s.an[f.index])
        if i == 1
            avl = seg_lengths[i]/2; # Minus end
        elseif i == length(s.an[f.index])
            avl = seg_lengths[end]/2; # Plus end
        else
            avl = (seg_lengths[i-1] + seg_lengths[i])/2; # Interior nodes
        end
        # px = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx; 
        # py = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy; # Dimensional positions
        px = s.an[f.index][i][1]*Lxx; 
        py = s.an[f.index][i][2]*Lyy; # Dimensional positions
        energy -= sqrt(2*parA.kb*parA.T*parA.lambda_a*avl/parN.dt)*( px*random[f.index][i][1] + py*random[f.index][i][2] ); # Thermal energy contribution
    end
    return energy
end

"Energy contribution of actin filament spring forces"
function energy_actin_spring(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    # Loop over segments
    for seg in f.segments
        Ls = get_segment_length(f, s, seg, Lxx, Lxy, Lyx, Lyy); # Physical, not translated
        energy += 0.5*parA.k*(Ls - seg.L_eq)^2; # Actin spring energy
    end
    return energy
end

"Energy contribution of drag between actin filaments and the cytoplasm"
function energy_actin_drag(f::Actin_Filament, s::State{T}, s_old::State{Float64}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    # Distribute drag along segments
    nPoints = 5; # Number of data points to sample per segment
    for i = 1:(length(s.an[f.index])-1) # Loop over segments
        Ls = get_segment_length(f, s, f.segments[i], Lxx, Lxy, Lyx, Lyy);
        avl = Ls/nPoints;
        for j = 1:nPoints
            # Extract data points
            pnx = s.an[f.index][i][1] + (j-0.5)/nPoints*(s.an[f.index][i+1][1]-s.an[f.index][i][1]);
            pny = s.an[f.index][i][2] + (j-0.5)/nPoints*(s.an[f.index][i+1][2]-s.an[f.index][i][2]); # Current dimensionless node positions
            pox = s_old.an[f.index][i][1] + (j-0.5)/nPoints*(s_old.an[f.index][i+1][1]-s_old.an[f.index][i][1]);
            poy = s_old.an[f.index][i][2] + (j-0.5)/nPoints*(s_old.an[f.index][i+1][2]-s_old.an[f.index][i][2]); # Old dimensionless node positions
            # Add energy contirbution
            energy += parA.lambda_a*avl*( ((pnx-pox)*Lxx + (pny-poy)*Lyx)^2 + ((pnx-pox)*Lxy + (pny-poy)*Lyy)^2 )/(2*parN.dt);
        end
    end
    return energy
end

"Energy contribution of actin filament bending"
function energy_actin_bending(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    # Apply bending energy at interior nodes
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
        L01 = get_segment_length(f, s, f.segments[i], Lxx, Lxy, Lyx, Lyy); # Left segment length
        L12 = get_segment_length(f, s, f.segments[i+1], Lxx, Lxy, Lyx, Lyy); # Right segment length
        avl = (L01 + L12)/2; # Average length
        ddfx = ((xp-xc)/L12 - (xc-xm)/L01)/avl; # Numerical second derivative (x)
        ddfy = ((yp-yc)/L12 - (yc-ym)/L01)/avl; # Numerical second derivative (y)
        # Add contribution to energy
        energy += 0.5*parA.kappa*avl*(ddfx^2 + ddfy^2); # Actin bending energy
    end
    return energy
end

"Energy contribution of cross-linker drag"
function energy_cross_link(l::Cross_Link, s::State{T}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T);
    # Extract cross-linker position (dimensionless, un-translated)
    x1, y1, x2, y2 = get_xl_pos(l, s); 
    # Convert to dimensional, translated positions
    lx1 = (x1 - l.t1[1])*Lxx + (y1 - l.t1[2])*Lyx;
    ly1 = (x1 - l.t1[1])*Lxy + (y1 - l.t1[2])*Lyy;
    lx2 = (x2 - l.t2[1])*Lxx + (y2 - l.t2[2])*Lyx;
    ly2 = (x2 - l.t2[1])*Lxy + (y2 - l.t2[2])*Lyy;
    # Energy
    energy += parA.lambda_pf*((lx1 - lx2)^2 + (ly1 - ly2)^2)/(2*parN.dt); # Cross-linker drag energy
    return energy
end

"Energy contribution of myosin spring forces"
function energy_myosin_spring(m::Myosin_Motor, s::State{T}, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T);
    # Extract motor binding sites (dimensionless, un-translated)
    x1, y1, x2, y2 = get_motor_pos(m, s, Lxx, Lxy, Lyx, Lyy); 
    # Convert to dimensional, translated positions
    mx1 = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx;
    my1 = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy;
    mx2 = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx;
    my2 = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy;
    # Energy contribution
    energy += 0.5*parM.k*((mx1 - mx2)^2 + (my1 - my2)^2); # Myosin spring energy
    return energy
end

"Energy contribution of actin--myosin interactions for a single motor"
function energy_myosin_actin(m::Myosin_Motor, s::State{T}, s_old::State{Float64}, parN, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    # Filament lengths, used to scale velocity
    L1 = zero(T); L2 = zero(T); # Pre-allocate filament lengths
    for seg in m.f1.segments
        L1 += get_segment_length(m.f1, s, seg, Lxx, Lxy, Lyx, Lyy)
    end
    for seg in m.f1.segments
        L2 += get_segment_length(m.f2, s, seg, Lxx, Lxy, Lyx, Lyy)
    end
    # Linear component
    energy += parM.Fs/parM.Vm*( L1*(s.mp[m.index][1] - s_old.mp[m.index][1] ))^2/( 2*parN.dt );
    energy += parM.Fs/parM.Vm*( L2*(s.mp[m.index][2] - s_old.mp[m.index][2] ))^2/( 2*parN.dt );
    # Constant component
    energy -= parM.Fs*L1*(s.mp[m.index][1] - s_old.mp[m.index][1]);
    energy -= parM.Fs*L2*(s.mp[m.index][2] - s_old.mp[m.index][2]);
    return energy
end