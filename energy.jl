# Actomyosin network energy functional
# Alex Tam, 12/10/2020

"Energy functional"
function energy_functional(x::Vector{T}, s_old::State{Float64}, af, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy) where{T}
    s = build_state(x, af, mm); # Rebuild state from vector input
    energy = zero(T); # Pre-allocate energy
    for f in af
        energy += energy_actin_spring(f, s, parA, Lxx, Lxy, Lyx, Lyy);
        energy += energy_actin_drag(f, s, s_old, parN, parA, Lxx, Lxy, Lyx, Lyy);
        energy += energy_actin_bending(f, s, parA, Lxx, Lxy, Lyx, Lyy);
    end
    for l in xl
        energy += energy_cross_link(l, s, parN, parA, Lxx, Lxy, Lyx, Lyy)
    end
    for m in mm
        energy += energy_myosin_spring(m, s, parN, parM, Lxx, Lxy, Lyx, Lyy);
        energy += energy_myosin_actin(m, s, s_old, parN, parM, Lxx, Lxy, Lyx, Lyy);
    end
    return energy
end

"Energy contribution of actin filament spring forces"
function energy_actin_spring(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    for seg in f.segments
        mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy); # Physical, not translated
        energy += 0.5*parA.k*(sqrt((px - mx)^2 + (py - my)^2) - seg.L_eq)^2; # Actin spring energy
    end
    return energy
end

"Energy contribution of drag between actin filaments and the cytoplasm"
function energy_actin_drag(f::Actin_Filament, s::State{T}, s_old::State{Float64}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    seg_lengths = get_segment_lengths(f, s, Lxx, Lxy, Lyx, Lyy); # Dimensional
    # Apply drag at nodes, scaled by segment lengths
    for i = 1:length(s.an[f.index]) # Loop over nodes
        pnx = s.an[f.index][i][1]; pny = s.an[f.index][i][2]; # Current dimensionless node positions
        pox = (s_old.an[f.index][i][1])*(Lxx/parN.lxx) + (s_old.an[f.index][i][2])*(Lyx/parN.lxx); # Old dimensionless node position in new co-ordinates (x)
        poy = (s_old.an[f.index][i][1])*(Lxy/parN.lyy) + (s_old.an[f.index][i][2])*(Lyy/parN.lyy); # Old dimensionless node position in new co-ordinates (y)
        # Obtain average length of segments adjacent to node
        if (i == 1)
            avl = seg_lengths[i]/2; # Minus end
        elseif (i == length(s.an[f.index]))
            avl = seg_lengths[end]/2; # Plus end
        else
            avl = (seg_lengths[i-1] + seg_lengths[i])/2; # Interior nodes
        end
        energy += parA.lambda_a*avl*( ((pnx-pox)*Lxx + (pny-poy)*Lyx)^2 + ((pnx-pox)*Lxy + (pny-poy)*Lyy)^2 )/(2*parN.dt); # Add scaled contribution to overall drag
    end
    return energy
end

"Energy contribution of actin filament bending"
function energy_actin_bending(f::Actin_Filament, s::State{T}, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    seg_lengths = get_segment_lengths(f, s, Lxx, Lxy, Lyx, Lyy); # Dimensional
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
        L01 = seg_lengths[i]; L12 = seg_lengths[i+1]; # Extract segment lengths
        avl = (L01 + L12)/2; # Average length
        ddfx = (xp-xc)/L12 - (xc-xm)/L01; # Numerical second derivative (x)
        ddfy = (yp-yc)/L12 - (yc-ym)/L01; # Numerical second derivative (y)
        # Add contribution to energy
        energy += 0.5*parA.kappa*(ddfx^2 + ddfy^2)/avl; # Actin bending energy
    end
    return energy
end

"Energy contribution of cross-linker drag"
function energy_cross_link(l::Cross_Link, s::State{T}, parN, parA, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T);
    # Extract cross-linker position (dimensionless, un-translated)
    x1, y1, x2, y2 = get_xl_pos(l, s); 
    # Convert to dimensional, translated positions
    lx1 = (x1 - l.t1[1]*Lxx/parN.lxx - l.t1[2]*Lyx/parN.lxx)*Lxx + (y1 - l.t1[1]*Lxy/parN.lyy - l.t1[2]*Lyy/parN.lyy)*Lyx;
    ly1 = (y1 - l.t1[1]*Lxy/parN.lyy - l.t1[2]*Lyy/parN.lyy)*Lyy + (x1 - l.t1[1]*Lxx/parN.lxx - l.t1[2]*Lyx/parN.lxx)*Lxy;
    lx2 = (x2 - l.t2[1]*Lxx/parN.lxx - l.t2[2]*Lyx/parN.lxx)*Lxx + (y2 - l.t2[1]*Lxy/parN.lyy - l.t2[2]*Lyy/parN.lyy)*Lyx;
    ly2 = (y2 - l.t2[1]*Lxy/parN.lyy - l.t2[2]*Lyy/parN.lyy)*Lyy + (x2 - l.t2[1]*Lxx/parN.lxx - l.t2[2]*Lyx/parN.lxx)*Lxy;
    # Energy
    energy += parA.lambda_xl*((lx1 - lx2)^2 + (ly1 - ly2)^2)/(2*parN.dt); # Cross-linker drag energy
    return energy
end

"Energy contribution of myosin spring forces"
function energy_myosin_spring(m::Myosin_Motor, s::State{T}, parN, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T);
    # Extract motor binding sites (dimensionless, un-translated)
    x1, y1, x2, y2 = get_motor_pos(m, s, Lxx, Lxy, Lyx, Lyy); 
    # Convert to dimensional, translated positions
    mx1 = (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxx + (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyx;
    my1 = (y1 - m.t1[1]*Lxy/parN.lyy - m.t1[2]*Lyy/parN.lyy)*Lyy + (x1 - m.t1[1]*Lxx/parN.lxx - m.t1[2]*Lyx/parN.lxx)*Lxy;
    mx2 = (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxx + (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyx;
    my2 = (y2 - m.t2[1]*Lxy/parN.lyy - m.t2[2]*Lyy/parN.lyy)*Lyy + (x2 - m.t2[1]*Lxx/parN.lxx - m.t2[2]*Lyx/parN.lxx)*Lxy;
    # Energy contribution
    energy += 0.5*parM.k*((mx1 - mx2)^2 + (my1 - my2)^2); # Myosin spring energy
    return energy
end

"Energy contribution of actin--myosin interactions for a single motor"
function energy_myosin_actin(m::Myosin_Motor, s::State{T}, s_old::State{Float64}, parN, parM, Lxx, Lxy, Lyx, Lyy) where {T}
    energy = zero(T); # Pre-allocate energy
    # Filament lengths, used to scale velocity
    L1 = sum(get_segment_lengths(m.f1, s, Lxx, Lxy, Lyx, Lyy));
    L2 = sum(get_segment_lengths(m.f2, s, Lxx, Lxy, Lyx, Lyy));
    # Linear component
    energy += parM.Fs/parM.Vm*( L1*(s.mp[m.index][1] - s_old.mp[m.index][1] ))^2/( 2*parN.dt );
    energy += parM.Fs/parM.Vm*( L2*(s.mp[m.index][2] - s_old.mp[m.index][2] ))^2/( 2*parN.dt );
    # Constant component
    energy -= parM.Fs*L1*(s.mp[m.index][1] - s_old.mp[m.index][1]);
    energy -= parM.Fs*L2*(s.mp[m.index][2] - s_old.mp[m.index][2]);
    return energy
end