# Data structure and methods for myosin motors

"Mutable data structure for a myosin motor"
mutable struct Myosin_Motor
    index::Int # Integer motor index
    f1::Actin_Filament # 1st actin filament
    f2::Actin_Filament # 2nd actin filament
    t1::Vector{Int} # 1st filament translation
    t2::Vector{Int} # 2nd filament translation
    to_remove::Bool # Binary variable: true if motor is to be removed (filament turnover)
    "Inner constructor to generate new myosin motors and update State"
    function Myosin_Motor(s::State, l::Cross_Link, ind::Int, Lxx, Lxy, Lyx, Lyy)
        push!(s.mp, get_motor_relative_pos_filament(s, l, Lxx, Lxy, Lyx, Lyy)); # Write relative positions to State
        return new(ind, l.f1, l.f2, l.t1, l.t2, false) # Generate motor
    end
end

"Convert relative position along segments to relative position along filaments"
function get_motor_relative_pos_filament(s::State, l::Cross_Link, Lxx, Lxy, Lyx, Lyy)
    # Compute dimensional segment lengths
    Ls1 = get_segment_lengths(l.f1, s, Lxx, Lxy, Lyx, Lyy); # Filament 1
    Ls2 = get_segment_lengths(l.f2, s, Lxx, Lxy, Lyx, Lyy); # Filament 2
    # Compute relative positions
    mp = Vector{}(undef, 2); # Pre-allocate vector of relative positions
    mp[1] = (sum(Ls1[1:l.seg1_index]) - (1-l.relative_pos1)*Ls1[l.seg1_index])/sum(Ls1);
    mp[2] = (sum(Ls2[1:l.seg2_index]) - (1-l.relative_pos2)*Ls2[l.seg2_index])/sum(Ls2);
    return mp
end

"Convert relative position along filaments to relative position along segments"
function get_motor_relative_pos_segment(m::Myosin_Motor, s::State, Lxx, Lxy, Lyx, Lyy)
    # Steps 1-3 must be dimensional for non-square domains
    # 1. Calculate dimensional segment lengths
    seg_lengths_1 = get_segment_lengths(m.f1, s, Lxx, Lxy, Lyx, Lyy);
    seg_lengths_2 = get_segment_lengths(m.f2, s, Lxx, Lxy, Lyx, Lyy);
    # 2. Determine index of relevant segments
    cumulative_length_1 = cumsum(seg_lengths_1);
    cumulative_length_2 = cumsum(seg_lengths_2);
    scaled_mp_1 = s.mp[m.index][1]*cumulative_length_1[end];
    scaled_mp_2 = s.mp[m.index][2]*cumulative_length_2[end];
    myosin_seg_1::Int = 1; myosin_seg_2::Int = 1;
    for i = 1:(length(cumulative_length_1)-1)
        if scaled_mp_1 > cumulative_length_1[i]
            myosin_seg_1 += 1;
        end
    end
    for i = 1:(length(cumulative_length_2)-1)
        if scaled_mp_2 > cumulative_length_2[i]
            myosin_seg_2 += 1;
        end
    end
    # 3. Determine relative position on segment
    if myosin_seg_1 == 1
        myosin_seg_rel_pos_1 = scaled_mp_1/(cumulative_length_1[myosin_seg_1])
    else
        myosin_seg_rel_pos_1 = (scaled_mp_1-cumulative_length_1[myosin_seg_1-1])/(cumulative_length_1[myosin_seg_1]-cumulative_length_1[myosin_seg_1-1])
    end
    if myosin_seg_2 == 1
        myosin_seg_rel_pos_2 = scaled_mp_2/(cumulative_length_2[myosin_seg_2])
    else
        myosin_seg_rel_pos_2 = (scaled_mp_2-cumulative_length_2[myosin_seg_2-1])/(cumulative_length_2[myosin_seg_2]-cumulative_length_2[myosin_seg_2-1])
    end
    return myosin_seg_1, myosin_seg_2, myosin_seg_rel_pos_1, myosin_seg_rel_pos_2
end

"Obtain dimensional, translated motor position"
function get_motor_pos(m::Myosin_Motor, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    # Find relative position along segment
    seg1, seg2, relative_pos1, relative_pos2 = get_motor_relative_pos_segment(m::Myosin_Motor, s::State, Lxx, Lxy, Lyx, Lyy)
    # Obtain segment nodes position
    m1x, m1y, p1x, p1y = get_segment_nodes(m.f1, m.f1.segments[seg1], s);
    m2x, m2y, p2x, p2y = get_segment_nodes(m.f2, m.f2.segments[seg2], s);
    # Compute motor positions
    x1 = m1x + relative_pos1*(p1x - m1x);
    y1 = m1y + relative_pos1*(p1y - m1y);
    x2 = m2x + relative_pos2*(p2x - m2x);
    y2 = m2y + relative_pos2*(p2y - m2y);
    return x1, y1, x2, y2
end