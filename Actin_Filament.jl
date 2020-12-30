# Data structures and methods for actin filaments
# Alex Tam, 12/10/2020

"Data structure for actin filament segments"
mutable struct Segment
    index::Int # Integer segment index
    L_eq::Float64 # [μm] Segment equilibrium length
    t::Vector{Vector{Int}} # Translations required by periodic BCs
end

"Mutable data structure for actin filaments"
mutable struct Actin_Filament
    index::Int # Integer filament index
    segments::Vector{Segment} # List of segments on actin filament
    "Inner constructor to generate actin filaments and update State"
    function Actin_Filament(s, ind, parA, Lxx, Lyy)
        # 1. Generate actin filament segments
        segments = Vector{Segment}(undef, parA.nSeg); # Pre-allocate
        for i = 1:parA.nSeg
            segments[i] = Segment(i, parA.LSeg, []); # Create segments
        end
        # 2. Generate random centre position and orientation
        cx = rand(); cy = rand(); # Dimensionless actin filament centre position
        angle = 2*pi*rand(); # Generate random actin filament orientation
        # 3. Calculate minus end position
        nodes = Vector{}(undef, parA.nSeg+1); # Pre-allocate
        Lf = parA.nSeg*parA.LSeg; # Total filament length
        nx = cx - 0.5*cos(angle)*Lf/Lxx; # Compute dimensionless x-position of minus end
        ny = cy - 0.5*sin(angle)*Lf/Lyy; # Compute dimensionless y-position of minus end
        nodes[1] = Vector{}(undef,2); # Pre-allocate
        nodes[1][1] = nx; nodes[1][2] = ny; # Store minus end position
        # 4. Calculate other node positions
        for i = 1:parA.nSeg
            nx += segments[i].L_eq*cos(angle)/Lxx; # Compute node x-position
            ny += segments[i].L_eq*sin(angle)/Lyy; # Compute node y-position
            nodes[i+1] = Vector{}(undef,2); # Pre-allocate
            nodes[i+1][1] = nx; nodes[i+1][2] = ny; # Store node positions
        end
        push!(s.an, nodes); # Store nodes in state
        # 5. Apply periodic BCs to each segment
        for i = 1:parA.nSeg
            segments[i].t = periodic(s.an[ind][i+1], s.an[ind][i]); # Store required translations
        end
        return new(ind, segments)
    end
end

"Obtain dimensional segment lengths"
function get_segment_lengths(f::Actin_Filament, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    Ls = Vector{T}(undef, length(f.segments)); # Pre-allocate segment lengths
    for seg in f.segments
        mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy); # Obtain segment nodes (dimensional, un-translated)
        Ls[seg.index] = sqrt((px - mx)^2 + (py - my)^2); # Compute segment length
    end
    return Ls
end

"Obtain dimensionless, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}) where {T}
    mx = s.an[f.index][seg.index][1];
    my = s.an[f.index][seg.index][2];
    px = s.an[f.index][seg.index+1][1];
    py = s.an[f.index][seg.index+1][2];
    return mx, my, px, py
end

"Obtain dimensionless, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}) where {T}
    mx = s.an[f.index][seg.index][1] - t[1];
    my = s.an[f.index][seg.index][2] - t[2];
    px = s.an[f.index][seg.index+1][1] - t[1];
    py = s.an[f.index][seg.index+1][2] - t[2];
    return mx, my, px, py
end

"Obtain dimensional, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = s.an[f.index][seg.index][1]*Lxx + s.an[f.index][seg.index][2]*Lyx;
    my = s.an[f.index][seg.index][1]*Lxy + s.an[f.index][seg.index][2]*Lyy;
    px = s.an[f.index][seg.index+1][1]*Lxx + s.an[f.index][seg.index+1][2]*Lyx;
    py = s.an[f.index][seg.index+1][1]*Lxy + s.an[f.index][seg.index+1][2]*Lyy;
    return mx, my, px, py
end

"Obtain dimensional, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}, parN, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = (s.an[f.index][seg.index][1] - t[1]*Lxx/parN.lxx - t[2]*Lyx/parN.lxx)*Lxx + (s.an[f.index][seg.index][2] - t[1]*Lxy/parN.lyy - t[2]*Lyy/parN.lyy)*Lyx;
    my = (s.an[f.index][seg.index][2] - t[1]*Lxy/parN.lyy - t[2]*Lyy/parN.lyy)*Lyy + (s.an[f.index][seg.index][1] - t[1]*Lxx/parN.lxx - t[2]*Lyx/parN.lxx)*Lxy;
    px = (s.an[f.index][seg.index+1][1] - t[1]*Lxx/parN.lxx - t[2]*Lyx/parN.lxx)*Lxx + (s.an[f.index][seg.index+1][2] - t[1]*Lxy/parN.lyy - t[2]*Lyy/parN.lyy)*Lyx;
    py = (s.an[f.index][seg.index+1][2] - t[1]*Lxy/parN.lyy - t[2]*Lyy/parN.lyy)*Lyy + (s.an[f.index][seg.index+1][1] - t[1]*Lxx/parN.lxx - t[2]*Lyx/parN.lxx)*Lxy;
    return mx, my, px, py
end