# Locate intersections between actin filaments
# Alex Tam, 12/10/2020

"Intersection search algorithm"
function intersection_search(s, af, mm)
    xl = Vector{Cross_Link}(); # Pre-allocate empty vector of cross-links
    # Loop over filament pairs
    for f1 in af
        for f2 in af
            if f1.index < f2.index # Avoid double couting
                # Loop over segments
                for seg1 in f1.segments
                    for seg2 in f2.segments
                        # Loop over translations
                        for i = 1:length(seg1.t)
                            for j = 1:length(seg2.t)
                                unoccupied = true; # Boolean for whether link contains a motor
                                # Get dimensionless, translated node positions
                                m1x, m1y, p1x, p1y = get_segment_nodes(f1, seg1, s, seg1.t[i]);
                                m2x, m2y, p2x, p2y = get_segment_nodes(f2, seg2, s, seg2.t[j]);
                                # Run function to check whether intersection is possible
                                int = intersection(m1x, m1y, p1x, p1y, m2x, m2y, p2x, p2y); 
                                # If so, compute intersection between segments
                                if int == true
                                    rel = [ (p1x - m1x) (m2x - p2x) ; (p1y - m1y) (m2y - p2y) ]\[ m2x - m1x ; m2y - m1y ]; # Solve linear system for relative positions
                                    ix = m1x + rel[1]*(p1x - m1x); # Intersection x-position
                                    iy = m1y + rel[2]*(p1y - m1y); # Intersection y-position
                                    # Check for a motor
                                    for m in mm
                                        if ((m.f1.index == f1.index) && (m.f2.index == f2.index)) || ((m.f1.index == f2.index) && (m.f2.index == f1.index))
                                            unoccupied = false;
                                        end
                                    end
                                    # Store intersection as a cross-link if unoccupied and in domain
                                    if (unoccupied == true) && (ix >= 0) && (ix <= 1) && (iy >= 0) && (iy <= 1) && all(rel .>= 0) && all(rel .<= 1)
                                        push!(xl, Cross_Link(f1, f2, seg1.index, seg2.index, rel[1], rel[2], seg1.t[i], seg2.t[j]));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return xl
end

"Check whether intersection of two line segments will occur"
function intersection(m1x, m1y, p1x, p1y, m2x, m2y, p2x, p2y)
    # Check whether segments overlap in x and y-directions
    x_overlap = ((m2x >= m1x) && (m2x <= p1x)) || ((m2x >= p1x) && (m2x <= m1x)) || ((p2x >= m1x) && (p2x <= p1x)) || ((p2x >= p1x) && (p2x <= m1x)) || ((m1x >= m2x) && (m1x <= p2x)) || ((m1x >= p2x) && (m1x <= m2x)) || ((p1x >= m2x) && (p1x <= p2x)) || ((p1x >= p2x) && (p1x <= m2x));
    y_overlap = ((m2y >= m1y) && (m2y <= p1y)) || ((m2y >= p1y) && (m2y <= m1y)) || ((p2y >= m1y) && (p2y <= p1y)) || ((p2y >= p1y) && (p2y <= m1y)) || ((m1y >= m2y) && (m1y <= p2y)) || ((m1y >= p2y) && (m1y <= m2y)) || ((p1y >= m2y) && (p1y <= p2y)) || ((p1y >= p2y) && (p1y <= m2y));
    if (x_overlap == true) && (y_overlap == true)
        int = true; # Return true if intersection is possible
    else
        int = false; # Return false if no intersection can occur
    end
    return int
end