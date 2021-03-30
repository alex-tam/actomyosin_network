# Implement random turnover of network components
# Alex Tam, 12/10/2020

"Random turnover of actin filaments"
function actin_turnover(s, af, mm, parN, parA, Lxx, Lyy)
    N = length(af); j = 0; # Indexing pre-allocation
    Pt = 1 - exp(-parA.k_off*parN.dt); # Probability that a filament will turn over in the time step
    # Remove actin filaments
    for i = 1:N
        j += 1; # Index (dynamic) of current filament
        random = rand(); # Random number [0, 1]
        if random <= Pt
            # Remove data from state and filament list
            deleteat!(s.an, j)
            deleteat!(af, j)
            # Update motors
            for m in mm
                # Mark motors for removal if filament turns over
                if m.f1.index == j
                    m.to_remove = true;
                end
                if m.f2.index == j
                    m.to_remove = true;
                end
                # Update index of filaments associated with motors
                if m.f1.index > j
                    m.f1.index -= 1; # Reduce index by 1 if necessary
                end
                if m.f2.index > j
                    m.f2.index -= 1; # Reduce index by 1 if necessary
                end
            end
            # Update indices in actin filament list
            for a = 1:length(af)
                af[a].index = a;
            end
            j -= 1; # Ensure correct index for next filament
        end
    end
    # Add new actin filaments
    nNew = parN.nA - length(af); # Immediately replace filaments
    # nNew = rand(Poisson(parN.nA*parA.k_off*parN.dt)); # Draw number of filaments to add from Poisson distribution
    nOld = length(af); # Initial number of filaments for indexing
    for i = 1:nNew
        push!(af, Actin_Filament(s, nOld+i, parA, Lxx, Lyy)); # Create new filament and state data
    end
    return af, mm, s
end

"Random turnover of myosin motors"
function myosin_turnover(s, xl, mm, parN, parM, Lxx, Lxy, Lyx, Lyy)
    N = length(mm); j = 0;
    for i = 1:N
        j += 1; # Index (dynamic) of current motor
        random = 1; Pt = 0; # Pre-allocate
        if mm[j].to_remove == false
            random = rand(); # Generate random number [0, 1] for current motor
            x1, y1, x2, y2 = get_motor_pos(mm[j], s, Lxx, Lxy, Lyx, Lyy); # Extract dimensional motor positions
            Lm = sqrt((x1-x2)^2 + (y1-y2)^2); # Compute motor length
            off_rate = parM.k_off*exp(parM.k*Lm/parM.F_ref); # Compute off rate of current motor using Bell's Law
            Pt = 1 - exp(-off_rate*parN.dt); # Probability that random turnover will occur at least once in the time step
        end
        # Replace motors
        if any(s.mp[j] .>= 1) || (random <= Pt) || (mm[j].to_remove == true)
            xl_index = rand(1:length(xl));
            # Delete data from State and motor list
            deleteat!(s.mp, j);
            deleteat!(mm, j);
            # Generate replacement motor
            push!(mm, Myosin_Motor(s, xl[xl_index], j, Lxx, Lxy, Lyx, Lyy));
            deleteat!(xl, xl_index);
            # Update order of motor indices
            for m = 1:N
                mm[m].index = m;
            end
            j -= 1; # Ensure correct index for next motor
        end
    end
    return mm, xl, s
end