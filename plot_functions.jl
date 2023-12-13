"""
    ground_trace(orbp::OrbitPropagator, Δt::Number; kwargs...)
    Compute the ground trace of the object with orbit defined by `orbp` by
    propagating the orbit by `Δt` [s] from the orbit epoch.
    By default, it considers that the orbit elements on the propagator are
    represented in the True Equator, Mean Equinox (TEME) reference frame and the
    ground trace will be computed in the Pseudo-Earth Fixed (PEF) reference frame.
    Hence, no EOP data is needed. However, this can be changed by the keywords
    presented as follows.
    # Keywords
    * `eop_data`: EOP data that will be used to convert the ECI reference frame to
                the ECEF reference frame. If `nothing`, then it will not be used
                (see `r_eci_to_ecef`). (**Default** = `nothing`)
    * `ECI`: ECI frame in which the orbit elements in `orbp` are represented.
            (**Default** = `TEME()`)
    * `ECEF`: ECEF frame that will be used to compute the ground trace.
            (**Default** = `PEF()`)
    * `dt`: Time interval between two samples [s]. (**Default** = 10.0)
    # Returns
    A vector of tuples with the pairs `(latitude,longitude)` of the ground trace.
"""
function ground_trace(orbp::OrbitPropagator, Δt::Number;
                      eop_data::Union{Nothing, EOPData_IAU1980, EOPData_IAU2000A} = nothing,
                      ECI = TEME(),
                      ECEF = PEF(),
                      dt::Number = 10)

    # Copy orbit structure so that it is not modified by `propagate`.
    orbp_c = deepcopy(orbp)

    # Timespan of the analysis.
    t = 0:dt:Δt

    # Compute the points represented in the inertial reference frame.
    r_i, ~ = propagate!(orbp_c, t)

    # Get the epochs in Julian Day of each instant.
    JD = get_epoch(orbp) .+ t./86400

    # Convert from the ECI to the ECEF frame.
    if eop_data == nothing
        r_e = map( (t,v_i)->r_eci_to_ecef(ECI, ECEF, t)*v_i, JD, r_i )
    else
        r_e = map( (t,v_i)->r_eci_to_ecef(ECI, ECEF, t, eop_data)*v_i, JD, r_i )
    end

    # Convert to Geodetic coordinates.
    geod = ecef_to_geodetic.(r_e)
    return map(x->(x[1],x[2]), geod)
end

function plotground(coord::Array)
    gcf();
    
    # Obter os valores de Latitude e Longitude para ser plotado em graus.
    y = map(x->x[1][1], coord);
    x = map(x->x[2][1], coord);
    
    y = rad2deg.(y);
    x = rad2deg.(x);
    
    # Plotar.
    f = plot(x,y, linestyle ="none", marker =".", markersize=5)
    title("Ground Trace")
    xlabel("Longitude")
    xlim(-180, 180)
    xticks([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180])
    ylabel("Latitude")
    ylim(-90, 90)
    return f
end

function plotpcd(location)
    plot(rad2deg(location[2]), rad2deg(location[1]), linestyle ="none", color="r", marker="d");
    return nothing
end

function plotrec(location)
    plot(rad2deg(location[2]), rad2deg(location[1]), linestyle ="none", color="b", marker="p");
    return nothing
end

function plot_all(constellation, rec_list, pcd_list, t)
    clf();

    ### Plot satellites ground traces 
    for satellite in constellation
            initorbit = init_orbit_propagator(Val(:J2), satellite);
            coordenadas = ground_trace(initorbit, t/10) # /10 para melhorar visual
            plotground(coordenadas)
    end
    
    ### Plot PCD positions
    for pcd in pcd_list
        plotpcd(pcd);
    end
    
    ### Plot Receiving Stations positions
    for station in rec_list
        plotrec(station)
    end
    
    gcf()

end

# Plots a point P from the Pareto Frontier
function plot_p_chart(P, nsat::Number; nvars::Number = 4, tprop::Number=24*3600*2)
    
    constellation = Vector{KeplerianElements}();
    
    for i=1:nsat
        push!(constellation, KeplerianElements(date_to_jd(epoca.y, epoca.m, epoca.d, epoca.h, epoca.min, epoca.s),
                                                    P.vars[(i-1)*nvars+1]*1e3,       # semiaxis_major
                                                    0,              # eccentricity
                                                    P.vars[(i-1)*nvars+2],   # inclination 
                                                    P.vars[(i-1)*nvars+3],   # raan
                                                    0,              # true_anom
                                                    P.vars[(i-1)*nvars+4])); # arg_perigee
    end
    
    #temp_fig = plot_all(constellation, rec_list, pcd_list, tprop);
    plot_all(constellation, rec_list, pcd_list, tprop);

end

function print_data(pf)
    for point in pf
        println(point.f[1],"        ",point.f[2])
    end
end

function print_vars(pf)
    for point in pf
        println(point.vars)
    end
end