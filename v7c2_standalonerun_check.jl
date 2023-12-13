using SatelliteToolbox
using Dates
using PyPlot
using Statistics
using MGEO
using JLD2
using LinearAlgebra

include("plot_functions.jl")


# function gettimepos(v::Number, T::Array)
#     pos = findfirst(isequal(v), T)
#     return pos
# end

 function get_numbits(sup::Number, inf::Number, res::Number)
     nmin = log((sup-inf)/res+1)/log(2);
     n = ceil(Int, nmin);
     return n;
 end

gs_reception_accessesm(orbp, rec_wgs84::Tuple, gs_wgs84::Tuple, vargs...; kwargs...) =
    gs_reception_accessesm(orbp, [rec_wgs84], [gs_wgs84], vargs...; kwargs...)

function gs_reception_accessesm(orbp, vrec_wgs84::AbstractVector{T}, vgs_wgs84::AbstractVector{T}, vargs...;
                             kwargs...) where T<:Tuple

    vrs_e = [geodetic_to_ecef(gs_wgs84...) for gs_wgs84 in vgs_wgs84]
    vrec_e = [geodetic_to_ecef(rec_wgs84...) for rec_wgs84 in vrec_wgs84]
    return gs_reception_accessesm(orbp, vrec_e, vrs_e, vargs...; kwargs...)
end

gs_reception_accessesm(orbp, rec_e::AbstractVector{T}, rs_e::AbstractVector{T}, vargs...; kwargs...) where
    T<:Number = gs_reception_accessesm(orbp, [rec_e], [rs_e], vargs...; kwargs...)

function gs_reception_accessesm(orbp, vrec_e::AbstractVector{T}, vrs_e::AbstractVector{T}, Δt::Number,
                             ECI::Union{T_ECIs, T_ECIs_IAU_2006},
                             ECEF::Union{T_ECEFs, T_ECEFs_IAU_2006},
                             position_vector::Array{Any},
                             alfa::Number, vargs...;
                             θ1::Number = 10*pi/180,
                             θ2::Number = 20*pi/180,
                             reduction::Function = v->|(v...),
                             step::Number = 60) where T<:AbstractVector

    # Time vector of the analysis.
    t = 0:step:Δt;

    # Get the epoch of the propagator.
    JD₀ = get_epoch(orbp);

    # Matrix that will contain the accesses.
    accesses = Matrix{DateTime}(undef,0,2);

    # State to help the computation.
    state = :initial

    # Lambda function to check if the ground station and
    # any reception station are visible.
    f(t)::Bool = begin
        #r_i,v_i  = propagate!(orbp, t);
        #r_e      = r_eci_to_ecef(DCM, ECI, ECEF, JD₀+t/86400, vargs...)*r_i;  -> calculado na f propagate_all.
        ntimer = 1+div(t, step);
        r_e = position_vector[ntimer, alfa];

        visibility_rec = [ground_station_visible(r_e, vec_e, θ1) for vec_e in vrec_e]
        visibility_rec = reduction(visibility_rec);
        
        if visibility_rec
            visibility_gs = [ground_station_visible(r_e, rs_e, θ2) for rs_e in vrs_e]
            visibility_gs = reduction(visibility_gs);
        else   
            visibility_gs = 0;
        end

        return  visibility_rec && visibility_gs
    end

    f2(t)::Bool = begin
        r_i,v_i  = propagate!(orbp, t);
        r_e      = r_eci_to_ecef(DCM, ECI, ECEF, JD₀+t/86400, vargs...)*r_i; # -> calculado na f propagate_all.
        #ntimer = 1+div(t, step);
        #r_e = position_vector[ntimer, alfa];

        visibility_rec = [ground_station_visible(r_e, vec_e, θ1) for vec_e in vrec_e]
        visibility_rec = reduction(visibility_rec);
        
        if visibility_rec
            visibility_gs = [ground_station_visible(r_e, rs_e, θ2) for rs_e in vrs_e]
            visibility_gs = reduction(visibility_gs);
        else   
            visibility_gs = 0;
        end

        return  visibility_rec && visibility_gs
    end



    access_beg = DateTime(now())
    access_end = DateTime(now())

    for k in t
        # Check if the ground station and reception station are visible.
        gs_visible = f(k)

        # Handle the initial case.
        if state == :initial
            if gs_visible
                access_beg = jd_to_date(DateTime, JD₀)
                state = :visible
            else
                state = :not_visible
            end
        # Handle transitions.
        elseif (state == :not_visible) && gs_visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing2(f2, k₀, k₁, false, true; Δ = 1e-4, max = 1000)
            #kc = (k - step)/2;

            state = :visible
            access_beg = jd_to_date(DateTime, JD₀ + kc/86400)

        elseif (state == :visible) && !gs_visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing2(f2, k₀, k₁, true, false; Δ = 1e-4, max = 1000)
            #kc = (k - step)/2;

            state = :not_visible
            access_end = jd_to_date(DateTime, JD₀ + kc/86400)

            accesses = vcat(accesses, [access_beg access_end])
        end
    end
    
    # If the analysis finished during an access, then just add the end of the
    # interval as the end of the access.
    if state == :visible
        access_end = jd_to_date(DateTime, JD₀ + Δt/86400)
        accesses   = vcat(accesses, [access_beg access_end])
    end

    return accesses
end


function propagate_all(constellation::Vector{KeplerianElements},
                        t::Number,
                        ECI::Union{T_ECIs, T_ECIs_IAU_2006},
                        ECEF::Union{T_ECEFs, T_ECEFs_IAU_2006};
                        step::Number = 60)

    tvector = 0:step:t;
    temp = Matrix(undef, length(tvector), length(constellation));

    for i=1:length(constellation)
        initorbit = init_orbit_propagator(Val(:J2), constellation[i]);
        r_i, v_i  = propagate!(initorbit, tvector);
        JD₀ = get_epoch(initorbit);
        r_e = map((x, y) -> r_eci_to_ecef(DCM, ECI, ECEF, JD₀+y/86400)*x, r_i, tvector);
        temp[:,i] = r_e;
    end 

    return temp;
end


function find_crossing2(f::Function, t₀::Number, t₁::Number, s₀, s₁, vargs...;
                       Δ = 1e-3, max = 100)
    it = 0

    T = typeof( (t₁ + t₀)/2 )

    while it <= max
        # Call the function at the middle of the interval.
        ti = (t₁ + t₀)/2
        si = f(ti, vargs...)

        # Compute the new interval.
        if si == s₀
            t₀ = ti
        elseif si == s₁
            t₁ = ti
        else
            error("The function `f` returned an unexpected state.")
        end

        # If the interval is small enough, then just return.
        (t₁ - t₀) < Δ && break

        it += 1
    end

    return T(t₁)
end

function start_data(numsat::Number)
    ### List of the Receiving Stations [Lat, Lon, H]
    rec_list = [(deg2rad(-02.338541), deg2rad(-44.404102), 058.81),   # Alcantara, RN
                (deg2rad(-15.555008), deg2rad(-56.069569), 237.00)];  # Cuiaba, MT


   
    # pcd_list = [(deg2rad(-07.256111), deg2rad(-49.172778), 265.0),  # S. Fe do Araguaia - TO
    #             (deg2rad(-22.508111), deg2rad(-45.023583), 626.0),  # Cruzeiro - SP
    #             (deg2rad(-22.412222), deg2rad(-42.798889), 871.0),  # Teresopolis - RJ
    #             (deg2rad(-23.464611), deg2rad(-50.342694), 743.0),  # Ribeirao do Pinhal - PR
    #             (deg2rad(-02.812111), deg2rad(-54.298694), 084.0),  # Curua Mantante - PA
    #             (deg2rad(-07.879389), deg2rad(-40.091889), 000.0),  # Ouricuri - PE
    #             (deg2rad(-30.552778), deg2rad(-52.406667), 420.0),  # Encruzilhada do Sul - RS
    #             (deg2rad(-02.593861), deg2rad(-44.211222), 062.0),  # Sao Luis - MA
    #             (deg2rad(-14.886111), deg2rad(-40.800833), 762.0),  # Vitoria da Conquista - BA
    #             (deg2rad(-10.789167), deg2rad(-65.280556), 150.0),  # Guajara Mirim - RO
    #             (deg2rad(-10.505833), deg2rad(-37.054444), 163.0),  # Capela - SE
    #             (deg2rad(-06.121944), deg2rad(-38.302222), 000.0),  # Lucrecia - RN
    #             (deg2rad(-05.421250), deg2rad(-40.481361), 343.0),  # Independencia - CE
    #             (deg2rad(-16.315361), deg2rad(-43.711250), 000.0)]; # Faz Analina - MG
    
    ### List of the Data Collection Platform - PCD 
    pcd_list = [(deg2rad(005.271944), deg2rad(-60.2125), 1465.0),
                (deg2rad(-33.751900), deg2rad(-53.3972), 0000.0),
                (deg2rad(-03.854600), deg2rad(-32.4234), 0000.0),
                (deg2rad(-07.155000), deg2rad(-34.7928), 0000.0),
                (deg2rad(-07.535830), deg2rad(-73.9906), 0193.0)];
    

    ### Design Restrictions
    a = (max = 8380, min = 6880);       # [km]  semi-major axis;
    #e = (max = 0.01, min = 0.0);       # [-]   eccentricity;
    i = (max = deg2rad(40), min = 0);  # [deg] inclination;
    Ω = (max = deg2rad(360), min = 0);          # [deg] RAAN;
    #ω = (max = 360, min = 0);          # [deg] argument of perigee;
    𝑓𝑜 = (max = deg2rad(360), min = 0); # [deg] initial true anomaly

    #      AbstractVector{Number}([a.min, e.min, i.min, Ω.min, ω.min, 𝑓𝑜.min]);
    vmin = AbstractVector{Number}([a.min, i.min, Ω.min, 𝑓𝑜.min]);
    vmax = AbstractVector{Number}([a.max, i.max, Ω.max, 𝑓𝑜.max]);

    #= bit resolution
    a = 8;      # 10km resolution for semi-major between 8380-6880 km
    e = 0;      # constant = 0
    i = 7;      # resolution 1° 
    Ω = 8;      # resolution 2°
    ω = 0;      # constant = 0      
    𝑓𝑜 = 8;     # resolution 2°  =#

    nbits = AbstractVector{Integer}([8, 7, 8, 8]);
    nsat = convert(Int64, numsat);   # number of satellites in the constellation

    var_min, var_max, var_names, var_bits = config_sat_vars(nsat, vmin, vmax, nbits);

    #theta1 = deg2rad(10);   # [rad]
    #theta2 = deg2rad(20);   # [rad]

    epoca = (y=2020,m=1,d=1,h=00,min=00,s=00);

    t = 20*24*3600; # [s] seconds in 20 days.
    tstep = 60; # [s]

    pECI = TEME();
    pECEF = PEF();

    #return t, epoca, theta1, theta2, pcd_list, rec_list, var_min, var_max, var_names, var_bits, nsat, pECI, pECEF, tstep;
    return t, epoca, pcd_list, rec_list, var_min, var_max, var_names, var_bits, nsat, pECI, pECEF, tstep;
end

# must be adjusted in more keplerian elements to be used as variables
function config_sat_vars(nsat::Integer,
                        vmin::AbstractVector{Number},
                        vmax::AbstractVector{Number},
                        nbits::AbstractVector{Integer})
    
    var_min = AbstractVector{Number}[];
    var_max = AbstractVector{Number}[];
    var_names = AbstractVector{String}[];
    var_bits = AbstractVector{Integer}[];
    for i=1:nsat
        var_min = vcat(var_min, vmin);
        var_max = vcat(var_max, vmax);
        var_names = vcat(var_names, ["Sat_$i a",  
                                     #"Sat_$i e",
                                     "Sat_$i i", 
                                     "Sat_$i Ω",
                                     #"Sat_$i ω", 
                                     "Sat_$i 𝑓𝑜"]);
        var_bits = vcat(var_bits, nbits);
    end

    var_min = convert(AbstractVector{Number}, var_min);
    var_max = convert(AbstractVector{Number}, var_max);
    var_names = convert(AbstractVector{String}, var_names);
    var_bits = convert(AbstractVector{Integer}, var_bits);

    return var_min, var_max, var_names, var_bits;
end

# function gen_sats_cons(vars, epoca::NamedTuple, num_sats::Integer)
#     #num_sats = convert(Integer, length(vars)/3); # sempre verificar length(vars)/NUM e vars[(i-1)*NUM+1ou2]
#     nvars = convert(Int64, convert(Integer, length(vars))/num_sats);
#     constellation = Vector{KeplerianElements}(undef, num_sats);
#     for i=1:num_sats
#         constellation[i] = KeplerianElements(date_to_jd(epoca.y, epoca.m, epoca.d, epoca.h, epoca.min, epoca.s),
#                                             vars[(i-1)*nvars+1]*1e3,    # a
#                                             0,                          # e
#                                             vars[(i-1)*nvars+2],        # i
#                                             vars[(i-1)*nvars+3],        # Ω
#                                             0,                          # ω
#                                             vars[(i-1)*nvars+4]);       # 𝑓𝑜   
#     end
#     return constellation;
# end

function standalone_const(epoca::NamedTuple, altitude::Number) ## standalone constellation for test cases
    
    constellation = Vector{KeplerianElements}();
    semiaxis = 6380+altitude;
    inclination = 30;
    #raan = [0 120 240];
    raan = [0 120 240]
    f_oi = [0 180];
    for point1 in raan
        for point2 in f_oi
            push!(constellation, KeplerianElements(date_to_jd(epoca.y, epoca.m, epoca.d, epoca.h, epoca.min, epoca.s),
                                                semiaxis*1e3,       # semiaxis_major
                                                0,                  # eccentricity
                                                inclination*pi/180, # inclination 
                                                deg2rad(point1),    # raan
                                                0,                  # arg_perigee
                                                deg2rad(point2)));  # true_anom
    
        end
    end
    return constellation;
end

function standalone_const2(epoca::NamedTuple, altitude::Number) ## standalone constellation for test cases
    
    constellation = Vector{KeplerianElements}();
    fo = pi/180 .*[0, 90, 180, 270];
    raan = pi/180 .* [0, 120, 240]
    inclination = 30*pi/180;
    #altitude = 750; #km
    for plane in raan
        for start in fo
            push!(constellation, KeplerianElements(date_to_jd(epoca.y, epoca.m, epoca.d, epoca.h, epoca.min, epoca.s),
                                                (altitude+6380)*1e3,       # semiaxis_major
                                                0,              # eccentricity
                                                inclination,    # inclination 
                                                plane,          # raan
                                                0,              # arg_perigee
                                                start));        # true_anom

        end
    end
    return constellation;
end


function satcon_optm(constellation::Vector{KeplerianElements},
                    t::Number,
                    pcd_list,
                    rec_list,
                    pECI::Union{T_ECIs, T_ECIs_IAU_2006},
                    pECEF::Union{T_ECEFs, T_ECEFs_IAU_2006};
                    theta1::Number = deg2rad(5),
                    theta2::Number = deg2rad(5),
                    tstep::Number = 60)

    valid = true;
    JD₀ = 0;
    tgap = 0;
    missed = false;
    pos_vector = propagate_all(constellation, t, pECI, pECEF);
    for pcd in pcd_list
        access_vector = Array{DateTime}(undef,0,2);
        gap_vector = Array{DateTime}(undef, 0, 2);
        for i=1:length(constellation)
            initorbit = init_orbit_propagator(Val(:J2), constellation[i]);
            JD₀ = get_epoch(initorbit);
            temp_gap = gs_reception_accessesm(initorbit, rec_list, [pcd], t, pECI, pECEF, pos_vector, i; θ1 = theta1, θ2 = theta2, step = tstep);
            if !isnothing(temp_gap) && !isempty(temp_gap)
                access_vector = vcat(access_vector, temp_gap);
            end
        end
        if !isnothing(access_vector) && !isempty(access_vector)
            access_vector = mergeIntervals(access_vector);
            gap_vector = gs_gaps(access_vector, JD₀, t);
            # @save "gap_vector.jld2" gap_vector
            tgap = max(gap_calc2(gap_vector), tgap);           
        else
            missed = true;
        end
             
    end
    if missed
        valid = false
    end

    return tgap, valid
end

function gap_calc2(gaps::Array{DateTime})
    max_gap = 0;
    ts = 60; # timescale in minutes
    gap_time = map((b,e)->(e-b).value/1000/ts, @view(gaps[:,1]), @view(gaps[:,2]));
    try
        max_gap, id_max = findmax(gap_time);
    catch
    end
    # if !isempty(gap_time)
    #     max_gap, id_max = findmax(gap_time);
    # end

    return max_gap;
end

function mean_a(constellation)
    sum = 0;
    for sat in constellation
        sum = sum + sat.a;
    end
    mean_a = sum/length(constellation);
    return mean_a;
end

function mergeIntervals(arr::AbstractArray)
    arr = sortslices(arr, dims=1);
    index = 1;
    narr = convert(Int64, size(arr,1));

    for i = 2:narr
        if arr[index, 2] >= arr[i, 1]
            arr[index, 2] = max(arr[index, 2], arr[i, 2]);
        else
            index = index+1;
            arr[index,:] = arr[i,:];
        end
    end

    return arr[1:index, 1:end]
end

function gs_gaps(accesses::Array{DateTime}, JD₀, Δt::Number)
    # Get the epoch of the propagator.
    DT₀ = jd_to_date(DateTime, JD₀)

    if isempty(accesses)
        return nothing;
    else
        # Compute the last propagation instant.
        JD₁ = JD₀ + Δt/86400
        DT₁ = jd_to_date(DateTime, JD₁)

        # Compute the gaps between accesses.
        gaps = Matrix{DateTime}(undef,0,2)

        # Check if the simulation did not start under the visibility of a ground
        # station.
        try
            if accesses[1,1] != DT₀
                gaps = vcat(gaps, [jd_to_date(DateTime, JD₀) accesses[1,1]])
            end
        
        # Check if the simulation did not end under the visibility of a ground
        # station.
            if accesses[end,2] != DT₁
                aux  = jd_to_date(DateTime, JD₁)
                gaps = vcat(gaps, [accesses[:,2] vcat(accesses[2:end,1], aux) ])
            else
                gaps = vcat(gaps, [accesses[1:end-1,2] accesses[2:end,1] ])
            end
            
            catch
        end
        return gaps
    end
end

function ground_station_visible(r_e::AbstractVector, rs_e::AbstractVector,θ::Number)
    # Check if the satellite is within the visibility circle of the station.
    dr_e = r_e - rs_e;
    cos_beta = dot(dr_e / norm(dr_e), rs_e / norm(rs_e));

    return cos_beta > cos(π / 2 - θ)
end







### iniciar dados
t, epoca, pcd_list, rec_list, var_min, var_max, var_names, var_bits, nsat, pECI, pECEF, passo = start_data(12);


# Execução
function run_opt()
    mapa_alt = [500 550 600 625 650 675 700 750 800 850 900 950 1000 1200 1400 1600 1800 2000];
    println("Altitude [km];    tgap[min]")
    for altitude in mapa_alt
        constellation = standalone_const2(epoca, altitude);
        tgap, valid = satcon_optm(constellation, t, pcd_list, rec_list, pECI, pECEF; tstep = passo);
        println("$altitude;    $tgap")
    end
end
