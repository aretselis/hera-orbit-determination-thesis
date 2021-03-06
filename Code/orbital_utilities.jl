function eccentric_anomaly_calculator(mean_anomaly, e)
    #=
    Compute E using Newton-Raphson
    Input:
    mean_anomaly, [rad]
    e, eccentricity []
    Output
    E, eccentric anomaly at desired time [rad]
    =#

    # Define initial value for E_0)
    if mean_anomaly < pi
        E_0 = mean_anomaly - e
    else
        E_0 = mean_anomaly + e
    end
    # Define f and f dot
    f(E) = mean_anomaly - E + e*sin(E)
    fdot(E) =  -1 + e*cos(E)
    # Stopping criteria
    N = 15  # Number of significant digits to be computed
    max_repetitions = 1000000
    es = 0.5 * 10^(2.0 - N)  # Scarborough Criterion
    ea = 100
    E_prev = E_0
    repetitions = 0
    # Main Newton-Raphson loop
    while ea > es
        repetitions = repetitions + 1
        E_next = E_prev - (f(E_prev) / fdot(E_prev))
        if E_next == 0
            return E_next
        end
        ea = abs((E_next - E_prev) * 100 / E_next)
        E_prev = E_next
        if repetitions > max_repetitions
            error("Max repetitions reached without achieving desired accuracy for E!")
        end
    end
    E = E_prev
    return E
end


function orbital_elements_to_cartesian(a, e, i, Omega, omega, ni, mu)
    #=
    Computes cartesian position and velocity vector given some orbital elements
    Input:
    a [m]
    e []
    i [deg]
    Omega [deg]
    omega [deg]
    ni [deg]
    mu [m^3/(kg*s^2)] (GM)
    Output:
    r_vector, v_vector
    =#

    if e < 10^-5
        if i < 10^-1
            # Circular equatorial orbit
            Omega = 0.0
            omega = 0.0 
        else
            # Circular inclined orbit
            omega = 0.0
        end
        # In this case, p=a
        p = a
    else
        if i < 10^-1
            # Elliptical equatorial orbit
            Omega = 0.0
        end
        p = a*(1-e^2)
    end
    r_pqw = [(p*cosd(ni))/(1+e*cosd(ni)); (p*sind(ni))/(1+e*cosd(ni)); 0.0]
    v_pqw = [-sqrt(mu/p)*sind(ni); sqrt(mu/p)*(e+cosd(ni)); 0.0]
    transformation = [cosd(Omega)*cosd(omega)-sind(Omega)*sind(omega)*cosd(i) -cosd(Omega)*sind(omega)-sind(Omega)*cosd(omega)*cosd(i) sind(Omega)*sind(i);
                      sind(Omega)*cosd(omega)+cosd(Omega)*sind(omega)*cosd(i) -sind(Omega)*sind(omega)+cosd(Omega)*cosd(omega)*cosd(i) -cosd(Omega)*sind(i);
                      sind(omega)*sind(i) cosd(omega)*sind(i) cosd(i)]
    r_vector = transformation*r_pqw
    v_vector = transformation*v_pqw
    return r_vector, v_vector
end 


function cartesian_to_orbital_elements(r_vector, v_vector, mu)
    #=
    Computes cartesian position and velocity vector given some orbital elements
    Input:
    r_vector [m] , v_vector [m/s]
    mu [m^3/(kg*s^2)] (GM)
    Output:
    a [m] (Semi-major axis)
    e [] (Eccentricity)
    i [deg] (Inclination)
    Omega [deg] (Longitude of the ascending node)
    omega [deg] (Argument of pericenter)
    ni [deg] (True anomaly)
    M [deg] (Mean anomaly)
    =#

    # Compute magnitude and angular momentum
    r = sqrt(r_vector[1]^2 + r_vector[2]^2 + r_vector[3]^2)
    v = sqrt(v_vector[1]^2 + v_vector[2]^2 + v_vector[3]^2)
    h_vector = cross(r_vector, v_vector)
    h = sqrt(h_vector[1]^2 + h_vector[2]^2 + h_vector[3]^2)
    n_vector = cross([0.0,0.0,1.0], h_vector)
    n = sqrt(n_vector[1]^2 + n_vector[2]^2 + n_vector[3]^2)
    e_vector = (((v^2 - (mu/r))*r_vector) - (dot(r_vector,v_vector)*v_vector))/mu
    e = sqrt(e_vector[1]^2 + e_vector[2]^2 + e_vector[3]^2)
    ksi = ((v^2)/2) - (mu/r)
    a = - mu/(2*ksi)
    i = acosd(h_vector[3]/h)
    current_orbit = dimorphos_orbit_type
    Omega = acosd(n_vector[1]/n)
    if n_vector[2] < 0
        Omega = 360 - Omega
    end
    if current_orbit == "EEO"
        omega = acosd(e_vector[1]/e)
        if e_vector[2] < 0
            omega = 360 - omega
        end
    else
        omega = acosd(dot(n_vector, e_vector)/(n*e))
        if e_vector[3] < 0
            omega = 360 - omega
        end
    end
    if current_orbit == "CIO"
        u = acosd(dot(n_vector, r_vector)/(n*r))
        if r_vector[3] < 0
            u = 360 - u
        end
        return a, e, i, Omega, omega, u
    elseif current_orbit == "CEO"
        lambda_true = acosd(r_vector[1]/r)
        if r_vector[2] < 0
            lambda_true = 360 - lambda_true
        end
        return a, e, i, Omega, omega, lambda_true
    else
        ni = acosd(dot(e_vector, r_vector)/(e*r))
        if dot(r_vector, v_vector) < 0
            ni = 360 - ni
        end
        M = deg2rad(ni) - 2*e*sind(ni) + ((3*e^2/4)+(e^4/8))*sind(2*ni) - (e^3/3)*sind(3*ni) + (5*e^4/32)*sind(4*ni) 
        #M = atan((sqrt(1-e^2)*sind(ni)/(1+e*cosd(ni))), (e+cosd(ni))/(1+e*cosd(ni))) - (e*(sqrt(1-e^2)*sind(ni)/(1+e*cosd(ni))))  
        M = rad2deg(M)
        return a, e, i, Omega, omega, M
    end
end


function sight_algorithm(r1_vector, r2_vector)
    #=
    Computes if there is sight between two points in space (assuming Didymos to be the obstuction)
    INPUT:
    r1_vector, Vector(Float64, 3), position vector of the first point in [m]
    r2_vector, Vector(Float64, 3), position vector of the second point in [m]
    OUTPUT:
    0, if there is no line of sight between the two points
    1, if there is line of sight between the two points
    =#
    radius_didymos = 325 # [m]
    r1_magnitude = sqrt(r1_vector[1]^2 + r1_vector[2]^2 + r1_vector[3]^2)
    r2_magnitude = sqrt(r2_vector[1]^2 + r2_vector[2]^2 + r2_vector[3]^2)
    tau_min = (r1_magnitude^2 - dot(r1_vector, r2_vector))/(r1_magnitude^2 + r2_magnitude^2 - (2*dot(r1_vector, r2_vector)))
    line_of_sight = false
    if tau_min < 0 || tau_min > 1
        line_of_sight = true
    else
        parametric_tau_min = ((1-tau_min)*r1_magnitude^2) + (dot(r1_vector, r2_vector)*tau_min)
        if parametric_tau_min >= radius_didymos
            line_of_sight = true
        end
    end
    if line_of_sight == false
        println("Shadow!")
    end
    return Int64(line_of_sight)
end 


function geometrical_shadow_check(x_dimorphos, y_dimorphos, z_dimorphos, x_sun, y_sun, z_sun)
    radius_didymos = 325 # [m]
    radius_sun = 696340000 # [m]
    distance = sqrt(x_sun^2 + y_sun^2 + z_sun^2)
    a_umbra = asin((radius_sun - radius_didymos)/distance)
    a_penumbra = asin((radius_sun + radius_didymos)/distance)
    # Initially assume we are in sunlight, compute magnitudes and dot product
    shadow = 1
    r_dimorphos_magnitude = sqrt(x_dimorphos^2 + y_dimorphos^2 + z_dimorphos^2)
    r_sun_magnitude = sqrt(x_sun^2 + y_sun^2 + z_sun^2)
    dot_product = x_dimorphos*x_sun + y_dimorphos*y_sun + z_dimorphos*z_sun
    # Decide if we should evaluate umbra/penumbra possibility
    if dot_product < 0
        dot_product = -x_dimorphos * x_sun - y_dimorphos * y_sun - z_dimorphos * z_sun
        angle_s = acos(dot_product/(r_dimorphos_magnitude*r_sun_magnitude))
        sat_horizontal = r_dimorphos_magnitude*cos(angle_s)
        sat_vertical = r_dimorphos_magnitude*sin(angle_s)
        x = radius_didymos/sin(a_penumbra)
        pen_vertical = tan(a_penumbra)*(x+sat_horizontal)
        if sat_vertical <= pen_vertical
            y = radius_didymos/sin(a_umbra)
            umb_vertical = tan(a_umbra)*(y-sat_horizontal)
            if sat_vertical <= umb_vertical
                shadow = 0
            else
                # Assume that Solar Intensity decreases linearly, calculate shadow value
                slope = (1-0)/(pen_vertical-umb_vertical)
                sat_vertical = sat_vertical - umb_vertical
                shadow = slope*sat_vertical
            end
        end
    end
    return shadow
end


function propagate_and_compute_dimorphos_pixel_points(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, start_time, end_time, step_size, spice_start_time, mu)
    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu, start_time, end_time, total_photos, enable_perturbation)
    # Select every xth element to match the photos taken
    photo_selector = LinRange(1, length(x_dimorphos), total_photos)
    index_vector = zeros(Int64, total_photos)
    for i = 1:length(index_vector)
        index_vector[i] = floor(Int64, photo_selector[i])
    end
    x_dimorphos = x_dimorphos[index_vector]
    y_dimorphos = y_dimorphos[index_vector]
    z_dimorphos = z_dimorphos[index_vector]
    vx_dimorphos = vx_dimorphos[index_vector]
    vy_dimorphos = vy_dimorphos[index_vector]
    vz_dimorphos = vz_dimorphos[index_vector]
    t_vector = t_vector[index_vector]
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Orbit start and end time (for SPICE)
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, total_photos, 3)
    barycenter_coordinates = zeros(Float64, total_photos, 3)
    barycenter_offset = zeros(Float64, 3)
    barycenter_initial_coordinates = zeros(Float64, 3)
    camera_position_didymos = zeros(Float64, total_photos, 3)
    camera_position_dimorphos = zeros(Float64, total_photos, 3)
    didymos_pixel_coordinates = zeros(Float64, total_photos, 2)
    dimorphos_pixel_coordinates = zeros(Float64, total_photos, 2)

    # HERA_AFC-1 coordinate system unit vectors
    e_x = [1.0, 0.0, 0.0]
    e_y = [0.0, 1.0, 0.0]
    e_z = [0.0, 0.0, 1.0]

    # Main SPICE calculations loop
    focal_length = 10.6*10^-5 # [km]
    iteration = 1
    for current_time in LinRange(start_time, end_time, total_photos)
        i = spice_start_time + current_time
        position_didymos, lt = spkpos("-658030", i, "J2000", "None", "-999")
        position_system_barycenter, lt = spkpos("2065803", i, "J2000", "None", "-999")          
        position_dimorphos = dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        for j in 1:3
            barycenter_offset[j] = 0.0
        end
        # Apply offset to all coordinates
        for j in 1:3
            barycenter_coordinates[iteration, j] = position_system_barycenter[j] + barycenter_offset[j]
            didymos_coordinates[iteration, j] = position_didymos[j] + barycenter_offset[j]
            dimorphos_coordinates[iteration, j] = barycenter_coordinates[iteration, j] + position_dimorphos[j]
        end
        # Rotate HERA_AFC-1 camera reference frame assuming it is always tracking the barycenter
        rotation_matrix = compute_rotation_matrix(barycenter_coordinates[iteration, :], e_z)
        barycenter_coordinates[iteration, :] = rotation_matrix * barycenter_coordinates[iteration, :]
        didymos_coordinates[iteration, :] =  rotation_matrix * didymos_coordinates[iteration, :]
        dimorphos_coordinates[iteration, :] =  rotation_matrix * dimorphos_coordinates[iteration, :]
        # Add pointing error 
        pointing_error_matrix = error_rotation_matrix(random_error[iteration], random_axis[iteration])
        camera_position_didymos[iteration, :] = pointing_error_matrix * didymos_coordinates[iteration, :]
        camera_position_dimorphos[iteration, :] = pointing_error_matrix * dimorphos_coordinates[iteration, :]
        # Store pixel coordinates for Didymos and Dimorphos
        for j in 1:2
            didymos_pixel_coordinates[iteration, j] = focal_length*camera_position_didymos[iteration, j]/camera_position_didymos[iteration, 3]
            dimorphos_pixel_coordinates[iteration, j] = focal_length*camera_position_dimorphos[iteration, j]/camera_position_dimorphos[iteration, 3]
        end
        iteration += 1
    end

    local_x_pixel_didymos, local_y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
    local_x_pixel_dimorphos, local_y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)

    # Reject images where we cannot see one object 
    for i in 1:length(local_x_pixel_didymos)
        if ismissing(local_x_pixel_didymos[i]) || ismissing(local_x_pixel_dimorphos[i])
            local_x_pixel_didymos[i] = missing
            local_y_pixel_didymos[i] = missing
            local_x_pixel_dimorphos[i] = missing
            local_y_pixel_dimorphos[i] = missing
        end
    end
    return local_x_pixel_dimorphos, local_y_pixel_dimorphos
end


function propagate_and_compute_dimorphos_3D_points(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, start_time, end_time, step_size, spice_start_time, mu)
    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu, start_time, end_time, total_photos, enable_perturbation)
    # Select every xth element to match the photos taken
    photo_selector = LinRange(1, length(x_dimorphos), total_photos)
    index_vector = zeros(Int64, total_photos)
    for i = 1:length(index_vector)
        index_vector[i] = floor(Int64, photo_selector[i])
    end
    x_dimorphos = x_dimorphos[index_vector]
    y_dimorphos = y_dimorphos[index_vector]
    z_dimorphos = z_dimorphos[index_vector]
    vx_dimorphos = vx_dimorphos[index_vector]
    vy_dimorphos = vy_dimorphos[index_vector]
    vz_dimorphos = vz_dimorphos[index_vector]
    t_vector = t_vector[index_vector]
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Orbit start and end time (for SPICE)
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, total_photos, 3)
    barycenter_coordinates = zeros(Float64, total_photos, 3)
    barycenter_offset = zeros(Float64, 3)
    barycenter_initial_coordinates = zeros(Float64, 3)
    camera_position_didymos = zeros(Float64, total_photos, 3)
    camera_position_dimorphos = zeros(Float64, total_photos, 3)
    didymos_pixel_coordinates = zeros(Float64, total_photos, 2)
    dimorphos_pixel_coordinates = zeros(Float64, total_photos, 2)

    # HERA_AFC-1 coordinate system unit vectors
    e_x = [1.0, 0.0, 0.0]
    e_y = [0.0, 1.0, 0.0]
    e_z = [0.0, 0.0, 1.0]

    # Main SPICE calculations loop
    focal_length = 10.6*10^-5 # [km]
    iteration = 1
    for current_time in LinRange(start_time, end_time, total_photos)
        i = spice_start_time + current_time
        position_didymos, lt = spkpos("-658030", i, "J2000", "None", "-999")
        position_system_barycenter, lt = spkpos("2065803", i, "J2000", "None", "-999")          
        position_dimorphos = dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        for j in 1:3
            barycenter_offset[j] = 0.0
        end
        # Apply offset to all coordinates
        for j in 1:3
            barycenter_coordinates[iteration, j] = position_system_barycenter[j] + barycenter_offset[j]
            didymos_coordinates[iteration, j] = position_didymos[j] + barycenter_offset[j]
            dimorphos_coordinates[iteration, j] = barycenter_coordinates[iteration, j] + position_dimorphos[j]
        end
        # Rotate HERA_AFC-1 camera reference frame assuming it is always tracking the barycenter
        rotation_matrix = compute_rotation_matrix(barycenter_coordinates[iteration, :], e_z)
        barycenter_coordinates[iteration, :] = rotation_matrix * barycenter_coordinates[iteration, :]
        didymos_coordinates[iteration, :] =  rotation_matrix * didymos_coordinates[iteration, :]
        dimorphos_coordinates[iteration, :] =  rotation_matrix * dimorphos_coordinates[iteration, :]
        # Add pointing error 
        pointing_error_matrix = error_rotation_matrix(random_error[iteration], random_axis[iteration])
        camera_position_didymos[iteration, :] = pointing_error_matrix * didymos_coordinates[iteration, :]
        camera_position_dimorphos[iteration, :] = pointing_error_matrix * dimorphos_coordinates[iteration, :]
        # Store pixel coordinates for Didymos and Dimorphos
        for j in 1:2
            didymos_pixel_coordinates[iteration, j] = focal_length*camera_position_didymos[iteration, j]/camera_position_didymos[iteration, 3]
            dimorphos_pixel_coordinates[iteration, j] = focal_length*camera_position_dimorphos[iteration, j]/camera_position_dimorphos[iteration, 3]
        end
        iteration += 1
    end
    return dimorphos_coordinates
end


function type_of_orbit(e, i)
    # Returns orbit type given eccentricity and inclination
    if e < 10^-5
        if i < 10^-1
            # Circular equatorial orbit
            answer = "CEO"
        else
            # Circular inclined orbit
            answer = "CIO"
        end
    else
        if i < 10^-1 
            # Elliptical equatorial orbit
            answer = "EEO"
        else 
            # Normal orbit
            answer = "Normal"
        end
    end
    return answer
end