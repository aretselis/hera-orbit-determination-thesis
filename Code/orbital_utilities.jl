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


function orbital_elements_to_cartesian(a, e, i, Omega, omega, M, mu)
    #=
    Computes cartesian position and velocity vector given some orbital elements
    Input:
    a [m]
    e []
    i [deg]
    Omega [deg]
    omega [deg]
    M [deg]
    mu [m^3/(kg*s^2)] (GM)
    Output:
    r_vector, v_vector
    =#

    # Convert M to radians
    M = deg2rad(M)
    # Compute E using Newton-Raphson
    E = eccentric_anomaly_calculator(M, e)
    # Compute x, xdot, y, ydot on the orbital plane
    x = a * (cos(E) - e)
    y = a * sqrt(1 - e^2) * sin(E)
    r = sqrt(x^2 + y^2)
    n = sqrt(mu / (a^3))  # Mean motion
    x_dot = -(n * a^2 / r) * sin(E)
    y_dot = (n * a^2 / r) * sqrt(1 - e^2) * cos(E)
    # Rotation Matrices definition
    Omega = deg2rad(Omega)
    omega = deg2rad(omega)
    i = deg2rad(i)
    P1 = [cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1]
    P2 = [1 0 0;0 cos(i) sin(i);0 sin(i) cos(i)]
    P3 = [cos(Omega) -sin(Omega) 0;sin(Omega) cos(Omega) 0;0 0 1]
    # Compute cartesian coordinates
    x_y_vector = [x, y, 0]
    x_y_dot_vector = [x_dot, y_dot, 0]
    r_vector = *(*(*(P3, P2), P1), x_y_vector)
    v_vector = *(*(*(P3, P2), P1), x_y_dot_vector)
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
    n_vector = cross([0,0,1], h_vector)
    n = sqrt(n_vector[1]^2 + n_vector[2]^2 + n_vector[3]^2)
    e_vector = (((v^2 - (mu/r))*r_vector) - (dot(r_vector,v_vector)*v_vector))/mu
    e = sqrt(e_vector[1]^2 + e_vector[2]^2 + e_vector[3]^2)
    ksi = ((v^2)/2) - (mu/r)
    a = - mu/(2*ksi)
    i = acos(h_vector[3]/h)
    Omega = acos(n_vector[1]/n)
    if n_vector[2] < 0
        Omega = 2*pi - Omega
    end
    omega = acos(dot(n_vector, e_vector)/(n*e))
    if e_vector[3] < 0
        omega = 2*pi - omega
    end
    ni = acos(dot(e_vector, r_vector)/(e*r))
    if dot(r_vector, v_vector) < 0
        ni = 2*pi - ni
    end
    M = atan((sqrt(1-e^2)*sin(ni)/(1+e*cos(ni))), (e+cos(ni))/(1+e*cos(ni))) - (e*(sqrt(1-e^2)*sin(ni)/(1+e*cos(ni))))  
    # Convert to degrees
    i = rad2deg(i)
    Omega = rad2deg(Omega)
    ni = rad2deg(ni)
    omega = rad2deg(omega)
    M = rad2deg(M)
    return a, e, i, Omega, omega, M
end


function orbit_evolution_and_cartesian_transform(a, e, i, n, po, tau, GM, time)
    #=
    Computes the orbit at the specified time step_size
    INPUT:
    a, e, i, n, po, tau: Orbital elements
    GM: GM of the System
    time: desired time step
    =#
    # Calculate perihelion distance
    q = a*abs(1.0-e)
    # Calculate longitude of pericenter
    p = n + po
    # Calculate mean motion
    nu = sqrt(GM)*a^(-1.5)
    # Calculate mean anomaly
    L = nu*(time-tau)
    # Compute orbit vectors
    x, y, z, v_x, v_y, v_z = mco_el2x
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


function propagate_and_compute_dimorphos_pixel_points(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, start_time, end_time, step_size, spice_start_time)
    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    propagation_steps = 10000
    propagation_step_size = (end_time-start_time)/propagation_steps
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_step_size, enable_perturbation)
    # Select every 10th element to match the photos
    x_dimorphos = x_dimorphos[1:10:end]
    y_dimorphos = y_dimorphos[1:10:end]
    z_dimorphos = z_dimorphos[1:10:end]
    vx_dimorphos = vx_dimorphos[1:10:end]
    vy_dimorphos = vy_dimorphos[1:10:end]
    vz_dimorphos = vz_dimorphos[1:10:end]
    t_vector = t_vector[1:10:end]
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Orbit start and end time (for SPICE)
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_offset = zeros(Float64, 3)
    barycenter_initial_coordinates = zeros(Float64, 3)
    camera_position_didymos = zeros(Float64, number_of_steps + 1, 3)
    camera_position_dimorphos = zeros(Float64, number_of_steps + 1, 3)
    didymos_pixel_coordinates = zeros(Float64, number_of_steps + 1, 2)
    dimorphos_pixel_coordinates = zeros(Float64, number_of_steps + 1, 2)

    # HERA_AFC-1 coordinate system unit vectors
    e_x = [1.0, 0.0, 0.0]
    e_y = [0.0, 1.0, 0.0]
    e_z = [0.0, 0.0, 1.0]

    # Main SPICE calculations loop
    focal_length = 10.6*10^-5 # [km]
    iteration = 1
    #println("\nSPICE calculations:")
    for i in (spice_start_time:step_size:spice_end_time)
        position_didymos, lt = spkpos("-658030", i, "HERA_AFC-1", "None", "-999")
        position_system_barycenter, lt = spkpos("2065803", i, "HERA_AFC-1", "None", "-999")          
        rotation_frame = pxform("J2000", "HERA_SPACECRAFT", i)
        position_dimorphos = rotation_frame * dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        if iteration == 1
            for j in 1:3
                barycenter_initial_coordinates[j] = position_system_barycenter[j]
            end
        else
            for j in 1:3
                barycenter_offset[j] = barycenter_initial_coordinates[j] - position_system_barycenter[j]
            end
        end
        # Apply offset to all coordinates
        for j in 1:3
            barycenter_coordinates[iteration, j] =  position_system_barycenter[j] + barycenter_offset[j]
            didymos_coordinates[iteration, j] =  position_didymos[j] + barycenter_offset[j]
            dimorphos_coordinates[iteration, j] =  position_dimorphos[j] + barycenter_initial_coordinates[j]
        end
        # Rotate HERA_AFC-1 camera reference frame assuming it is always tracking the barycenter
        rotation_matrix = compute_rotation_matrix(barycenter_coordinates[iteration, :], e_z)
        barycenter_coordinates[iteration, :] = rotation_matrix * barycenter_coordinates[iteration, :]
        didymos_coordinates[iteration, :] =  rotation_matrix * didymos_coordinates[iteration, :]
        dimorphos_coordinates[iteration, :] =  rotation_matrix * dimorphos_coordinates[iteration, :]
        # Rotate Didymos-Dimorphos points
        rotation_matrix = compute_rotation_matrix(e_z, barycenter_coordinates[iteration, :])
        camera_position_didymos[iteration, :] = rotation_matrix * didymos_coordinates[iteration, :]
        camera_position_dimorphos[iteration, :] = rotation_matrix * dimorphos_coordinates[iteration, :]
        # Store pixel coordinates for Didymos and Dimorphos
        for j in 1:2
            didymos_pixel_coordinates[iteration, j] = focal_length*camera_position_didymos[iteration, j]/camera_position_didymos[iteration, 3]
            dimorphos_pixel_coordinates[iteration, j] = focal_length*camera_position_dimorphos[iteration, j]/camera_position_dimorphos[iteration, 3]
        end
        iteration += 1
    end

    x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)

    # Add centroid pixel error
    x_pixel_dimorphos = x_pixel_dimorphos .+ x_centroid_pixel_error
    y_pixel_dimorphos = y_pixel_dimorphos .+ y_centroid_pixel_error

    return x_pixel_dimorphos, y_pixel_dimorphos
end


function propagate_and_compute_dimorphos_3D_points(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, start_time, end_time, step_size, spice_start_time)
    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    propagation_steps = 10000
    propagation_step_size = (end_time-start_time)/propagation_steps
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_step_size, enable_perturbation)
    # Select every 10th element to match the photos
    x_dimorphos = x_dimorphos[1:10:end]
    y_dimorphos = y_dimorphos[1:10:end]
    z_dimorphos = z_dimorphos[1:10:end]
    vx_dimorphos = vx_dimorphos[1:10:end]
    vy_dimorphos = vy_dimorphos[1:10:end]
    vz_dimorphos = vz_dimorphos[1:10:end]
    t_vector = t_vector[1:10:end]
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Orbit start and end time (for SPICE)
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_offset = zeros(Float64, 3)
    barycenter_initial_coordinates = zeros(Float64, 3)
    camera_position_didymos = zeros(Float64, number_of_steps + 1, 3)
    camera_position_dimorphos = zeros(Float64, number_of_steps + 1, 3)
    didymos_pixel_coordinates = zeros(Float64, number_of_steps + 1, 2)
    dimorphos_pixel_coordinates = zeros(Float64, number_of_steps + 1, 2)

    # HERA_AFC-1 coordinate system unit vectors
    e_x = [1.0, 0.0, 0.0]
    e_y = [0.0, 1.0, 0.0]
    e_z = [0.0, 0.0, 1.0]

    # Main SPICE calculations loop
    focal_length = 10.6*10^-5 # [km]
    iteration = 1
    #println("\nSPICE calculations:")
    for i in (spice_start_time:step_size:spice_end_time)
        position_didymos, lt = spkpos("-658030", i, "HERA_AFC-1", "None", "-999")
        position_system_barycenter, lt = spkpos("2065803", i, "HERA_AFC-1", "None", "-999")          
        rotation_frame = pxform("J2000", "HERA_SPACECRAFT", i)
        position_dimorphos = rotation_frame * dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        if iteration == 1
            for j in 1:3
                barycenter_initial_coordinates[j] = position_system_barycenter[j]
            end
        else
            for j in 1:3
                barycenter_offset[j] = barycenter_initial_coordinates[j] - position_system_barycenter[j]
            end
        end
        # Apply offset to all coordinates
        for j in 1:3
            barycenter_coordinates[iteration, j] =  position_system_barycenter[j] + barycenter_offset[j]
            didymos_coordinates[iteration, j] =  position_didymos[j] + barycenter_offset[j]
            dimorphos_coordinates[iteration, j] =  position_dimorphos[j] + barycenter_initial_coordinates[j]
        end
        # Rotate HERA_AFC-1 camera reference frame assuming it is always tracking the barycenter
        rotation_matrix = compute_rotation_matrix(barycenter_coordinates[iteration, :], e_z)
        barycenter_coordinates[iteration, :] = rotation_matrix * barycenter_coordinates[iteration, :]
        didymos_coordinates[iteration, :] =  rotation_matrix * didymos_coordinates[iteration, :]
        dimorphos_coordinates[iteration, :] =  rotation_matrix * dimorphos_coordinates[iteration, :]
        # Rotate Didymos-Dimorphos points
        rotation_matrix = compute_rotation_matrix(e_z, barycenter_coordinates[iteration, :])
        camera_position_didymos[iteration, :] = rotation_matrix * didymos_coordinates[iteration, :]
        camera_position_dimorphos[iteration, :] = rotation_matrix * dimorphos_coordinates[iteration, :]
        # Store pixel coordinates for Didymos and Dimorphos
        for j in 1:2
            didymos_pixel_coordinates[iteration, j] = focal_length*camera_position_didymos[iteration, j]/camera_position_didymos[iteration, 3]
            dimorphos_pixel_coordinates[iteration, j] = focal_length*camera_position_dimorphos[iteration, j]/camera_position_dimorphos[iteration, 3]
        end
        iteration += 1
    end
    return dimorphos_coordinates
end