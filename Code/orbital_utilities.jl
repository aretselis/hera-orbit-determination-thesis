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
    # Compute a, e and i 
    a = mu/((2*mu/r) - v^2)
    e = sqrt(1 - ((h^2)/(mu*a)))
    i = acos(h_vector[3]/h)
    # Compute Omega (Longitude of the ascending node)
    if h_vector[3] > 0
        sin_Omega =   h_vector[1] / (h * sin(i))
        cos_Omega = - h_vector[2] / (h * sin(i))
    else
        sin_Omega = - h_vector[1] / (h * sin(i))
        cos_Omega =   h_vector[2] / (h * sin(i))
    end
    if sin_Omega >= 0
        Omega = acos(cos_Omega)
    else
        Omega = 2*pi - acos(cos_Omega)
    end
    # Compute E 
    cos_E = (a-r)/(a*e)
    sin_E = dot(r_vector, v_vector) / (e*sqrt(mu*a))
    # Compute true anomaly
    cos_ni = (cos_E - e) / (1 - (e*cos_E))
    sin_ni = (sqrt(1 - e^2) * sin_E) / (1 - (e * cos_E))
    if sin_ni >=0
        ni = acos(cos_ni)
    else
        ni = 2*pi - acos(cos_ni)
    end
    # Compute mean anomaly
    if sin_E >=0
        E = acos(cos_E)
    else
        E = 2*pi - acos(cos_E)
    end
    M = E - e*sin_E
    # Compute omega (argument of pericenter)
    sin_omega_plus_ni = r_vector[3] / (r*sin(i))
    cos_omega_plus_ni = ((r_vector[1]*cos_Omega) + (r_vector[2]*sin_Omega))/r
    if sin_omega_plus_ni >=0
        omega = acos(cos_omega_plus_ni) - ni
    else
        omega = 2*pi - acos(cos_omega_plus_ni) - ni
    end
    # Convert to degrees
    i = rad2deg(i)
    Omega = rad2deg(Omega)
    ni = rad2deg(ni)
    M = rad2deg(M)
    omega = rad2deg(omega)
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
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_step_size)
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
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, step_size)
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