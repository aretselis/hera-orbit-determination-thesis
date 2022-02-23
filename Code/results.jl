function compute_percentage_error(actual_elements, predicted_elements)
    #=
    Computes percentage error given two sets of orbital elements
    Input:
    actual_elements, Vector(Float64, 6) containing actual orbital elements
    predicted_elements, Vector(Float64, 6) containing predicted orbital elements 
    Output:
    percentage_errors, Vector(Float64, 6) containing percentage (%) errors
    =#

    percentage_errors = zeros(Float64, 6)
    for i in 1:6
        if i>=3 && i<=6
            # Ensure conversion to radians for angles
            actual_elements[i] = deg2rad(actual_elements[i])
            predicted_elements[i] = deg2rad(predicted_elements[i])
        end
        percentage_errors[i] = 100 * (abs(predicted_elements[i] - actual_elements[i])/actual_elements[i])
        println("Percentage error for element " * string(i) * " is " * string(percentage_errors[i]))
    end
    return Nothing
end


function compute_mean_absolute_percentage_error(actual_elements, predicted_elements)
    #=
    Computes mean absolute percentage errors for two orbits, provided with the initial orbital elements 
    Input:
    actual_elements, Vector(Float64, 6) containing actual orbital elements
    predicted_elements, Vector(Float64, 6) containing predicted orbital elements 
    Output:
    percentage_errors, Vector(Float64, 6) mean absolute percentage errors (%) errors
    =#

    # Dimorphos orbit
    a_dimorphos = actual_elements[1] # Semi-major axis
    e_dimorphos = actual_elements[2] # Eccentricity
    i_dimorphos = actual_elements[3] # Inclination
    Omega_dimorphos = actual_elements[4] # Argument of periapsis
    omega_dimorphos = actual_elements[5] # Longitude of the ascending node
    M_dimorphos = actual_elements[6] # Mean anomaly

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, total_photos, enable_perturbation)
    # Select every xth element to match the photos taken
    photo_selector = Int64(round((size(x_dimorphos)[1]-1)/total_photos))
    x_dimorphos = x_dimorphos[1:photo_selector:end]
    y_dimorphos = y_dimorphos[1:photo_selector:end]
    z_dimorphos = z_dimorphos[1:photo_selector:end]
    vx_dimorphos = vx_dimorphos[1:photo_selector:end]
    vy_dimorphos = vy_dimorphos[1:photo_selector:end]
    vz_dimorphos = vz_dimorphos[1:photo_selector:end]
    t_vector = t_vector[1:photo_selector:end]
 
    # Compute orbital elements vs time for the fitted orbit
    vector_size = size(x_dimorphos)[1]
    a_vector = zeros(Float64, vector_size)
    e_vector = zeros(Float64, vector_size)
    i_vector = zeros(Float64, vector_size)
    Omega_vector = zeros(Float64, vector_size)
    omega_vector = zeros(Float64, vector_size)
    M_vector = zeros(Float64, vector_size)
    for i in 1:vector_size
        a_vector[i], e_vector[i], i_vector[i], Omega_vector[i], omega_vector[i], M_vector[i] =  cartesian_to_orbital_elements([x_dimorphos[i], y_dimorphos[i], z_dimorphos[i]], [vx_dimorphos[i], vy_dimorphos[i], vz_dimorphos[i]], mu_system)
        i_vector[i] = deg2rad(i_vector[i])
        Omega_vector[i] = deg2rad(Omega_vector[i])
        omega_vector[i] = deg2rad(omega_vector[i])
        M_vector[i] = deg2rad(M_vector[i])
    end

    observed_values = [a_vector e_vector i_vector Omega_vector omega_vector M_vector]

    # Extract final guess from optimization
    a_final = predicted_elements[1] 
    e_final = predicted_elements[2]
    i_final = predicted_elements[3]
    Omega_final = predicted_elements[4]
    omega_final = predicted_elements[5]
    M_final = predicted_elements[6]

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_final, e_final, i_final, Omega_final, omega_final,M_final, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, total_photos, enable_perturbation)
    # Select every xth element to match the photos taken
    photo_selector = Int64(round((size(x_dimorphos)[1]-1)/total_photos))
    x_dimorphos = x_dimorphos[1:photo_selector:end]
    y_dimorphos = y_dimorphos[1:photo_selector:end]
    z_dimorphos = z_dimorphos[1:photo_selector:end]
    vx_dimorphos = vx_dimorphos[1:photo_selector:end]
    vy_dimorphos = vy_dimorphos[1:photo_selector:end]
    vz_dimorphos = vz_dimorphos[1:photo_selector:end]
    t_vector = t_vector[1:photo_selector:end]

    # Compute orbital elements vs time for the fitted orbit
    vector_size = size(x_dimorphos)[1]
    a_vector = zeros(Float64, vector_size)
    e_vector = zeros(Float64, vector_size)
    i_vector = zeros(Float64, vector_size)
    Omega_vector = zeros(Float64, vector_size)
    omega_vector = zeros(Float64, vector_size)
    M_vector = zeros(Float64, vector_size)
    for i in 1:vector_size
        a_vector[i], e_vector[i], i_vector[i], Omega_vector[i], omega_vector[i], M_vector[i] = cartesian_to_orbital_elements([x_dimorphos[i], y_dimorphos[i], z_dimorphos[i]], [vx_dimorphos[i],vy_dimorphos[i], vz_dimorphos[i]], mu_system)
        i_vector[i] = deg2rad(i_vector[i])
        Omega_vector[i] = deg2rad(Omega_vector[i])
        omega_vector[i] = deg2rad(omega_vector[i])
        M_vector[i] = deg2rad(M_vector[i])
    end

    predicted_values = [a_vector e_vector i_vector Omega_vector omega_vector M_vector]

    percentage_errors = zeros(Float64, 6)
    for i in 1:6
        sum = 0.0
        for j in 1:vector_size
            sum += abs((observed_values[j, i] - predicted_values[j, i])/observed_values[j, i])
        end
        percentage_errors[i] = 100 * sum / (vector_size)
        println("Mean absolute percentage error for element " * string(i) * " is " * string(percentage_errors[i]) * " %")
    end
    return Nothing
end