function propagate_and_plot_orbital_elements_vs_time(OE_original, OE_fitted)
    # Dimorphos orbit
    a_dimorphos = OE_original[1] # Semi-major axis
    e_dimorphos = OE_original[2] # Eccentricity
    i_dimorphos = OE_original[3] # Inclination
    Omega_dimorphos = OE_original[4] # Argument of periapsis
    omega_dimorphos = OE_original[5] # Longitude of the ascending node
    M_dimorphos = OE_original[6] # Mean anomaly

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    propagation_steps = 100000
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_steps, enable_perturbation)
 
    # Compute orbital elements vs time for the original orbit
    vector_size = size(x_dimorphos)[1]
    a_vector = zeros(Float64, vector_size)
    e_vector = zeros(Float64, vector_size)
    i_vector = zeros(Float64, vector_size)
    Omega_vector = zeros(Float64, vector_size)
    omega_vector = zeros(Float64, vector_size)
    M_vector = zeros(Float64, vector_size)
    for i in 1:vector_size
        a_vector[i], e_vector[i], i_vector[i], Omega_vector[i], omega_vector[i], M_vector[i] =  cartesian_to_orbital_elements([x_dimorphos[i], y_dimorphos[i], z_dimorphos[i]], [vx_dimorphos[i], vy_dimorphos[i], vz_dimorphos[i]], mu_system)
    end

    plt_a_vs_time = plot(t_vector, a_vector, label= "a(t) original", xlabel = "Time, t, [seconds]", ylabel = "Semi-major axis, a, [m]")
    plt_e_vs_time = plot(t_vector, e_vector, label= "e(t)", xlabel = "Time, t, [seconds]", ylabel = "Eccentricity, e, []")
    plt_i_vs_time = plot(t_vector, i_vector, label= "i(t)", xlabel = "Time, t, [seconds]", ylabel = "Inclination, i, [deg]")
    plt_Omega_vs_time = plot(t_vector, Omega_vector, label= "Omega(t)", xlabel = "Time, t, [seconds]", ylabel = "Longitude of the ascending node, Omega, [deg]")
    plt_periapsis_vs_time = plot(t_vector, omega_vector, label= "omega(t)", xlabel = "Time, t, [seconds]", ylabel = "Argument of periapsis, omega, [deg]")
    plt_M_vs_time = plot(t_vector, M_vector, label= "M(t)", xlabel = "Time, t, [seconds]", ylabel = "Mean anomaly, M, [deg]")

    # Extract final guess from optimization
    a_final = OE_fitted[1] 
    e_final = OE_fitted[2]
    i_final = OE_fitted[3]
    Omega_final = OE_fitted[4]
    omega_final = OE_fitted[5]
    M_final = OE_fitted[6]

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_final, e_final, i_final, Omega_final, omega_final, M_final, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    propagation_steps = 100000
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_steps, enable_perturbation)
 
    # Compute orbital elements vs time for the original orbit
    vector_size = size(x_dimorphos)[1]
    a_vector = zeros(Float64, vector_size)
    e_vector = zeros(Float64, vector_size)
    i_vector = zeros(Float64, vector_size)
    Omega_vector = zeros(Float64, vector_size)
    omega_vector = zeros(Float64, vector_size)
    M_vector = zeros(Float64, vector_size)
    for i in 1:(vector_size)
        a_vector[i], e_vector[i], i_vector[i], Omega_vector[i], omega_vector[i], M_vector[i] =  cartesian_to_orbital_elements([x_dimorphos[i], y_dimorphos[i], z_dimorphos[i]], [vx_dimorphos[i], vy_dimorphos[i], vz_dimorphos[i]], mu_system)
    end

    plot!(plt_a_vs_time, t_vector, a_vector, label= "a(t) Prediction")
    plot!(plt_e_vs_time, t_vector, e_vector, label= "e(t) Prediction")
    plot!(plt_i_vs_time, t_vector, i_vector, label= "i(t) Prediction")
    plot!(plt_Omega_vs_time, t_vector, Omega_vector, label= "Omega(t) Prediction")
    plot!(plt_periapsis_vs_time, t_vector, omega_vector, label= "omega(t) Prediction")
    plot!(plt_M_vs_time, t_vector, M_vector, label= "M(t) Prediction")

    savefig(plt_a_vs_time, ".\\Results\\a_vs_time.pdf")
    savefig(plt_e_vs_time, ".\\Results\\e_vs_time.pdf")
    savefig(plt_i_vs_time, ".\\Results\\i_vs_time.pdf")
    savefig(plt_Omega_vs_time, ".\\Results\\Omega_vs_time.pdf")
    savefig(plt_periapsis_vs_time, ".\\Results\\periapsis_vs_time.pdf")
    savefig(plt_M_vs_time, ".\\Results\\M_vs_time.pdf")
    return Nothing
end


function propagate_and_plot_xyz(OE_original)
    # Dimorphos orbit
    a_dimorphos = OE_original[1] # Semi-major axis
    e_dimorphos = OE_original[2] # Eccentricity
    i_dimorphos = OE_original[3] # Inclination
    Omega_dimorphos = OE_original[4] # Argument of periapsis
    omega_dimorphos = OE_original[5] # Longitude of the ascending node
    M_dimorphos = OE_original[6] # Mean anomaly

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    propagation_steps = 100000
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_steps, enable_perturbation)

    plt_position = plot(t_vector, x_dimorphos, label= "x(t) original", xlabel = "Time, t, [seconds]", ylabel = "[m]")
    plot!(plt_position, t_vector, y_dimorphos, label= "y(t) original", xlabel = "Time, t, [seconds]", ylabel = "[m]")
    plot!(plt_position, t_vector, z_dimorphos, label= "z(t) original", xlabel = "Time, t, [seconds]", ylabel = " [m]")
    display(plt_position)
end