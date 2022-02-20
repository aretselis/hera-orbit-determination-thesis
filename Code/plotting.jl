function propagate_and_plot_orbital_elements_vs_time(OE_original, OE_fitted)
    # Dimorphos orbit
    a_dimorphos = OE_original[1] # Semi-major axis
    e_dimorphos = OE_original[2] # Eccentricity
    i_dimorphos = OE_original[3] # Inclination
    Omega_dimorphos = OE_original[4] # Argument of periapsis
    omega_dimorphos = OE_original[5] # Longitude of the ascending node
    M_dimorphos = OE_original[6] # Mean anomaly

    # Dimorphos properties
    area_dimorphos = 2*pi*85.0^2
    c_p = 2.2

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
     x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_step_size, enable_perturbation, area_dimorphos, c_p)
     # Select every 10th element to match the photos
     x_dimorphos = x_dimorphos[1:10:end]
     y_dimorphos = y_dimorphos[1:10:end]
     z_dimorphos = z_dimorphos[1:10:end]
     vx_dimorphos = vx_dimorphos[1:10:end]
     vy_dimorphos = vy_dimorphos[1:10:end]
     vz_dimorphos = vz_dimorphos[1:10:end]
     t_vector = t_vector[1:10:end]
     dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)
 
     # Compute orbital elements vs time for the original orbit
     a_vector = zeros(Float64, number_of_steps + 1)
     e_vector = zeros(Float64, number_of_steps + 1)
     i_vector = zeros(Float64, number_of_steps + 1)
     Omega_vector = zeros(Float64, number_of_steps + 1)
     omega_vector = zeros(Float64, number_of_steps + 1)
     M_vector = zeros(Float64, number_of_steps + 1)
    for i in 1:(number_of_steps+1)
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
    propagation_steps = 10000
    propagation_step_size = (end_time-start_time)/propagation_steps
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, propagation_step_size, enable_perturbation, area_dimorphos, c_p)
    # Select every 10th element to match the photos
    x_dimorphos = x_dimorphos[1:10:end]
    y_dimorphos = y_dimorphos[1:10:end]
    z_dimorphos = z_dimorphos[1:10:end]
    vx_dimorphos = vx_dimorphos[1:10:end]
    vy_dimorphos = vy_dimorphos[1:10:end]
    vz_dimorphos = vz_dimorphos[1:10:end]
    t_vector = t_vector[1:10:end]
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Compute orbital elements vs time for the fitted orbit
    a_vector = zeros(Float64, number_of_steps + 1)
    e_vector = zeros(Float64, number_of_steps + 1)
    i_vector = zeros(Float64, number_of_steps + 1)
    Omega_vector = zeros(Float64, number_of_steps + 1)
    omega_vector = zeros(Float64, number_of_steps + 1)
    M_vector = zeros(Float64, number_of_steps + 1)
    for i in 1:(number_of_steps+1)
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
