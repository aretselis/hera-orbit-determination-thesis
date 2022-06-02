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

    if dimorphos_orbit_type == "Normal"
        plt_a_vs_time = plot(t_vector, a_vector, label= L"a(t)" * " original", xlabel = "Time, " * L"t,\ [s]", ylabel = "Semi-major axis, " * L"a,\ [m]", widen=true, formatter=:plain, legend = :topright)
        plt_e_vs_time = plot(t_vector, e_vector, label= L"e(t)" * " original", xlabel = "Time, " * L"t,\ [s]", ylabel = "Eccentricity, " * L"e,\ [\ ]", widen=true, formatter=:plain, legend = :topright)
        plt_i_vs_time = plot(t_vector, i_vector, label= L"i(t)" * " original", xlabel = "Time, " * L"t,\ [s]", ylabel = "Inclination, " * L"i,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_Omega_vs_time = plot(t_vector, Omega_vector, label= L"\Omega(t)" * " original", xlabel = "Time, " * L"t,\ [s]", ylabel = "Longitude of the ascending node, " * L"\Omega,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_periapsis_vs_time = plot(t_vector, omega_vector, label= L"\omega(t)" * " original", xlabel = "Time, " * L"t,\ [s]", ylabel = "Argument of periapsis, " * L"\omega,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_M_vs_time = plot(t_vector, M_vector, label= L"M(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Mean anomaly, " * L"M,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
    elseif dimorphos_orbit_type == "EEO"
        plt_a_vs_time = plot(t_vector, a_vector, label= L"a(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Semi-major axis, " * L"a,\ [m]", widen=true, formatter=:plain, legend = :topright)
        plt_e_vs_time = plot(t_vector, e_vector, label= L"e(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Eccentricity, " * L"e,\ [\ ]", widen=true, formatter=:plain, legend = :topright)
        plt_i_vs_time = plot(t_vector, i_vector, label= L"i(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Inclination, " * L"i,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_Omega_vs_time = plot(t_vector, Omega_vector, label= L"\tilde{\omega}_{true}(t)\ \textnormal{original}", xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True longitude of periapsis},\  \tilde{\omega}_{true},\ [deg]", widen=true, formatter=:plain, legend = :topright)
        plt_periapsis_vs_time = plot(t_vector, omega_vector, label= L"\tilde{\omega}_{true}(t)\ \textnormal{original}", xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True longitude of periapsis},\  \tilde{\omega}_{true},\ [deg]", widen=true, formatter=:plain, legend = :topright)
        plt_M_vs_time = plot(t_vector, M_vector, label= L"M(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = ylabel = "Mean anomaly, " * L"M,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
    elseif dimorphos_orbit_type == "CIO"
        plt_a_vs_time = plot(t_vector, a_vector, label= L"a(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Semi-major axis, " * L"a,\ [m]", widen=true, formatter=:plain, legend = :topright)
        plt_e_vs_time = plot(t_vector, e_vector, label= L"e(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Eccentricity, " * L"e,\ [\ ]", widen=true, formatter=:plain, legend = :topright)
        plt_i_vs_time = plot(t_vector, i_vector, label= L"i(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Inclination, " * L"i,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_Omega_vs_time = plot(t_vector, Omega_vector, label= L"\Omega(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Longitude of the ascending node, " * L"\Omega,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_periapsis_vs_time = plot(t_vector, omega_vector, label= L"u(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Argument of latitude, "* L"u,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_M_vs_time = plot(t_vector, M_vector, label= L"u(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Argument of latitude, "* L"u,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
    elseif dimorphos_orbit_type == "CEO"
        plt_a_vs_time = plot(t_vector, a_vector, label= L"a(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Semi-major axis, " * L"a,\ [m]", widen=true, formatter=:plain, legend = :topright)
        plt_e_vs_time = plot(t_vector, e_vector, label= L"e(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Eccentricity, " * L"e,\ [\ ]", widen=true, formatter=:plain, legend = :topright)
        plt_i_vs_time = plot(t_vector, i_vector, label= L"i(t)" * " original", xlabel = "Time, "* L"t,\ [s]", ylabel = "Inclination, " * L"i,\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_Omega_vs_time = plot(t_vector, Omega_vector, label= L"\lambda_{true}\ \textnormal{original}", xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True Longitude},\ \lambda_{true},\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_periapsis_vs_time = plot(t_vector, omega_vector, label=L"\lambda_{true}\ \textnormal{original}", xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True Longitude},\ \lambda_{true},\ [\deg]", widen=true, formatter=:plain, legend = :topright)
        plt_M_vs_time = plot(t_vector, M_vector, label= L"\lambda_{true}(t)\ \textnormal{original}", xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True Longitude},\ \lambda_{true},\ [\deg]", widen=true, formatter=:plain, legend = :topright)
    end

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

    if dimorphos_orbit_type == "Normal"
        plot!(plt_a_vs_time, t_vector, a_vector, label= L"a(t)" * " prediction")
        plot!(plt_e_vs_time, t_vector, e_vector, label= L"e(t)" * " prediction")
        plot!(plt_i_vs_time, t_vector, i_vector, label= L"i(t)" * " prediction")
        plot!(plt_Omega_vs_time, t_vector, Omega_vector, label= L"\Omega(t)" * " prediction")
        plot!(plt_periapsis_vs_time, t_vector, omega_vector, label= L"\omega(t)" * " prediction")
        plot!(plt_M_vs_time, t_vector, M_vector, label= L"M(t)" * " prediction")
    elseif dimorphos_orbit_type == "EEO"
        plot!(plt_a_vs_time, t_vector, a_vector, label= L"a(t)" * " prediction")
        plot!(plt_e_vs_time, t_vector, e_vector, label= L"e(t)" * " prediction")
        plot!(plt_i_vs_time, t_vector, i_vector, label= L"i(t)" * " prediction")
        plot!(plt_Omega_vs_time, t_vector, Omega_vector, label= L"\tilde{\omega}_{true}(t)\ \textnormal{prediction}")
        plot!(plt_periapsis_vs_time, t_vector, omega_vector, label= L"\tilde{\omega}_{true}(t)\ \textnormal{prediction}")
        plot!(plt_M_vs_time, t_vector, M_vector, label= L"M(t)" * " prediction")
    elseif dimorphos_orbit_type == "CIO"
        plot!(plt_a_vs_time, t_vector, a_vector, label= L"a(t)" * " prediction")
        plot!(plt_e_vs_time, t_vector, e_vector, label= L"e(t)" * " prediction")
        plot!(plt_i_vs_time, t_vector, i_vector, label= L"i(t)" * " prediction")
        plot!(plt_Omega_vs_time, t_vector, Omega_vector, label= L"\Omega(t)" * " prediction")
        plot!(plt_periapsis_vs_time, t_vector, omega_vector, label= L"u(t)" * " prediction")
        plot!(plt_M_vs_time, t_vector, M_vector, label= L"u(t)" * " prediction")
    elseif dimorphos_orbit_type == "CEO"
        plot!(plt_a_vs_time, t_vector, a_vector, label= L"a(t)" * " prediction")
        plot!(plt_e_vs_time, t_vector, e_vector, label= L"e(t)" * " prediction")
        plot!(plt_i_vs_time, t_vector, i_vector, label= L"i(t)" * " prediction")
        plot!(plt_Omega_vs_time, t_vector, Omega_vector, label= L"\lambda_{true}\ \textnormal{prediction}")
        plot!(plt_periapsis_vs_time, t_vector, omega_vector, label= L"\lambda_{true}(t)\ \textnormal{prediction}")
        plot!(plt_M_vs_time, t_vector, M_vector, label= L"\lambda_{true}(t)\ \textnormal{prediction}")
    end

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

    plt_position = plot(t_vector, x_dimorphos, label= "x(t) original", xlabel = "Time, "* L"t,\ [s]", ylabel = "[m]")
    plot!(plt_position, t_vector, y_dimorphos, label= "y(t) original", xlabel = "Time, "* L"t,\ [s]", ylabel = "[m]")
    plot!(plt_position, t_vector, z_dimorphos, label= "z(t) original", xlabel = "Time, "* L"t,\ [s]", ylabel = " [m]")
    display(plt_position)
end


function plot_accuracy_vs_photo_analysis(photo_vector)
    
    accuracy_array =  load("C:\\Users\\retse\\repos\\hera-orbit-determination\\Results\\accuracy_array.jld")["accuracy_array"]
    number_of_runs = size(accuracy_array)[1] - 1
    a_error_vector = zeros(Float64, number_of_runs)
    e_error_vector = zeros(Float64, number_of_runs)
    i_error_vector = zeros(Float64, number_of_runs)
    M_error_vector = zeros(Float64, number_of_runs)
    for i=1:length(photo_vector)
        a_error_vector[i] = accuracy_array[i, 1]
        e_error_vector[i] = accuracy_array[i, 2]
        i_error_vector[i] = accuracy_array[i, 3]
        M_error_vector[i] = accuracy_array[i, 6]
    end
    pgfplotsx()
    accuracy_plot = plot(photo_vector, a_error_vector, xaxis=:log, label = L"a", xlabel = "Number of photos used", ylabel = "Mean Absolute Performance Error [%]", widen=true, formatter=:plain, legend = :topright)
    plot!(photo_vector, e_error_vector, xaxis=:log, label = L"e")
    #plot!(photo_vector, i_error_vector, xaxis=:log, label = L"i")
    plot!(photo_vector, M_error_vector, xaxis=:log, label = L"\lambda_{true}")
    savefig(".\\Results\\accuracy_vs_photos_plot.pdf")
end

