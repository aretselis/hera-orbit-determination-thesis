using SPICE, ProgressBars, Plots, LinearAlgebra, Optim
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")
include("optimization.jl")
include("results.jl")


function load_hera_spice_kernels()
    # Extract Hera orbit from SPICE
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)
    # Load all planets, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
    didymos_barycenter_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\didymos_hor_000101_500101_v01.bsp"
    hera_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_PO_v01.bsp"
    hera_ck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\ck\\hera_sc_PO_EMA_20270209_20270727_f20181203_v03.bc"
    hera_equipment = "C:\\Users\\retse\\repos\\hera-data\\kernels\\fk\\hera_v07.tf"
    hera_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\hera_didymos_v04.tpc"
    bodies_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\de-403-masses.tpc"
    hera_sclk = "C:\\Users\\retse\\repos\\hera-data\\kernels\\sclk\\hera_fict_20181203.tsc"
    hera_instruments = "C:\\Users\\retse\\repos\\hera-data\\kernels\\ik\\hera_afc_v03.ti"
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP1_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3VCF_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_ECP_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP1_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3VCF_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_ECP_v01.bsp")
    furnsh(planets_kernel)
    furnsh(didymos_kernel)
    furnsh(dimorphos_kernel)
    furnsh(didymos_barycenter_kernel)
    furnsh(hera_kernel)
    furnsh(hera_ck)
    furnsh(hera_equipment)
    furnsh(hera_pck)
    furnsh(bodies_pck)
    furnsh(hera_sclk)
    furnsh(hera_instruments)
end


function main()
    # Dimorphos orbit
    a_dimorphos = 1183.725 # Semi-major axis
    e_dimorphos = 0.00039 # Eccentricity
    i_dimorphos = 6.9 # Inclination
    Omega_dimorphos = 7.2 # Argument of periapsis
    omega_dimorphos = 5.8 # Longitude of the ascending node
    M_dimorphos = 3.85 # Mean anomaly

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

    plotlyjs()
    plt_oe_vs_time = plot(t_vector, a_vector, label= "a(t)")
    #scatter!(t_vector, e_vector, label= "e(t)")
    #scatter!(t_vector, i_vector, label= "i(t)")
    #scatter!(t_vector, Omega_vector, label= "Omega(t)")
    #scatter!(t_vector, omega_vector, label= "omega(t)")
    #scatter!(t_vector, M_vector, label= "M(t)")

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
    println("\nSPICE calculations:")
    for i in ProgressBar(spice_start_time:step_size:spice_end_time)
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

    # Compute boundary vectors/points for the camera
    boresight_vector = barycenter_coordinates[iteration-1, :]
    rotation_matrix = compute_rotation_matrix(boresight_vector, e_z)
    camera_id = -999110
    boundary_vectors_matrix = get_boundary_vectors(camera_id)
    line_generator = LinRange(0, 30, 10)

    # Rotate camera boundary vectors to assume camera is always tracking barycenter
    for i in 1:4
        result = rotation_matrix * boundary_vectors_matrix[i, :]
        for j in 1:3
            boundary_vectors_matrix[i,j] = result[j]
        end
    end

    # Plot 3D image for reference    
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    #plt_3d = scatter3d(x_didymos, y_didymos, z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel = "z [km]", label="Didymos", markersize = 1)
    #scatter3d!(x_dimorphos, y_dimorphos, z_dimorphos, color = "orange", markersize = 1, label = "Dimorphos (RK4)")

    # Compute and plot camera boundary lines
    sensor_boundary_points = zeros(Float64, 4, 3)
    global x_boundaries = zeros(Float64, 4)
    global y_boundaries = zeros(Float64, 4)
    for i in 1:4
        x_line, y_line, z_line = boundary_vectors_matrix[i, 1] * line_generator, boundary_vectors_matrix[i, 2] * line_generator, boundary_vectors_matrix[i, 3] * line_generator
        sensor_boundary_points[i, 1] = x_line[length(x_line)]
        sensor_boundary_points[i, 2] = y_line[length(y_line)]
        sensor_boundary_points[i, 3] = z_line[length(z_line)]
        x_boundaries[i] = focal_length*x_line[length(x_line)]/z_line[length(z_line)]
        y_boundaries[i] = focal_length*y_line[length(y_line)]/z_line[length(z_line)]
        camera_label = "Camera boundary line " * string(i)
        #plot3d!(x_line, y_line, z_line, color="green", label=camera_label) 
    end

    # Compute coordinates in pixels
    x_pixel_didymos, y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)
    global x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)
    x_pixel_boundaries, y_pixel_boundaries = convert_to_pixels(x_boundaries,y_boundaries, x_boundaries, y_boundaries)

    # Bounds for optimization
    lower = [1190.0-30, 0.0, 0.0, 0.0, 1.0, 0.0]
    upper = [1190.0+30, 0.001, 10.0, 10.0, 10.0, 10.0]
    initial_guess = [1200.0, 0.0001, 4.0, 5.0, 3.0, 7.0]

    # Minimize square mean error to find best orbital elements
    global res = optimize(residuals, initial_guess, ParticleSwarm(; lower, upper, n_particles = 50), Optim.Options(iterations = 150, g_tol = 1e-10, show_trace = true))
    print(res)
    println(Optim.minimizer(res))

    # Extract final guess from optimization
    a_final = Optim.minimizer(res)[1]
    e_final = Optim.minimizer(res)[2]
    i_final = Optim.minimizer(res)[3]
    Omega_final = Optim.minimizer(res)[4]
    omega_final = Optim.minimizer(res)[5]
    M_final = Optim.minimizer(res)[6]

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

    # plotlyjs()
    #plt_oe_vs_time = plot(t_vector, a_vector, label= "a(t)")
    plot!(t_vector, a_vector, label= "a(t) (fit)")
    #scatter!(t_vector, i_vector, label= "i(t)")
    #scatter!(t_vector, Omega_vector, label= "Omega(t)")
    #scatter!(t_vector, omega_vector, label= "omega(t)")
    #scatter!(t_vector, M_vector, label= "M(t)")

    # Compute percentage error of the predicted orbit relative to the initial orbit 
    initial_elements = [a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos]
    final_elements = [a_final, e_final, i_final, Omega_final, omega_final, M_final]
    compute_percentage_error(initial_elements, final_elements)

    # Compute pixel and 3D points from dimorphos as a result of the optimization
    x_pixel_fit_dimorphos, y_pixel_fit_dimorphos = propagate_and_compute_dimorphos_pixel_points(a_final, e_final, i_final, Omega_final, omega_final, M_final, start_time, end_time, step_size, spice_start_time)
    dimorphos_fit_coordinates = propagate_and_compute_dimorphos_3D_points(a_final, e_final, i_final, Omega_final, omega_final, M_final, start_time, end_time, step_size, spice_start_time)
    x_fit_dimorphos = dimorphos_fit_coordinates[:, 1]
    y_fit_dimorphos = dimorphos_fit_coordinates[:, 2]
    z_fit_dimorphos = dimorphos_fit_coordinates[:, 3]
    #scatter3d!(x_fit_dimorphos, y_fit_dimorphos, z_fit_dimorphos, color = "purple", markersize = 1, label = "Dimorphos (Final guess)")


    # Plot 2D image simulation (camera plane)  
    plt_2d = scatter(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], label= "Didymos")
    scatter!(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], label = "Dimorphos")
    scatter!(x_boundaries, y_boundaries, label = "Camera boundaries")

    # Plot image using pixels 
    plt_pixel = scatter(x_pixel_didymos, y_pixel_didymos, label = "Didymos")
    scatter!(x_pixel_dimorphos, y_pixel_dimorphos, label = "Dimorphos")
    scatter!(x_pixel_fit_dimorphos, y_pixel_fit_dimorphos, label = "Dimorphos best fit")
    scatter!(x_pixel_boundaries, y_pixel_boundaries, label= "Pixel boundaries")

    display(plt_oe_vs_time)
end

function cb(os)
    println(os.metadata)
    return false    
end

load_hera_spice_kernels()

# System properties
G = 6.67430*10^-11
mass_didymos = 5.32*10^11
global mass_dimorphos = 4.94*10^11
global mu_system = G*(mass_didymos+mass_dimorphos)
global c = 3.0*10^8
global enable_perturbation = true

# Time properties
hour = 3600.0
number_of_steps = 1000
global start_time = 0.0
global end_time = 13*hour
global step_size = (end_time-start_time)/number_of_steps
global spice_start_time = utc2et("2027-02-25T08:14:58")
global flux = 1358

main()

# Free used kernels
kclear()