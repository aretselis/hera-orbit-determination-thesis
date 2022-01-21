using SPICE, ProgressBars, Plots, LinearAlgebra
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")


function load_hera_spice_kernels()
    # Extract Hera orbit from SPICE
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)
    # Load all planets, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
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
    furnsh(hera_kernel)
    furnsh(hera_ck)
    furnsh(hera_equipment)
    furnsh(hera_pck)
    furnsh(bodies_pck)
    furnsh(hera_sclk)
    furnsh(hera_instruments)
end


function main()
    # System properties
    G = 6.67430*10^-11
    mass_didymos = 5.32*10^11
    mass_dimorphos = 4.94*10^11
    mu_system = G*(mass_didymos+mass_dimorphos)

    # Time properties
    hour = 3600.0
    number_of_steps = 1000
    start_time = 0.0
    end_time = 200*hour
    step_size = (end_time-start_time)/number_of_steps

    # Dimorphos orbit
    a_dimorphos = 1183.0 # Semi-major axis
    e_dimorphos = 0.001 # Eccentricity
    i_dimorphos = 0.005 # Inclination
    Omega_dimorphos = 0.0 # Argument of periapsis
    omega_dimorphos = 0.0 # Longitude of the ascending node
    M_dimorphos = 0.0 # Time of periapsis passage

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

    load_hera_spice_kernels()

    # Orbit start and end time (for SPICE)
    spice_start_time = utc2et("2027-02-25T08:14:58")
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

    # Free used kernels
    kclear()

    # Plot 3D image for reference    
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    plotlyjs()
    plt_3d = scatter3d(x_didymos, y_didymos, z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel = "z [km]", label="Didymos", markersize = 1)
    scatter3d!(x_dimorphos, y_dimorphos, z_dimorphos, color = "orange", markersize = 1, label = "Dimorphos (RK4)")

    # Compute and plot camera boundary lines
    sensor_boundary_points = zeros(Float64, 4, 3)
    x_boundaries = zeros(Float64, 4)
    y_boundaries = zeros(Float64, 4)
    for i in 1:4
        x_line, y_line, z_line = boundary_vectors_matrix[i, 1] * line_generator, boundary_vectors_matrix[i, 2] * line_generator, boundary_vectors_matrix[i, 3] * line_generator
        sensor_boundary_points[i, 1] = x_line[length(x_line)]
        sensor_boundary_points[i, 2] = y_line[length(y_line)]
        sensor_boundary_points[i, 3] = z_line[length(z_line)]
        x_boundaries[i] = focal_length*x_line[length(x_line)]/z_line[length(z_line)]
        y_boundaries[i] = focal_length*y_line[length(y_line)]/z_line[length(z_line)]
        camera_label = "Camera boundary line " * string(i)
        plot3d!(x_line, y_line, z_line, color="green", label=camera_label) 
    end

    # Compute coordinates in pixels
    x_pixel_didymos, y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)
    x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2],x_boundaries, y_boundaries)
    x_pixel_boundaries, y_pixel_boundaries = convert_to_pixels(x_boundaries,y_boundaries, x_boundaries, y_boundaries)

    # Plot 2D image simulation (camera plane)  
    plt_2d = scatter(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], label= "Didymos")
    scatter!(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], label = "Dimorphos")
    scatter!(x_boundaries, y_boundaries, label = "Camera boundaries")

    # Plot image using pixels 
    plt_pixel = scatter(x_pixel_didymos, y_pixel_didymos, label = "Didymos")
    scatter!(x_pixel_dimorphos, y_pixel_dimorphos, label = "Dimorphos")
    scatter!(x_pixel_boundaries, y_pixel_boundaries, label= "Pixel boundaries")

    display(plt_pixel)
end

main()
