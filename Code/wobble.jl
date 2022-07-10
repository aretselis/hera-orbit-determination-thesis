using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics, LaTeXStrings
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")
include("optimization.jl")
include("results.jl")
include("plotting.jl")


function main()
    # Dimorphos orbit
    global a_dimorphos = rand(Uniform(1160, 1220)) # Semi-major axis
    global e_dimorphos = rand(Uniform(0.0000001, 0.00001)) # Eccentricity
    global i_dimorphos = rand(Uniform(0.00000001, 0.00001)) # Inclination
    global Omega_dimorphos = 0.0 # Longitude of the ascending node
    global omega_dimorphos = 0.0 # Argument of periapsis
    global M_dimorphos = rand(Uniform(0.0, 359.99)) # True anomaly 

    # Determine what kind of orbit we are dealing with
    global dimorphos_orbit_type = type_of_orbit(e_dimorphos, i_dimorphos)
    println("Orbit type is "*dimorphos_orbit_type)

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
    # Select every xth element to match the photos taken and convert to km
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

    x_didymos = (G*mass_dimorphos/mu_system_original)*x_dimorphos
    y_didymos = (G*mass_dimorphos/mu_system_original)*y_dimorphos
    z_didymos = (G*mass_dimorphos/mu_system_original)*z_dimorphos
    didymos_coordinates = hcat(x_didymos/1000, y_didymos/1000, z_didymos/1000)

    r_didymos_original = zeros(Float64, length(x_didymos))

    for kati=1:length(r_didymos_original)
        r_didymos_original[kati] = sqrt(x_didymos[kati]^2 + y_didymos[kati]^2 + z_didymos[kati]^2)
    end
    pgfplotsx()   
    plt_wobble = plot(t_vector, r_didymos_original, label = "Didymos barycentric orbit (Original)", widen=true, formatter=:plain, legend = :topright, xlabel = "Time, t, [s]", ylabel = "Position, r, [m]")

    # Orbit start and end time (for SPICE)
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    #didymos_coordinates = zeros(Float64, total_photos, 3)
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
    println("\nSPICE calculations:")
    for current_time in ProgressBar(LinRange(start_time, end_time, total_photos))
        i = spice_start_time + current_time
        position_didymos = didymos_coordinates[iteration, :]
        position_system_barycenter, lt = spkpos("2065803", i, "J2000", "None", "-999")          
        position_dimorphos = dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        for j in 1:3
            barycenter_offset[j] = rand(barycenter_error_distribution) + rand(hera_error_distribution)
        end
        # Apply offset to all coordinates
        for j in 1:3
            barycenter_coordinates[iteration, j] = position_system_barycenter[j] + barycenter_offset[j]
            didymos_coordinates[iteration, j] = barycenter_coordinates[iteration, j] + position_didymos[j]
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
    plt_3d = scatter3d(x_didymos, y_didymos, z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel = "z [km]", label="Didymos", markersize = 1.5)
    scatter3d!(plt_3d, x_dimorphos, y_dimorphos, z_dimorphos, color = "orange", markersize = 1.5, label = "Dimorphos (RK4)")

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
        #plot3d!(plt_3d, x_line, y_line, z_line, color="green", label=camera_label) 
    end

    # Compute coordinates in pixels
    global x_pixel_didymos, y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
    x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
    x_pixel_boundaries, y_pixel_boundaries = convert_to_pixels(x_boundaries,y_boundaries, x_boundaries, y_boundaries)

    # Add centroid pixel error
    global x_pixel_didymos = x_pixel_didymos
    global y_pixel_didymos = y_pixel_didymos 

        # Circular equatorial orbit, OD with 4 orbital elements
        # Bounds for optimization
        bounds = [0.6*mass_dimorphos 0.0;
                  1.4*mass_dimorphos 0.0]
        result = optimize(wobble, bounds, ECA(options = Options(debug=false, iterations=50, store_convergence = true)))
        # Extract final guess from optimization
        m_dimorphos_guess = minimizer(result)[1]

    #compute_percentage_error(initial_elements, final_elements)
    dimorphos_percentage_error = 100*abs(mass_dimorphos-m_dimorphos_guess)/mass_dimorphos
    println("Initial mass " * string(mass_dimorphos) * ", Final mass: " *string(m_dimorphos_guess))
    println("Error on M_dimorphos is " * string(dimorphos_percentage_error) * " %")
    # Plot image using pixels 
    plt_pixel = scatter(x_pixel_didymos, y_pixel_didymos, label = "Didymos", widen=true, formatter=:plain, legend = :topright, xlabel = "X Pixel Coordinate", ylabel = "Y Pixel Coordinate")

     # Compute initial position and velocity vector from orbital elements
     r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, G*(mass_didymos+m_dimorphos_guess))
     x = r_vector[1]
     y = r_vector[2]
     z = r_vector[3]
     vx = v_vector[1]
     vy = v_vector[2]
     vz = v_vector[3]
 
     # Propagate Dimorphos
     x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, G*(mass_didymos+m_dimorphos_guess), start_time, end_time, total_photos, enable_perturbation)
     # Select every xth element to match the photos taken and convert to km
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
 
     x_didymos = (G*mass_dimorphos/mu_system_original)*x_dimorphos
     y_didymos = (G*mass_dimorphos/mu_system_original)*y_dimorphos
     z_didymos = (G*mass_dimorphos/mu_system_original)*z_dimorphos
     didymos_coordinates = hcat(x_didymos/1000, y_didymos/1000, z_didymos/1000)

    r_didymos_original = zeros(Float64, length(x_didymos))
    for kati=1:length(r_didymos_original)
        r_didymos_original[kati] = sqrt(x_didymos[kati]^2 + y_didymos[kati]^2 + z_didymos[kati]^2)
    end
   
    plot!(plt_wobble, t_vector, r_didymos_original, label = "Didymos barycentric orbit (Predicted)", widen=true, formatter=:plain, legend = :topright, xlabel = "Time, t, [s]", ylabel = "Position, r, [m]")
    display(plt_wobble)
    
 
     # Orbit start and end time (for SPICE)
     spice_end_time = spice_start_time + end_time
 
     # Initialize position arrays (x, y, z)
     #didymos_coordinates = zeros(Float64, total_photos, 3)
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
     #println("\nSPICE calculations:")
     for current_time in LinRange(start_time, end_time, total_photos)
         i = spice_start_time + current_time
         position_didymos = didymos_coordinates[iteration, :]
         position_system_barycenter, lt = spkpos("2065803", i, "J2000", "None", "-999")          
         position_dimorphos = dimorphos_coordinates[iteration, :]
         # Get didymos reference position or translate hera_coordinates
         for j in 1:3
             barycenter_offset[j] = rand(barycenter_error_distribution) + rand(hera_error_distribution)
         end
         # Apply offset to all coordinates
         for j in 1:3
             barycenter_coordinates[iteration, j] = position_system_barycenter[j] + barycenter_offset[j]
             didymos_coordinates[iteration, j] = barycenter_coordinates[iteration, j] + position_didymos[j]
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
     pgfplotsx()   
     x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
     x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
     plt_3d = scatter3d(x_didymos, y_didymos, z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel = "z [km]", label="Didymos", markersize = 1.5)
     scatter3d!(plt_3d, x_dimorphos, y_dimorphos, z_dimorphos, color = "orange", markersize = 1.5, label = "Dimorphos (RK4)")
 
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
         #plot3d!(plt_3d, x_line, y_line, z_line, color="green", label=camera_label) 
     end
 
     # Compute coordinates in pixels
     global x_pixel_didymos, y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
     x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
     x_pixel_boundaries, y_pixel_boundaries = convert_to_pixels(x_boundaries,y_boundaries, x_boundaries, y_boundaries)
 
     # Add centroid pixel error
     global x_pixel_didymos = x_pixel_didymos
     global y_pixel_didymos = y_pixel_didymos
     scatter!(plt_pixel, x_pixel_didymos, y_pixel_didymos, label = "Didymos best fit")
     savefig(plt_wobble, ".\\Results\\wobble_position.pdf")
end

load_hera_spice_kernels()


#for random_counter in 1:20
#    println("=====RUN " *string(random_counter)* "======")
    # System properties
    global G = 6.67430*10^-11
    mass_error = 0.02 # in percentage
    mass_error_distribution = Uniform(-mass_error, mass_error)
    global mass_didymos = 5.32*10^11 + rand(mass_error_distribution)*5.32*10^11
    global radius_didymos = 390 # [m]
    global J2_didymos = 0.011432722
    global mass_dimorphos = 4.94*10^9 + rand(mass_error_distribution)*4.94*10^9
    global mu_system = G*(mass_didymos+mass_dimorphos)
    global mu_system_original = mu_system
    println(mu_system_original)
    global c = 3.0*10^8
    global enable_perturbation = true

    # Time properties
    hour = 3600.0
    number_of_orbits = 25
    # Photos to be taken 
    photos_per_orbit = 40
    global total_photos = photos_per_orbit * number_of_orbits
    global start_time = 0.0
    global end_time = 300*hour
    global step_size = floor((end_time-start_time)/total_photos)
    global spice_start_time = utc2et("2027-02-28T08:14:58")
    global flux = 1358

    # Errors definition
    # Centroid pixel error
    pixel_error = 4
    pixel_error_distribution = Uniform(-pixel_error, pixel_error)
    global x_centroid_pixel_error = Int64(round(rand(pixel_error_distribution)))
    global y_centroid_pixel_error = Int64(round(rand(pixel_error_distribution)))
    # Barycenter position error
    global error_barycenter = 10/1000 # [km]
    global barycenter_error_distribution = Normal(0.0, error_barycenter)
    # HERA position error
    global error_hera_position = 10/1000 # [km]
    global hera_error_distribution = Normal(0.0, error_hera_position)
    # HERA pointing error
    global pointing_error = deg2rad(1)
    global pointing_error_distribution = Normal(0.0, pointing_error)
    global random_error = rand(pointing_error_distribution, total_photos)
    global random_axis = rand(1:2, total_photos)
    # Randomly dropped images
    number_of_dropped_images = Int64(round(0.04*total_photos))
    dropped_image_distribution = Uniform(1, total_photos)
    global dropped_images_index = zeros(Int64, number_of_dropped_images) 
    for p=1:number_of_dropped_images
        dropped_images_index[p] = Int64(round(rand(dropped_image_distribution)))
    end
    main()
#end

# Free used kernels
kclear()