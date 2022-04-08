using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")
include("optimization.jl")
include("results.jl")
include("plotting.jl")


function main()
    # Dimorphos orbit
    global a_dimorphos = 1183.593 # Semi-major axis
    global e_dimorphos = 0.00000098 # Eccentricity
    global i_dimorphos = 0.035 # Inclination
    global Omega_dimorphos = 0.0 # Longitude of the ascending node
    global omega_dimorphos = 0.0 # Argument of periapsis
    global M_dimorphos = 172.18 # True anomaly 

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
    photo_selector = Int64(floor((length(x_dimorphos)-1)/total_photos))
    x_dimorphos = x_dimorphos[1:photo_selector:end]
    y_dimorphos = y_dimorphos[1:photo_selector:end]
    z_dimorphos = z_dimorphos[1:photo_selector:end]
    vx_dimorphos = vx_dimorphos[1:photo_selector:end]
    vy_dimorphos = vy_dimorphos[1:photo_selector:end]
    vz_dimorphos = vz_dimorphos[1:photo_selector:end]
    t_vector = t_vector[1:photo_selector:end]
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
    println("\nSPICE calculations:")
    for current_time in ProgressBar(LinRange(start_time, end_time, total_photos))
        i = spice_start_time + current_time
        position_didymos, lt = spkpos("-658030", i, "J2000", "None", "-999")
        position_system_barycenter, lt = spkpos("2065803", i, "J2000", "None", "-999")          
        position_dimorphos = dimorphos_coordinates[iteration, :]
        # Get didymos reference position or translate hera_coordinates
        for j in 1:3
            barycenter_offset[j] = rand(barycenter_error_distribution) + rand(hera_error_distribution)
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
        # Rotate Didymos-Dimorphos points
        rotation_matrix = compute_rotation_matrix(e_z, barycenter_coordinates[iteration, :])
        camera_position_didymos[iteration, :] = rotation_matrix * didymos_coordinates[iteration, :]
        camera_position_dimorphos[iteration, :] = rotation_matrix * dimorphos_coordinates[iteration, :]
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
    plotlyjs()    
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    plt_3d = scatter3d(x_didymos, y_didymos, z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel = "z [km]", label="Didymos", markersize = 1)
    scatter3d!(plt_3d, x_dimorphos, y_dimorphos, z_dimorphos, color = "orange", markersize = 1, label = "Dimorphos (RK4)")

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
        plot3d!(plt_3d, x_line, y_line, z_line, color="green", label=camera_label) 
    end

    # Compute coordinates in pixels
    x_pixel_didymos, y_pixel_didymos = convert_to_pixels(didymos_pixel_coordinates[:, 1], didymos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
    global x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
    x_pixel_boundaries, y_pixel_boundaries = convert_to_pixels(x_boundaries,y_boundaries, x_boundaries, y_boundaries)

    # Reject images where we cannot see one object 
    for i in 1:length(x_pixel_didymos)
        if ismissing(x_pixel_didymos[i]) || ismissing(x_pixel_dimorphos[i])
            x_pixel_didymos[i] = missing
            y_pixel_didymos[i] = missing
            x_pixel_dimorphos[i] = missing
            y_pixel_dimorphos[i] = missing
        end
    end
    usable_photos = length(collect(skipmissing(x_pixel_dimorphos)))

    # Add centroid pixel error
    x_pixel_dimorphos = x_pixel_dimorphos .+ x_centroid_pixel_error
    y_pixel_dimorphos = y_pixel_dimorphos .+ y_centroid_pixel_error

    # Perform orbit determination
    if dimorphos_orbit_type == "Normal"
        # Full OD problem, 6 orbital elements
        bounds = [1190.0-30 0.0001 0.0 0.0 0.0 0.0; 
                  1190.0+30 0.1 90.0 359.99 359.99 359.99]
        # Minimize square mean error to find best orbital elements
        result = optimize(circular_residuals, bounds, ECA(options = Options(debug=true, iterations=100)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = minimizer(result)[4]
        omega_final = minimizer(result)[5]
        M_final = minimizer(result)[6]
    elseif dimorphos_orbit_type == "EEO"
        # Elliptical Equatorial orbit, OD with 5 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.000001 0.0 0.0;
                  1190.0+30 0.00001 0.0001 359.9 359.99]
        result = optimize(residuals, bounds, PSO(;
        N  = 46,
        C1 = 2.0,
        C2 = 2.0,
        ω  = 0.8,
        information = Information(),options = Options(debug=true, iterations=50)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = Omega_dimorphos
        omega_final = minimizer(result)[4]
        M_final = minimizer(result)[5]
    elseif dimorphos_orbit_type == "CIO"
        # Circular inclined orbit, OD with 5 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.0 0.0 0.0;
                  1190.0+30 0.00001 90.0 359.99 359.99]
        result = optimize(residuals, bounds, PSO(;
        N  = 46,
        C1 = 2.0,
        C2 = 2.0,
        ω  = 0.8,
        information = Information(),options = Options(debug=true, iterations=50)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = minimizer(result)[4]
        omega_final = omega_dimorphos
        M_final = minimizer(result)[5]
    elseif dimorphos_orbit_type == "CEO"
        # Circular equatorial orbit, OD with 4 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.001 0.0;
                  1190.0+30 0.00001 0.1 359.99]
        result = optimize(residuals, bounds, ECA(options = Options(debug=true, iterations=50)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = Omega_dimorphos
        omega_final = omega_dimorphos
        M_final = minimizer(result)[4]
    end
    
    # Save orbital elements plots vs time
    original_orbital_elements = [a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos]
    fitted_orbital_elements = [a_final, e_final, i_final, Omega_final, omega_final, M_final]
    propagate_and_plot_orbital_elements_vs_time(original_orbital_elements, fitted_orbital_elements)

    # Compute percentage error of the predicted orbit relative to the initial orbit 
    initial_elements = [a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos]
    final_elements = [a_final, e_final, i_final, Omega_final, omega_final, M_final]
    println(final_elements)
    compute_percentage_error(initial_elements, final_elements)
    compute_mean_absolute_percentage_error(initial_elements, final_elements)

    # Print number of usable photos
    println("Total photos: " * string(total_photos))
    println("Usable photos: " * string(usable_photos))
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

    display(plt_3d)
end

function cb(os)
    println(os.metadata)
    return false    
end

load_hera_spice_kernels()

# System properties
G = 6.67430*10^-11
mass_didymos = 5.32*10^11
global radius_didymos = 390 # [m]
global J2_didymos = 0.011432722
global mass_dimorphos = 4.94*10^9
global mu_system = G*(mass_didymos+mass_dimorphos)
global c = 3.0*10^8
global enable_perturbation = true

# Time properties
hour = 3600.0
number_of_orbits = 10
# Photos to be taken 
photos_per_orbit = 50
global total_photos = photos_per_orbit * number_of_orbits
global start_time = 0.0
global end_time = 120*hour
global step_size = floor((end_time-start_time)/total_photos)
global spice_start_time = utc2et("2027-01-28T08:14:58")
global flux = 1358

# Errors definition
pixel_error = 4
pixel_error_distribution = Uniform(-pixel_error, pixel_error)
global x_centroid_pixel_error = Int64(round(rand(pixel_error_distribution)))
global y_centroid_pixel_error = Int64(round(rand(pixel_error_distribution)))
global error_barycenter = 10/1000 # [km]
global barycenter_error_distribution = Normal(0.0, error_barycenter)
global error_hera_position = 10/1000 # [km]
global hera_error_distribution = Normal(0.0, error_hera_position)
global pointing_error = deg2rad(1)
global pointing_error_distribution = Normal(0.0, pointing_error)
global random_error = zeros(total_photos) #rand(pointing_error_distribution, total_photos)
global random_axis = rand(1:3, total_photos)

main()

# Free used kernels
kclear()