using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics, LaTeXStrings
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")
include("optimization.jl")
include("results.jl")
include("plotting.jl")

function get_all_pictures(a_value, e_value, i_value, Omega_value, omega_value, M_value)
    # Determine what kind of orbit we are dealing with
    global dimorphos_orbit_type = "CEO"#type_of_orbit(e_dimorphos, i_dimorphos)
    println("Orbit type is "*dimorphos_orbit_type)
    global enable_perturbation = true

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_value, e_value, i_value, Omega_value, omega_value, M_value, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, total_photos, enable_perturbation)
    # Select every xth element to match the photos taken and convert to km
    println(length(x_dimorphos))
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
    global time_vector = t_vector[index_vector]
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
    x_pixel_dimorphos, y_pixel_dimorphos = convert_to_pixels(dimorphos_pixel_coordinates[:, 1], dimorphos_pixel_coordinates[:, 2], x_boundaries, y_boundaries)
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

    # Add centroid pixel error
    original_x_pixel_dimorphos = x_pixel_dimorphos .+ x_centroid_pixel_error
    original_y_pixel_dimorphos = y_pixel_dimorphos .+ y_centroid_pixel_error
    return original_x_pixel_dimorphos, original_y_pixel_dimorphos
end


function main()
    # Determine what kind of orbit we are dealing with
    global dimorphos_orbit_type = "CEO"#type_of_orbit(e_dimorphos, i_dimorphos)
    println("Orbit type is "*dimorphos_orbit_type)
    usable_photos = length(collect(skipmissing(x_pixel_dimorphos)))
    println("Usable photos: " *string(usable_photos))
    global enable_perturbation = false
    # Perform orbit determination
    if dimorphos_orbit_type == "Normal"
        minimum_error = residuals([a_dimorphos e_dimorphos i_dimorphos Omega_dimorphos omega_dimorphos M_dimorphos])
        println("Minimum error is " *string(minimum_error)* " pixels")
        # Full OD problem, 6 orbital elements
        bounds = [1190.0-30 0.0001 0.0 0.0 0.0 0.0; 
                  1190.0+30 0.9 90.0 359.99 359.99 359.99]
        # Minimize square mean error to find best orbital elements
        result = optimize(residuals, bounds, ECA(information = Information(f_optimum = minimum_error), options = Options(debug=true, iterations=100, store_convergence = true)))
        #result = optimize(residuals, bounds, PSO(information = Information(f_optimum = minimum_error),options = Options(debug=true, iterations=100, store_convergence = true)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = minimizer(result)[4]
        omega_final = minimizer(result)[5]
        M_final = minimizer(result)[6]
    elseif dimorphos_orbit_type == "EEO"
        minimum_error = residuals([a_dimorphos e_dimorphos i_dimorphos omega_dimorphos M_dimorphos])
        println("Minimum error is " *string(minimum_error)* " pixels")
        # Elliptical Equatorial orbit, OD with 5 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.001 0.0 0.0;
                  1190.0+30 0.9 0.1 359.9 359.99]
        result = optimize(residuals, bounds, ECA(information = Information(f_optimum = minimum_error), options = Options(debug=true, iterations=100, store_convergence = true)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = Omega_dimorphos
        omega_final = minimizer(result)[4]
        M_final = minimizer(result)[5]
    elseif dimorphos_orbit_type == "CIO"
        minimum_error = residuals([a_dimorphos e_dimorphos i_dimorphos Omega_dimorphos M_dimorphos])
        println("Minimum error is " *string(minimum_error)* " pixels")
        # Circular inclined orbit, OD with 5 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.0 0.0 0.0;
                  1190.0+30 0.00001 90.0 359.99 359.99]
        result = optimize(residuals, bounds, ECA(information = Information(f_optimum = minimum_error), options = Options(debug=true, iterations=100, store_convergence = true)))
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = minimizer(result)[3]
        Omega_final = minimizer(result)[4]
        omega_final = omega_dimorphos
        M_final = minimizer(result)[5]
    elseif dimorphos_orbit_type == "CEO"
        minimum_error = residuals([a_dimorphos e_dimorphos i_dimorphos M_dimorphos])
        println("Minimum error is " *string(minimum_error)* " pixels")
        # Circular equatorial orbit, OD with 4 orbital elements
        # Bounds for optimization
        bounds = [1190.0-30 0.00000001 0.0 0.8*mu_system_original;
                  1190.0+30 0.004 359.99 1.2*mu_system_original]
        metaheuristic_success = false
        #while metaheuristic_success == false
        #    try
                result = optimize(residuals, bounds, ECA(information = Information(f_optimum = minimum_error), options = Options(debug=true, iterations=100, store_convergence = true)))
        #        metaheuristic_success = true
        #    catch
        #        println("Run failed, trying again...")
        #    end
        #end
        # Extract final guess from optimization
        a_final = minimizer(result)[1]
        e_final = minimizer(result)[2]
        i_final = i_dimorphos
        Omega_final = Omega_dimorphos
        omega_final = omega_dimorphos
        M_final = minimizer(result)[3]
        mu_final = minimizer(result)[4]
    end
    return a_final, e_final, i_final, Omega_final, omega_final, M_final, mu_final
end

function start_plot(orbital_elements)
    global enable_perturbation = true
    # Dimorphos orbit
    a_dimorphos = orbital_elements[1] # Semi-major axis
    e_dimorphos = orbital_elements[2] # Eccentricity
    i_dimorphos = orbital_elements[3] # Inclination
    Omega_dimorphos = orbital_elements[4] # Argument of periapsis
    omega_dimorphos = orbital_elements[5] # Longitude of the ascending node
    M_dimorphos = orbital_elements[6] # Mean anomaly

    # Compute initial position and velocity vector from orbital elements
    r_vector, v_vector = orbital_elements_to_cartesian(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, mu_system)
    x = r_vector[1]
    y = r_vector[2]
    z = r_vector[3]
    vx = v_vector[1]
    vy = v_vector[2]
    vz = v_vector[3]
    #global enable_perturbation = true
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
    global plt_a_vs_time = plot(t_vector, a_vector, xlabel = "Time, "* L"t,\ [s]", ylabel = "Semi-major axis, " * L"a,\ [m]", legend=false)
    global plt_e_vs_time = plot(t_vector, e_vector, xlabel = "Time, " * L"t,\ [s]", ylabel = "Eccentricity, " * L"e,\ [\ ]", legend=false)
    global plt_M_vs_time = plot(t_vector, M_vector, xlabel = "Time, "* L"t,\ [s]", ylabel = L"\textnormal{True Longitude},\ \lambda_{true},\ [\deg]", legend=false)
end

function add_to_plots(orbital_elements)
    global enable_perturbation = false
    # Dimorphos orbit
    a_dimorphos = orbital_elements[1] # Semi-major axis
    e_dimorphos = orbital_elements[2] # Eccentricity
    i_dimorphos = orbital_elements[3] # Inclination
    Omega_dimorphos = orbital_elements[4] # Argument of periapsis
    omega_dimorphos = orbital_elements[5] # Longitude of the ascending node
    M_dimorphos = orbital_elements[6] # Mean anomaly

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

    # Compute orbital elements vs time
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
    crop_start_index = findall(x->x==current_timestamp[1], t_vector)[1]
    crop_end_index = findall(x->x==current_timestamp[length(current_timestamp)], t_vector)[1]
    a_vector = a_vector[crop_start_index:crop_end_index]
    e_vector = e_vector[crop_start_index:crop_end_index]
    M_vector = M_vector[crop_start_index:crop_end_index]
    t_vector = t_vector[crop_start_index:crop_end_index]
    plot!(plt_a_vs_time, t_vector, a_vector)
    plot!(plt_e_vs_time, t_vector, e_vector)
    plot!(plt_M_vs_time, t_vector, M_vector)
end

load_hera_spice_kernels()

# System properties
G = 6.67430*10^-11
mass_error = 0.02 # in percentage
mass_error_distribution = Uniform(-mass_error, mass_error)
mass_didymos = 5.32*10^11 + rand(mass_error_distribution)*5.32*10^11
global radius_didymos = 390 # [m]
global J2_didymos = 0.011432722
global mass_dimorphos = 4.94*10^9 + rand(mass_error_distribution)*4.94*10^9
global mu_system = G*(mass_didymos+mass_dimorphos)
print(mu_system)
# Dimorphos initial orbit definition
a_value = 1178.523 # Semi-major axis
e_value = 0.0000078 # Eccentricity
i_value = 0.0078 # Inclination
Omega_value = 0.0 # Longitude of the ascending node
omega_value = 0.0 # Argument of periapsis
M_value = 42.236 # True anomaly 
start_plot([a_value, e_value, i_value, Omega_value, omega_value, M_value])
orbital_period = 2*pi*sqrt(a_value^3/mu_system)
global c = 3.0*10^8
global enable_perturbation = true

# Time properties
hour = 3600.0
number_of_orbits = 25
# Photos to be taken 
photos_per_orbit = 10
global total_photos = 3000#photos_per_orbit * number_of_orbits
global start_time = 0
global end_time = orbital_period
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
global all_x_pixel_dimorphos, all_y_pixel_dimorphos = get_all_pictures(a_value, e_value, i_value, Omega_value, omega_value, M_value)

global random_counter = 1
global start_index = 1
global end_index = 100
pgfplotsx()
global debug = false
global guess_array = zeros(Float64, 30, 4)
for random_counter=1:30
    global enable_perturbation = false
    global x_pixel_dimorphos = all_x_pixel_dimorphos[start_index:end_index]
    global y_pixel_dimorphos = all_y_pixel_dimorphos[start_index:end_index]
    global current_timestamp = time_vector[start_index:end_index]
    a_guess, e_guess, i_guess, Omega_guess, omega_guess, M_guess, mu_guess = main()
    add_to_plots([a_guess, e_guess, i_guess, Omega_guess, omega_guess, M_guess])
    guess_array[random_counter, 1] = a_guess
    guess_array[random_counter, 2] = e_guess
    guess_array[random_counter, 3] = M_guess
    guess_array[random_counter, 4] = mu_guess
    global start_index = start_index + 100
    global end_index = end_index + 100
end

display(plt_e_vs_time)

# Free used kernels
kclear()

global start_index = 1
global end_index = 100
for random_counter=1:30
    output = compute_mean_absolute_percentage_error([a_value, e_value, i_value, Omega_value, omega_value, M_value], [guess_array[random_counter, 1], guess_array[random_counter, 2], i_value, Omega_value, omega_value, guess_array[random_counter, 3]])
    println(output) 
end