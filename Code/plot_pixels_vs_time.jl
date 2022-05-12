using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics
include("orbital_utilities.jl")
include("propagators.jl")
include("spice_utilities.jl")
include("camera_utilities.jl")
include("optimization.jl")
include("results.jl")
include("plotting.jl")

function main(original, fitted)
    # Extract original orbit
    a_original = original[1]
    e_original = original[2] 
    i_original = original[3]
    Omega_original = original[4] 
    omega_original = original[5]
    M_original = original[6] 

    # Extract fitted orbit
    a_fitted = fitted[1]
    e_fitted = fitted[2] 
    i_fitted = fitted[3]
    Omega_fitted = fitted[4] 
    omega_fitted = fitted[5]
    M_fitted = fitted[6] 

    # Get pixel points with no errors (original orbit)
    current_error = false
    x_didymos_no_error, y_didymos_no_error, x_dimorphos_no_error, y_dimorphos_no_error, time = propagate_and_compute_pixel_points(a_original, e_original, i_original, Omega_original, omega_original, M_original, current_error)

    # Get pixel points with errors (original orbit)
    current_error = true
    x_didymos_error, y_didymos_error, x_dimorphos_error, y_dimorphos_error, time = propagate_and_compute_pixel_points(a_original, e_original, i_original, Omega_original, omega_original, M_original, current_error)

    # Get pixel points with errors (fitted orbit)
    current_error = true
    x_didymos_fitted_error, y_didymos_fitted_error, x_dimorphos_fitted_error, y_dimorphos_fitted_error, time = propagate_and_compute_pixel_points(a_fitted, e_fitted, i_fitted, Omega_fitted, omega_fitted, M_fitted, current_error)

    # Plot image 
    plotlyjs()
    plot_x_pixels = scatter(time, x_didymos_no_error, label = "Didymos (without error)")
    scatter!(plot_x_pixels, time, x_didymos_error, label = "Didymos (with error)")
    #scatter!(plot_x_pixels, time, x_fitted_error, label = "Fitted (with error)")

    display(plot_x_pixels)
end

function propagate_and_compute_pixel_points(a_dimorphos, e_dimorphos, i_dimorphos, Omega_dimorphos, omega_dimorphos, M_dimorphos, error)
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
        t_vector = LinRange(start_time, end_time, total_photos)
    
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
                if error == true
                    barycenter_offset[j] = rand(barycenter_error_distribution) + rand(hera_error_distribution)
                else 
                    barycenter_offset[j] = 0 
                end
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
            if error == true
                pointing_error_matrix = error_rotation_matrix(random_error[iteration], random_axis[iteration])
                camera_position_didymos[iteration, :] = pointing_error_matrix * didymos_coordinates[iteration, :]
                camera_position_dimorphos[iteration, :] = pointing_error_matrix * dimorphos_coordinates[iteration, :]
            else
                camera_position_didymos[iteration, :] = didymos_coordinates[iteration, :]
                camera_position_dimorphos[iteration, :] = dimorphos_coordinates[iteration, :]
            end
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
        if error == true 
            x_pixel_dimorphos = x_pixel_dimorphos .+ x_centroid_pixel_error
            y_pixel_dimorphos = y_pixel_dimorphos .+ y_centroid_pixel_error
        end
        return x_pixel_didymos, y_pixel_didymos, x_pixel_dimorphos, y_pixel_dimorphos, t_vector
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
number_of_orbits = 25
# Photos to be taken 
photos_per_orbit = 40
global total_photos = 50
global start_time = 0.0
global end_time = 300*hour
global step_size = floor((end_time-start_time)/total_photos)
global spice_start_time = utc2et("2027-02-28T08:14:58")
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
global random_error = rand(pointing_error_distribution, total_photos)
global random_axis = rand(1:2, total_photos)

original_elements = [1183.593 0.00000098 0.035 0.0 0.0 172.18]
fitted_elements = [1183.5823172254468 3.212248590433731e-6 0.04662985538588903 0.0 0.0 172.12056152365003]

main(original_elements, fitted_elements)

# Free used kernels
kclear()