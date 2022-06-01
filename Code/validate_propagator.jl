using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics, DataFrames, CSV
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

    println("Initial State Vector")
    println(x)
    println(y)
    println(z)
    println(vx)
    println(vy)
    println(vz)

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, total_photos, enable_perturbation)

    vallado_x = CSV.read("x_vector.csv", DataFrame, header=0)
    vallado_y = CSV.read("y_vector.csv", DataFrame, header=0)
    vallado_z = CSV.read("z_vector.csv", DataFrame, header=0)
    vallado_vx = CSV.read("vx_vector.csv", DataFrame, header=0)
    vallado_vy = CSV.read("vy_vector.csv", DataFrame, header=0)
    vallado_vz = CSV.read("vz_vector.csv", DataFrame, header=0)

    pgfplotsx()
    x_plot = plot(t_vector, x_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "X, [m]")
    plot!(x_plot, t_vector, vallado_x, label = "PKEPLER routine")
    y_plot = plot(t_vector, y_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "Y, [m]")
    plot!(y_plot, t_vector, vallado_y, label = "PKEPLER routine")
    z_plot = plot(t_vector, z_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "Z, [m]")
    plot!(z_plot, t_vector, vallado_z, label = "PKEPLER routine")
    vx_plot = plot(t_vector, vx_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "Vx, [m/s]")
    plot!(vx_plot, t_vector, vallado_vx, label = "PKEPLER routine")
    vy_plot = plot(t_vector, vy_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "Vy, [m/s]")
    plot!(vy_plot, t_vector, vallado_vy, label = "PKEPLER routine")
    vz_plot = plot(t_vector, vz_dimorphos, label = "RK4 Two-Body with J2", xlabel = "Time, t, [seconds]", ylabel = "Vz, [m/s]")
    plot!(vz_plot, t_vector, vallado_vz, label = "PKEPLER routine")
    savefig(x_plot, ".\\Results\\validate_x_data.pdf")
    savefig(y_plot, ".\\Results\\validate_y_data.pdf")
    savefig(z_plot, ".\\Results\\validate_z_data.pdf")
    savefig(vx_plot, ".\\Results\\validate_vx_data.pdf")
    savefig(vy_plot, ".\\Results\\validate_vy_data.pdf")
    savefig(vz_plot, ".\\Results\\validate_vz_data.pdf")
end

# System properties
G = 6.67430*10^-11
mass_error = 0.02 # in percentage
mass_error_distribution = Uniform(-mass_error, mass_error)
mass_didymos = 5.32*10^11 + rand(mass_error_distribution)*5.32*10^11
global radius_didymos = 390 # [m]
global J2_didymos = 0.011432722
global mass_dimorphos = 4.94*10^9  + rand(mass_error_distribution)*4.94*10^9
global mu_system = G*(mass_didymos+mass_dimorphos)
global c = 3.0*10^8
global enable_perturbation = true

# Time properties
hour = 3600.0
number_of_orbits = 25
# Photos to be taken 
photos_per_orbit = 40
global total_photos = photos_per_orbit * number_of_orbits
global start_time = 0.0
global end_time = 120*hour
global step_size = floor((end_time-start_time)/total_photos)

main()