using SPICE, ProgressBars, Plots
include("orbital_utilities.jl")
include("Propagators.jl")

function main()
    # System properties
    G = 6.67430*10^-11
    mass_didymos = 5.32*10^11
    mass_dimorphos = 4.94*10^11
    mu_system = G*(mass_didymos+mass_dimorphos)
    
    # Dimorphos orbit
    a_dimorphos = 1183.0 # Semi-major axis
    e_dimorphos = 0.01 # Eccentricity
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

    # Time properties
    hour = 3600.0
    number_of_steps = 10000
    start_time = 0.0
    end_time = 200*hour
    time = LinRange(start_time, end_time, number_of_steps)
    step_size = (end_time-start_time)/number_of_steps
    
    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, step_size)

    # Extract Hera orbit from SPICE
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)
    # Load all planets, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
    hera_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp"
    furnsh(planets_kernel)
    furnsh(didymos_kernel)
    furnsh(dimorphos_kernel)
    furnsh(hera_kernel)

    # Orbit start and end time
    spice_start_time = utc2et("2027-04-25T08:14:58")
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    hera_coordinates = zeros(Float64, number_of_steps, 3)

    # Main spice calculations loop
    iteration = 1
    println("\nSPICE calculations:")
    for i in ProgressBar(spice_start_time:step_size:spice_end_time-step_size)
        hera_position = spkpos("-999", i, "J2000", "none", "2065803")[1]
        for j in 1:3
            hera_coordinates[iteration, j] = hera_position[j]
        end
        iteration += 1
    end

    # Extract the values we need from the matrices 
    x_hera, y_hera, z_hera = hera_coordinates[:, 1], hera_coordinates[:, 2], hera_coordinates[:, 3]

    # Plot Results
    
    plotlyjs()
    plot3d(x_dimorphos/1000,y_dimorphos/1000,z_dimorphos/1000,title="System Trajectories", xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Dimorphos")
    plot3d!(x_hera,y_hera,z_hera, xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Hera")
end

main()
