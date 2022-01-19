using SPICE, ProgressBars, Plots
include("orbital_utilities.jl")
include("Propagators.jl")
include("errors.jl")

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

    # Orbit start and end time
    spice_start_time = utc2et("2027-02-25T08:14:58")
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    hera_coordinates = zeros(Float64, number_of_steps, 3)
    spatial_resolution = zeros(Float64, number_of_steps)

    # Main spice calculations loop
    iteration = 1
    println("\nSPICE calculations:")
    for i in ProgressBar(spice_start_time:step_size:spice_end_time-step_size)
        hera_position = spkpos("-999", i, "J2000", "none", "2065803")[1]
        for j in 1:3
            hera_coordinates[iteration, j] = hera_position[j]
        end
        distance = sqrt(hera_coordinates[iteration, 1]^2+hera_coordinates[iteration, 2]^2+hera_coordinates[iteration, 3]^2)
        spatial_resolution[iteration] = spatial_resolution_calculator(distance)
        iteration += 1
    end

    # Extract the values we need from the matrices 
    x_hera, y_hera, z_hera = hera_coordinates[:, 1], hera_coordinates[:, 2], hera_coordinates[:, 3]

    # Plot Results
    
    plotlyjs()
    plt_3d = plot3d(x_dimorphos/1000,y_dimorphos/1000,z_dimorphos/1000,title="System Trajectories", xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Dimorphos")
    plot3d!(x_hera,y_hera,z_hera, xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Hera")

    plt_distance = plot(t_vector[1:10000], spatial_resolution, label="Spatial resolution", xlabel="Time", ylabel="Spatial resolution, [m]", title="Spatial resolution vs Time")
    display(plt_distance)
end

main()
