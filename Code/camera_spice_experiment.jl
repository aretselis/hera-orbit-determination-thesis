using SPICE, ProgressBars, Plots
# using FileIO, GLMakie
include("orbital_utilities.jl")
include("Propagators.jl")


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

    # Propagate Dimorphos
    x_dimorphos, y_dimorphos, z_dimorphos, vx_dimorphos, vy_dimorphos, vz_dimorphos, t_vector = runge_kutta_4(x, y, z, vx, vy, vz, mu_system, start_time, end_time, step_size)
    dimorphos_coordinates = hcat(x_dimorphos/1000, y_dimorphos/1000, z_dimorphos/1000)

    # Extract Hera orbit from SPICE
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)
    # Load all planets, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
    hera_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp"
    hera_ck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\ck\\hera_sc_PO_EMA_20270209_20270727_f20181203_v03.bc"
    hera_equipment = "C:\\Users\\retse\\repos\\hera-data\\kernels\\fk\\hera_v06.tf"
    hera_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\hera_didymos_v03.tpc"
    bodies_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\de-403-masses.tpc"
    hera_sclk = "C:\\Users\\retse\\repos\\hera-data\\kernels\\sclk\\hera_fict_20181203.tsc"
    hera_instruments = "C:\\Users\\retse\\repos\\hera-data\\kernels\\ik\\hera_afc_v03.ti"
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

    # Orbit start and end time (for SPICE)
    spice_start_time = utc2et("2027-04-25T08:14:58")
    spice_end_time = spice_start_time + end_time
    println(spice_end_time)

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, number_of_steps + 1, 3)
    #dimorphos_coordinates = zeros(Float64, number_of_steps, 3)
    spice_dimorphos_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_coordinates = zeros(Float64, number_of_steps + 1, 3)
    barycenter_offset = zeros(Float64, 1, 3)
    barycenter_initial_coordinates = zeros(Float64, 1, 3)

    # Main spice calculations loop
    iteration = 1
    println("\nSPICE calculations:")
    for i in ProgressBar(spice_start_time:step_size:spice_end_time)
        position_didymos, lt = spkpos("-658030", i, "HERA_SPACECRAFT", "None", "-999")
        spice_position_dimorphos, lt = spkpos("-658031", i, "HERA_SPACECRAFT", "None", "-999")
        rotation_frame = pxform("J2000", "HERA_SPACECRAFT", i)
        position_dimorphos = rotation_frame * dimorphos_coordinates[iteration, :]
        position_system_barycenter, lt = spkpos("2065803", i, "HERA_SPACECRAFT", "None", "-999")
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
        for j in 1:3
            didymos_coordinates[iteration, j] = barycenter_offset[j] + position_didymos[j]
            spice_dimorphos_coordinates[iteration, j] = barycenter_offset[j] + spice_position_dimorphos[j]
            barycenter_coordinates[iteration, j] = barycenter_offset[j] + position_system_barycenter[j]
            dimorphos_coordinates[iteration, j] =  position_dimorphos[j] + barycenter_initial_coordinates[j]
        end
        iteration += 1
    end
    kclear()

    # Plot
    
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    x_spice_dimorphos, y_spice_dimorphos, z_spice_dimorphos = spice_dimorphos_coordinates[:, 1], spice_dimorphos_coordinates[:, 2], spice_dimorphos_coordinates[:, 3]
    plotlyjs()
    scatter3d(x_didymos, y_didymos,z_didymos, color = "blue", xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Didymos", markersize = 1)
    scatter3d!(x_dimorphos, y_dimorphos,z_dimorphos, color = "orange", markersize = 1, label = "Dimorphos")
    scatter3d!(x_spice_dimorphos, y_spice_dimorphos, z_spice_dimorphos, color = "red", markersize = 1, label = "SPICE Dimorphos")
    #=
    # Save image as .gif figure (warning, might take a lot of time)
    plt = scatter3d(1, title="Hera Trajectory", xaxis=("x",(-20,20)), yaxis=("y",(-20,20)), zaxis=("z",(-25,5)), markersize=1)
    anim = @animate for i in ProgressBar(1:1:length(x_didymos))
        push!(plt, x_didymos[i], y_didymos[i], z_didymos[i])
        push!(plt, x_dimorphos[i], y_dimorphos[i], z_dimorphos[i])
    end
    gif(anim, "hera_orbit_new_rf.gif", fps=15)
    
    s = Scene()
    didymos = load("C:\\Users\\retse\\repos\\hera-data\\kernels\\dsk\\g_50677mm_rad_obj_dida_0000n00000_v001.obj")
    dimorphos = load("C:\\Users\\retse\\repos\\hera-data\\kernels\\dsk\\g_08438mm_lgt_obj_didb_0000n00000_v002.obj")
    mesh!(s, didymos, color = "gray")
    #translate!(s, 2.0)
    mesh!(s, dimorphos, color = "blue")
    s
    =#
 
end

main()
