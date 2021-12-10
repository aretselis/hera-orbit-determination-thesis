using SPICE, ProgressBars, FileIO, GLMakie
include("orbital_utilities.jl")
include("Propagators.jl")


function main()
    # Time properties
    hour = 3600.0
    number_of_steps = 1000
    start_time = 0.0
    end_time = 2*hour
    time = LinRange(start_time, end_time, number_of_steps)
    step_size = (end_time-start_time)/number_of_steps
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

    hera_ck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\ck\\hera_sc_PO_EMA_20270209_20270727_f20181203_v03.bc"
    hera_equipment = "C:\\Users\\retse\\repos\\hera-data\\kernels\\fk\\hera_v06.tf"
    hera_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\hera_didymos_v03.tpc"
    bodies_pck = "C:\\Users\\retse\\repos\\hera-data\\kernels\\pck\\de-403-masses.tpc"
    hera_sclk = "C:\\Users\\retse\\repos\\hera-data\\kernels\\sclk\\hera_fict_20181203.tsc"
    furnsh(hera_ck)
    furnsh(hera_equipment)
    furnsh(hera_pck)
    furnsh(bodies_pck)
    furnsh(hera_sclk)
    
    time = utc2et("2027-04-25T08:14:58")
    position, lt = spkpos("-658030", time, "HERA_AFC-1", "None", "-999")

    # Orbit start and end time
    spice_start_time = utc2et("2027-04-25T08:14:58")
    spice_end_time = spice_start_time + end_time

    # Initialize position arrays (x, y, z)
    didymos_coordinates = zeros(Float64, number_of_steps, 3)
    dimorphos_coordinates = zeros(Float64, number_of_steps, 3)
    didymos_initial_coordinates = zeros(Float64, 1, 3)
    didymos_offset = zeros(Float64, 1, 3)
    didymos_initial_coordinates = zeros(Float64, 1, 3)
    didymos_offset = zeros(Float64, 1, 3)

    # Main spice calculations loop
    iteration = 1
    println("\nSPICE calculations:")
    for i in ProgressBar(spice_start_time:step_size:spice_end_time-step_size)
        position_didymos, lt = spkpos("-658030", i, "HERA_SPACECRAFT", "None", "-999")
        position_dimorphos, lt = spkpos("-658031", i, "HERA_SPACECRAFT", "None", "-999")
        # Get didymos reference position or translate hera_coordinates
        if iteration == 1
            for j in 1:3
                didymos_initial_coordinates[j] = position_didymos[j]
            end
        else
            for j in 1:3
                didymos_offset[j] = didymos_initial_coordinates[j] - position_didymos[j]
            end
        end
        for j in 1:3
            didymos_coordinates[iteration, j] = position_didymos[j] + didymos_offset[j]
            dimorphos_coordinates[iteration, j] = position_dimorphos[j] + didymos_offset[j]
        end
        iteration += 1
    end
    kclear()

    # Plot
    
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos = dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    scatter(x_didymos, y_didymos,z_didymos, color = "blue")
    scatter(x_dimorphos, y_dimorphos,z_dimorphos, color = "orange")
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
