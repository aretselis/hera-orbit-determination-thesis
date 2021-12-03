using SPICE, ProgressBars, Plots

function main()
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)

    # Load a planetars, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
    hera_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp"
    furnsh(planets_kernel)
    furnsh(didymos_kernel)
    furnsh(dimorphos_kernel)
    furnsh(hera_kernel)

    # Define 1 et_minute
    et_start = utc2et("2027-04-25T09:47:00")
    et_end = utc2et("2027-04-25T09:48:00")
    et_minute = et_end-et_start

    # Orbit start and end time
    start_time = utc2et("2027-04-25T08:14:58")
    end_time = utc2et("2027-05-05T08:14:58")
    step_size = 100*et_minute
    steps = convert(Int64, ceil((end_time-start_time)/step_size))

    # Initialize position arrays (x, y, z)
    hera_coordinates = zeros(Float64, steps, 3)
    didymos_coordinates = zeros(Float64, steps, 3)
    dimorphos_coordinates = zeros(Float64, steps, 3)

    # Main spice calculations loop
    iteration = 1
    for i in ProgressBar(start_time:step_size:end_time)
        hera_position = spkpos("-999", i, "J2000", "none", "2065803")[1]
        didymos_position = spkpos("-658030", i, "J2000", "none", "2065803")[1]
        dimorphos_position = spkpos("-658031", i, "J2000", "none", "2065803")[1]
        for j in 1:3
            hera_coordinates[iteration, j] = hera_position[j]
            didymos_coordinates[iteration, j] = didymos_position[j]
            dimorphos_coordinates[iteration, j] = dimorphos_position[j]
        end
        iteration += 1
    end

    # Extract the values we need from the matrices 

    x_hera, y_hera, z_hera = hera_coordinates[:, 1], hera_coordinates[:, 2], hera_coordinates[:, 3]
    x_didymos,y_didymos,z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]
    x_dimorphos, y_dimorphos, z_dimorphos,= dimorphos_coordinates[:, 1], dimorphos_coordinates[:, 2], dimorphos_coordinates[:, 3]
    
    # Plot Hera Static img
    pyplot()
    plot3d(x_hera,y_hera,z_hera,title="Hera Trajectory", xlabel="x", ylabel="y", zlabel="z", label="Hera")
    plot3d!(x_didymos,y_didymos,z_didymos, label="Didymos")
    plot3d!(x_dimorphos, y_dimorphos, z_dimorphos, label="Dimorphos")

    # Save image as .gif figure (warning, might take a lot of time)
    plt = scatter3d(1, title="Hera Trajectory", xaxis=("x",(-20,20)), yaxis=("y",(-20,20)), zaxis=("z",(-25,5)), markersize=1)
    anim = @animate for i in ProgressBar(1:1:length(x_hera))
        push!(plt, x_hera[i], y_hera[i], z_hera[i])
        push!(plt, x_didymos[i], y_didymos[i], z_didymos[i])
        push!(plt, x_dimorphos[i], y_dimorphos[i], z_dimorphos[i])
    end
    gif(anim, "hera_orbit.gif", fps=15)

end

main()