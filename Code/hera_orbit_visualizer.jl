using SPICE, ProgressBars, Plots

function main()
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernerl)

    # Load a planetars, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymain_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    hera_kernerl = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp"
    furnsh(planets_kernel)
    furnsh(didymain_kernel)
    furnsh(hera_kernerl)

    # Define 1 et_minute
    et_start = utc2et("2027-04-25T09:47:00")
    et_end = utc2et("2027-04-25T09:48:00")
    et_minute = et_end-et_start

    x = []
    y = []
    z = []

    # Orbit start and end time
    start_time = utc2et("2027-04-25T08:14:58")
    end_time = utc2et("2027-06-24T08:14:58")

    for i in ProgressBar(start_time:100*et_minute:end_time)
        position = spkpos("-999", i, "J2000", "none", "2065803")[1]
        push!(x, position[1])
        push!(y, position[2])
        push!(z, position[3])
    end

    # Plot Hera Static img
    pyplot()
    plot3d(x,y,z,title="Hera Trajectory", xlabel="x", ylabel="y", zlabel="z")

    # Save image as .gif figure (warning, might take a lot of time)
    plt = scatter3d(1, title="Hera Trajectory", xaxis=("x",(-20,20)), yaxis=("y",(-20,20)), zaxis=("z",(-25,5)), markersize=1)
    anim = @animate for i in ProgressBar(1:1:length(x))
        push!(plt, x[i], y[i], z[i])
    end
    gif(anim, "hera_orbit.gif", fps=15)

end

main()