using SPICE, ProgressBars, Plots

function load_hera_spice_kernels()
    # Extract Hera orbit from SPICE
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
    furnsh(leap_kernel)
    # Load all planets, Didymos-Dimorphos and HERA kernels
    planets_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp"
    didymos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp"
    dimorphos_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_DCP3_v01.bsp"
    didymos_barycenter_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\didymos_hor_000101_500101_v01.bsp"
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
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_ECP_v01.bsp")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\didymos_hor_000101_500101_v01.bsp")
    furnsh(planets_kernel)
    furnsh(didymos_kernel)
    furnsh(dimorphos_kernel)
    furnsh(didymos_barycenter_kernel)
    furnsh(hera_kernel)
    furnsh(hera_ck)
    furnsh(hera_equipment)
    furnsh(hera_pck)
    furnsh(bodies_pck)
    furnsh(hera_sclk)
    furnsh(hera_instruments)
end


function main()
    # Define 1 et_minute
    et_start = utc2et("2027-04-25T09:47:00")
    et_end = utc2et("2027-04-25T09:48:00")
    et_minute = et_end-et_start

    # Orbit start and end time
    start_time = utc2et("2027-02-25T08:14:58")
    #start_time = utc2et("2027-01-29T08:14:59")
    #end_time = utc2et("2027-03-21T08:14:58")
    end_time = start_time + 120*60*60
    step_size = et_minute
    steps = convert(Int64, ceil((end_time-start_time)/step_size))

    # Initialize position arrays (x, y, z)
    hera_coordinates = zeros(Float64, steps, 3)
    didymos_coordinates = zeros(Float64, steps, 3)
    dimorphos_coordinates = zeros(Float64, steps, 3)

    # Main spice calculations loop
    iteration = 1
    for i in ProgressBar(start_time:step_size:end_time-1)
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
    pgfplotsx()
    hera_plot = plot3d(x_hera,y_hera,z_hera, xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", label="Hera", aspect_ratio=:1)
    plot3d!(hera_plot, x_didymos,y_didymos,z_didymos, label="Didymos", aspect_ratio=:1)
    plot3d!(hera_plot, x_dimorphos, y_dimorphos, z_dimorphos, label="Dimorphos", aspect_ratio=:1)
    #display(hera_plot)
    savefig(".\\Results\\ECP_thesis_orbit_plot.pdf")
    #=
    # Save image as .gif figure (warning, might take a lot of time)
    plt = scatter3d(1, title="Hera Trajectory", xaxis=("x",(-20,20)), yaxis=("y",(-20,20)), zaxis=("z",(-25,5)), markersize=1)
    anim = @animate for i in ProgressBar(1:1:length(x_hera))
        push!(plt, x_hera[i], y_hera[i], z_hera[i])
        push!(plt, x_didymos[i], y_didymos[i], z_didymos[i])
        push!(plt, x_dimorphos[i], y_dimorphos[i], z_dimorphos[i])
    end
    gif(anim, "hera_orbit.gif", fps=15)
    =#
end

load_hera_spice_kernels()
main()
kclear()