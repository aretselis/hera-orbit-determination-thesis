using SPICE, ProgressBars, Plots, LinearAlgebra, Distributions, Random,  DifferentialEquations, Metaheuristics, LaTeXStrings

function load_hera_spice_kernels()
    #=
    Loads the necessary kernels from SPICE to run thesis related computed

    NOTE: Make sure to clear the kernels after you are done by using kclear()
    =#
    # Load leap seconds kernel
    leap_kernel = "C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls"
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
    ECP_didymos = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_ECP_v01.bsp"
    ECP_dimorpos = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymoon_ECP_v01.bsp"
    ECP_hera = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_ECP_v01.bsp"
    barycenter_spk = "C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\didymos_hor_000101_500101_v01.bsp"
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\fk\\hera_v07.tf")
    furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\fk\\hera_ops_v02.tf")
    furnsh(leap_kernel)
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
    furnsh(ECP_didymos)
    furnsh(ECP_dimorpos)
    furnsh(ECP_hera)
    furnsh(barycenter_spk)
end

function spatial_resolution_calculator(distance)
    #=
    Computes spatial resolution for the Hera Asteroid Framing Camera, based on its angular_resolution
    
    Input:
    # (Float64) Distance to the object from the camera in [km]
    Output:
    # (Float64) Spatial resolution in m/pixel
    =#
    angular_resolution = 94.1*10^-6 # [rad/pixel]
    distance = distance*10^5 # [cm]
    spatial_resolution = angular_resolution*distance / 100
    return spatial_resolution
end

function main()
    # Define 1 et_minute
    et_start = utc2et("2027-04-25T09:47:00")
    et_end = utc2et("2027-04-25T09:48:00")
    et_minute = et_end-et_start

    # Orbit start and end time
    start_time = utc2et("2027-02-28T08:14:58")
    #start_time = utc2et("2027-01-29T08:14:59")
    #end_time = utc2et("2027-03-21T08:14:58")
    end_time = start_time + 330*60*60
    step_size = et_minute
    steps = convert(Int64, ceil((end_time-start_time)/step_size))

    # Initialize position arrays (x, y, z)
    hera_coordinates = zeros(Float64, steps, 3)
    distance = zeros(Float64, steps)
    spatial_resolution = zeros(Float64, steps)
    time = zeros(Float64, steps)
    didymos_coordinates = zeros(Float64, steps, 3)

    # Main spice calculations loop
    iteration = 1
    for i in ProgressBar(start_time:step_size:end_time-1)
        time[iteration] = iteration*60
        hera_position = spkpos("-999", i, "J2000", "none", "2065803")[1]
        didymos_position = spkpos("-658030", i, "J2000", "none", "-999")[1]
        for j in 1:3
            hera_coordinates[iteration, j] = hera_position[j]
            didymos_coordinates[iteration, j] = didymos_position[j]
        end
        iteration += 1
    end

    # Extract the values we need from the matrices 

    x_hera, y_hera, z_hera = hera_coordinates[:, 1], hera_coordinates[:, 2], hera_coordinates[:, 3]
    x_didymos, y_didymos, z_didymos = didymos_coordinates[:, 1], didymos_coordinates[:, 2], didymos_coordinates[:, 3]

    for i=1:length(x_hera)
        distance[i] = sqrt((x_didymos[i])^2 + (y_didymos[i])^2 + (z_didymos[i])^2)
        spatial_resolution[i] = spatial_resolution_calculator(distance[i])
    end
            
    # Plot spatial resolution
    pgfplotsx()
    #plotlyjs()
    spatial_resolution_plot = plot(time, spatial_resolution, xlabel="Time, t, [sec]", ylabel="Spatial resolution, [m]", label="Spatial resolution", widen=true, formatter=:plain, legend = :bottomright, minorgrid=true)
    display(spatial_resolution_plot)
    savefig(".\\Results\\spatial_resolution_plot.pdf")
end

load_hera_spice_kernels()
main()
kclear()