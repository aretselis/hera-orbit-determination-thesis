using SPICE


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


function get_boundary_vectors(instrument)
    # Get boundary vectors of a specific instrument
    # Input:
    # (Int) Instrument ID in SPICE 
    # Output:
    # 4x3 matrix containing the 4 three-dimensional boundary vectors 
    shape, name, boresight, boundary_vectors = getfov(instrument)
    boundary_vectors_matrix = zeros(Float64, 4, 3)
    for i in 1:4
        boundary_vectors_matrix[i, :] = boundary_vectors[i]
    end
    return boundary_vectors_matrix
end
