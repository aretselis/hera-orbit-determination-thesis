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