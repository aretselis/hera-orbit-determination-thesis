function compute_rotation_matrix(a, b)
    # Make vectors unit vectors
    a = a/sqrt(dot(a, a))
    b = b/sqrt(dot(b, b))
    v = cross(a, b)
    c = dot(a, b)
    I = [1.0 0.0 0.0; 
         0.0 1.0 0.0; 
         0.0 0.0 1.0]
    skew_symmetric_matrix = [0.0 -v[3] v[2]; 
                             v[3] 0.0 -v[1];
                            -v[2] v[1] 0.0]
    rotation_matrix = I + skew_symmetric_matrix + skew_symmetric_matrix*skew_symmetric_matrix*(1/(1+c))
    return rotation_matrix
end


function convert_to_pixels(x_data, y_data, x_boundaries, y_boundaries)
    #=
    Computes corresponding pixels for two vectors of x and y observations, given the corresponding x and y boundaries
    
    Input:
    # (Vector, Float64) x_data, coordinates of x observations
    # (Vector, Float64) y_data, coordinates of y observations
    # (Vector, Float64) x_boundaries, coordinates of x boundaries
    # (Vector, Float64) y_boundaries, coordinates of y boundaries
    Output:
    # (Vector, Float64) x_pixel, vector containing x pixel coordinate of observations
    # (Vector, Float64) y_pixel, vector containing y pixel coordinate of observations
    =#
    
    # Define pixels (assuming rectangular image)
    pixels = 1020
    image_size = pixels + 1 

    # Bottom left pixel as [0, 0]
    x_bottom_left = minimum(x_boundaries)
    y_bottom_left = minimum(y_boundaries)
    # Top left pixel as [0, image_size]
    y_top_left = maximum(y_boundaries)
    # Bottom right pixel as [image_size, 0]
    x_bottom_right = maximum(x_boundaries)

    # Generate pixels 
    x_pixel_range = LinRange(x_bottom_left, x_bottom_right, image_size)
    y_pixel_range = LinRange(y_bottom_left, y_top_left, image_size)
    x_pixel = zeros(Int64, length(x_data))
    y_pixel = zeros(Int64, length(x_data))
    for i=1:length(x_data)
        # Assign x_pixel coordinate
        for j=1:length(x_pixel_range)-1
            if x_data[i] >= x_pixel_range[j] && x_data[i] <= x_pixel_range[j+1]
                x_pixel[i] = j
                break
            end
        end
        # Assign y_pixel coordinate
        for j=1:length(y_pixel_range)-1
            if y_data[i] >= y_pixel_range[j] && y_data[i] <= y_pixel_range[j+1]
                y_pixel[i] = j
                break
            end
        end
    end
    return x_pixel, y_pixel
end