using SPICE

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
