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