function compute_percentage_error(actual_elements, predicted_elements)
    #=
    Computes percentage error given two sets of orbital elements
    Input:
    actual_elements, Vector(Float64, 6) containing actual orbital elements
    predicted_elements, Vector(Float64, 6) containing predicted orbital elements 
    Output:
    percentage_errors, Vector(Float64, 6) containing percentage (%) errors
    =#

    percentage_errors = zeros(Float64, 6)
    for i in 1:6
        if i>=3 && i<=6
            # Ensure conversion to radians for angles
            actual_elements[i] = deg2rad(actual_elements[i])
            predicted_elements[i] = deg2rad(predicted_elements[i])
        end
        percentage_errors[i] = 100 * (abs(predicted_elements[i] - actual_elements[i])/actual_elements[i])
        println("Percentage error for element " * string(i) * " is " * string(percentage_errors[i]))
    end
    return Nothing
end