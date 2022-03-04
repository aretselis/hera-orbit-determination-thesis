function residuals(orbit)
    # Select elements
    a = orbit[1]
    e = orbit[2]
    i = orbit[3]
    Omega = orbit[4]
    omega = orbit[5]
    M = orbit[6]
    # Compute dimorphos pixel points for given orbit and compute xhi2
    observed_x, observed_y = x_pixel_dimorphos, y_pixel_dimorphos
    predicted_x, predicted_y = propagate_and_compute_dimorphos_pixel_points(a, e, i, Omega, omega, M, start_time, end_time, step_size, spice_start_time)
    xhi2 = sum_of_squared_residuals(observed_x, observed_y, predicted_x, predicted_y)
    return Float64(xhi2)
end


function sum_of_squared_residuals(observed_x, observed_y, predicted_x, predicted_y)
    x_residuals = observed_x - predicted_x
    y_residuals = observed_y - predicted_y
    # Clean up missing values
    x_residuals = collect(skipmissing(x_residuals))
    y_residuals = collect(skipmissing(y_residuals))
    squared_sum = sum(x_residuals.^2) + sum(y_residuals.^2)
    return squared_sum
end


