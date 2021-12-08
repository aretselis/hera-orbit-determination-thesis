function eccentric_anomaly_calculator(mean_anomaly, e)
    #=
    Compute E using Newton-Raphson
    Input:
    mean_anomaly, [rad]
    e, eccentricity []
    Output
    E, eccentric anomaly at desired time [rad]
    =#

    # Define initial value for E_0)
    if mean_anomaly < pi
        E_0 = mean_anomaly - e
    else
        E_0 = mean_anomaly + e
    end
    # Define f and f dot
    f(E) = mean_anomaly - E + e*sin(E)
    fdot(E) =  -1 + e*cos(E)
    # Stopping criteria
    N = 15  # Number of significant digits to be computed
    max_repetitions = 1000000
    es = 0.5 * 10^(2.0 - N)  # Scarborough Criterion
    ea = 100
    E_prev = E_0
    repetitions = 0
    # Main Newton-Raphson loop
    while ea > es
        repetitions = repetitions + 1
        E_next = E_prev - (f(E_prev) / fdot(E_prev))
        if E_next == 0
            return E_next
        end
        ea = abs((E_next - E_prev) * 100 / E_next)
        E_prev = E_next
        if repetitions > max_repetitions
            error("Max repetitions reached without achieving desired accuracy for E!")
        end
    end
    E = E_next
    return E
end


function orbital_elements_to_cartesian(a, e, i, Omega, omega, M, mu)
    #=
    Computes cartesian position and velocity vector given some orbital elements
    Input:
    a [m]
    e []
    i [deg]
    Omega [deg]
    omega [deg]
    M [deg]
    mu [m^3/(kg*s^2)] (GM)
    Output:
    r_vector, v_vector
    =#

    # Convert M to radians
    M = deg2rad(M)
    # Compute E using Newton-Raphson
    E = eccentric_anomaly_calculator(M, e)
    # Compute x, xdot, y, ydot on the orbital plane
    x = a * (cos(E) - e)
    y = a * sqrt(1 - e^2) * sin(E)
    r = sqrt(x^2 + y^2)
    n = sqrt(mu / (a^3))  # Mean motion
    x_dot = -(n * a^2 / r) * sin(E)
    y_dot = (n * a^2 / r) * sqrt(1 - e^2) * cos(E)
    # Rotation Matrices definition
    Omega = deg2rad(Omega)
    omega = deg2rad(omega)
    i = deg2rad(i)
    P1 = [cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1]
    P2 = [1 0 0;0 cos(i) sin(i);0 sin(i) cos(i)]
    P3 = [cos(Omega) -sin(Omega) 0;sin(Omega) cos(Omega) 0;0 0 1]
    # Compute cartesian coordinates
    x_y_vector = [x, y, 0]
    x_y_dot_vector = [x_dot, y_dot, 0]
    r_vector = *(*(*(P3, P2), P1), x_y_vector)
    v_vector = *(*(*(P3, P2), P1), x_y_dot_vector)
    return r_vector, v_vector
end 


function orbit_evolution_and_cartesian_transform(a, e, i, n, po, tau, GM, time)
    #=
    Computes the orbit at the specified time step_size
    INPUT:
    a, e, i, n, po, tau: Orbital elements
    GM: GM of the System
    time: desired time step
    =#
    # Calculate perihelion distance
    q = a*abs(1.0-e)
    # Calculate longitude of pericenter
    p = n + po
    # Calculate mean motion
    nu = sqrt(GM)*a^(-1.5)
    # Calculate mean anomaly
    L = nu*(time-tau)
    # Compute orbit vectors
    x, y, z, v_x, v_y, v_z = mco_el2x
end