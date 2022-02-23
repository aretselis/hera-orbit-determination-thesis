function runge_kutta_4(x_0, y_0, z_0, vx_0, vy_0, vz_0, mu, t_start, t_end, step, perturbation)
    function two_body_perturbed!(du, u, p, t)
        x, y, z, v_x, v_y, v_z = u
        du[1] = v_x
        du[2] = v_y
        du[3] = v_z
        du[4] = -mu*x/(sqrt(x^2 + y^2 + z^2)^3) - ((3*J2_didymos*mu*(radius_didymos^2)*x)/(2*(sqrt(x^2 + y^2 + z^2)^5)) * (1 - ((5*z^2)/(sqrt(x^2 + y^2 + z^2)^2))))
        du[5] = -mu*y/(sqrt(x^2 + y^2 + z^2)^3) - ((3*J2_didymos*mu*(radius_didymos^2)*y)/(2*(sqrt(x^2 + y^2 + z^2)^5)) * (1 - ((5*z^2)/(sqrt(x^2 + y^2 + z^2)^2))))
        du[6] = -mu*z/(sqrt(x^2 + y^2 + z^2)^3) - ((3*J2_didymos*mu*(radius_didymos^2)*z)/(2*(sqrt(x^2 + y^2 + z^2)^5)) * (3 - ((5*z^2)/(sqrt(x^2 + y^2 + z^2)^2))))
    end

    rv_initial = [x_0, y_0, z_0, vx_0, vy_0, vz_0]
    t_span = (t_start, t_end)

    prob = ODEProblem(two_body_perturbed!,rv_initial,t_span)
    integrator = init(prob, Tsit5(), dt=10, adaptive=false)
    sol = solve!(integrator)
    x_vector = sol[1, :]
    y_vector = sol[2, :]
    z_vector = sol[3, :]
    vx_vector = sol[4, :]
    vy_vector = sol[5, :]
    vz_vector = sol[6, :]
    t_vector = sol.t[:]
    return x_vector, y_vector, z_vector, vx_vector, vy_vector, vz_vector, t_vector
end