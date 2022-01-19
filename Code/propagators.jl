function runge_kutta_4(x_0, y_0, z_0, vx_0, vy_0, vz_0, mu, t_start, t_end, step)
    #=
    # 4th order Runge-Kutta orbit propagator assuming n no perturbations
    =#

    # Preallocate required vectors
    size = Int64(ceil((t_end-t_start)/step)) + 1
    tn = Vector{Float64}(undef, size)
    tn[1] = t_start
    xn = Vector{Float64}(undef, size)
    xn[1] = x_0
    yn = Vector{Float64}(undef, size)
    yn[1] = y_0
    zn = Vector{Float64}(undef, size)
    zn[1] = z_0 
    vxn = Vector{Float64}(undef, size)
    vxn[1] = vx_0
    vyn = Vector{Float64}(undef, size)
    vyn[1] = vy_0
    vzn = Vector{Float64}(undef, size)
    vzn[1] = vz_0
    counter = 1

    println("\nPropagating:")
    iter = ProgressBar(t_start:step:t_end-step)
    for i in iter
        # Calculate k1 values
        k1_x = fx(vx_0)
        k1_y = fy(vy_0)
        k1_z = fz(vz_0)
        k1_vx = fv_x(x_0, y_0, z_0, mu)
        k1_vy = fv_y(x_0, y_0, z_0, mu)
        k1_vz = fv_z(x_0, y_0, z_0, mu)
        # Calculate midpoint values
        mid_x = x_0 + k1_x * step/2
        mid_y = y_0 + k1_y * step/2
        mid_z = z_0 + k1_z * step/2
        mid_vx = vx_0 + k1_vx * step/2
        mid_vy = vy_0 + k1_vy * step/2
        mid_vz = vz_0 + k1_vz * step/2
        # Calculate k2 values
        k2_x = fx(mid_vx)
        k2_y = fy(mid_vy)
        k2_z = fz(mid_vz)
        k2_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k2_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k2_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Calculate next midpoint values
        mid_x = x_0 + k2_x * step / 2
        mid_y = y_0 + k2_y * step / 2
        mid_z = z_0 + k2_z * step / 2
        mid_vx = vx_0 + k2_vx * step / 2
        mid_vy = vy_0 + k2_vy * step / 2
        mid_vz = vz_0 + k2_vz * step / 2
        # Calculate k3 values
        k3_x = fx(mid_vx)
        k3_y = fy(mid_vy)
        k3_z = fz(mid_vz)
        k3_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k3_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k3_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Calculate next midpoint values
        mid_x = x_0 + k3_x * step
        mid_y = y_0 + k3_y * step
        mid_z = z_0 + k3_z * step
        mid_vx = vx_0 + k3_vx * step
        mid_vy = vy_0 + k3_vy * step
        mid_vz = vz_0 + k3_vz * step
        # Calculate k4 values
        k4_x = fx(mid_vx)
        k4_y = fy(mid_vy)
        k4_z = fz(mid_vz)
        k4_vx = fv_x(mid_x, mid_y, mid_z, mu)
        k4_vy = fv_y(mid_x, mid_y, mid_z, mu)
        k4_vz = fv_z(mid_x, mid_y, mid_z, mu)
        # Compute r, v values and append to list
        xn[counter + 1] = xn[counter] + (step / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        yn[counter + 1] = yn[counter] + (step / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y)
        zn[counter + 1] = zn[counter] + (step / 6) * (k1_z + 2 * k2_z + 2 * k3_z + k4_z)
        vxn[counter + 1] = vxn[counter] + (step / 6) * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx)
        vyn[counter + 1] = vyn[counter] + (step / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy)
        vzn[counter + 1] = vzn[counter] + (step / 6) * (k1_vz + 2 * k2_vz + 2 * k3_vz + k4_vz)
        tn[counter + 1] = tn[counter] + step
        # Prepare for the next iteration and reset values
        counter += 1
        x_0 = xn[counter]
        y_0 = yn[counter]
        z_0 = zn[counter]
        vx_0 = vxn[counter]
        vy_0 = vyn[counter]
        vz_0 = vzn[counter]
    end
    return xn, yn, zn, vxn, vyn, vzn, tn 
end


function fx(v_x)
    # Assuming x'(t)=v_x
    v_x
end


function fy(v_y)
    # Assuming y'(t)=v_y
    v_y
end


function fz(v_z)
    # Assuming z'(t)=v_z
    v_z
end


function fv_x(x, y, z, mu)
    # Assuming v_x'(t)=-mu*x/(sqrt(x^2+y^2+z^2))^3
    -mu*x/(^(^(^(x, 2) + ^(y, 2) + ^(z, 2), 1/2), 3))
end


function fv_y(x, y, z, mu)
    # Assuming v_y'(t)=-mu*y/(sqrt(x^2+y^2+z^2))^3
    -mu*y/(^(^(^(x, 2) + ^(y, 2) + ^(z, 2), 1/2), 3))
end


function fv_z(x, y, z, mu)
    # Assuming v_z'(t)=-mu*z/(sqrt(x^2+y^2+z^2))^3
    -mu*z/(^(^(^(x, 2) + ^(y, 2) + ^(z, 2), 1/2), 3))
end