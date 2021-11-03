using SPICE, ProgressBars, Plots

# Load leap seconds kernel
furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\lsk\\naif0012.tls")

# Convert the calendar date to ephemeris seconds past J2000
et = utc2et("2027-04-25T09:47:00")
print(et)

# Load a planetary ephemeris kernel
furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\de432s.bsp")
furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_didymain_DCP3_v01.bsp")
furnsh("C:\\Users\\retse\\repos\\hera-data\\kernels\\spk\\HERA_sc_DCP3_v01.bsp")

# Get the position of Mars at `et` w.r.t. Earth

et_start = utc2et("2027-04-25T09:47:00")
et_end = utc2et("2027-04-25T09:48:00")
et_minute = et_end-et_start

x = []
y = []
z = []

start_time = utc2et("2027-04-25T09:47:00")
end_time = utc2et("2027-05-10T19:47:00")

for i in ProgressBar(start_time:et_minute:end_time)
    position = spkpos("-999", i, "J2000", "none", "2065803")[1]
    push!(x, position[1])
    push!(y, position[2])
    push!(z, position[3])
end

# Plot Hera Static img

#pyplot()
#plot3d(x,y,z,title="Hera Trajectory", xaxis=("x",(-10,20)), yaxis=("x",(-15,10)), zaxis=("x",(-15,10)))

# Save image as figure (warning, might take a lot of time)

plt = scatter3d(1, title="Hera Trajectory", xaxis=("x",(-15,20)), yaxis=("y",(-15,10)), zaxis=("z",(-20,5)), markersize=1)

anim = @animate for i in ProgressBar(1:100:length(x))
    push!(plt, x[i], y[i], z[i])
end
gif(anim, "hera_orbit.gif", fps=15)
