mass_didymos = 5.32*10^11 # [kg]

a = 0.416194 # km
b = 0.418765 # km
c = 0.39309 # km

I_x = 31436196699045452.00000000000000000000  # mass_didymos/5 * (b^2 + c^2)
I_y = 32009618023928268.00000000000000000000 # mass_didymos/5 * (a^2 + c^2)
I_z = 32882223794965416.00000000000000000000 # mass_didymos/5 * (a^2 + b^2)

R = (a + b)/2

C22 = (I_x-I_y) / (-4*mass_didymos*R^2)
C20 = (I_y-I_z)/(mass_didymos*R^2) - 2 * C22

J2 = -C20

println(J2/10^6)