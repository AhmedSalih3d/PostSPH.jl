## Given Parameters

using Plots

H     = 0.3;        # meter, wave height
T     = 2;          # seconds, period
delta = 0;          # radians,initial phase

## Derived Parameters
# not sure how S0 calcuated
# not sure how L  is calculated
S0    = 0.230913;   # meter, piston stroke
L     = 5.49035;    # meter, wave length
omega = (2*pi)/T;   # radians, angular frequency
k     = (2*pi)/L;   # radians, wave number
d     = 1.2;        # meter, water depth
g     = 9.81;       # gravity, m/s^2
rho   = 1000;       #kg/m^3, density


## Linear Theory - Regular Waves

# Far-field free surface elevation
t = collect(0:0.05:40); #Fixed t just for visual
#etaL = @(x) H/2 * cos(omega*t - k*x + delta);
#fplot(etaL)

# Piston Movement

e1 = (S0/2) * sin.(omega.*t .+ delta);
plot(e1)
## Force
# NOTE sigma = omega, so I use omega from now on
# NOTE in their formula h is water depth, I use d
# NOTE evanescent modes are neglected for now
# Since I could not get 2nd order to work properly,
# I use S0. They are pretty much equal anyways.

# Due to l assumed infinite (I guess domain is taken
# to be very large), second term in bracket in A is
# approximated as zero
A = ((2*omega*S0)/(k*(sinh(2*k*d)+2*k*d))) * sinh(k*d);

# Wave maker force
F = (((rho*omega*A*sinh(k*d)))/k) .* cos.(omega.*t) .+ (rho*g*d^2)/2;

# Work

Work = F .* e1;

plot(F)
