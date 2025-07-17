
# Interactive
lambd = 2.0

# Values for array in delta
delta_min = 1e-2 # If you go at lower delta, the numerical computation starts beeing problematic
delta_max = 0.1 # Delta = 0.0 has a divergent problem, solve this
delta_step = 2.5e-3 # The step can be set to whatever

# Number of modes: Not necessary to modify
N_modes = 1000

# Range of plot for Free Energy
x_min = 1e-5
x_max = 1e6

# for making the grid in mu_tilde
mu_tilde_min = 0.0
mu_tilde_max = 3.0
step_mu_tilde = 1e-3 # Precision depends on lambda, for lambda>2 use 1e-3
step_mu_tilde_homogeneous = step_mu_tilde * 1e-2

# Used Primarily for finding the chemical potential at which the Areas are equal
N_samples = 100