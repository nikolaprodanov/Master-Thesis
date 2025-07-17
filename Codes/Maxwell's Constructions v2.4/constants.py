
# Interactive
lambd = 2.0

# Values for array in delta
delta_step = 2e-3 # The step can be set to whatever
delta_min = 0.008 # If you go at lower delta, the numerical computation starts beeing problematic
delta_max = 0.1 # Delta = 0.0 has a divergent problem, solve this

# Number of modes: Not necessary to modify
N_modes = 1000

# Range of plot for Free Energy
x_min = 1e-5
x_max = 1e6

# for making the grid in mu_tilde
mu_tilde_min = 0.0
mu_tilde_max = 3.0
step_mu_tilde = 1e-4 # Precision depends on lambda, for lambda>2 use 1e-3
step_mu_tilde_homogeneous = step_mu_tilde * 1e-1

# Used Primarily for finding the chemical potential at which the Areas are equal
N_samples = 100

# Used in only n(delta = 0)
epsilon = 1e-50