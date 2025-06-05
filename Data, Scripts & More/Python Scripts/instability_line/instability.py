# ****************************************************
#    PEIERLS INSTABILITY LINE FOR THE COUPLED KTV
# ****************************************************
#        What do you want the program to do?
# ****************************************************

make_plot_susceptibility = False
make_plot_delta_vs_lambda = True
make_plot_delta_vs_electron_density = False

# ****************************************************
#                       Constants
# ****************************************************

from constants import N_modes, lambd

# ****************************************************
#                       Libraries
# ****************************************************

import functions
import figures 
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os

def run_parallel(mu_values):
    print(f"Number of CPU cores available: {os.cpu_count()}")
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(functions.find_instability_point, mu_values))
    return results

if __name__ == "__main__":

    # ****************************************************
    #                     Playground
    # ****************************************************

    # Plot the Susceptibility for desired Delta and mu
    figures.plot_susceptibility(make_plot_susceptibility, delta = 0.2, mu = 1.99, N_modes = N_modes)

    # Plot the Instability line for Delta vs Lamba with desired mu
    mu_values = [0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,1.80,1.90,1.95,1.99]
    #mu_values = [2.001, 2.005, 2.010, 2.050, 2.100, 2.200, 2.300, 2.400, 3.000]
    figures.plot_instability_delta_VS_lambda(make_plot_delta_vs_lambda, N_modes, mu_values)

    # ****************************************************
    #   Parallelization Algorithm for finding solutions
    # ****************************************************

    print(f'Making Instability Line for the Phase Diagram: {make_plot_delta_vs_electron_density}')
    if make_plot_delta_vs_electron_density == False:
        exit()

    # Prepare mu_values
    mu_values = functions.construct_mu_values()

    # Find the solutions
    results = run_parallel(mu_values)

    # Sorting Algorithm
    functions.sorting_algorithm(results)
    
    # Save the results in a file
    functions.save_solutions_in_file(lambd, results)

    # Plot the results in a diagram Delta vs electron density
    figures.plot_superior_instability(results)
    

