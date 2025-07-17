# ****************************************************
#                   Definition of variables
# ****************************************************
# The variables are defined as following:
# delta == \Delta / t
# mu == \mu / t                 
# x == ( g * X ) / t
# y == ( g * X_0 ) / t = - lambd * n    , where n is the electron density
# mu_tilde == \tilde{\mu} / t = ( \mu - g * X_0) / t = mu + lambd * n
# The non topological phase is formed for mu_tilde >= 2.0 * t
# ****************************************************
#                    Libraries
# ****************************************************

from constants import N_modes, lambd, delta_min, delta_max, delta_step, x_min, x_max
import functions as fun 
import plotting
import algorithms 

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import time

def run_maxwell_construction(delta):
    start_time = time.time()
    result = algorithms.Make_Maxwells_Constructions(delta, N_modes)
    elapsed = time.time() - start_time
    print(f'Finalized RUN for Delta = {delta:.4f} | process time = {elapsed:.1f} s.')
    return result

if __name__ == "__main__":

    playground_mode = True

    # ****************************************************
    #                      PLAYGROUND
    # ****************************************************

    if playground_mode == True:

        delta, mu_tilde, x = 7e-2, 1e-3, 0.0

        #plotting.plot_spectrum(delta=delta,mu_tilde=mu_tilde,x=x, N_modes = N_modes)
 
        #plotting.plot_Free_Energy_vs_x(delta=delta,mu_tilde=mu_tilde,N_modes=N_modes) 
 
        plotting.mu_vs_n_with_Maxwells_Constructions(delta, N_modes)

        # Instability lines of homogeneous and dimer
        #delta_values = fun.assing_delta_list(delta_min, delta_max, delta_step)
        #algorithms.find_homogeneous_spinodals_in_dimer(delta_values, N_modes)
        #algorithms.find_dimer_spinodals_in_homogeneous(delta_values, N_modes)

        # Exit the program
        exit()

    # ****************************************************
    #                   MAIN PROGRAM
    # ****************************************************

    # Assign Values for Delta
    delta_values = fun.assing_delta_list(delta_min, delta_max, delta_step)

    # Prepare array for storing all data (12 inputs): see in algorithms >> Make_Maxwells_Constructions
    DATA = []

    # Use all available CPUs
    num_processes = multiprocessing.cpu_count()  
    print(f'Initializing Parallel Runs with {num_processes} CPUs...Please wait...')

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # Submit all tasks
        futures = [executor.submit(run_maxwell_construction, delta) for delta in delta_values]

        # Collect results as they complete (unordered)
        for future in as_completed(futures):
            result = future.result()
            DATA.append(result)

    # Sort the data with increasing Delta
    DATA = sorted(DATA, key=lambda x: x[0])

    # Make the phase diagram: delta vs n
    plotting.make_phase_diagram(DATA)


    # ****************************************************
    #                   EXIT PROGRAM
    # ****************************************************
    

