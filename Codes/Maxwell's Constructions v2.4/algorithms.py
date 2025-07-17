import numpy as np
from constants import lambd, x_max, x_min, mu_tilde_max, mu_tilde_min, step_mu_tilde, N_samples, step_mu_tilde_homogeneous
import functions as fun

# not exactly necessary but for testing
import matplotlib.pyplot as plt

def Make_Maxwells_Constructions(delta, N_modes):
    #print(f'for {delta} looking for the solutions')
    # Get the data of densities, chemical potentials, mu_tilde, x
    n_solutions, mu_solutions, mu_tilde_solutions, x_solutions = construct_data_mutilde_n_x(delta, N_modes)
    #print(f'for {delta} looking for extremas')
    # find extremas (n_extrema_1, mu_extrema_1, n_extrema_2, mu_extrema_2, n_phase_change, which phase transition)
    extremas = find_extrema_points(n_solutions, mu_solutions, mu_tilde_solutions, x_solutions)

    # Assign the array for solutions
    maxwell_points = []

    # Initiate Algorithm and store data
    for extrema in extremas:
        
        # Find all the points of interest
        n_phase_separation_1, n_phase_separation_2, mu_star = Minimization_of_Area(mu_solutions, n_solutions, extrema)
        n_extrema_1 = extrema[0]
        mu_extrema_1 = extrema[1]
        n_extrema_2 = extrema[2]
        mu_extrema_2 = extrema[3]
        n_phase_change = extrema[4]
        phase_separation = extrema[5]

        # Store data
        maxwell_points.append( (
            delta, phase_separation,                                    # Information
            n_phase_separation_1, n_phase_separation_2, mu_star,        # Phase Separation Points
            n_extrema_1, mu_extrema_1, n_extrema_2, mu_extrema_2,   # extrema points
            n_phase_change,                                             # For plotting colors
            n_solutions, mu_solutions                                   # for plotting mu vs n
        )
        )

    return maxwell_points

def log_decade_sampling(x_lower, x_higher):
    """
    Generate logarithmically structured sampling values between x_lower and x_higher
    without duplicating values like 10.0, 100.0, etc.
    """
    #decade_steps = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    decade_steps = np.linspace(1, 10, 100)
    log_min = int(np.floor(np.log10(x_lower)))
    log_max = int(np.ceil(np.log10(x_higher)))

    x = []
    for exp in range(log_min, log_max + 1):
        # Skip last step (10) if it's not the final exponent to avoid duplication
        steps = decade_steps if exp == log_max else decade_steps[:-1]
        values = steps * (10 ** exp) 
        x.extend(values)

    x = np.array(x)
    x = x[(x >= x_lower) & (x <= x_higher)]

    # insert the zero value
    x = np.insert(x, 0, 0.0)

    return x

def find_global_minima_free_energy(delta, mu_tilde, N_modes):

    # generate a grid
    x_grid = log_decade_sampling(x_min*lambd, x_max*lambd)

    # calculate the free energy
    free_energy = []
    for x in x_grid:
        free_energy.append(fun.Free_Helmholtz_Energy(delta, mu_tilde, x, N_modes))
    
    # make the plot nicer
    min_value = min(free_energy)

    # sometimes the solution is really small which would lead the algorithm to think that the phase is dimerized when it actually isn't
    # make an algorithm that solves this issue?
    # Not necessary if x is 'big enough'

    return x_grid[free_energy.index(min_value)]

def assign_array_mu_tilde(step):
    mu_tilde = mu_tilde_min
    array = []
    j = 0
    while mu_tilde < mu_tilde_max:
        mu_tilde = mu_tilde_min + j * step
        if mu_tilde == 2.0:
            j = j + 1
            continue
        array.append(mu_tilde)
        j = j + 1
    return array

def construct_data_mutilde_n_x(delta, N_modes):
    '''
    Returns lists of n, mu, mu_tilde, x
    '''

    # assign mu_tilde array
    mu_tilde_array = assign_array_mu_tilde(step=step_mu_tilde)

    # Arrays of solutions to store data
    n_solutions, mu_solutions, mu_tilde_solutions, x_solutions = [],[],[],[]

    # Find where is the last dimer
    for mu_tilde in mu_tilde_array:

        # Find the solution
        x_global_minima = find_global_minima_free_energy(delta, mu_tilde, N_modes)

        # Stop if the dimerized phase is not formed anymore
        if fun.determine_phase(x_global_minima, mu_tilde) != 'dimer':
            break
        
        # Assign data for dimer
        n = fun.density(delta, mu_tilde, x_global_minima, N_modes)
        n_solutions.append(n)
        mu_solutions.append(mu_tilde - lambd*n)
        mu_tilde_solutions.append(mu_tilde)
        x_solutions.append(x_global_minima)

        #print(f'n = {n} | mu_tilde = {mu_tilde} | x = {x_global_minima}')

    # Last density of the dimerized phase before the homogeoneous starts
    if n_solutions != []:
        n_last_dimer = n_solutions[-1]
    else:
        n_last_dimer = -1 # means the dimerized phase never formed

    # Assign data for the homogeneous phase
    mu_tilde_array = assign_array_mu_tilde(step=step_mu_tilde_homogeneous) # When the phase is homogeneous you need a higher step in mu_tilde
    for mu_tilde in mu_tilde_array:
        n = fun.density(delta, mu_tilde, 0.0, N_modes)
        if n >= n_last_dimer:
            #print(f'{mu_tilde:.8f}')
            n_solutions.append(n)
            mu_solutions.append(mu_tilde - lambd*n)
            mu_tilde_solutions.append(mu_tilde)
            x_solutions.append(0.0) 

    return n_solutions, mu_solutions, mu_tilde_solutions, x_solutions

def find_extrema_points(n_solutions, mu_solutions, mu_tilde_solutions, x_solutions):
    '''
    Finds the pair points of extrema points in the corresponding phase change
    '''
    extrema_points = []

    phase_change_densities = find_phase_change_densities(n_solutions, x_solutions, mu_tilde_solutions)

    for i in range(0, len(phase_change_densities)):
        n_phase_change = phase_change_densities[i][0]
        phase_change = phase_change_densities[i][1]
    
        # Find the First extrema Point
        filtered = [(mu, n) for mu, n in zip(mu_solutions, n_solutions) if n <= n_phase_change]
        if filtered:
            mu_extrema_1, n_extrema_1 = max(filtered, key=lambda x: x[0])
        else:
            mu_extrema_1, n_extrema_1 = None, None

        # Find the Second extrema Point
        filtered = [(mu, n) for mu, n in zip(mu_solutions, n_solutions) if n >= n_phase_change and mu <= mu_extrema_1]

        if filtered:
            mu_extrema_2, n_extrema_2 = min(filtered, key=lambda x: x[0])
        else:
            mu_extrema_2, n_extrema_2 = None, None

        extrema_points.append((n_extrema_1, mu_extrema_1, n_extrema_2, mu_extrema_2, n_phase_change, phase_change))
        

    return extrema_points

def find_phase_change_densities(n_solutions, x_solutions, mu_tilde_solutions):

    phase_change_densities = []

    old_phase = fun.determine_phase(x_solutions[0], mu_tilde_solutions[0])

    for i in range(0, len(x_solutions)):
        current_phase = fun.determine_phase(x_solutions[i], mu_tilde_solutions[i])
        if old_phase != current_phase:
            phase_change_densities.append((n_solutions[i], f'{old_phase}'+'-to-'+f'{current_phase}'))
            old_phase = current_phase
    
    return phase_change_densities

def assign_mu_star_array(mu_extrema_1, mu_extrema_2):
    '''
    Assigns values to mu_star to search the solution where the areas are equal bellow and above the line
    Returns, array 
    '''
    array = []
    step = ( mu_extrema_1 - mu_extrema_2)/(N_samples-1)
    for i in range(0,N_samples):
        array.append(mu_extrema_2 + step*i)
    return array

def given_mu_star_find_integration_limits(mu_solutions, n_solutions, extrema, mu_star):
    '''
    Finds the integration points n1 and n2 corresponding to mu_star for the area study
    '''
    n_extrema_1 = extrema[0]
    n_extrema_2 = extrema[2]

    # Find the lower integration limit n1-------------------------------------------
    # Filter for indices where n[i] <= n_extrema
    valid_indices = [i for i in range(len(n_solutions)) if n_solutions[i] <= n_extrema_1]

    if not valid_indices:
        return None  # No valid values found

    # Find the index with the closest mu[i] to mu_star among valid indices
    closest_index = min(valid_indices, key=lambda i: abs(mu_solutions[i] - mu_star))

    n1 = n_solutions[closest_index]

    # Find the upper integration limit n2--------------------------------------------
    # Filter for indices where n[i] >= n_extrema_2
    valid_indices = [i for i in range(len(n_solutions)) if n_solutions[i] >= n_extrema_2]

    if not valid_indices:
        return None  # No valid values found

    # Find the index with the closest mu[i] to mu_star among valid indices
    closest_index = min(valid_indices, key=lambda i: abs(mu_solutions[i] - mu_star))

    n2 = n_solutions[closest_index]

    return n1, n2

def Minimization_of_Area(mu_solutions, n_solutions, extrema):

    # Get the extrema Points
    mu_extrema_1 = extrema[1]
    mu_extrema_2 = extrema[3]

    # Assign an array for searching for the mu_star
    mu_star_array = assign_mu_star_array(mu_extrema_1, mu_extrema_2)

    min_area = float('inf')
    best_mu_star = None
    n_phase_separation_1, n_phase_separation_2 = None, None

    for mu_star in mu_star_array:
        n1, n2 = given_mu_star_find_integration_limits(mu_solutions, n_solutions, extrema, mu_star)

        if n1 is None or n2 is None:
            continue  # skip invalid entries

        area = abs(fun.Maxwell_Area_Calculator(n_solutions, mu_solutions, n1, n2, mu_star))

        if area < min_area:
            min_area = area
            best_mu_star = mu_star
            n_phase_separation_1, n_phase_separation_2 = n1, n2

    return n_phase_separation_1, n_phase_separation_2, best_mu_star


# ***************************************************
#               ADDITIONAL ALGORITHMS
# ***************************************************

def array_derivative_F(delta, mu_tilde, N_modes):
    '''
    Calculates the derivative of F in each point of the grid x, except the first and last point
    '''
    x = log_decade_sampling(x_min*lambd, x_max*lambd)
    array_derivative = []
    for i in range(1, len(x) - 1):
        array_derivative.append((fun.Free_Helmholtz_Energy(delta, mu_tilde, x[i+1], N_modes)-fun.Free_Helmholtz_Energy(delta, mu_tilde, x[i-1], N_modes))/(x[i+1]-x[i-1]))

    return array_derivative, x

def find_extremas_of_F(delta, mu_tilde, N_modes):
    '''
    Given a F(x,delta,mu_tilde) with only x variable, it finds if there exists a local minima in F(x) for x>0
    if there is none, it returs None
    '''
    dF_dx, x = array_derivative_F(delta, mu_tilde, N_modes)

    x_local_maxima, x_local_minima = None, None

    for i in range(1, len(dF_dx)):
        if dF_dx[i] < 0 and dF_dx[i-1] >= 0 and x_local_maxima == None:
            x_local_maxima = x[i]
        if dF_dx[i] > 0 and dF_dx[i-1] <= 0 and x_local_minima == None:
            x_local_minima = x[i]
        # exit loop strategy:
        if x_local_maxima != None and x_local_minima != None:
            break

    return x_local_minima, x_local_maxima

def find_local_minimum_of_F(x_grid, delta, mu_tilde, N_modes):
    y_values = []
    for i, x in enumerate(x_grid):
        y_values.append(fun.Free_Helmholtz_Energy(delta, mu_tilde, x, N_modes))

    x_local_minimum = None

    for i in range(1, len(x_grid) - 1):
        if y_values[i] < y_values[i - 1] and y_values[i] < y_values[i + 1]:
            x_local_minimum = x_grid[i]
    return x_local_minimum

def find_local_maximum_of_F(x_grid, delta, mu_tilde, N_modes):
    y_values = []
    for i, x in enumerate(x_grid):
        y_values.append(fun.Free_Helmholtz_Energy(delta, mu_tilde, x, N_modes))

    x_local_maximum = None

    for i in range(1, len(x_grid) - 1):
        if y_values[i] > y_values[i - 1] and y_values[i] > y_values[i + 1]:
            x_local_maximum = x_grid[i]
    return x_local_maximum

def find_homogeneous_spinodals_in_dimer(delta, N_modes):

    mu_tilde = assign_array_mu_tilde(step=step_mu_tilde)
    x_grid = log_decade_sampling(x_min*lambd, x_max*lambd)

    delta_spinodals_upper, delta_spinodals_lower = [], []
    n_spinodals_upper, n_spinodals_lower = [], []
    mu_spinodals_lower, mu_spinodals_upper = [], []

    delta_upper_last = float('inf')
    delta_lower_last = float('-inf')

    for i in range(0,len(mu_tilde)):

        dimer_exists = False
        spinodal_found = False

        # search for the upper homogeneous spinodal point
        print('\n\nWe are in the upper part of the algorithm')
        for j in range(len(delta) -1, -1, -1):

            if delta[j] > delta_upper_last:
                continue

            current_phase = fun.determine_phase( find_global_minima_free_energy(delta[j], mu_tilde[i], N_modes), mu_tilde[i])
            #print(f'mu_tilde = {mu_tilde[i]} | delta = {delta[j]} | current_phase = {current_phase}')

            if current_phase == 'dimer':
                # Check if a dimer ever formed
                if dimer_exists == False:
                    dimer_exists = True
                # Check if it is a spinodal point
                x_local_maximum = find_local_maximum_of_F(x_grid, delta[j], mu_tilde[i], N_modes)
                if x_local_maximum == None:
                    delta_spinodals_upper.append(delta[j])
                    n_spinodals_upper.append(fun.density(delta[j], mu_tilde[i], 0.0, N_modes))
                    mu_spinodals_upper.append(mu_tilde[i])
                    print('I found a spinodal upper!\n')
                    spinodal_found = True
                    # for making the algorithm faster
                    delta_upper_last = delta[j]
                    break
        
        if dimer_exists == False or spinodal_found == False:
            break

        # search for the lower homogeneous spinodal point
        print('We are in the lower part')
        for j in range(0, len(delta)):

            if delta[j] < delta_lower_last:
                continue

            current_phase = fun.determine_phase( find_global_minima_free_energy(delta[j], mu_tilde[i], N_modes), mu_tilde[i])
            #print(f'mu_tilde = {mu_tilde[i]} | delta = {delta[j]} | current_phase = {current_phase}')

            if current_phase == 'dimer':
                x_local_maximum = find_local_maximum_of_F(x_grid, delta[j], mu_tilde[i], N_modes)
                if x_local_maximum == None:
                    delta_spinodals_lower.append(delta[j])
                    n_spinodals_lower.append(fun.density(delta[j], mu_tilde[i], 0.0, N_modes))
                    mu_spinodals_lower.append(mu_tilde[i])
                    print('I found a spinodal lower!\n')
                    # for making the algorithm faster
                    delta_lower_last = delta[j]
                    break

    # Pairing
    paired_spinodals_upper = list(zip(delta_spinodals_upper, n_spinodals_upper, mu_spinodals_upper))
    paired_spinodals_lower = list(zip(delta_spinodals_lower, n_spinodals_lower, mu_spinodals_lower))

    # Sorting
    paired_spinodals_upper = sorted(paired_spinodals_upper, key=lambda x: (-x[0], x[1]))
    paired_spinodals_lower = sorted(paired_spinodals_lower, key=lambda x: (-x[0], -x[1]))

    # Unpairing
    delta_spinodals_upper, n_spinodals_upper, mu_spinodals_upper = zip(*paired_spinodals_upper)
    delta_spinodals_lower, n_spinodals_lower, mu_spinodals_lower = zip(*paired_spinodals_lower)

    # Assuming these are lists
    combined_delta_spinodals = delta_spinodals_upper + delta_spinodals_lower
    combined_n_spinodals = n_spinodals_upper + n_spinodals_lower
    combined_mu_spinodals = mu_spinodals_upper + mu_spinodals_lower

    with open(f'spinodal_homogeneous_lambda={lambd}.dat', 'w') as f:
        for delta, n, mu in zip(combined_delta_spinodals, combined_n_spinodals, combined_mu_spinodals):
            f.write(f"{n} {delta} {mu}\n")

    plt.plot(combined_n_spinodals, combined_delta_spinodals, marker = 'x', linestyle = '-')
    #plt.plot(n_spinodals_upper, delta_spinodals_upper, marker = 'x', linestyle='-')
    #plt.plot(n_spinodals_lower, delta_spinodals_lower, marker = 'x', linestyle = '-')
    plt.show()

    plt.close()

    plt.plot(combined_mu_spinodals, combined_delta_spinodals, marker = 'x', linestyle = '-')
    plt.show()

    return combined_delta_spinodals, combined_n_spinodals

def find_dimer_spinodals_in_homogeneous(delta, N_modes):

    mu_tilde = assign_array_mu_tilde(step=step_mu_tilde)
    x_grid = log_decade_sampling(x_min*lambd, x_max*lambd)

    delta_spinodals_upper, delta_spinodals_lower = [], []
    n_spinodals_upper, n_spinodals_lower = [], []
    mu_spinodals_lower, mu_spinodals_upper = [], []

    delta_upper_last = float('inf')
    delta_lower_last = float('-inf')

    for i in range(0,len(mu_tilde)):

        homogeneous_exists = False
        spinodal_found = False

        # search for the upper dimer spinodal point
        print('\n\nWe are in the upper part of the algorithm')
        for j in range(len(delta) -1, -1, -1):

            if delta[j] > delta_upper_last:
                continue

            current_phase = fun.determine_phase( find_global_minima_free_energy(delta[j], mu_tilde[i], N_modes), mu_tilde[i])
            #print(f'mu_tilde = {mu_tilde[i]} | delta = {delta[j]} | current_phase = {current_phase}')

            if current_phase == 'topo':
                # Check if a dimer ever formed
                if homogeneous_exists == False:
                    homogeneous_exists = True
                # Check if it is a spinodal point
                x_local_minimum = find_local_minimum_of_F(x_grid, delta[j], mu_tilde[i], N_modes)
                if x_local_minimum != None:
                    delta_spinodals_upper.append(delta[j])
                    n_spinodals_upper.append(fun.density(delta[j], mu_tilde[i], 0.0, N_modes))
                    mu_spinodals_upper.append(mu_tilde[i])
                    print('I found a spinodal upper!\n')
                    spinodal_found = True
                    # for making the algorithm faster
                    delta_upper_last = delta[j]
                    break
        
        if homogeneous_exists == False or spinodal_found == False:
            break

        # search for the lower homogeneous spinodal point
        #print('We are in the lower part')
        for j in range(0, len(delta)):

            break

            if delta[j] < float('inf'):
                continue

            current_phase = fun.determine_phase( find_global_minima_free_energy(delta[j], mu_tilde[i], N_modes), mu_tilde[i])
            #print(f'mu_tilde = {mu_tilde[i]} | delta = {delta[j]} | current_phase = {current_phase}')

            if current_phase == 'topo':
                x_local_minimum = find_local_minimum_of_F(x_grid, delta[j], mu_tilde[i], N_modes)
                if x_local_minimum != None:
                    delta_spinodals_lower.append(delta[j])
                    n_spinodals_lower.append(fun.density(delta[j], mu_tilde[i], x_last_minimum, N_modes))
                    mu_spinodals_lower.append(mu_tilde[i])
                    print('I found a spinodal lower!\n')
                    # for making the algorithm faster
                    delta_lower_last = delta[j]
                    break
                x_last_minimum = x_local_minimum

    # Assuming these are lists
    combined_delta_spinodals = delta_spinodals_upper + delta_spinodals_lower
    combined_n_spinodals = n_spinodals_upper + n_spinodals_lower
    combined_mu_spinodals = mu_spinodals_upper + mu_spinodals_lower

# Pair them
    paired_spinodals = list(zip(combined_delta_spinodals, combined_n_spinodals, combined_mu_spinodals))

# Sort by decreasing delta, then decreasing n
    sorted_spinodals = sorted(paired_spinodals, key=lambda x: (-x[0], x[1]))

# Unzip if needed
    sorted_delta_spinodals, sorted_n_spinodals, sorted_mu_spinodals = zip(*sorted_spinodals)

    with open(f'spinodal_dimer_lambda={lambd}.dat', 'w') as f:
        for delta, n, mu in zip(sorted_delta_spinodals, sorted_n_spinodals, sorted_mu_spinodals):
            f.write(f"{n} {delta} {mu}\n")

    plt.plot(sorted_n_spinodals, sorted_delta_spinodals, marker = 'x', linestyle = '-')
    plt.show()

    plt.close()

    plt.plot(sorted_mu_spinodals, sorted_delta_spinodals, marker = 'x', linestyle = '-')
    plt.show()

    return combined_delta_spinodals, combined_n_spinodals

            








