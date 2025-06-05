# ****************************************************
#    Functions for the instability_line.py
# ****************************************************
#   a = 1
#   t = 1
#   Volume = (N-1) * a

from constants import lambd, N_modes

g_lambda = lambd
g_N_modes = N_modes

import math
import numpy as np
from scipy.integrate import simpson

def energy(k, delta, mu):
    return np.sqrt( (np.cos(k) - mu/2) ** 2 + (delta*np.sin(k)) ** 2 )

def f_q(q, k, delta, mu):
    
    E_plus = energy(k + q/2, delta, mu)
    E_minus = energy(k - q/2, delta, mu)

    ksi_plus = np.cos(k + q/2) - mu/2
    ksi_minus = np.cos(k - q/2) - mu/2

    delta_plus = delta * np.sin(k + q/2)
    delta_minus = delta * np.sin(k - q/2)
    
    if q != 0.0:
        return 0.5 * ( E_plus*E_minus - ksi_plus*ksi_minus + delta_plus*delta_minus )/( E_plus*E_minus * ( E_plus + E_minus ) )
    else:
        return 0.5 * ( (delta_plus ** 2)/(E_plus ** 3) )

def chi(q, delta, mu, N_modes):
    k = np.linspace(0, 2 * np.pi, N_modes)
    return np.sum(f_q(q, k, delta, mu)) * 1/N_modes

def chi_RPA(chi, lambd):
    return chi/(1-2*lambd*chi)

def electron_density(delta, mu, N_modes):
    k = np.linspace(0, 2 * np.pi, N_modes)
    return 0.5 + (1/(2*N_modes)) * np.sum( (mu/2-np.cos(k))/(energy(k, delta, mu)) )

def maximum_chi_q(delta, mu, N_modes):
    q = np.linspace(0, 2 * np.pi, N_modes)
    chi_q = np.array([chi(q_i, delta, mu, N_modes) for q_i in q])
    return np.max(chi_q)

def converge_to_solution(delta_current, delta_previous, number_of_conv, mu, N_modes):
    j = 0
    delta_c = delta_current
    delta_p = delta_previous
    while j <= number_of_conv:

        j = j + 1
        delta_list = [delta_p + i * (delta_c - delta_p) / 10 for i in range(11)]
        #print(delta_list)
        sign_old = 1 - g_lambda * maximum_chi_q(delta_list[0], mu, N_modes) * 0.5

        for i in range(len(delta_list)):
            sign_new = 1 - g_lambda * maximum_chi_q(delta_list[i], mu, N_modes) * 0.5

            if sign_new * sign_old <= 0:
                delta_c = delta_list[i]
                delta_p = delta_list[i-1]
                #print(f'{delta_c} {delta_p}')
                break

            sign_old = sign_new
    return ( delta_c + delta_p ) * 0.5


def find_instability_point(mu):

    # This Algorithm finds the solutions and saves them in a list
    print(f'Processing {mu}...')
    # [(n1, delta1, mu1), (n2, delta2, mu2), ...]

    delta_step = 0.1
    delta, i = 0.0, 0
    solutions = []

    # iterate through delta
    while delta < 2.0:

        # Calculate the equation (here the sign is only important)
        sign_new = 1 - g_lambda * maximum_chi_q(delta, mu, g_N_modes) * 0.5

        # Assign the sign of the first iteration
        if i == 0:
            sign_old = sign_new
        
        # Checking if there is a instability point
        if sign_old*sign_new <= 0:
            
            # The solution is in between delta_current and delta_previous
            delta = converge_to_solution(delta, delta - delta_step, 2, mu, g_N_modes)

            # Append the solution
            solutions.append((electron_density(delta, mu, g_N_modes), delta, mu))
            #print(f'n = {electron_density(delta, mu, g_N_modes)}, delta = {delta}, mu = {mu}')


        # Upgrade delta and old sign
        sign_old = sign_new
        i = i + 1
        delta =  i * delta_step
    #print(solutions)
    return solutions

def save_solutions_in_file(lambd, results):
    with open(f"instability_line_lambda={lambd}.txt", "w") as f:
        f.write(f'#n delta mu\n')
        for solutions in results:
            # Extract the tuple inside the single-item list
            for solution in solutions:
                val1, val2, val3 = solution
                # Write formatted line with 3 columns: n delta mu
                f.write(f"{val1:15.8f} {val2:15.8f} {val3:15.8f}\n")

def construct_mu_values():
    mu, mu_step, mu_step_precision = 0.0, 0.1, 0.001
    i = 0
    mu_values = []
    while mu < 3.0:
        
        # Skip mu == 2 
        if mu == 2.0:
            continue
        
        # High steps away from 2
        if abs(mu - 1.99) >= 4*mu_step:
            mu_values.append(mu)
            mu = mu + mu_step
        # Precision steps for mu close to 2
        else:
            mu_values.append(mu)
            mu = mu + mu_step_precision
    return mu_values

def sorting_algorithm(results):
    array_mu_less_2 = []
    array_mu_bigger_2 = []
    for solutions in results:
        # Extract the tuple inside the single-item list
        for solution in solutions:
            val1, val2, val3 = solution
            if val3 <= 2.0:
                array_mu_less_2.append((val1, val2, val3))
            else:
                array_mu_bigger_2.append((val1, val2, val3))
    
    # Sort with increasing mu the array with mu < 2
    array_mu_less_2.sort(key=lambda x: x[2])
    # Sort with increasing electron density the array with mu > 2
    array_mu_bigger_2.sort(key=lambda x: x[0])

    return array_mu_less_2 + array_mu_bigger_2



