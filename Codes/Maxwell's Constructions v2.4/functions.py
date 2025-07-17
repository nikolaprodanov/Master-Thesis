# ****************************************************
#           Definitions of Functions
# ****************************************************
#   t = 1

from constants import lambd, x_min, x_max
import math
import numpy as np
from scipy.integrate import simpson

def my_sign(x):
    print(x)
    if x == 0:
        return 0
    else:
        if x > 0:
            return 1
        else:
            return -1

def energy(k, nu, delta, mu_tilde, x):
    '''
    Calculates E_k^nu spectrum of the Kitaev Model for different cases:
    case 1: delta = 0, mu_tilde = 0  Half - Filling Specific Case
    case 3: delta = 0, mu_tilde > 0  
    case 2: delta > 0, mu_tilde = 0  Half - Filling General Case
    case 4: delta > 0, mu_tilde > 0
    Variables:
        k is the mode
        nu takes values + 1 or -1
    '''

    ksi_k = np.cos(k)

    if delta == 0:
        if mu_tilde == 0:
            return np.sqrt( ksi_k**2 + (x/2)**2 ) # it is independent of nu but this does not create an issue in the sums
        else:
            return abs( np.sqrt( ksi_k**2 + (x/2)**2 ) + nu * (mu_tilde/2) )
    if delta != 0:

        delta_k = delta * np.sin(k)

        if mu_tilde == 0:
            return np.sqrt( ksi_k**2 + (delta_k + nu * (x/2))**2 )
        else:
            return np.sqrt( ksi_k**2 + delta_k**2 + (mu_tilde/2)**2 + (x/2)**2 + nu * np.sqrt( (ksi_k*mu_tilde)**2 + (x*delta_k)**2 + (x*mu_tilde/2)**2 ) )
  
def density(delta, mu_tilde, x, N_modes):
    '''
    Calculates electron density in the Kitaev model as the ratio of number of electrons and number of sites
    case 1: mu_tilde = 0      Half - Filling Case
    case 2: delta = 0, mu_tilde > 0  
    case 3: delta > 0, mu_tilde > 0
    Variables:
        k is the mode
        nu takes values + 1 or -1
    '''
    if mu_tilde == 0: # add also mu_tilde = 2
        return 0.5 
    else:
        # define the reciprocal space
        k = np.linspace(0, 2 * np.pi * (N_modes/2 - 1) / N_modes, N_modes//2)

        if delta == 0:
            sum_nu_pos = (+1) * np.sum( my_sign( np.sqrt(np.cos(k)**2 + (x/2)**2) + (+1)*mu_tilde/2 ) ) 
            sum_nu_neg = (-1) * np.sum( my_sign( np.sqrt(np.cos(k)**2 + (x/2)**2) + (-1)*mu_tilde/2 ) )
            return 0.5 + (sum_nu_pos + sum_nu_neg)/(2*N_modes)
        else:
            sum_nu_pos = np.sum( ( 0.5 + ((np.cos(k))**2+(x/2)**2)/np.sqrt((np.cos(k)*mu_tilde)**2+(x*delta * np.sin(k))**2+(x*mu_tilde/2)**2))/energy(k, +1, delta, mu_tilde,x) )
            sum_nu_neg = np.sum( ( 0.5 - ((np.cos(k))**2+(x/2)**2)/np.sqrt((np.cos(k)*mu_tilde)**2+(x*delta * np.sin(k))**2+(x*mu_tilde/2)**2))/energy(k, -1, delta, mu_tilde,x) )
            return 0.5 + (mu_tilde/2) * (sum_nu_pos + sum_nu_neg)/N_modes

def Free_Helmholtz_Energy(delta, mu_tilde, x, N_modes):
    '''
    Calculates the Free Helmholtz Energy, here it is defined as F/N_sites
    Variables:
        k is the mode
        nu takes values + 1 or -1
    '''

    # define the reciprocal space
    k = np.linspace(0, 2 * np.pi * (N_modes/2 - 1) / N_modes, N_modes//2)

    # calculate the electron density
    n = density(delta, mu_tilde, x, N_modes)

    # calculate the grand canonical potential
    omega = -np.sum( energy(k, +1, delta, mu_tilde, x) + energy(k, -1, delta, mu_tilde, x))/N_modes + 0.5* (x**2)/lambd + 0.5* lambd * n**2 - 0.5*mu_tilde
    
    # calculate the free energy
    return omega + mu_tilde * n - lambd*n**2

def derivative_Free_Helmholtz_Energy(delta, mu_tilde, x, N_modes): # this one is actually wrong !!!!
    '''
    Calculates the Free Helmholtz Energy, here it is defined as F/N_sites
    case 1: delta = 0, mu_tilde = 0  Half - Filling Specific Case
    case 3: delta = 0, mu_tilde > 0  
    case 2: delta > 0, mu_tilde = 0  Half - Filling General Case
    case 4: delta > 0, mu_tilde > 0
    Variables:
        k is the mode
        nu takes values + 1 or -1
    '''

    # define the reciprocal space
    k = np.linspace(0, 2 * np.pi * (N_modes/2 - 1) / N_modes, N_modes//2)

    if delta == 0:
        if mu_tilde == 0:
            return x/lambd - x/(2*N_modes)*np.sum(1/energy(k,+1, delta, mu_tilde, x)) # there is no sum over positive and negative nu, it has been taken into account!
        else: 
            sum_nu_pos = np.sum(my_sign(np.sqrt((np.cos(k))**2+(x/2)**2) + (mu/2)) * (x/2)/(np.sqrt((np.cos(k))**2+(x/2)**2))) 
            sum_nu_neg = np.sum(my_sign(np.sqrt((np.cos(k))**2+(x/2)**2) - (mu/2)) * (x/2)/(np.sqrt((np.cos(k))**2+(x/2)**2))) 
            return x/lambd - x/(2*N_modes) * (sum_nu_pos + sum_nu_neg)
    else:
        if mu_tilde == 0:
            sum_nu_pos = np.sum((x+delta*np.sin(k))/energy(k,+1,delta,mu_tilde,x))
            sum_nu_neg = np.sum((x-delta*np.sin(k))/energy(k,-1,delta,mu_tilde,x))
            return x/lambd - 1/(2*N_modes) * (sum_nu_pos + sum_nu_neg)
        else:
            sum_nu_pos = np.sum(((x/2 + x*((np.cos(k))**2+(x/2)**2))/np.sqrt((np.cos(k)*mu_tilde)**2+(x*delta * np.sin(k))**2+(x*mu_tilde/2)**2))/energy(k, +1, delta, mu_tilde,x))
            sum_nu_neg = np.sum(((x/2 - x*((np.cos(k))**2+(x/2)**2))/np.sqrt((np.cos(k)*mu_tilde)**2+(x*delta * np.sin(k))**2+(x*mu_tilde/2)**2))/energy(k, -1, delta, mu_tilde,x))
            return x/lambd - ( sum_nu_pos + sum_nu_neg )/(2*N_modes)

def determine_phase(x, mu_tilde):
    if x == 0:
        if mu_tilde < 2.0:
            return 'topo'
        else:
            return 'nontopo'
    else:
        return 'dimer'

def Maxwell_Area_Calculator(n_solutions, mu_solutions, n1, n2, mu_star):
    """
    Computes the integral of (mu(n) - mu_star) dn from n1 to n2
    using the trapezoidal rule.
    Returns:
        float: The computed integral value.
    """
    # Sanity check
    if len(n_solutions) != len(mu_solutions):
        raise ValueError("n_solutions and mu_solutions must have the same length")

    # Sort data by n
    sorted_data = sorted(zip(n_solutions, mu_solutions))
    n_sorted, mu_sorted = zip(*sorted_data)
    n_arr = np.array(n_sorted)
    mu_arr = np.array(mu_sorted)

    # Select data in the range [n1, n2]
    mask = (n_arr >= n1) & (n_arr <= n2)
    n_segment = n_arr[mask]
    mu_segment = mu_arr[mask]

    # If range is empty or has fewer than 2 points, return 0
    if len(n_segment) < 2:
        return 0.0

    # Shift mu by mu_star
    f_segment = mu_segment - mu_star

    # Integrate
    integral = np.trapz(f_segment, x=n_segment)
    return integral

def assing_delta_list(delta_min, delta_max, delta_step):
    array = []
    delta = delta_min
    j = 1
    while delta <= delta_max:
        array.append(delta)
        delta = delta_min + delta_step * j
        j = j + 1
    return array