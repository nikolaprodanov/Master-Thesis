# ****************************************************
#    Figures for the instability_line.py
# ****************************************************

import functions
import matplotlib.pyplot as plt
import numpy as np
from constants import lambd, N_modes

def plot_susceptibility(logical, delta, mu, N_modes):
    
    if logical == False:
        print(f'Making Susceptibility Figure: {logical}')
        return
    else:
        print(f'Making Susceptibility Figure: {logical}')

    q = np.linspace(0, 2 * np.pi, N_modes)
    chi_q = np.array([functions.chi(q_i, delta, mu, N_modes) for q_i in q]) 

    plt.plot(q, chi_q, marker = None, linestyle = '-', label = r'$\Delta =$ '+f'{delta}, ' + r'$\mu =$ '+f'{mu}')
    #plt.xlim(0, 2 * np.pi)
    plt.ylim(0)
    plt.xlabel('q (rad)')
    plt.ylabel(r'$\chi(q)$')
    plt.title('Susceptibility without RPA')
    plt.legend()
    plt.savefig(f'susceptibility_Delta={delta}_mu={mu}.png')
    plt.close()

def plot_instability_delta_VS_lambda(logical, N_modes, mu_values):

    if logical == False:
        print(f'Making Diagram Delta VS lambda: {logical}')
        return
    else:
        print(f'Making Diagram Delta VS lambda: {logical}')

    delta_step = 0.01

    for mu in mu_values:
        
        array_delta = np.array([])
        array_lambd = np.array([])
        delta = delta_step
        i = 1

        print(f'Calculating for '+'mu =' +f' {mu}')

        while delta <= 1.5:

            # Find the Maximum Value of the susceptibility
            q = np.linspace(0, 2 * np.pi, N_modes)
            chi_q = np.array([functions.chi(q_i, delta, mu, N_modes) for q_i in q])

            max_chi = np.max(chi_q)

            # Calculate Lambda
            lambd = 2/max_chi # The real relation is 1 - lambd * chi * 0.5 = 0

            # Store the data
            array_delta = np.append(array_delta, delta)
            array_lambd = np.append(array_lambd, lambd)

            # Set the new delta
            delta = i * delta_step
            i = i + 1

        #plot that line
        plt.plot(array_lambd, array_delta, marker = None, linestyle = '-', label = r'$\mu = $'+f'{mu}')    

    # Plot the lines where lambda is constant
    plt.axvline(x=2, color='r', linestyle='--', label=r'$\lambda = 2.0$')
    plt.axvline(x=4, color='b', linestyle='--', label=r'$\lambda = 4.0$')

    plt.title('Instability Line '+r'$\Delta$'+' VS '+r'$\lambda$')
    plt.xlabel(r'$\lambda(\chi(\Delta, \mu))$')
    plt.ylabel(r'$\Delta$')
    #plt.xscale('log')
    plt.legend()
    plt.savefig('Diagram_Delta_VS_lambda.png')
    plt.close()

def plot_instability_line_phase_diagram(logical, lambd, N_modes):


    if logical == False:
        print(f'Making Instability Line: {logical}')
        return
    else:
        print(f'Making Instability Line: {logical}')
    
    # **************************************************************
    # Algorithm for finding the transition line PART 1: mu < 2
    # **************************************************************
    # start from delta zero and increase
    mu_step, delta_step = 0.001, 0.001
    i,j = 1,0
    mu = 0.0

    array_delta = np.array([])
    array_n = np.array([])
    q = np.linspace(0, 2 * np.pi, N_modes)

    while mu < 2.0:

        print(f'Status(PART1): mu = {mu}')
        instability_point_found = False
        delta = delta_step

        while instability_point_found == False:
            # Assign Delta each iteration
            delta = i * delta_step

            # Calculate the susceptibility and its maximum value
            chi_q = np.array([functions.chi(q_i, delta, mu, N_modes) for q_i in q])
            max_chi = np.max(chi_q)

            # Check if the Susceptibility is bellow the divergence line (1-lambda/2 * chi = 0)
            if max_chi <= 2/lambd:
                instability_point_found = True

                # Storing Data
                array_delta = np.append(array_delta, delta)
                array_n = np.append(array_n, functions.electron_density(delta, mu, N_modes))
                print(f'density = {functions.electron_density(delta, mu, N_modes)} and delta = {delta}')
            else:
                # If the instability point is not found, continue searching
                i = i + 1
        
        # Reset Delta for next iteration
        i = 1

        # Increase mu
        j = j + 1
        mu =  j * mu_step
    
    # **************************************************************
    # Algorithm for finding the transition line PART 2: mu > 2
    # **************************************************************
    # NOTICE: for lambda = 4 the program did not find convergences for mu>2
    # **************************************************************
    # When mu > 2 the max susceptibility is at q = 0 or q = 2pi and is always smaller than 2/lambda.
    # STRATEGY: here small delta will give a small susceptibility --> we increase delta untill the max susceptibility > 2/lambda

    mu = 2.0 + mu_step
    i,j = 1,1
    max_delta_iterations = 300

    while mu < 2.0:

        print(f'Status(PART 2): mu = {mu}')
        instability_point_found = False
        delta = delta_step

        while instability_point_found == False:
            # Assign Delta each iteration
            delta = i * delta_step
            #print(f'Current delta search: {delta}')

            # Calculate the susceptibility and its maximum value
            max_chi = functions.chi(0.0, delta, mu, N_modes)

            # Check if the Susceptibility is bellow the divergence line (1-lambda/2 * chi = 0)
            if max_chi >= 2/lambd:
                instability_point_found = True

                # Storing Data
                array_delta = np.append(array_delta, delta)
                array_n = np.append(array_n, functions.electron_density(delta, mu, N_modes))
                print(f'I found a {delta}')
            else:
                # If the instability point is not found, continue searching
                i = i + 1
            if i > max_delta_iterations:
                break
        
        # Reset Delta for next iteration
        i = 1

        # Increase mu
        j = j + 1
        mu = 2.0 + j * mu_step


    # Plot the instability line for the phase diagram
    plt.plot(array_n, array_delta, marker= 'x', linestyle ='-', label = r'$\lambda =$'+f' {lambd}')
    plt.title('Instability Line for '+r'$\lambda =$'+f' {lambd}')
    plt.xlabel('Electron Density')
    plt.ylabel(r'$\Delta$')
    plt.ylim(0,2)
    plt.xlim(0.5,1)
    plt.savefig(f'Instability_lambda={lambd}.png')
    plt.close()

    with open('instability_line.txt', 'w') as file:
        for a, b in zip(array_n, array_delta):
            file.write(f"{a}\t{b}\n")

def plot_superior_instability(results):

    n = []
    delta = []

    for solutions in results:
        # Extract the tuple inside the single-item list
        for solution in solutions:
            val1, val2, val3 = solution
            # Write formatted line with 3 columns
            n.append(val1)
            delta.append(val2)


    plt.scatter(n, delta, color='black', linestyle='-', marker=None)
    plt.title('Instability Line for '+r'$\lambda =$'+f' {lambd}')
    plt.xlabel('Electron Density')
    plt.ylabel(r'$\Delta$')
    plt.ylim(0,2)
    plt.xlim(0.5,1)
    plt.savefig(f'Instability_line_lambda={lambd}.png')
    plt.close()