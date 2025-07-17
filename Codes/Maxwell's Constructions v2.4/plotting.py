
from constants import lambd, x_min, x_max, N_samples, N_modes
import functions as fun 
import algorithms

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes,mark_inset

def plot_spectrum(delta,mu_tilde, x, N_modes):
    # Plot the Free Energy vs x
    k = np.linspace(0, 2 * np.pi * (N_modes/2 - 1) / N_modes, N_modes//2)
    spectrum_1 = []
    spectrum_2 = []
    spectrum_3 = []
    spectrum_4 = []
    for i in range(len(k)):
        spectrum_1.append(fun.energy(k[i],+1,delta,mu_tilde,x))
        spectrum_2.append(fun.energy(k[i],-1,delta,mu_tilde,x))
        spectrum_3.append(-spectrum_1[i])
        spectrum_4.append(-spectrum_2[i])

    # plot
    plt.plot(k, spectrum_1, marker='x', linestyle = None, color = 'b')
    plt.plot(k, spectrum_2, marker='x', linestyle = None, color = 'b')
    plt.plot(k, spectrum_3, marker='x', linestyle = None, color = 'b')
    plt.plot(k, spectrum_4, marker='x', linestyle = None, color = 'b')

    plt.xlabel('k')
    plt.ylabel('Spectrum E_k')
    plt.title(r'$\lambda=$'+f'{lambd}, '+r'$\Delta=$'+f'{delta}, '+r'$\tilde{\mu}=$'+f'{mu_tilde}, '+r'$x=$'+f'{x}')
    plt.show()
    plt.close()
    return

def plot_Free_Energy_vs_x(delta,mu_tilde,N_modes): 

    # Plot the Free Energy vs x
    x = algorithms.log_decade_sampling(x_min*lambd,x_max*lambd)
    free_energy = []
    empty_list = []
    F_of_zero = fun.Free_Helmholtz_Energy(delta,mu_tilde, 0.0, N_modes)
    for i in range(len(x)):
        free_energy.append(fun.Free_Helmholtz_Energy(delta,mu_tilde,x[i], N_modes))
        empty_list.append(F_of_zero)
        
    # Plot settings
    plt.xlabel(r'$x=\frac{gX_0}{t}$')
    plt.ylabel('Free Energy')
    plt.title('Free Energy(x; '+r'$\lambda=$'+f'{lambd}, '+r'$\Delta=$'+f'{delta}, '+r'$\tilde{\mu}=$'+f'{mu_tilde})')

    # plot
    plt.plot(x, free_energy, marker = None, linestyle = '-', label=r'$F(x)$', color = 'black')
    plt.plot(x, empty_list, linestyle = '--', marker=None, label=r'$F(0) = $'+f'{F_of_zero}',color = 'red')

    #plt.yscale('symlog', linthresh=1e-5)  # near 0 behaves linear

    # Add other points
    # Global minima:
    min_value = min(free_energy)
    min_index = free_energy.index(min_value)
    plt.plot(x[min_index], min_value, marker = 'x', linestyle = None, label =r'F(x(global min)='+f'{x[min_index]}'+r'$) = $'+f'{min_value}', color = 'black')
    # Local minimum
    x_local_minimum = algorithms.find_local_minimum_of_F(x, delta, mu_tilde, N_modes)
    if x_local_minimum != None and x_local_minimum != x[min_index]:
        plt.plot(x_local_minimum, fun.Free_Helmholtz_Energy(delta, mu_tilde, x_local_minimum, N_modes), marker = 'x', linestyle = None,label = f'x_local_minimum = {x_local_minimum}', color = 'black')
    if x_local_minimum == None:
        plt.plot([],[], marker = 'x', linestyle = None ,label = 'x(local min) = None', color = 'black')
    if x_local_minimum == x[min_index]:
        plt.plot([],[], marker = 'x',linestyle = None,label = 'x(local min) = x(global min)', color = 'black')

    # Local maximum
    x_local_maximum = algorithms.find_local_maximum_of_F(x, delta, mu_tilde, N_modes)
    if x_local_maximum != None:
        plt.plot(x_local_maximum, fun.Free_Helmholtz_Energy(delta, mu_tilde, x_local_maximum, N_modes), marker = 'x', linestyle = None,label = f'x(local max) = {x_local_maximum}', color = 'black')
    if x_local_maximum == None and x[min_index] != 0:
        plt.plot(0.0, fun.Free_Helmholtz_Energy(delta, mu_tilde, 0.0, N_modes), marker = 'x', linestyle = None,label = f'x(local max) = 0.0', color = 'black')
        x_local_maximum = 0.0
    if x_local_maximum == None and x[min_index] == 0:
        plt.plot([] , [], marker = 'x', linestyle = None, label = 'x(local max) = None', color = 'black')

    # Adjust the scale of the plot
    plt.xscale('symlog', linthresh=x[1])
    plt.ylim(min_value-abs(min_value)*0.5, min_value+abs(min_value)*0.5)

    # Get the current axes (the main plot)
    main_ax = plt.gca()

    # Create an inset manually in figure coordinates [left, bottom, width, height]
    if x_local_minimum != None and x_local_maximum != None:
        inset_ax = inset_axes(main_ax, width="30%", height="30%", loc='upper right')

        f_maximum = fun.Free_Helmholtz_Energy(delta, mu_tilde, x_local_maximum, N_modes)
        f_minimum = fun.Free_Helmholtz_Energy(delta, mu_tilde, x_local_minimum, N_modes)

        inset_ax.plot(x, free_energy, marker = None, linestyle = '-', color = 'black')
        inset_ax.plot(x, empty_list, linestyle = '--', marker=None, color = 'red')
        inset_ax.plot(x_local_minimum, f_minimum, marker = 'x', color = 'black')
        inset_ax.plot([],[])
        inset_ax.plot(x_local_maximum, f_maximum, N_modes, marker = 'x', color = 'black')

        if x_local_maximum == 0.0:
            x_local_maximum = -x[1]
        
        inset_ax.set_xlim(x_local_maximum - abs(x_local_maximum)*0.5, x_local_minimum + abs(x_local_minimum)*0.5)
        inset_ax.set_ylim(f_minimum - abs(f_minimum)*0.01, f_maximum + abs(f_maximum)*0.01)

        inset_ax.set_xticks([])
        inset_ax.set_yticks([])

    main_ax.legend()
    plt.savefig(f'Free_Energy_vs_x_for_lambda={lambd}_delta={delta}_mu_tilde={mu_tilde}.png')
    plt.show()
    plt.close()
    return

def plot_derivative_Free_Energy_vs_x(delta,mu_tilde,N_modes):

    x = algorithms.log_decade_sampling(x_min,x_max)

    free_energy = []
    empty_list = []
    for i in range(len(x)):
        free_energy.append(fun.derivative_Free_Helmholtz_Energy(delta,mu_tilde,x[i], N_modes))
        empty_list.append(0.0)

    # plot
    plt.plot(x, free_energy, marker='x', linestyle = None)
    plt.plot(x, empty_list, linestyle='--')
    plt.xlabel('x')
    plt.ylabel('Derivative Free Energy')
    plt.title(r'$\lambda=$'+f'{lambd}, '+r'$\Delta=$'+f'{delta}, '+r'$\tilde{\mu}=$'+f'{mu_tilde} ')
    plt.xscale('log')
    plt.yscale('symlog', linthresh=1e-2)  # near 0 behaves linear
    plt.show()
    return

def plot_topology_change_line():

    delta = np.linspace(0,2.0, N_samples)
    densities = []
    for d in delta:
        densities.append(fun.density(d,2.0,0,N_modes))

    plt.plot(densities, delta)
    plt.show()
    return

def make_phase_diagram(DATA):

    print('Producing The Phase Diagram')

    # See the data
    #for i in range(0,len(DATA)):
    #    print(DATA[i])

    # Sorting Data and Preparing for plots-------------------------------

    # Arrays for the dimer-to-topo phase separation
    delta_list_1, n_ps_dimer_topo_1, n_ps_dimer_topo_2, n_extrema_dimer_topo_1, n_extrema_dimer_topo_2 = [], [], [], [], []

    # Arrays for the topo-to-nontopo phase separation
    delta_list_2, n_ps_topo_nontopo_1, n_ps_topo_nontopo_2, n_extrema_topo_nontopo_1, n_extrema_topo_nontopo_2 = [], [], [], [], []

    for i in range(0, len(DATA)):

        for j in range(0, len(DATA[i])):

            # What is the current phase separation?
            ps_j = DATA[i][j][1]
            
            if ps_j == 'dimer-to-topo': # Assign Data to dimer-to-topo phase separation
                delta_list_1.append(DATA[i][j][0])
                n_ps_dimer_topo_1.append(DATA[i][j][2])
                n_ps_dimer_topo_2.append(DATA[i][j][3])
                n_extrema_dimer_topo_1.append(DATA[i][j][5])
                n_extrema_dimer_topo_2.append(DATA[i][j][7])
            else: # Assign Data to topo-to-nontopo phase separation
                delta_list_2.append(DATA[i][j][0])
                n_ps_topo_nontopo_1.append(DATA[i][j][2])
                n_ps_topo_nontopo_2.append(DATA[i][j][3])
                n_extrema_topo_nontopo_1.append(DATA[i][j][5])
                n_extrema_topo_nontopo_2.append(DATA[i][j][7])
    
    # End Sorting Data and Preparing for plots---------------------------
    
    # Make the plots phase diagrams
    if delta_list_1 != []:
        plt.plot(n_ps_dimer_topo_1, delta_list_1, marker='x', linestyle='-', color='green', label='PS Dimer')
        plt.plot(n_ps_dimer_topo_2, delta_list_1, marker='x', linestyle='-', color='red')
        #plt.plot(n_extrema_dimer_topo_1, delta_list_1, marker='x', linestyle='--', color='grey')
        #plt.plot(n_extrema_dimer_topo_2, delta_list_1, marker='x', linestyle='--', color='grey')

        with open('PS_dimer_topo.dat', 'w') as f:
            f.write('#delta n_PS_dimer n_PS_topo\n')
            for item1, item2, item3 in zip(delta_list_1, n_ps_dimer_topo_1, n_ps_dimer_topo_2):
                f.write(f'{item1} {item2} {item3}\n')

    if delta_list_2 != []:
        plt.plot(n_ps_topo_nontopo_1, delta_list_2, marker= None, linestyle='-', color='red', label='PS topo')
        plt.plot(n_ps_topo_nontopo_2, delta_list_2, marker= None, linestyle='-', color='blue', label='PS non-topo')
        plt.plot(n_extrema_topo_nontopo_1, delta_list_2, marker= None, linestyle='--', color='grey', label='extremas')
        plt.plot(n_extrema_topo_nontopo_2, delta_list_2, marker= None, linestyle='--', color='grey')

    
    plt.title(f'Phase Diagram lambda = {lambd}')
    plt.xlabel('n (electron density)')
    plt.ylabel('Delta')
    plt.legend()
    plt.xlim(0.5, 1.0)
    plt.ylim(0.0)
    #plt.ylim(0.0, 2.0)

    plt.savefig(f'Phase Diagram lambda = {lambd:.1f}.png')
    plt.show()
    plt.close()

    return

# How DATA[i][j] looks like:
# Store data
#        maxwell_points.append( (
#            delta, phase_separation,                                    # Information
#            n_phase_separation_1, n_phase_separation_2, mu_star,        # Phase Separation Points
#            n_extrema_1, mu_extrema_1, n_extrema_2, mu_extrema_2,   # extrema points
#            n_phase_change#,                                             # For plotting colors
#           n_solutions, mu_solutions                                   # for plotting mu vs n
# ))

def plot_color(phase):
    if phase == 'dimer':
        return 'green'
    if phase == 'topo':
        return 'red'
    if phase == 'nontopo':
        return 'blue'
    return

def mu_vs_n_with_Maxwells_Constructions(delta, N_modes):
    # Get data
    data_mu_vs_n = algorithms.Make_Maxwells_Constructions(delta, N_modes)

    fig, ax = plt.subplots()
    dimer_inset_drawn = False
    all_phase_lines = []

    for i in range(len(data_mu_vs_n)):
        delta = data_mu_vs_n[i][0]
        phase_separation = data_mu_vs_n[i][1]
        n_phase_change = data_mu_vs_n[i][9]
        mu_solutions, n_solutions = data_mu_vs_n[i][11], data_mu_vs_n[i][10]

        if i == 0:
            n_lower_limit = 0.5

        filtered_n = [ni for ni in n_solutions if n_lower_limit <= ni < n_phase_change]
        filtered_mu = [mu_solutions[i] for i, ni in enumerate(n_solutions) if n_lower_limit <= ni < n_phase_change]

        phase_label = phase_separation.split('-to-')[0]
        phase_color = plot_color(phase_label)

        ax.plot(filtered_n, filtered_mu, linestyle='-', marker='x', color=phase_color, label=phase_label)
        all_phase_lines.append((filtered_n, filtered_mu, phase_color))

        # Check for dimer phase and create inset
        if 'dimer' in phase_separation and not dimer_inset_drawn:
            n_phase_separation_2 = data_mu_vs_n[i][3]
            mu_star_dimer = data_mu_vs_n[i][4]

            # Get μ values in dimer range
            dimer_mu_values = [mu_solutions[j] for j, nj in enumerate(n_solutions)
                               if 0.5 <= nj <= n_phase_separation_2]
            mu_max_dimer = max(dimer_mu_values)
            mu_min_dimer = min(dimer_mu_values)
            mu_margin = 0.05 * (mu_max_dimer - mu_min_dimer)
            y_min_inset = mu_min_dimer - mu_margin
            y_max_inset = mu_max_dimer + mu_margin

            # Create inset axis, placed slightly lower
            axins = inset_axes(ax, width="40%", height="30%", 
                   #bbox_to_anchor=(0.55, 0.35, 0.4, 0.3),  # x0, y0, width, height
                   bbox_transform=ax.transAxes,
                   loc='upper center')
            dimer_x_range = (0.5, (0.5+n_phase_separation_2)/2)
            dimer_y_range = (y_min_inset, y_max_inset)

            dimer_inset_drawn = True

        n_lower_limit = n_phase_change

    # Plot last phase
    filtered_n = [ni for ni in n_solutions if n_lower_limit <= ni]
    filtered_mu = [mu_solutions[i] for i, ni in enumerate(n_solutions) if n_lower_limit <= ni]
    phase_label = phase_separation.split('-to-')[1]
    phase_color = plot_color(phase_label)
    ax.plot(filtered_n, filtered_mu, linestyle='-', marker='x', color=phase_color, label=phase_label)
    all_phase_lines.append((filtered_n, filtered_mu, phase_color))

    # Plot Maxwell constructions (main and inset)
    for i in range(len(data_mu_vs_n)):
        n1, n2, mu_star = data_mu_vs_n[i][2], data_mu_vs_n[i][3], data_mu_vs_n[i][4]
        ax.hlines(y=mu_star, xmin=n1, xmax=n2, color='black', linestyle='-')
        ax.plot([n1, n2], [mu_star, mu_star], 'ko')

        if dimer_inset_drawn:
            axins.hlines(y=mu_star, xmin=n1, xmax=n2, color='black', linestyle='-')
            axins.plot([n1, n2], [mu_star, mu_star], 'ko')

    # Plot inset content
    if dimer_inset_drawn:
        for n_vals, mu_vals, color in all_phase_lines:
            axins.plot(n_vals, mu_vals, linestyle='-', marker='x', color=color)

        axins.set_xlim(*dimer_x_range)
        axins.set_ylim(*dimer_y_range)
        #axins.set_xscale('log')
        #axins.set_title("Zoom on Dimer Phase", fontsize=8)

    # Final main plot settings
    ax.set_xlabel('electron density (n)')
    ax.set_ylabel(r'$\mu = \tilde{\mu} - \lambda \cdot n$')
    ax.set_xlim(0.5, 1)
    ax.set_title(r'$mu$'+' VS '+ r'$n$'+' for '+r'$(\lambda = $'+f'{lambd}'r', $\Delta = $'+f'{delta})')
    ax.legend()
    plt.show()
    plt.close()

    return
