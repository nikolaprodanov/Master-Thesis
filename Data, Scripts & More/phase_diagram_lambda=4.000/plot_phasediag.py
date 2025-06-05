import pandas as pd
import matplotlib.pyplot as plt
import os

lambd = 4.0

# === Configuration: Set your file info here ===
file_config = {
    'maxwell_dimmer_homogeneous_1.txt': {'label': 'PS dimer', 'color': 'green'},
    'maxwell_dimmer_homogeneous_2.txt': {'label': 'PS topo', 'color': 'red'},
    'maxwell_nontopo_topo_1.txt': {'color': 'red'},  # No label
    'maxwell_nontopo_topo_2.txt': {'label': 'PS non-topo', 'color': 'blue'},
    'spinodal_dimmer_homogeneous_1.txt': {'label': 'Spinodal', 'color': 'grey'},
    'spinodal_dimmer_homogeneous_2.txt': {'color': 'grey'},  # No label
    'spinodal_nontopo_topo_1.txt': {'color': 'grey'},  # No label
    'spinodal_nontopo_topo_2.txt': {'color': 'grey'},  # No label
    'instability_line_lambda=4.0.txt': {'label':'Instability Line', 'color':'black'}
}

# === Plotting ===
#plt.figure(figsize=(10, 6))

for filename, props in file_config.items():
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        continue

    try:
        # Read file (assuming whitespace-separated, no header, 2 columns)
        data = pd.read_csv(filename, delim_whitespace=True, header=None)
        x, y = data[0], data[1]

        # Get label and color if provided
        label = props.get('label', None)
        color = props.get('color', None)

        # Plot with optional label and color
        if label == 'Spinodal' or label == None or label == 'Instability Line':
            if label != 'Instability Line':
                plt.plot(x, y, label=label, color=color, linewidth=2)
            else:
                plt.plot(x, y, linestyle= 'None', marker= '.', label=label, color=color)
        else:
            plt.plot(x, y, label=label, color=color)
    except Exception as e:
        print(f"Error reading '{filename}': {e}")

# Plot settings
plt.title(f"Phase Diagram for "+r'$\lambda =$' + f" {lambd}")
plt.xlabel("Electron Density")
plt.ylabel(r"$\Delta [t]$")
plt.xlim(0.5,1)
plt.ylim(0,4)
plt.legend()
plt.tight_layout()
plt.savefig(f'phase_diagram_lambda={lambd}.png')

# Show plot
plt.show()

