import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_first_two_columns_skip_comments():
    current_folder = os.getcwd()  # Folder where the script is located
    fig, ax = plt.subplots(figsize=(10, 6))

    # Store all data in a list to reuse for inset
    data_list = []

    for filename in os.listdir(current_folder):
        if filename.endswith(".dat") or filename.endswith(".txt"):
            filepath = os.path.join(current_folder, filename)

            try:
                # Try reading with whitespace delimiter and skip comment lines starting with #
                data = pd.read_csv(filepath, delim_whitespace=True, header=None, comment='#')

                # If not enough columns, try comma separator with comment skipping
                if data.shape[1] < 2:
                    data = pd.read_csv(filepath, sep=',', header=None, comment='#')

                # Take only the first two columns
                data = data.iloc[:, :2]

                x = data.iloc[:, 0]
                y = data.iloc[:, 1]

                ax.plot(x, y, label=filename)
                data_list.append((x, y))
            except Exception as e:
                print(f"Skipping {filename}: {e}")

    ax.set_xlabel("n")
    ax.set_ylabel("Delta")
    ax.set_title("Spinodals of Dimer-Homogeneous for lambda = 4.0")
    ax.set_xlim(0.5, 1)  # Main plot x-axis limit
    ax.set_ylim(0)
    ax.legend()

    # Create inset axes [x0, y0, width, height] in figure coordinates
    inset_ax = fig.add_axes([0.1, 0.4, 0.35, 0.35])

    # Plot the same data in the inset zoomed between 0.5 and 0.52
    for (x, y) in data_list:
        inset_ax.plot(x, y)

    inset_ax.set_xlim(0.5, 0.505)
    # Set y limits dynamically based on data in the zoomed range 0.5 to 0.52
    y_min = min(
        y[(x >= 0.5) & (x <= 0.52)].min()
        for x, y in data_list
        if ((x >= 0.5) & (x <= 0.52)).any()
    )
    y_max = max(
        y[(x >= 0.5) & (x <= 0.52)].max()
        for x, y in data_list
        if ((x >= 0.5) & (x <= 0.52)).any()
    )
    inset_ax.set_ylim(0, 0.6)

    inset_ax.set_title("Zoom: 0.5 to 0.52")
    inset_ax.tick_params(axis='both', which='major', labelsize=8)

    plt.tight_layout()
    plt.savefig('results_lambda=4.0.png')
    plt.show()
    plt.close()

# Run it
plot_first_two_columns_skip_comments()

