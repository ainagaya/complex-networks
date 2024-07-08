import matplotlib.pyplot as plt
import argparse

# Setting the global style for plots
plt.style.use('seaborn-whitegrid')

# Customizing the plot appearance
def set_plot_style():
    plt.rcParams['figure.figsize'] = (10, 6)  # Figure size
    plt.rcParams['axes.titlesize'] = 18       # Title font size
    plt.rcParams['axes.labelsize'] = 14       # Axes labels font size
    plt.rcParams['axes.grid'] = True          # Show grid by default
    plt.rcParams['grid.alpha'] = 0.3          # Grid transparency
    plt.rcParams['grid.linestyle'] = '--'     # Grid line style
    plt.rcParams['xtick.labelsize'] = 12      # X-tick labels font size
    plt.rcParams['ytick.labelsize'] = 12      # Y-tick labels font size
    plt.rcParams['legend.fontsize'] = 12      # Legend font size
    plt.rcParams['legend.frameon'] = True     # Show legend frame
    plt.rcParams['legend.loc'] = 'best'       # Legend location
    plt.rcParams['scatter.marker'] = 'o'      # Scatter plot marker
    plt.rcParams['lines.markersize'] = 8      # Marker size
    plt.rcParams['font.family'] = 'serif'     # Font family
    plt.rcParams['font.serif'] = ['Times New Roman']  # Font type

set_plot_style()  

parser = argparse.ArgumentParser()
parser.add_argument('--network_file', type=str, help='Path to the edge list file')
args = parser.parse_args()

knn_list = []
degrees = []

# Read the calculated degree distribution from the file
with open(f'output/knn_{args.network_file}.dat', 'r') as f:
    lines = f.readlines()
    for line in lines:
        i, knn = line.strip().split()
        if float(knn) != 0:
            degrees.append(float(i))
            knn_list.append(float(knn))


# Plot the degree distribution
plt.plot(degrees, knn_list, '.-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$k$')
plt.ylabel('$K_{nn}/\kappa$')
plt.title('$K_{nn}$ Distribution')
plt.savefig(f'figures/knn_distribution_{args.network_file}.png')

# Plot the clustering coefficient distribution

cc_list = []
degrees = []

plt.clf()

with open(f'output/cc_{args.network_file}.dat', 'r') as f:
    lines = f.readlines()
    for line in lines:
        i, cc = line.strip().split()
        # zero values are not considered since in the log-log
        # plot they would be -inf
        if float(cc) > 10e-20:
            degrees.append(float(i))
            cc_list.append(float(cc))

plt.plot(degrees, cc_list, '.-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$k$')
plt.ylabel('$\\bar{c}(k)$')
#plt.title('cc Distribution')
plt.savefig(f'figures/cc_distribution_{args.network_file}.png')

