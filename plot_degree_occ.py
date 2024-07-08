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

degree_occurrences = []
acc_degree_occurrences = []
degrees = []

set_plot_style()

parser = argparse.ArgumentParser()
parser.add_argument('--network_file', type=str, help='Path to the edge list file')
args = parser.parse_args()

# Read the calculated degree distribution from the file
with open(f'output/pk_{args.network_file}.dat', 'r') as f:
    lines = f.readlines()
    for line in lines:
        i, degree_occurrence, acc_degree_occurrence = line.strip().split()
        degrees.append(float(i))
        degree_occurrences.append(float(degree_occurrence))
        acc_degree_occurrences.append(float(acc_degree_occurrence))

#print(degrees)
#print(degree_occurrences)

# Plot the degree distribution
plt.plot(degrees, degree_occurrences, '-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('Occurrences')
plt.title('Degree Distribution')
plt.savefig(f'figures/degree_distribution_{args.network_file}.png')

plt.clf()
# Plot the degree accumulated distribution
plt.plot(degrees, acc_degree_occurrences, '-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('Accumulated degree occurrences')
plt.title('Accumulated Degree Distribution')
plt.savefig(f'figures/acc_degree_distribution_{args.network_file}.png')

plt.clf()
# Plot the complemetntary degree accumulated distribution
comp_acc_degree_occurrences = [1 - x for x in acc_degree_occurrences]
plt.plot(degrees, comp_acc_degree_occurrences, '-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('Accumulated degree occurrences')
plt.title('Complementary Accumulated Degree Distribution')
plt.savefig(f'figures/comp_acc_degree_distribution_{args.network_file}.png')