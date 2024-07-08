# Script to read a network from an edge list file and visualize it
# Usage: python read_net.py --network_file <path_to_edge_list_file> --draw <circular or spring>
# It also calculates the average degree and the number of nodes in the network
# Author: Aina Gaya-Àvila

import networkx as nx
import argparse

import numpy as np

import matplotlib.pyplot as plt

from collections import defaultdict

import matplotlib as mpl

parser = argparse.ArgumentParser()
parser.add_argument('--network_file', type=str, help='Path to the edge list file')
parser.add_argument('--draw', type=str, help='Circular or spring')
args = parser.parse_args()

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

# Read the network from the edge list file
G = nx.read_edgelist(args.network_file)

# Visualize and save the network as an image
if args.draw == 'spring':
    pos = nx.spring_layout(G, seed=3068, iterations=3)
    nx.draw_networkx(G, pos=pos, with_labels=False, node_size=1, width=0.3)   # default spring_layout
elif args.draw == 'circular':
    nx.draw_circular(G, with_labels=False, node_size=1, width=0.1)
else:
    print('Invalid draw option. Use either spring or circular')
    exit()

plt.savefig(f'figures/network_{args.network_file}_{args.draw}.png')


# Get the list of degrees for each node
degrees = [degree for _, degree in G.degree()]
# Calculate the average degree
average_degree = sum(degrees) / G.number_of_nodes()
# Get the number of nodes
num_nodes = G.number_of_nodes()
# Print the results
print("Average degree:", average_degree)
print("Number of nodes:", num_nodes)

# Make a histogram of the degrees
plt.clf()
plt.hist(degrees, bins=200, color='green', edgecolor='black', histtype='stepfilled')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Degree k')
plt.ylabel('Frequency')
plt.title('Degree Distribution P(k)')
plt.savefig(f'figures/degree_distribution_{args.network_file}_networkx.png')


average_neighbor_degree = nx.average_neighbor_degree(G)

# Compute the average neighbor degree of each degree
degree_avg_neighbor_degree = np.zeros(num_nodes)
print(degree_avg_neighbor_degree)
for node in average_neighbor_degree.keys():
    degree = degrees[int(node) - 1]
    degree_avg_neighbor_degree[degree] =+ 1

# Plot the average neighbor degree
plt.clf()
plt.plot(degrees, degree_avg_neighbor_degree, '.')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Degree')
plt.ylabel('Average Neighbor Degree')
plt.title('Average Neighbor Degree vs Degree')
plt.savefig(f'figures/average_neighbor_degree_{args.network_file}_networkx.png')

# Calculate the average clustering coefficient
average_cluster_coefficient = nx.average_clustering(G)

print("Average cluster coefficient:", average_cluster_coefficient)

plt.clf()
# Plot and save Kcores
# build a dictionary of k-level with the list of nodes
kcores = defaultdict(list)
for n, k in nx.core_number(G).items():
    kcores[k].append(n)

print(kcores)
# compute position of each node with shell layout
pos = nx.layout.shell_layout(G, list(kcores.values()))

n_lines = 5
cmap = mpl.colormaps['rainbow']
# Take colors at regular intervals spanning the colormap.
colors = cmap(np.linspace(0, 1, n_lines))
colors = colors.tolist()
print(colors)


# draw nodes, edges and labels
for kcore, nodes in kcores.items():
    print(f'Kcore {kcore}: {nodes}')
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors[kcore])
nx.draw_networkx_edges(G, pos, width=0.05)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.savefig(f'figures/Kcores_{args.network_file}.png')