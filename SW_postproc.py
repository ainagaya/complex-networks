import networkx as nx
import matplotlib.pyplot as plt
import glob
import argparse
import numpy as np

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

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-k', type=int, choices=[2, 4], help='Value of k (2 or 4)')
args = parser.parse_args()

# Append k2 or k4 to the file pattern
if args.k == 2:
    file_pattern = "output/edge_list_SW_p_*_k2_s_*.dat"
elif args.k == 4:
    file_pattern = "output/edge_list_SW_p_*_k4_s_*.dat"
else:
    print("Invalid value of k. Please choose either 2 or 4.")

# Read the file with p=0
file_name_p0 = f"output/edge_list_SW_p_0_k{args.k}.dat"
graph_p0 = nx.read_edgelist(file_name_p0)

# Calculate the cluster coefficient and path length for p=0
if nx.is_connected(graph_p0):
    cluster_coefficient_p0 = nx.average_clustering(graph_p0)
    path_length_p0 = nx.average_shortest_path_length(graph_p0)
else:
    print("Graph is not connected.")
    exit()

print(f"Cluster coefficient for p=0: {cluster_coefficient_p0}")
print(f"Path length for p=0: {path_length_p0}")
    


# Initialize lists to store the cluster coefficients and path lengths
cluster_coefficients = []
path_lengths = []
probabilities = []

print("File pattern", file_pattern)

# Iterate over the files matching the pattern
for file_name in sorted(glob.glob(file_pattern)):

    print(file_name)
    # Extract the value of p from the file name
    p = float(file_name.split("_")[4].split("_s_*.dat")[0])
    print(p)
    # Append the value of p to the array of probabilities
    probabilities.append(p)
    
    # Read the edge list from the file
    graph = nx.read_edgelist(file_name)

    # print(graph.edges())

    # Plot the graph
    # plt.figure()
    # nx.draw_networkx(graph, with_labels=True,  node_size=1, width=0.1)
    # plt.title(file_name)
    # plt.show()

    # Calculate the cluster coefficient and path length
    if nx.is_connected(graph):
        cluster_coefficient = nx.average_clustering(graph)
        path_length = nx.average_shortest_path_length(graph)
    else:
        print("Graph is not connected. Adding connections between connected components.")
        Gcc = sorted(nx.connected_components(graph), key=len, reverse=True)
        G0 = graph.subgraph(Gcc[0])

        
        # Add links between connected components until the network is connected
        for i in range(1, len(Gcc)):
            Gcc_i = graph.subgraph(Gcc[i])
            nodes_G0 = list(G0.nodes())
            nodes_Gcc_i = list(Gcc_i.nodes())
            node_G0 = nodes_G0[0]
            node_Gcc_i = nodes_Gcc_i[0]
            graph.add_edge(node_G0, node_Gcc_i)
            G0 = graph.subgraph(list(G0.nodes()) + list(Gcc_i.nodes()))


        cluster_coefficient = nx.average_clustering(G0)
        path_length = nx.average_shortest_path_length(G0)

    # Plot the graph
    plt.figure()
    nx.draw_circular(graph, with_labels=False, node_size=1, width=0.1)
    # Remove the "output/" part from the file name
    file_name = file_name.replace("output/", "")
    plt.title(file_name)
    plt.savefig(f'figures/network_{file_name}.png')
    plt.clf()
        
    # Append the results to the lists
    cluster_coefficients.append(cluster_coefficient)
    path_lengths.append(path_length)

    print(f"Cluster coefficient: {cluster_coefficient}")
    print(f"Path length: {path_length}")

print(probabilities)
print(cluster_coefficients)
print(path_lengths)

# Convert the lists to numpy arrays
probabilities = np.array(probabilities)
cluster_coefficients = np.array(cluster_coefficients)
path_lengths = np.array(path_lengths)

# Find unique probabilities
unique_probabilities = np.unique(probabilities)

# Initialize lists to store the averaged cluster coefficients and path lengths
averaged_cluster_coefficients = []
averaged_path_lengths = []
error_bars_cluster_coefficients = []
error_bars_path_lengths = []

# Iterate over the unique probabilities
for p in unique_probabilities:
    # Find the indices where the probability matches
    indices = np.where(probabilities == p)
    # Extract the cluster coefficients and path lengths for the matching indices
    cc_values = cluster_coefficients[indices]/cluster_coefficient_p0
    pl_values = path_lengths[indices]/path_length_p0
    # Calculate the mean and standard deviation
    mean_cc = np.mean(cc_values)
    mean_pl = np.mean(pl_values)
    std_cc = np.std(cc_values)
    std_pl = np.std(pl_values)
    # Append the averaged values and error bars to the lists
    averaged_cluster_coefficients.append(mean_cc)
    averaged_path_lengths.append(mean_pl)
    error_bars_cluster_coefficients.append(std_cc)
    error_bars_path_lengths.append(std_pl)

# Plot the cluster coefficients and path lengths with error bars
plt.errorbar(unique_probabilities, averaged_cluster_coefficients, yerr=error_bars_cluster_coefficients, fmt='.', label='Cluster Coefficient', color = 'red')
plt.errorbar(unique_probabilities, averaged_path_lengths, yerr=error_bars_path_lengths, fmt='.', label='Average Path Length', color = 'green')
plt.xscale('log')
plt.xlabel('Probability')
plt.ylabel('')
plt.legend()
plt.savefig(f'figures/cluster_coefficient_path_length_k{args.k}.png')