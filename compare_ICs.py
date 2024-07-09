import matplotlib.pyplot as plt
import numpy as np

# List of file names
random = 'output/GIL_out.ego-facebook_lambda_0.50_IC_random-infected.dat' 
hubs = 'output/GIL_out.ego-facebook_lambda_0.50_IC_hubs-infected.dat'
marginal = 'output/GIL_out.ego-facebook_lambda_0.50_IC_margin-infected.dat'

file_names = [random, hubs, marginal]
colors = ['magenta', 'tomato', 'brown']

# Read data from each file
for file_name in file_names:
    data = np.loadtxt(file_name)
    time = data[:, 0]
    rho = data[:, 1]
    label = "GIL" + file_name.split('/')[-1].split('_')[-1].split('-')[0]
    plt.plot(time, rho, label=label, color=colors.pop(0))

random = 'output/RK_out.ego-facebook_lambda_0.50_IC_random-infected.dat' 
hubs = 'output/RK_out.ego-facebook_lambda_0.50_IC_hubs-infected.dat'
marginal = 'output/RK_out.ego-facebook_lambda_0.50_IC_margin-infected.dat'

file_names = [random, hubs, marginal]
colors = ["cyan", "green", "blue"]

# Read data from each file
for file_name in file_names:
    data = np.loadtxt(file_name)
    time = data[:, 0]
    rho = data[:, 1]
    label = "RK" + file_name.split('/')[-1].split('_')[-1].split('-')[0]
    plt.plot(time, rho, label=label, color=colors.pop(0))

# Plot the data
plt.xlabel('t')
plt.ylabel(r'$\rho$')
plt.grid(True)
plt.legend()
plt.savefig('figures/compare_ICs.png')