import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib as mpl

plt.rcParams['text.usetex'] = True

def extract_steady_state(file_path):
    data = np.loadtxt(file_path)
    print(file_path)
    if data.size == 0:
        print(f'File {file_path} is empty')
        return 0
    time = data[:, 0]
    rho = data[:, 1]
    steady_state = np.mean(rho[-1000:])  # Assuming steady state is reached in the last 100 time steps
    strd_desv = np.std(rho[-1000:])
    print(f'Steady state: {steady_state}')
    return steady_state, strd_desv

def plot_average_rho_vs_lambda(file_pattern):
    files = glob.glob(file_pattern)
    lambdas = []
    average_rho_values = []
    files.sort()

    for file in files:
        print(file)
        lambda_value = file.split('_')[-1].split('.')[0:2]
        lambda_value = '.'.join(lambda_value)
        method = file.split('_')[0].split('/')[-1]
        print(lambda_value)
        steady_state, strd_desv = extract_steady_state(file)
        lambdas.append(float(lambda_value))
        average_rho_values.append(steady_state)

    if method == "RK":
        plt.errorbar(lambdas, average_rho_values, yerr=strd_desv, fmt='-', label='RK', color="peru")
    elif method == "GIL":
        plt.errorbar(lambdas, average_rho_values, yerr=strd_desv, fmt='.', markersize = 10 ,label='GIL', color="saddlebrown")
    plt.xlabel(r' $\lambda$ ')
    plt.ylabel(r'$\rho_{st} $ ')
#    plt.title('Average Rho vs Lambda')
    plt.grid(True)

def plot_rho_vs_t(file_pattern):
    files = glob.glob(file_pattern)
    files.sort()
    print(files)
    n_lines = files.__len__()
    cmap = mpl.colormaps['rainbow']

    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))
    colors = colors.tolist()
    print(f"Colors: {colors}")

    plt.show()
    for file in files:
        color = colors.pop(0)
        print(file)
        lambda_value = file.split('_')[-1].split('.')[0:2]
        lambda_value = '.'.join(lambda_value)
        data = np.loadtxt(file)
        if data.size > 0:
            time = data[:, 0]
            rho = data[:, 1]
            plt.plot(time, rho, label=f'$\lambda = $ {lambda_value}', color = color)
    plt.xlabel('$t$')
    plt.ylabel(r'$\rho$')
#    plt.title('Rho vs Time')
    plt.grid(True)
    plt.legend()

def plot_rho_vs_t_compare(file_pattern1, file_pattern2):
    files1 = glob.glob(file_pattern1)
    files2 = glob.glob(file_pattern2)
    files1.sort()
    files2.sort()
    colors_1 = ["lightcoral", "greenyellow", "lightskyblue", "plum"]
    colors_2 = ["darkred", "darkolivegreen", "deepskyblue", "indigo"]
    for file in files1:
        data1 = np.loadtxt(file)
        lambda_value = file.split('_')[-1].split('.')[0:2]
        lambda_value = '.'.join(lambda_value)
        if data1.size > 0:
            if lambda_value in list_to_plot:
                color = colors_1.pop(0)
                time1 = data1[:, 0]
                rho1 = data1[:, 1]
                plt.plot(time1, rho1, '-', label=f"RK $\lambda =$ {lambda_value}", color = color)
    for file in files2:
        data2 = np.loadtxt(file)
        lambda_value = file.split('_')[-1].split('.')[0:2]
        lambda_value = '.'.join(lambda_value)
        if data2.size > 0:
            if lambda_value in list_to_plot:
                color = colors_2.pop(0)
                time2 = data2[:, 0]
                rho2 = data2[:, 1]
                plt.plot(time2, rho2, '.', markersize = 0.5, label=f"GIL $\lambda =$ {lambda_value}", color = color)    
    plt.xlabel('Time (step)')
    plt.ylabel(r'$\rho$')
    plt.title('Rho vs Time')
    plt.grid(True)
    plt.legend() 

parser = argparse.ArgumentParser()
parser.add_argument('--network_file', type=str, help='Path to the edge list file')
args = parser.parse_args()

file_pattern1 = 'output/RK_*_*_*.dat' 

print("Plotting rho vs t for RK")
#plot_rho_vs_t(file_pattern1)
#plt.savefig(f'figures/average_rho_vs_t_{args.network_file}_RK.png')
plt.clf()

print("Plotting rho vs lambda for RK")
#plot_average_rho_vs_lambda(file_pattern1)
#plt.savefig(f'figures/average_rho_vs_lambda_{args.network_file}_RK.png')
#plt.clf()
#print("Plotting rho vs lambda for GIL")
file_pattern2 = 'output/GIL_*_*_*.dat' 
#plot_average_rho_vs_lambda(file_pattern2)
#plt.legend()
#plt.savefig(f'figures/compare_rho_vs_lambda_{args.network_file}_GIL.png')
#plt.clf()

print("Plotting rho vs t for GIL")
plot_rho_vs_t(file_pattern2)
plt.xscale('log')
plt.savefig(f'figures/average_rho_vs_t_{args.network_file}_GIL.png')
plt.clf()

list_to_plot = ["1.00", "2.00", "5.00", "9.00"]

plot_rho_vs_t_compare(file_pattern1, file_pattern2)
plt.savefig(f'figures/compare_rho_vs_t_{args.network_file}.png')
plt.clf()


file_GIL="output/GIL_out.ego-facebook_lambda_1.00.dat"
file_RK="output/RK_out.ego-facebook_lambda_1.00.dat"
lambda_value = file_GIL.split('_')[-1].split('.')[0:2]
lambda_value = '.'.join(lambda_value)
data = np.loadtxt(file_RK)
time = data[:, 0]
rho = data[:, 1]
plt.plot(time, rho, label=f"RK $\lambda =$ {lambda_value}", color = "peru", linewidth = 2)
data = np.loadtxt(file_GIL)
time = data[:, 0]
rho = data[:, 1]
plt.plot(time, rho, '.', label=f"GIL $\lambda =$ {lambda_value}", color = "saddlebrown", markersize = 1 )
plt.xlabel('$t$')
plt.ylabel(r'$\rho$')
#    plt.title('Rho vs Time')
plt.grid(True)
plt.legend()
plt.savefig(f'figures/lambda_1.00.png')


