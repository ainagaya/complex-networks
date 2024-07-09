# Network analysis and dynamics

This project aims to provide tools for analyzing and understanding network structures and dynamics. It consists in 3 parts:

- Network structure and analysis
- Small-world network building
- Dynamics: SIS model

## Table of Contents
- [Requirements](#requirements)
- [Usage](#usage)

## Requirements
- Software:
  - Fortran90
  - Python3 with the following libraries: NetworkX, Matplotlib, Numpy, argparse

- Network to analize

## Usage

All the output files will be placed in the `output` directory. All the figures generated will be placed in the `figures` directory.

### Network structure and analysis

1) In the Makefile, put the name of the network that you wish to analyze in `NETWORK=out.ego-facebook`. 
2) Run `make run_analize`. This will give the number of nodes, edges, the average degree of the network, the list of degrees (degree of each node), the maximum degree of the network, the average cluster coefficient, and will write k_nn(k), P(k) and P_cc(k) to the `output` directory, marked with the name of the network that its being analyzed.
3) Run `make plot` to visualize your results. The figures will be placed in the `figures` directory.

### Small-world network building

SW networks of average degree $\langle k \rangle = 2, 4$ can be created, by executing `make run_SW_k2` or `make run_SW_k4` respectively. This will iterate over the probability of rewiring, with probability values going from 0 to 1, equally spaced in logscale. If more stadistics is needed, the `seed` can be changed.


### Dynamics

1) Place the parameters of your simulation in `input.nml` file.
2) Generate the list of infected nodes by running `make genICs`.
3) Run `make run_dyn` to run Runge Kutta 4 integration, or `make run_dyn_g` to use Gillespie algorithm. If you want to iterate over values, tune and run `submit.sh` script.
4) Generate the corresponding plots by running `make plot_dyn`.

