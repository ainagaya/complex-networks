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

### Network structure and analysis

1) In the Makefile, put the name of the network that you wish to analyze in `NETWORK=out.ego-facebook`. 
2) Run `make run_analize`. This will give the number of nodes, edges, the average degree of the network, the list of degrees (degree of each node), the maximum degree of the network, the average cluster coefficient, and will write k_nn(k), P(k) and P_cc(k) to the `output` directory, marked with the name of the network that its being analyzed.
3) Run `make plot` to visualize your results. The figures will be placed in the `figures` directory.

### Small-world


### Dynamics

