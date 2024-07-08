# Fortran compiler
FC=gfortran

COMP_D_FLAGS=-Wall -Wextra -g -fcheck=all -Wall -Wextra -pedantic -std=f2008 -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan -finit-integer=-9999

#NETWORK=inf-openflights.edges
NETWORK=out.ego-facebook

net.o: read_net.f90 nets_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@

nets_subs.o: subroutines_XC.f90
	$(FC) $(COMP_D_FLAGS) -c $^ -o $@

SW_k4.o: SW_k4.f90 nets_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@

SW_k2.o: SW_k2.f90 nets_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@

dyn_subs.o: routines_dynamics.f90
	$(FC) $(COMP_D_FLAGS) -c $^ -o $@

dyn.o: dyn.f90 nets_subs.o dyn_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@ 

dyn_g.o: dyn_gillespie.f90 nets_subs.o dyn_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@

genICs.o: generate_IC.f90 nets_subs.o dyn_subs.o
	$(FC) $(COMP_D_FLAGS) $^ -o $@


.PHONY.:
run_analize: net.o output
	echo "${NETWORK}" > input
	./$< < input

output:
	mkdir -p output

figures:
	mkdir -p figures

.PHONY.: run_SW_k2 
run_SW_k2: SW_k2.o output
	./$^

.PHONY.: run_SW_k4 
run_SW_k4: SW_k4.o output
	./$^

.PHONY.: post_SW figures
post_SW:
	python3 SW_postproc.py -k 4

.PHONY.: run_dyn
run_dyn: dyn.o output
	echo "${NETWORK}" > input
	./$^ < input

.PHONY.: run_dyn_g
run_dyn_g: dyn_g.o output
	echo "${NETWORK}" > input
	./$^ < input

.PHONY.: genICs
genICs: genICs.o
	echo "${NETWORK}" > input
	./$^ < input	

.PHONY.: plot_dyn
plot_dyn: figures
	python3 plot_dynamics.py --network_file ${NETWORK}


.PHONY.: plot
plot: figures
	python3 read_net.py --network_file ${NETWORK} --draw spring
	python3 plot_degree_occ.py --network_file ${NETWORK}
	python3 plot_knn.py --network_file ${NETWORK}
