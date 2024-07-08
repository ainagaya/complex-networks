! This program simulates the dynamics of a disease spreading over a network
! using a Runge Kutta 4. The network is read from a file, and the
! initial conditions are read from another file. The parameters of the
! simulation are read from a namelist file. The output is written to a file
! with a name that includes the network filename and the value of lambda.
! Author: Aina Gaya Àvila
!         July 2024
! Usage: Modify inputs.nml to set the parameters of the simulation, then
!       (if desired) run genICs to generate initial conditions, and finally
!       `make run_dyn` to compile and run the program

program dynamicsRK

    use :: network_analysis
    use :: dynamics

    implicit none

    ! Number of initially infected nodes (defined in namelist)
    integer :: N_infected 
    ! Initial time, final time, time step, lambda, delta (defined in namelist)
    real :: ti, tf, dt, lambda, delta 
    ! Number of nodes, number of edges (defined in network)
    integer :: N, edges 
    ! Network filename
    character(len=100) :: network   
     ! List of infected nodes (defined in initial conditions)
    integer, allocatable :: list_infected(:)
    ! Loop variable
    integer :: i 
    ! Number of time steps
    integer :: Nsteps 
    ! Time
    real :: t 
    ! List of degrees, list of neighbors, pointers, pairs (defined in network)
    integer, allocatable :: list_of_degrees(:), V(:), pointers(:, :), pair_list(:, :)
    ! Density of infected nodes 
    real, allocatable :: rho(:), rho_new(:) 
    ! String for lambda
    character(len=10) :: lambda_str 
    ! File with initial conditions
    character(len=100) :: IC_file 
    ! Seed for random number generator
    integer, allocatable :: seed(:)
    integer :: nn, s
    common/parameters/delta, lambda

    ! Read values from namelist
    namelist /parameters/ N_infected, ti, tf, dt, lambda, delta, IC_file

    ! Set seed for random number generator
    s = 5
    call random_seed(size=nn)
    allocate(seed(nn))
    seed = s    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)

    print*, "***************************************"
    print*, "Starting dynamicsRK."
    print*, "***************************************"

    print*, "(When running interactively): Enter the network filename: "  
    read(*,*) network

    ! Open namelist file and read parameters
    open(unit=10, file='input.nml', status='old')
    read(10, parameters)
    close(10)

    print*, "RUNNING WITH PARAMETERS:"

    print*, "N_infected: ", N_infected
    print*, "ti: ", ti
    print*, "tf: ", tf
    print*, "dt: ", dt
    print*, "lambda: ", lambda
    print*, "delta: ", delta

    call read_network(network, pair_list, N, edges)
    print*, "Finishing reading net: ", network

    allocate(list_of_degrees(N))
    allocate(V(2*edges))
    allocate(pointers(N, 2))
    allocate(list_infected(N))

    call build_list_of_degrees(pair_list, edges, list_of_degrees)
    print*, "Finished building list of degrees"
    call build_pointers_from_file(list_of_degrees, N, edges, V, pointers)
    print*, "Finished building pointers"

    print*, "Reading initial conditions from file", IC_file

    open(unit=10, file=IC_file, status='old')
    do i = 1, N
        read(10, *) list_infected(i)
    end do

    ! Compute the number of timsteps from the parameters defined in the namelist
    Nsteps = int((tf - ti)/dt)
    ! Initialize time
    t = ti

    allocate(rho(N))
    call build_rho(N, list_infected, rho)

    allocate(rho_new(N))
    
    write(lambda_str, '(f4.2)') lambda

    open(12, file="output/RK_" // trim(network) // "_lambda_" // trim(adjustl(lambda_str)) // ".dat")

    ! Main loop: iterate over time steps
    do i = 1, Nsteps
        t = ti + i*dt
        if (mod(i, 1000) == 0) then
            print*, "Current t: ", t
        end if
        call RK4(t, dt, rho, rho_new, V, pointers, N)
        rho = rho_new
        write(12,*) t, sum(rho)/N
    end do

    close(12)

    print*, "Run complete."

end program