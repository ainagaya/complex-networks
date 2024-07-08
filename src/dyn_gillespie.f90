program dynamicsGIL

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
    ! Current recover probability
    real :: recover_probability 
    ! Number of current active edges
    integer :: N_active_edges
    ! List of current active edges
    integer, allocatable :: list_active_edges(:,:)
    ! Step counter
    integer :: step
    ! Random number, time increment
    real :: rnd, tau
    ! Count of nonzero elements in list_active_edges (contingency)
    integer :: count_nonzero
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
    print*, "Starting dynamicsGIL."
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

    ! Initial time
    t = ti
    step = 0

    allocate(rho(N))
    call build_rho(N, list_infected, rho)

    allocate(list_active_edges(2*edges, 2))

    write(lambda_str, '(f4.2)') lambda

    open(12, file="output/GIL_" // trim(network) // "_lambda_" // trim(adjustl(lambda_str)) // ".dat")

    ! Main loop: iterate until final time
    do while (t < tf)
        call build_active_edges(V, pointers, N, N_infected, edges, list_infected, N_active_edges, list_active_edges)     
        recover_probability = N_infected / (N_infected + lambda*N_active_edges)
        
        if (mod(step,1000) == 0) then
            print*, "Step ", step, " at time ", t, " of ", tf
            print*, "Recover probability: ", recover_probability
            print*, "N_infected: ", N_infected
            print*, "N_active_edges: ", N_active_edges
        end if 

        ! Time increment mechanism
        call random_number(rnd)
        tau = -log(1 - rnd) / (N_infected + lambda * N_active_edges)

        if (tau > 0.1*dt) then
            t = t + tau
        else
            t = t + dt
        end if
        
        call random_number(rnd)

        if (rnd < recover_probability) then
            ! Recover a node
            if (N_infected > 0) then
                call recover_gillespie(N, N_infected, list_infected)
            else
                ! No infected nodes to recover
                continue
            end if
        else
            ! Infect a node
            if (N_infected < N) then
                call infect_gillespie(N, N_infected, list_infected, N_active_edges, list_active_edges)
            else
                ! No susceptible nodes to infect
                continue
            end if
        end if

        step = step + 1

        ! Contingency mesure
        count_nonzero = count(list_infected /= 0)
        if (count_nonzero /= N_infected) then
            print*, "Error in infected list"
            print*, "Resetting N_infected to ", count_nonzero
            print*, "N_infected: ", N_infected
            print*, "lst_infected: ", list_infected
            N_infected = count_nonzero
        end if

        call build_rho(N, list_infected, rho)
        write(12,*) t, sum(rho)/N

    end do

    close(12)

    print*, "Run complete."

    print*, "Saving final infected list"
    open(unit=13, file="output/infected_list_" // trim(network) // "_lambda_" // trim(adjustl(lambda_str)) // ".dat")
    do i = 1, N_infected
        write(13, *) list_infected(i)
    end do
    close(13)
    
end program