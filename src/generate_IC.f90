program generate_ICs

    use :: network_analysis
    use :: dynamics

    implicit none

    integer :: N_infected
    real :: ti, tf, dt, lambda, delta
    integer :: N
    integer, allocatable :: list_infected(:)
    integer :: edges
    integer, allocatable :: pair_list(:, :)
    character(len=100) :: network
    integer :: i
    character(len=100) :: IC_file

    ! Read values from namelist
    namelist /parameters/ N_infected, ti, tf, dt, lambda, delta, IC_file


    print*, "Enter the network filename: "  
    read(*,*) network

    open(unit=10, file='input.nml', status='old')
    read(10, parameters)
    close(10)

    call random_seed()

    call read_network(network, pair_list, N, edges)

    print*, N

    allocate(list_infected(N))

    list_infected = 0

    print*, N, N_infected

    call infect_random(N, N_infected, list_infected)

    print*, "Randomly infected nodes", list_infected(:N_infected)

    open(unit=10, file=IC_file)

    do i = 1, N
        write(10, *) list_infected(i)
    end do

    close(10)



end program