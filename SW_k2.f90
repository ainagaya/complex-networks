program SW

    use :: network_analysis 

    implicit none

    real :: probability
    integer, allocatable :: edge_list_SW(:,:), list_of_degrees(:), pointers(:, :), V(:)
    real, allocatable :: c(:)
    integer :: N, average_degree, i, k, nn, s
    integer, allocatable :: seed(:)
    real, allocatable :: degree_occurrences(:), acc_degree_occurrences(:)
    integer :: npoints
    real :: loga, logb, logstep, a, b
    ! Network models
    ! Small world network

    character(len=50) :: file_name

    s = 5

    call random_seed(size=nn)
    allocate(seed(nn))
    seed = s    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)


    N = 2888

    allocate(edge_list_SW(N, 2))
    allocate(list_of_degrees(N))
    allocate(pointers(N, 2))
    allocate(V(2*size(edge_list_SW, 1)))

    k = 0
    probability = 0

    call SW_model(N, probability, edge_list_SW)

    write(file_name, "(A, I1, A)") "output/edge_list_SW_p_", k, "_k2.dat"
    open(14, file=file_name)
    do i = 1, N
        write(14, *) edge_list_SW(i, 1), edge_list_SW(i, 2)
    end do
    close(14)

    ! Define the range and number of points
    a = 1.0d-4       ! Start value (10^(-4))
    b = 1.0d0        ! End value (1)
    npoints = 10           ! Number of points

    ! Calculate logarithms of a and b
    loga = log10(a)
    logb = log10(b)

    ! Calculate the spacing in the log scale
    logstep = (logb - loga) / real(npoints-1)

    ! Generate the values in log scale

    do k = 0, npoints-1
    !k = 0
        probability = 10**(loga + k*logstep)
        print*, "Probability: ", probability
        call SW_model(N, probability, edge_list_SW)

        print*, size(edge_list_SW, 1)
        
        call build_list_of_degrees(edge_list_SW, size(edge_list_SW, 1), list_of_degrees)
        call build_pointers_from_array(edge_list_SW, list_of_degrees, N, size(edge_list_SW, 1), V, pointers)

        allocate(degree_occurrences(maxval(list_of_degrees)))
        allocate(acc_degree_occurrences(maxval(list_of_degrees)))
        call calculate_pk(list_of_degrees, degree_occurrences, acc_degree_occurrences)
        call cluster_coefficient(list_of_degrees, pointers, V, N, degree_occurrences, c)
        deallocate(degree_occurrences)
        deallocate(acc_degree_occurrences)
        
        ! Save edge_list_SW to a file
        
        write(file_name, "(A, F6.4, A, I1, A)") "output/edge_list_SW_p_", probability, "_k2_s_", s, ".dat"
        open(14, file=file_name)
        do i = 1, N
            write(14, *) edge_list_SW(i, 1), edge_list_SW(i, 2)
        end do
        close(14)
    end do



    contains 

    subroutine SW_model(N, probability, edge_list)
        implicit none
        integer, intent(in) :: N
        real, intent(in) :: probability
        integer, dimension(N):: first_array, second_array
        integer, dimension(:, :), intent(out) :: edge_list
        logical :: rewire
    
        integer                  :: i, node, j
        real                     :: ran, ran_node
    

        ! Create a ring lattice
        do i = 1, N
            first_array(i) = i
            second_array(i) = i + 1
        end do

        ! Close the ring
        second_array(N) = 1


       ! second_array(N) = 1

        ! Rewire the network
        do i = 1, N
            call random_number(ran)
            if (ran < probability) then
                rewire = .true.
                do while (rewire)
                    call random_number(ran_node)
                    node = int(ran_node*N) + 1
                    if (node == i) then ! Avoid self loops
                        rewire = .true.
                    else if (node == first_array(i) .or. node == second_array(i)) then ! Avoid multiple edges
                        rewire = .true.
                    else
                        rewire = .false.
                        ! Rewire the edge
                    end if
                end do
                second_array(i) = node
            end if
        end do

        ! Create the edge list
        do i = 1, N
            edge_list(i, 1) = first_array(i)
            edge_list(i, 2) = second_array(i)
        end do

        return
    
    
    end subroutine

end program