
program read_net
    ! This program reads a network from a file and calculates 
    ! structural properties of the netowrk such as the degree distribution,
    ! the clustering coefficient and the average degree.
    
    ! The network is read from a file with the following format:
    ! node1 node2

    ! The output is written to the output folder.
    
    ! Author: Aina Gaya-Àvila
    ! Date: 2024-06-20

    use :: network_analysis 

    implicit none

    ! Number of edges and nodes of the network
    integer :: edges, N
    ! iterator
    integer :: i
    ! Variables to store the network
    integer, allocatable :: pair_list(:, :), list_of_degrees(:), V(:), pointers(:,:)
    ! Variables to store the average degree
    real :: average_degree, average_degree_2
    ! Variables to store the degree distribution
    real, allocatable :: degree_occurrences(:), acc_degree_occurrences(:), k_nn(:)
    ! Name of the network file
    character(len=100) :: network
    ! Variables to store the clustering coefficient
    real, allocatable :: c(:)

    call random_seed()

    ! Read the network filename
    print*, "Enter the network filename: "
    read(*,*) network

    print*, "**********************************"
    print*, "Analyzing network: ", network
    print*, "**********************************"

    call read_network(network, pair_list, N, edges)

    allocate(V(2*edges))
    allocate(pointers(N, 2))
    allocate(list_of_degrees(N)) 
    
    call build_list_of_degrees(pair_list, edges, list_of_degrees)
    call build_pointers_from_file(list_of_degrees, N, edges, V, pointers)
    
    ! The average degree can be calculated in two ways: 
    ! 1) sum of all degrees divided by the number of nodes
    ! 2) 2*number of edges divided by the number of nodes
    average_degree_2 = real(2*edges)/real(N)

    call calculate_average_degree(list_of_degrees, average_degree)

    print*, "number of nodes: ", N
    print*, "number of edges:", edges
    print*, "average degree: ", average_degree, average_degree_2

    print*, "list of degrees: ", list_of_degrees

    print*, "***********************************************"

    close(11)

    ! Calculate P(k) and cc P(k)
    allocate(degree_occurrences(maxval(list_of_degrees)))
    allocate(acc_degree_occurrences(maxval(list_of_degrees)))

    print*, "max degree (out sub): ", maxval(list_of_degrees)

    open(11, file="output/pk_" // trim(network) // ".dat")

    call calculate_pk(list_of_degrees, degree_occurrences, acc_degree_occurrences)

    do i = 1, maxval(list_of_degrees)
        write(11, *) i, degree_occurrences(i), acc_degree_occurrences(i)
    end do

    close(11)

    ! calculate k_nn(k)

    open(12, file="output/knn_" // trim(network) // ".dat")

    call calculate_knn(list_of_degrees, degree_occurrences, pointers, V, N, k_nn)

    do i = 1, maxval(list_of_degrees)
        write(12, *) i, k_nn(i)
    end do

    close(12)

    open(13, file="output/cc_" // trim(network) // ".dat")

    ! calculate clustering coefficient
    call cluster_coefficient(list_of_degrees, pointers, V, N, degree_occurrences, c)

    do i = 1, N
        write(13, *) i, c(i)
    end do

    close(13)

    print*, "average clustering coefficient: ", sum(c)/N

    print*, "***********************************************"
    print*, "Analysis completed"
    print*, "k_nn(k), P(k) and P_cc(k) written to output folder"
    print*, "***********************************************"

end program

