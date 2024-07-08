! Module with subroutines to analyze networks
! Author: Aina Gaya Avila
module network_analysis
    implicit none

    contains

    subroutine read_network(network, pair_list, N, edges)
        ! Read network from file
        ! network: filename
        ! pair_list: list of pairs of nodes
        ! N: number of nodes
        ! edges: number of edges
        implicit none
        character(len=100), intent(in) :: network
        integer, allocatable, intent(out) :: pair_list(:,:)
        integer, intent(out) :: N, edges
        integer, allocatable :: max_list_of_nodes(:)
        integer :: current_node_1, current_node_2, i, k, stat
        
        open(11, file=network)

        print*, "Reading network from file: ", network
        print*, "**********************************"

        edges = 0

        ! Count number of edges (rows in the file)
        read_loop : do while ( .true. )
            read(11, *, iostat=stat) k, k
            if (stat .eq. 0) then
                edges = edges + 1   
            else if (stat .eq. -1) then ! end of file
                exit read_loop
            end if 
                
        end do read_loop

        rewind(11)

        allocate(max_list_of_nodes(edges))
        allocate(pair_list(edges, 2))

        max_list_of_nodes = 0
        N = 0

        ! Build the list of nodes and the list of pairs
        do i = 1, edges
            read(11, *) current_node_1, current_node_2
            pair_list(i,1) = current_node_1
            pair_list(i,2) = current_node_2
            if (.not. any(current_node_1 == max_list_of_nodes)) then
                N = N + 1
                max_list_of_nodes(N) = current_node_1
            end if

            if (.not. any(current_node_2 == max_list_of_nodes)) then
                N = N + 1
                max_list_of_nodes(N) = current_node_2
            end if
        end do
    
    end subroutine

    ! ******************************************************************
    subroutine build_list_of_degrees(pair_list, edges, list_of_degrees)
        ! Build list of degrees from pair list
        ! pair_list: list of pairs of nodes
        ! edges: number of edges
        ! list_of_degrees: list with the degree of each node
        implicit none
        integer, intent(in) :: pair_list(:,:), edges
        integer, intent(out) :: list_of_degrees(:)

        integer :: i

        list_of_degrees = 0
        
        do i = 1, edges
            list_of_degrees(pair_list(i,1)) = list_of_degrees(pair_list(i,1)) + 1
            list_of_degrees(pair_list(i,2)) = list_of_degrees(pair_list(i,2)) + 1
        end do

        return
    end subroutine build_list_of_degrees

    ! *************************************************************************
    subroutine build_pointers_from_file(list_of_degrees, N, edges, V, pointers)
        ! Build pointers from file 11
        ! list_of_degrees: list with the degree of each node
        ! N: number of nodes
        ! edges: number of edges
        ! V: list of edges
        ! pointers: pointers to the list of edges for each node i, 
        !           pointers(i, 1) is the initial position, 
        !           pointers(i, 2) is the final position
        implicit none
        integer, intent(in) :: list_of_degrees(:), N, edges
        integer, intent(out) :: V(:), pointers(:,:)

        integer :: i, current_node_1, current_node_2

        pointers = 0

        do i = 1, N
            ! initial position of the pointer
            pointers(i, 1) = sum(list_of_degrees(:i-1)) + 1
            ! initialize final position of the pointer
            pointers(i, 2) = pointers(i, 1) 
        end do 
    
        rewind(11)
        
        ! Read again the edges
        do i = 1, edges
            read(11, *) current_node_1, current_node_2
            V(pointers(current_node_1, 2)) = current_node_2
            pointers(current_node_1, 2) = pointers(current_node_1, 2) + 1 
            V(pointers(current_node_2, 2)) = current_node_1
            pointers(current_node_2, 2) = pointers(current_node_2, 2) + 1
        end do 
    
        ! Fix pointers
        pointers(:, 2) = pointers(:, 2) - 1
    
        return
    end subroutine build_pointers_from_file

    ! *************************************************************************
    subroutine build_pointers_from_array(list_of_nodes, list_of_degrees, N, edges, V, pointers)
        ! Build pointers from array
        ! list_of_nodes: list of pairs of nodes with the format node1 node2
        ! list_of_degrees: list with the degree of each node
        ! N: number of nodes
        ! edges: number of edges
        ! V: list of edges
        ! pointers: pointers to the list of edges for each node i
        
        implicit none
        integer, intent(in) :: list_of_nodes(:, :), N, edges, list_of_degrees(:)
        integer, intent(out) :: V(:), pointers(:,:)

        integer :: i, current_node_1, current_node_2

        pointers = 0

        do i = 1, N
            ! initial position of the pointer
            pointers(i, 1) = sum(list_of_degrees(:i-1)) + 1
            ! initialize final position of the pointer
            pointers(i, 2) = pointers(i, 1) 
        end do 
    
    
        do i = 1, edges
            current_node_1 = list_of_nodes(i, 1)
            current_node_2 = list_of_nodes(i, 2) 
            V(pointers((current_node_1), 2)) = current_node_2
            pointers(current_node_1, 2) = pointers(current_node_1, 2) + 1 
            V(pointers((current_node_2), 2)) = current_node_1
            pointers(current_node_2, 2) = pointers(current_node_2, 2) + 1
        end do 
    
        pointers(:, 2) = pointers(:, 2) - 1
    
        return

    end subroutine build_pointers_from_array

    ! *************************************************************************
    subroutine calculate_pk(list_of_degrees, degree_occurrences, acc_degree_occurrences)
        ! Calculate P(k) and cc P(k)
        ! list_of_degrees: list with the degree of each node
        ! degree_occurrences: P(k)
        ! acc_degree_occurrences: cc P(k)
        implicit none
    
        integer :: i
        integer, intent(in) :: list_of_degrees(:)
        real, intent(out) :: degree_occurrences(:), acc_degree_occurrences(:)
    
        degree_occurrences = 0
        acc_degree_occurrences = 0
    
        do i = 1, size(list_of_degrees)
            degree_occurrences(list_of_degrees(i)) = degree_occurrences(list_of_degrees(i)) + 1
        end do 
    
        ! Normalize
        degree_occurrences = degree_occurrences/size(list_of_degrees)
    
        do i = 1, maxval(list_of_degrees)
            acc_degree_occurrences(i) = sum(degree_occurrences(:i))
        end do 
    
        return
    end subroutine calculate_pk
        
    ! *************************************************************************  
    subroutine calculate_knn(list_of_degrees, degree_occurrences, pointers, V, N, k_nn)
        ! Calculate k_nn: average nearest neighbor degree
        ! list_of_degrees: list with the degree of each node
        ! degree_occurrences: P(k)
        ! pointers: pointers to the list of edges for each node i
        ! V: list of edges
        ! N: number of nodes
        
        implicit none
    
        integer :: i, z, neigh, degree_node, degree_neigh
        integer, intent(in) :: list_of_degrees(:), pointers(:, :), V(:), N
        real :: degree_occurrences(:)
        real, allocatable, intent(out) :: k_nn(:)

        real :: average_degree, average_degree_squared, kappa
    
        allocate(k_nn(maxval(list_of_degrees))) 
        
        k_nn = 0
    
        do i = 1, N
            degree_node = list_of_degrees(i)
            do z = pointers(i, 1), pointers(i, 2)
                neigh = V(z)
                degree_neigh = list_of_degrees(neigh)
                k_nn(degree_node) = k_nn(degree_node) + degree_neigh/(degree_node*degree_occurrences(degree_node))
            end do
        end do
    
        call calculate_average_degree(list_of_degrees, average_degree)
        call calculate_average_degree_squared(list_of_degrees, average_degree_squared)

        kappa = average_degree_squared/average_degree
        k_nn = k_nn/kappa
    
    ! print*, "k_nn", k_nn
    
        return
    end subroutine calculate_knn
 
    ! *************************************************************************  
    subroutine cluster_coefficient(list_of_degrees, pointers, V, N, degree_occurrences, c)
        ! average clustering coefficient:  for each
        ! node i (1 loop) increment the corresponding entry in (1
        ! vector) with the local contribution of all neighbors of i (number
        ! of triangles normalized by the maximum possible number of
        ! triangles given the degree of node i and by number of nodes in
        ! degree class k)
        ! To calculate the number of triangles attached to a node, check
        ! connections between all possible pairs of neighbors of the node
        ! (3 nested loops)

        implicit none

        integer :: i, j, k, triangles, current_degree, max_triangles
        integer, intent(in) :: list_of_degrees(:), pointers(:, :), V(:), N
        real, intent(in) :: degree_occurrences(:)
        real, intent(out), allocatable :: c(:)


        allocate(c(N))
        c = 0

        ! no haig de iterar per nodes sino per degrees!!!
        do i = 1, N
            triangles = 0
            current_degree = list_of_degrees(i)
            if (current_degree > 2) then
                max_triangles = current_degree*(current_degree-1)*(current_degree-2)/6
            else
                max_triangles = 0
            end if
            do j = V(pointers(i, 1)), V(pointers(i, 2))
                do k = V(pointers(j, 1)), V(pointers(j, 1))
                    triangles = triangles + count(V(pointers(k, 1):pointers(k, 2)) == i)
    !                print*, i, j, k
                end do
            end do
            if ((max_triangles /= 0).and.(degree_occurrences(current_degree) /= 0 ) ) then
                c(current_degree) = c(current_degree) + triangles/(max_triangles*degree_occurrences(current_degree))
            end if
        end do

    end subroutine cluster_coefficient

     ! *************************************************************************
    ! subroutine average_shortest_path_lenght(list_of_degrees, pointers, V, N, l)
    !     implicit none
    !     integer, intent(in) :: list_of_degrees(:), pointers(:, :), V(:), N
    !     real, intent(out) :: l(:)

    !     integer :: i, j

    !     ! average shortest path length: for each node i (1 loop) calculate
    !     ! the shortest path length to all other nodes (N-1 loops) and
    !     ! increment the corresponding entry in the vector (1 vector) 

    !     ! We need to link N*(N-1)/2 pairs of nodes

    !     do i = 1, N
    !         do j = i+1, N
    !             ! calculate shortest path length
    !             ! TODO
    !         end do
    !     end do

    ! end subroutine average_shortest_path_lenght

    subroutine calculate_average_degree(list_of_degrees, average_degree)
        ! Calculate the average degree of the network
        ! list_of_degrees: list with the degree of each node
        ! average_degree: average degree of the network
        implicit none
        integer, intent(in) :: list_of_degrees(:)
        real, intent(out) :: average_degree

        average_degree = real(sum(list_of_degrees))/size(list_of_degrees)

        return
    end subroutine calculate_average_degree

    subroutine calculate_average_degree_squared(list_of_degrees, average_degree_squared)
        ! Calculate the average degree squared of the network
        ! list_of_degrees: list with the degree of each node
        ! average_degree_squared: average degree squared of the network
        implicit none
        integer, intent(in) :: list_of_degrees(:)
        real, intent(out) :: average_degree_squared

        average_degree_squared = sum(list_of_degrees**2)/size(list_of_degrees)

        return
    end subroutine calculate_average_degree_squared

end module network_analysis