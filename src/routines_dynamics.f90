module dynamics
    implicit none

    contains

    subroutine infect_random(N, N_infected, list_infected)
        implicit none
        integer, intent(in) :: N
        integer, intent(inout) :: list_infected(:)
        integer, intent(in) :: N_infected
        integer :: i
        integer :: j
        real :: rnd

        do i = 1, N_infected
            print*, "i", i
            do while (.true.)
                call random_number(rnd)
                j = int(rnd*N) + 1
                if (.not.  any(list_infected == j) ) then
                    list_infected(i) = j
                    exit
                else
                    continue
                end if
            end do
        end do

    end subroutine infect_random

    subroutine active_edges()

    end subroutine active_edges

    ! *********************************************************************
!                       RUNGE-KUTTA 4
! *********************************************************************
!     Mètode que calcula 1 pas de Runge Kutta 4.
!     INPUTS: 
!           - t: variable independent on busquem el valor de yyout
!           - dt: pas entre dos t's consecutius
!           - yyin(nequs): vector amb les variables necessàries per
!                          calcular yyout
!           - nequs: nombre d'equacions amb les que treballem.
!     OUPUTS:
!           - yyout(nequs): vector amb les variables de sortida, és a
!                           dir, les imatges de t.
!     És necessària una subrutina "derivada" que calcula la part dreta 
!     de l'equació diferencial.

    Subroutine RK4(t, dt, yyin, yyout, V, pointers, N)
        Implicit none
        integer, intent(in) :: V(:), pointers(:, :), N
        real, intent(in) :: yyin(:)
        real, intent(out) :: yyout(:)
        real, dimension(N) :: k1, k2, k3, k4, vec
        real :: t, dt, tin
        
        tin = t
        Call derivada (yyin, k1, V, pointers, N)
        
        vec = yyin + dt/2.*k1
        tin = t + dt/2.
        Call derivada (vec, k2, V, pointers, N)
  
        vec = yyin + dt/2.*k2
        tin = t + dt/2.
        Call derivada (vec, k3, V, pointers, N)
       
        vec = yyin + dt*k3
  
        tin = t + dt
        Call derivada (vec, k4, V, pointers, N)
  
        yyout = yyin + dt/6. * (k1 + 2.*k2 + 2.*k3 + k4)
  
        Return
    End subroutine
  
    Subroutine derivada(yin, yout, V, pointers, N)
        ! Calcula les derivades que necessitem pel problema.
        ! Per aquest problema: 
        ! {dotrho} = {- delta*rho + lambda sum_j a_j rho_j (1-rho)}      
        ! INPUTS: - yin: vector amb nvars varibles del quan hem de fer les
                ! derivades. En aquest problema, yin = {phi, w}
                ! - h: interval de temps considerat
                ! - nvars: nombre d'EDO's que té el problema
        ! OUTPUTS: yout: vector amb nvars variables, les derivades. En 
        !               aquest problema, {dotphi}
        Implicit none
        integer :: i, j
        real :: delta, lambda, sum
        integer, intent(in) :: V(:), pointers(:, :), N
        real, intent(in) :: yin(:)
        real, intent(out) :: yout(:)
        integer :: neighbor
        common/parameters/delta, lambda

        do i = 1, N
            sum = 0
            ! iterate over neighbours
            do j = pointers(i, 1), pointers(i, 2)
                neighbor = V(j)
                sum = sum + yin(neighbor)*(1-yin(i))

            end do
            yout(i) = -delta*yin(i) + lambda*sum
        end do

  
        Return
    End Subroutine

    subroutine build_rho(N, list_infected, rho)
        implicit none
        integer, intent(in) :: N
        integer, intent(in) :: list_infected(:)
        real, intent(inout) :: rho(:)
        integer :: i

        do i = 1, N
            if (any(list_infected == i)) then
                rho(i) = 1.
            else
                rho(i) = 0
            end if
        end do

    end subroutine build_rho

    subroutine build_active_edges(V, pointers, N, N_infected, edges, list_infected, N_active_edges, list_active_edges)
        implicit none
        integer, intent(in) :: V(:), pointers(:, :), N
        integer, intent(in) :: N_infected
        integer, intent(in) :: list_infected(:)
        integer, intent(out) :: N_active_edges
        integer, intent(inout) :: list_active_edges(:,:)
        integer, intent(in) :: edges
        integer :: i, j
        integer :: infected_node
        integer :: neighbor


        N_active_edges = 0
        list_active_edges = 0

        if (N_infected == 0) then
            return
        end if

        do i = 1, N_infected
            infected_node = list_infected(i)
            if (infected_node /= 0) then
                do j = pointers(infected_node, 1), pointers(infected_node, 2)
                    neighbor = V(j)
                    if (.not. any(list_infected == neighbor)) then
                        N_active_edges = N_active_edges + 1
                        list_active_edges(N_active_edges, 1) = infected_node
                        list_active_edges(N_active_edges, 2) = neighbor
                    else
                        ! The neighbor is already infected
                        continue
                    end if
                end do
            end if
        end do



    end subroutine build_active_edges

    subroutine infect_gillespie(N, N_infected, list_infected, N_active_edges, list_active_edges)
        implicit none
        integer, intent(in) :: N
        integer, intent(inout) :: list_infected(:)
        integer, intent(inout) :: N_infected
        integer, intent(in) :: N_active_edges
        integer, intent(in) :: list_active_edges(:,:)
        integer :: infected_node, infect
        integer :: i
        real :: j

        ! infecting a random node choosing it from the active links list
        call random_number(j)
        infect = int(j*N_active_edges) + 1

        infected_node = list_active_edges(infect, 2)


        N_infected = N_infected + 1
        list_infected(N_infected) = infected_node
        

    end subroutine infect_gillespie

    subroutine recover_gillespie(N, N_infected, list_infected)
        implicit none
        integer, intent(in) :: N
        integer, intent(inout) :: list_infected(:)
        integer, intent(inout) :: N_infected
        integer :: recovered_node, recover
        integer :: i
        real :: j

        ! recovering a random node
        call random_number(j)
        recover = int(j*N_infected) + 1
  

        recovered_node = list_infected(recover)


        if (recover /= N_infected) then
            list_infected(recover) = list_infected(N_infected)
            list_infected(N_infected) = 0
            N_infected = N_infected - 1
        else if (recover == N_infected) then
            list_infected(recover) = 0
            N_infected = N_infected - 1
        else 
            print*, "Error in recovering node"
            stop
        end if
        

    end subroutine recover_gillespie

end module dynamics