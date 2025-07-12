! =============*FORTRAN*============ !
!   _ __ ___  _ __ (_) __| | ___| |  !
!  | '_ ` _ \| '_ \| |/ _` |/ __| |  !
!  | | | | | | |_) | | (_| | (__| |  !
!  |_| |_| |_| .__/|_|\__,_|\___|_|  !
!            |_|                     !
! ================================== !
! ================================================================= !
!  Copyright (C) 2025, Simon Kebinger
!
!  This file is part of the MPI decomposition library "mpidcl" for
!  structured multidmensional domains.
!
!  This library is distributed under the BSD 3-Clause License.
! ================================================================= !

submodule(MOD_MPI_decomposition) SMOD_info
    !! Submodule contains various subroutines to retrieve information about the
    !! decomposition, such as cartesian rank layout or the index decomposition
    !! summary.

    implicit none
contains

    module subroutine print_decomposition_summary(info, comm)
        !! Prints the domain decomposition summary for all MPI ranks.
        !!
        !! Gathers the local index ranges from all ranks and prints a formatted
        !! table on rank 0 showing the (i, j, k) bounds for each rank's block.
        !! Only active dimensions (m_xi, m_eta, m_tau) are shown.
        !!
        !! Output is dimension-aware and adjusts formatting accordingly.

        class(decomp_info) :: info          !! Object containing the decomposition info
        type(MPI_Comm), intent(in) :: comm  !! MPI communicator (usually MPI_COMM_WORLD)

        ! Local variables
        integer :: rank, size, ierr
        integer, allocatable :: all_info(:,:)
        integer :: i, ndim
        integer :: send_buf(6)

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, size, ierr)

        ! Determine number of active dimensions
        ndim = count(info%dims >= 1)

        ! Prepare buffer for all 3 dimensions (unused entries = 0)
        send_buf = 0
        send_buf(1) = info%ilow
        send_buf(2) = info%ihigh
        send_buf(3) = info%jlow
        send_buf(4) = info%jhigh
        send_buf(5) = info%klow
        send_buf(6) = info%khigh

        ! Allocate full buffer (6 entries per rank)
        allocate(all_info(6, size))

        ! Gather data to rank 0
        call MPI_Gather(send_buf, 6, MPI_INTEGER, all_info, 6, MPI_INTEGER, 0, comm, ierr)

        if (rank == 0) then
            select case (ndim)
              case (1)
                print *, 'Domain decomposition summary (1D block layout):'
                print *, '-------------------------------------------------------------'
                print *, ' Rank | i-range        '
                print *, '-------------------------------------------------------------'
                do i = 1, size
                    write(*,'(I5, 2X, "[", I6, ":", I6, "]")') i-1, all_info(1,i), all_info(2,i)
                end do
                print *, '-------------------------------------------------------------'
              case (2)
                print *, 'Domain decomposition summary (2D block layout):'
                print *, '-------------------------------------------------------------'
                print *, ' Rank | i-range        | j-range        '
                print *, '-------------------------------------------------------------'
                do i = 1, size
                    write(*,'(I5, 2X, "[", I6, ":", I6, "]", 3X, "[", I6, ":", I6, "]")') i-1, &
                        all_info(1,i), all_info(2,i), all_info(3,i), all_info(4,i)
                end do
                print *, '-------------------------------------------------------------'
              case (3)
                print *, 'Domain decomposition summary (3D block layout):'
                print *, '--------------------------------------------------------------------------'
                print *, ' Rank | i-range        | j-range        | k-range        '
                print *, '--------------------------------------------------------------------------'
                do i = 1, size
                    write(*,'(I5, 2X, "[", I6, ":", I6, "]", 3X, "[", I6, ":", I6, "]", 3X, "[", I6, ":", I6, "]")') i-1, &
                        all_info(1,i), all_info(2,i), all_info(3,i), all_info(4,i), all_info(5,i), all_info(6,i)
                end do
                print *, '--------------------------------------------------------------------------'
            end select
        end if

        deallocate(all_info)
    end subroutine print_decomposition_summary


    module subroutine print_cartesian_rank_layout(info)
        !! Prints a visual representation of the Cartesian rank layout.
        !! Supports 1D, 2D, or 3D decompositions.
        !!
        !! For 1D: prints a simple line of ranks.
        !! For 2D: prints a chessboard-style grid (existing logic).
        !! For 3D: prints a 2D layout for each k-layer (z-slice).

        class(decomp_info), intent(in) :: info  !! Object containing the decomposition info

        integer :: rank, size, ierr, ndim
        integer, allocatable :: coords(:)
        integer, allocatable :: local_info(:)
        integer, allocatable :: recvbuf(:)
        integer, allocatable :: all_ranks(:)
        integer, allocatable :: all_coords(:,:)
        integer :: i, j, k, idx, d, slice

        call MPI_Comm_rank(info%comm_cart, rank, ierr)
        call MPI_Comm_size(info%comm_cart, size, ierr)

        ! Determine number of dimensions
        ndim = count(info%dims >= 1)
        allocate(coords(ndim))
        allocate(local_info(ndim + 1))

        ! Get Cartesian coordinates of current rank
        call MPI_Cart_coords(info%comm_cart, rank, ndim, coords, ierr)

        ! Prepare send buffer: [rank, coord_1, coord_2, ..., coord_n]
        local_info(1) = rank
        local_info(2:) = coords

        if (rank == 0) allocate(recvbuf(size * (ndim + 1)))

        ! Gather all ranks + coordinates to root
        call MPI_Gather(local_info, ndim + 1, MPI_INTEGER, recvbuf, ndim + 1, MPI_INTEGER, 0, info%comm_cart, ierr)

        if (rank == 0) then
            allocate(all_ranks(size))
            allocate(all_coords(ndim, size))

            ! Unpack gathered info
            do i = 1, size
                idx = (ndim + 1) * (i - 1)
                all_ranks(i) = recvbuf(idx + 1)
                do d = 1, ndim
                    all_coords(d, i) = recvbuf(idx + 1 + d)
                end do
            end do

            select case (ndim)
              case (1)
                print *, "MPI Cartesian 1D layout:"
                do i = 0, info%dims(1) - 1
                    do idx = 1, size
                        if (all_coords(1, idx) == i) then
                            write(*,'(A,I2)', advance='no') " [", all_ranks(idx)
                            write(*,'(A)', advance='no') "]"
                            exit
                        end if
                    end do
                end do
                print *, ""

              case (2)
                print *, "MPI Cartesian 2D layout (dims = ", info%dims(1), ",", info%dims(2), "):"
                print *, "Coordinates (x,y), y=0 at bottom"
                call print_grid_2d(info%dims(1), info%dims(2), size, all_ranks, all_coords)

              case (3)
                print *, "MPI Cartesian 3D layout (dims = ", info%dims(1), ",", info%dims(2), ",", info%dims(3), "):"
                do slice = info%dims(3) - 1, 0, -1
                    print *, "Z-slice k =", slice
                    call print_grid_3d_slice(slice, info%dims(1), info%dims(2), size, all_ranks, all_coords)
                end do

            end select
        end if

    contains

        subroutine print_grid_2d(nx, ny, size, ranks, coords)
            integer, intent(in) :: nx, ny, size
            integer, intent(in) :: ranks(:)
            integer, intent(in) :: coords(:,:)
            integer :: i, j, idx, found_rank

            write(*,'(A)', advance='no') "  +"
            do i = 1, nx
                write(*,'(A)', advance='no') "---+"
            end do
            print *

            do j = ny - 1, 0, -1
                write(*,'(A)', advance='no') trim(itoa(j)) // " |"
                do i = 0, nx - 1
                    found_rank = -1
                    do idx = 1, size
                        if (coords(1, idx) == i .and. coords(2, idx) == j) then
                            found_rank = ranks(idx)
                            exit
                        end if
                    end do
                    if (found_rank >= 0) then
                        write(*,'(A,I2,A)', advance='no') " ", found_rank, " |"
                    else
                        write(*,'(A)', advance='no') "  - |"
                    end if
                end do
                print *
                write(*,'(A)', advance='no') "  +"
                do i = 1, nx
                    write(*,'(A)', advance='no') "---+"
                end do
                print *
            end do

            write(*,'(A)', advance='no') "    "
            do i = 0, nx - 1
                write(*,'(A,I2,A)', advance='no') " ", i, "  "
            end do
            print *
        end subroutine print_grid_2d

        subroutine print_grid_3d_slice(k, nx, ny, size, ranks, coords)
            integer, intent(in) :: k, nx, ny, size
            integer, intent(in) :: ranks(:)
            integer, intent(in) :: coords(:,:)
            integer :: i, j, idx, found_rank

            write(*,'(A)', advance='no') "  +"
            do i = 1, nx
                write(*,'(A)', advance='no') "---+"
            end do
            print *

            do j = ny - 1, 0, -1
                write(*,'(A)', advance='no') trim(itoa(j)) // " |"
                do i = 0, nx - 1
                    found_rank = -1
                    do idx = 1, size
                        if (coords(1, idx) == i .and. coords(2, idx) == j .and. coords(3, idx) == k) then
                            found_rank = ranks(idx)
                            exit
                        end if
                    end do
                    if (found_rank >= 0) then
                        write(*,'(A,I2,A)', advance='no') " ", found_rank, " |"
                    else
                        write(*,'(A)', advance='no') "  - |"
                    end if
                end do
                print *
                write(*,'(A)', advance='no') "  +"
                do i = 1, nx
                    write(*,'(A)', advance='no') "---+"
                end do
                print *
            end do

            write(*,'(A)', advance='no') "    "
            do i = 0, nx - 1
                write(*,'(A,I2,A)', advance='no') " ", i, "  "
            end do
            print *
        end subroutine print_grid_3d_slice

        pure function itoa(i) result(str)
            integer, intent(in) :: i
            character(len=12) :: str
            write(str, '(I0)') i
        end function itoa

    end subroutine print_cartesian_rank_layout

end submodule
