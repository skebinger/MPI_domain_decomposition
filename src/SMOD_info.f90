submodule(MOD_MPI_decomposition) SMOD_info
    !! Submodule contains various subroutines to retrieve information about the 
    !! decomposition, such as cartesian rank layout or the index decomposition 
    !! summary.

    implicit none
contains

    module subroutine print_decomposition_summary(info,comm)
        !! Prints the domain decomposition summary for all MPI ranks.
        !!
        !! Gathers the local index ranges from all ranks and prints a formatted
        !! table on rank 0 showing the (i,j) bounds for each rank's block.
        class(decomp_info) :: info !! Object containing the decomposition info
        integer, intent(in) :: comm !! MPI communicator (usually MPI_COMM_WORLD)

        ! local variables
        integer :: rank, size, ierr
        integer, allocatable :: all_info(:,:)
        integer :: i
        integer :: send_buf(4)

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, size, ierr)

        allocate(all_info(4, size))

        send_buf(1) = info%ilow
        send_buf(2) = info%ihigh
        send_buf(3) = info%jlow
        send_buf(4) = info%jhigh

        ! Collect the decomposition info on the rank 0 (master rank)
        call MPI_Gather(send_buf, 4, MPI_INTEGER, all_info, 4, MPI_INTEGER, 0, comm, ierr)

        if (rank == 0) then
            print *, 'Domain decomposition summary (2D block layout):'
            print *, '-------------------------------------------------------------'
            print *, ' Rank | i-range        | j-range        '
            print *, '-------------------------------------------------------------'
            do i = 1, size
                write(*,'(I5, 2X, "[", I6, ":", I6, "]", 3X, "[", I6, ":", I6, "]")') i-1, &
                    all_info(1,i), all_info(2,i), all_info(3,i), all_info(4,i)
            end do
            print *, '-------------------------------------------------------------'
        end if

        deallocate(all_info)
    end subroutine print_decomposition_summary

    module subroutine print_cartesian_rank_layout(info)
        !! Prints a visual representation of the rank layout in the form of a chessboard.
        use mpi
        class(decomp_info), intent(in) :: info !! Object containing the decomposition info

        ! local variables
        integer :: rank, ierr, size
        integer :: coords(2)
        integer :: i, j, idx
        integer, dimension(3) :: local_info
        integer, allocatable :: recvbuf(:)
        integer, allocatable :: all_ranks(:)
        integer, allocatable :: all_coords(:,:)
        integer :: found_rank

        call MPI_Comm_rank(info%comm_cart, rank, ierr)
        call MPI_Comm_size(info%comm_cart, size, ierr)

        call MPI_Cart_coords(info%comm_cart, rank, 2, coords, ierr)

        local_info = [rank, coords(1), coords(2)]

        if (rank == 0) then
            allocate(recvbuf(3*size))
        end if

        if (rank == 0) then
            call MPI_Gather(local_info, 3, MPI_INTEGER, recvbuf, 3, MPI_INTEGER, 0, info%comm_cart, ierr)
        else
            call MPI_Gather(local_info, 3, MPI_INTEGER, MPI_BOTTOM, 3, MPI_INTEGER, 0, info%comm_cart, ierr)
        end if

        ! only print on rank 0 (master rank)
        if (rank == 0) then
            allocate(all_ranks(size))
            allocate(all_coords(2,size))

            do i = 1, size
                idx = 3*(i-1)
                all_ranks(i) = recvbuf(idx+1)
                all_coords(1,i) = recvbuf(idx+2)
                all_coords(2,i) = recvbuf(idx+3)
            end do

            print *, "MPI Cartesian grid layout (dims = ", info%dims(1), ",", info%dims(2), "):"
            print *, "Coordinates (x,y), y=0 at bottom"

            ! Print top border of the grid
            write(*,'(A)', advance='no') "  +"
            do i = 1, info%dims(1)
                write(*,'(A)', advance='no') "---+"
            end do
            print *

            ! Loop over rows from top (y=dims(2)-1) to bottom (y=0)
            do j = info%dims(2)-1, 0, -1
                ! Print rank line
                write(*,'(A)', advance='no') trim(adjustl(itoa(j)))// " |"
                do i = 0, info%dims(1)-1
                    found_rank = -1
                    do idx = 1, size
                        if (all_coords(1,idx) == i .and. all_coords(2,idx) == j) then
                            found_rank = all_ranks(idx)
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

                ! Print separator line
                write(*,'(A)', advance='no') "  +"
                do i = 1, info%dims(1)
                    write(*,'(A)', advance='no') "---+"
                end do
                print *
            end do

            ! Print x-axis labels below the grid
            write(*,'(A)', advance='no') "    "
            do i = 0, info%dims(1)-1
                write(*,'(A,I2,A)', advance='no') " ", i, "  "
            end do
            print *, ""

        end if

    contains

        pure function itoa(i) result(str)
            integer, intent(in) :: i
            character(len=12) :: str
            write(str,'(I0)') i
        end function itoa

    end subroutine print_cartesian_rank_layout

end submodule
