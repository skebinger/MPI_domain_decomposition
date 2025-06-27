submodule(MOD_MPI_decomposition) SMOD_init_decomposition
    !! Deals with the initialization of the domain decomposition
    implicit none
contains

    module subroutine initialize_decomposition(info, m_xi, m_eta, comm)
        !! Initializes the domain decomposition for a 2D grid using MPI.
        !!
        !! The decomposition is performed only by the master rank (rank 0), which computes
        !! a 2D processor grid decomposition of the global domain into chunks.
        !! Each chunk corresponds to a rectangular subdomain of the global mesh.
        !!
        !! Rank 0 then sends the computed local index ranges (bounds) to all other ranks
        !! using point-to-point MPI communication. Each rank receives its assigned local
        !! ranges for both i (xi) and j (eta) directions.
        !!
        !! This approach avoids all ranks calling MPI_Dims_create, thus preventing
        !! potential segmentation faults and ensuring consistent decomposition.
        !!
        !! The subroutine uses a helper `compute_range` to divide each dimension evenly
        !! among the processors assigned along that dimension, distributing any remainder
        !! cells to the lowest ranks.
        !!
        !! Input:
        !! - m_xi: Total number of cells in the i-direction (xi-axis).
        !! - m_eta: Total number of cells in the j-direction (eta-axis).
        !! - comm: MPI communicator over which ranks are defined.
        !!
        !! Output (module private variables set):
        !! - ilow, ihigh: Local i-direction bounds (inclusive).
        !! - jlow, jhigh: Local j-direction bounds (inclusive).
        class(decomp_info), intent(out) :: info !! Object containing the decomposition info
        integer, intent(in) :: m_xi             !! Total number of cells in the i-direction (xi-axis)
        integer, intent(in) :: m_eta            !! Total number of cells in the j-direction (eta-axis)
        integer, intent(in) :: comm             !! MPI communicator (normally 'MPI_COMM_WORLD')


        ! Internal variables
        integer :: rank, size, ierr
        integer :: dims_send(2),dims_rec(2)             ! Processor grid in i and j directions; buffers for MPI communication
        integer :: i_procs, j_procs    ! Number of processors along i and j
        integer, allocatable :: ilows(:), ihighs(:), jlows(:), jhighs(:)
        integer :: i_rank, j_rank
        integer :: i
        integer :: sendbuf(4), recvbuf(4)
        integer :: il, ih, jl, jh
        integer :: r                    ! Flattened rank index

        character(len=40) :: rank_str

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, size, ierr)

        if (rank == 0) then
            ! Initialize dims to zero for MPI_Dims_create
            info%dims = 0

            ! Create a 2D processor grid decomposition
            ! let MPI decide the layout
            call MPI_Dims_create(size, 2, info%dims, ierr)
            i_procs = info%dims(1)
            j_procs = info%dims(2)

            ! Allocate arrays to store decomposition info for all ranks
            allocate(ilows(size), ihighs(size), jlows(size), jhighs(size))

            ! Compute local index bounds for each rank in the 2D grid
            do i_rank = 0, i_procs - 1
                do j_rank = 0, j_procs - 1
                    ! Flatten 2D grid coords into linear rank
                    r = i_rank * j_procs + j_rank

                    ! Compute i-range for processor in i direction
                    call compute_range(1, m_xi, i_procs, i_rank, il, ih)

                    ! Compute j-range for processor in j direction
                    call compute_range(1, m_eta, j_procs, j_rank, jl, jh)

                    ilows(r+1) = il
                    ihighs(r+1) = ih
                    jlows(r+1) = jl
                    jhighs(r+1) = jh
                end do
            end do
        end if

        ! communicate the dimensions of the grid to all ranks
        dims_send=info%dims
        if(rank==0)then
            ! Rank 0 sends the grid dimensions to all other ranks
            do i=1,size-1
                call MPI_Send(dims_send, 2, MPI_INTEGER, i, 0, comm, ierr)
            end do
        else
            ! Other ranks recieve the dimensions
            call MPI_Recv(dims_rec, 2, MPI_INTEGER, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
            info%dims=dims_rec
        end if

        ! Each rank prepares to receive or send its decomposition info
        if (rank == 0) then
            ! Rank 0 sends each rank its local bounds (except itself)
            do i = 0, size - 1
                sendbuf(1) = ilows(i+1)
                sendbuf(2) = ihighs(i+1)
                sendbuf(3) = jlows(i+1)
                sendbuf(4) = jhighs(i+1)
                if (i == 0) then
                    ! Rank 0 keeps its own bounds locally
                    recvbuf = sendbuf
                else
                    ! Send to other ranks
                    call MPI_Send(sendbuf, 4, MPI_INTEGER, i, 0, comm, ierr)
                end if
            end do
        else
            ! Other ranks receive their local decomposition bounds from rank 0
            call MPI_Recv(recvbuf, 4, MPI_INTEGER, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
        end if

        ! Assign local decomposition bounds from received buffer
        info%ilow = recvbuf(1)
        info%ihigh = recvbuf(2)
        info%jlow = recvbuf(3)
        info%jhigh = recvbuf(4)

        ! Cleanup allocated arrays on rank 0
        if (rank == 0) then
            deallocate(ilows, ihighs, jlows, jhighs)
        end if

        write(rank_str, '(I0)') rank
        ! Write decomposition to disk
        open(unit=1,file=adjustl('output/domain_rank_' // trim(rank_str)//'.dat'),status='unknown',form='formatted')
        write(1,*) rank
        write(1,*) info%ilow
        write(1,*) info%ihigh
        write(1,*) info%jlow
        write(1,*) info%jhigh
        close(unit=1)

    end subroutine initialize_decomposition

    module subroutine setup_cartesian_topology(info,comm_in)
        !! Create a 2D Cartesian communicator and get neighbor ranks.
        !!
        !! The number of processes in each dimension is taken from previous
        !! calculation executed in initialize_decomposition. Sets the global
        !! module variable 'comm_cart'.
        class(decomp_info) :: info !! Object containing the decomposition info
        integer, intent(in) :: comm_in !! Input communicator (usually MPI_COMM_WORLD)

        integer :: ierr, rank, size
        integer :: coords(2)
        logical :: periods(2)

        call MPI_Comm_rank(comm_in, rank, ierr)
        call MPI_Comm_size(comm_in, size, ierr)

        periods = (/ .false., .false. /)  ! Non-periodic boundaries for now

        ! Create Cartesian communicator
        call MPI_Cart_create(comm_in, 2, info%dims, periods, .true., info%comm_cart, ierr)

        ! Get coordinates of this rank
        call MPI_Cart_coords(info%comm_cart, rank, 2, coords, ierr)

        ! Find left/right neighbors (dimension 1: x)
        call MPI_Cart_shift(info%comm_cart, 0, 1, info%nbr_left, info%nbr_right, ierr)

        ! Find bottom/top neighbors (dimension 2: y)
        call MPI_Cart_shift(info%comm_cart, 1, 1, info%nbr_bottom, info%nbr_top, ierr)

    end subroutine setup_cartesian_topology

    subroutine compute_range(start_idx, n, nprocs, rank, lo, hi)
        !! Computes a 1D index range for a given process.
        !!
        !! Performs block decomposition in one dimension, including handling of
        !! remainders if the number of elements is not divisible by the number of ranks.
        integer :: start_idx, n, nprocs, rank
        integer :: lo, hi
        integer :: base, extra, offset

        base  = n / nprocs
        extra = mod(n, nprocs)

        offset = rank * base + min(rank, extra)
        lo = start_idx + offset
        hi = lo + base - 1
        if (rank < extra) hi = hi + 1
    end subroutine compute_range

end submodule
