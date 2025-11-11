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

submodule(MOD_MPI_decomposition) SMOD_init_decomposition
    !! Deals with the initialization of the domain decomposition
    implicit none
contains

    module subroutine initialize_decomposition(info, m_xi, m_eta, m_tau, comm)
        !! Initializes the domain decomposition for a N-dimensional grid using MPI.
        !!
        !! The decomposition supports 1D, 2D, and 3D domain partitioning depending on the
        !! non-zero inputs among m_xi, m_eta, and m_tau. A zero value for a dimension
        !! disables decomposition in that direction.
        !!
        !! The decomposition is performed only by the master rank (rank 0), which computes
        !! an N-D processor grid decomposition of the global domain into chunks.
        !! Each chunk corresponds to a rectangular subdomain of the global mesh.
        !!
        !! Rank 0 then sends the computed local index ranges (bounds) to all other ranks
        !! using point-to-point MPI communication. Each rank receives its assigned local
        !! ranges for all active directions.
        !!
        !! This approach avoids all ranks calling MPI_Dims_create, thus preventing
        !! potential segmentation faults and ensuring consistent decomposition.
        !!
        !! The subroutine uses a helper `compute_range` to divide each dimension evenly
        !! among the processors assigned along that dimension, distributing any remainder
        !! cells to the lowest ranks.
        !!
        !! Input:
        !! - m_xi: Total number of cells in the i-direction (xi-axis). 0 = inactive.
        !! - m_eta: Total number of cells in the j-direction (eta-axis). 0 = inactive.
        !! - m_tau: Total number of cells in the k-direction (tau-axis). 0 = inactive.
        !! - comm: MPI communicator over which ranks are defined.
        !!
        !! Output (module private variables set):
        !! - ilow, ihigh: Local i-direction bounds (inclusive).
        !! - jlow, jhigh: Local j-direction bounds (inclusive).
        !! - klow, khigh: Local k-direction bounds (inclusive), if active.
        class(decomp_info), intent(out) :: info        !! Object containing the decomposition info
        integer, intent(in) :: m_xi                    !! Total number of cells in the i-direction (xi-axis)
        integer, intent(in) :: m_eta                   !! Total number of cells in the j-direction (eta-axis)
        integer, intent(in) :: m_tau                   !! Total number of cells in the k-direction (tau-axis)
        type(MPI_Comm), intent(in) :: comm             !! MPI communicator (normally 'MPI_COMM_WORLD')

        ! Internal variables
        integer :: rank, size, ierr
        integer :: ndim                                ! Number of active dimensions
        integer :: i_procs, j_procs, k_procs           ! Number of processors in each direction
        integer :: i_rank, j_rank, k_rank              ! Processor indices in each direction
        integer :: il, ih, jl, jh, kl, kh              ! Computed index bounds for each direction
        integer :: r                                   ! Flattened rank index
        integer :: i                                   ! Loop index
        integer :: sendbuf(6), recvbuf(6)              ! Buffers for sending/receiving bounds
        integer :: dims_send(3), dims_rec(3)           ! Processor grid dimensions
        integer :: global_cells(3)                     ! Global grid sizes
        integer, allocatable :: bounds_lows(:,:), bounds_highs(:,:) ! Full decomposition map

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, size, ierr)

        ! initialize all index bounds to 1, so later only the index ranges of the active spatial
        ! dimension is overwritten
        info%ilow   = 1
        info%ihigh  = 1
        info%jlow   = 1
        info%jhigh  = 1
        info%klow   = 1
        info%khigh  = 1

        ! Determine active dimensions and count them
        ndim = 0
        if (m_xi >= 1) ndim = ndim + 1
        if (m_eta >= 1) ndim = ndim + 1
        if (m_tau >= 1) ndim = ndim + 1

        ! Initialize global grid size array
        global_cells = (/ m_xi, m_eta, m_tau /)

        if (rank == 0) then
            ! Initialize processor grid dimensions to zero
            info%dims = 0

            ! Compute processor grid layout using MPI_Dims_create for active dimensions
            call MPI_Dims_create(size, ndim, info%dims, ierr)

            ! Set number of processors in each direction (1 for inactive directions)
            i_procs = merge(info%dims(1), 1, m_xi >= 1)
            j_procs = merge(info%dims(2), 1, m_eta >= 1)
            k_procs = merge(info%dims(3), 1, m_tau >= 1)

            ! Allocate arrays to hold bounds for all ranks
            allocate(bounds_lows(3, size), bounds_highs(3, size))

            ! Loop over all processor grid coordinates and compute bounds
            do i_rank = 0, i_procs - 1
                do j_rank = 0, j_procs - 1
                    do k_rank = 0, k_procs - 1
                        r = i_rank * (j_procs * k_procs) + j_rank * k_procs + k_rank

                        if (m_xi >= 1) then
                            call compute_range(1, m_xi, i_procs, i_rank, il, ih)
                            bounds_lows(1, r+1)  = il
                            bounds_highs(1, r+1) = ih
                        end if
                        if (m_eta >= 1) then
                            call compute_range(1, m_eta, j_procs, j_rank, jl, jh)
                            bounds_lows(2, r+1)  = jl
                            bounds_highs(2, r+1) = jh
                        end if
                        if (m_tau >= 1) then
                            call compute_range(1, m_tau, k_procs, k_rank, kl, kh)
                            bounds_lows(3, r+1)  = kl
                            bounds_highs(3, r+1) = kh
                        end if
                    end do
                end do
            end do
        end if

        ! Share processor grid dimensions with all ranks
        dims_send = info%dims
        if (rank == 0) then
            do i = 1, size - 1
                call MPI_Send(dims_send, 3, MPI_INTEGER, i, 0, comm, ierr)
            end do
        else
            call MPI_Recv(dims_rec, 3, MPI_INTEGER, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
            info%dims = dims_rec
        end if

        ! Derive processor grid dimensions per rank
        i_procs = merge(info%dims(1), 1, m_xi >= 1)
        j_procs = merge(info%dims(2), 1, m_eta >= 1)
        k_procs = merge(info%dims(3), 1, m_tau >= 1)

        ! Compute rank's coordinates in the processor grid
        k_rank = mod(rank, k_procs)
        j_rank = mod(rank / k_procs, j_procs)
        i_rank = rank / (j_procs * k_procs)
        r = i_rank * (j_procs * k_procs) + j_rank * k_procs + k_rank

        if (rank == 0) then
            ! Rank 0 fills its own bounds directly
            sendbuf = 0
            if (m_xi >= 1) then
                sendbuf(1) = bounds_lows(1, r+1)
                sendbuf(2) = bounds_highs(1, r+1)
            end if
            if (m_eta >= 1) then
                sendbuf(3) = bounds_lows(2, r+1)
                sendbuf(4) = bounds_highs(2, r+1)
            end if
            if (m_tau >= 1) then
                sendbuf(5) = bounds_lows(3, r+1)
                sendbuf(6) = bounds_highs(3, r+1)
            end if

            recvbuf = sendbuf

            ! Send bounds to all other ranks
            do i = 1, size - 1
                sendbuf = 0
                if (m_xi >= 1) then
                    sendbuf(1) = bounds_lows(1, i+1)
                    sendbuf(2) = bounds_highs(1, i+1)
                end if
                if (m_eta >= 1) then
                    sendbuf(3) = bounds_lows(2, i+1)
                    sendbuf(4) = bounds_highs(2, i+1)
                end if
                if (m_tau >= 1) then
                    sendbuf(5) = bounds_lows(3, i+1)
                    sendbuf(6) = bounds_highs(3, i+1)
                end if

                call MPI_Send(sendbuf, 6, MPI_INTEGER, i, 0, comm, ierr)
            end do
        else
            ! Non-root ranks receive their assigned bounds
            call MPI_Recv(recvbuf, 6, MPI_INTEGER, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
        end if

        ! Assign local bounds from buffer
        if (m_xi >= 1) then
            info%ilow = recvbuf(1)
            info%ihigh = recvbuf(2)
        end if
        if (m_eta >= 1) then
            info%jlow = recvbuf(3)
            info%jhigh = recvbuf(4)
        end if
        if (m_tau >= 1) then
            info%klow = recvbuf(5)
            info%khigh = recvbuf(6)
        end if

        ! Deallocate only on root
        if (rank == 0) then
            deallocate(bounds_lows, bounds_highs)
        end if

    end subroutine initialize_decomposition


    module subroutine setup_cartesian_topology(info, comm_in)
        !! Create a 1D, 2D, or 3D Cartesian communicator and get neighbor ranks.
        !!
        !! The number of processes in each dimension is taken from previous
        !! calculation executed in initialize_decomposition. Sets the global
        !! module variable 'comm_cart' and neighbor ranks.
        !!
        !! Neighbor ranks include:
        !! - nbr_left/nbr_right (i-direction)
        !! - nbr_bottom/nbr_top (j-direction)
        !! - nbr_back/nbr_front (k-direction), if active

        class(decomp_info) :: info               !! Object containing the decomposition info
        type(MPI_Comm), intent(in) :: comm_in    !! Input communicator (usually MPI_COMM_WORLD)

        integer :: ierr, rank, size, ndim
        integer, allocatable :: coords(:)
        logical, allocatable :: periods(:)

        call MPI_Comm_rank(comm_in, rank, ierr)
        call MPI_Comm_size(comm_in, size, ierr)

        ! Determine active dimensions
        ndim = count(info%dims >= 1)

        ! Allocate arrays dynamically for ndim dimensions
        allocate(coords(ndim))
        allocate(periods(ndim))
        periods = .false.  ! Use non-periodic boundaries for now

        ! Create Cartesian communicator with correct dimensionality
        call MPI_Cart_create(comm_in, ndim, info%dims(1:ndim), periods, .true., info%comm_cart, ierr)

        ! Get coordinates of this rank in the Cartesian grid
        call MPI_Cart_coords(info%comm_cart, rank, ndim, coords, ierr)

        ! Find neighbors in each active direction
        if (ndim >= 1) then
            call MPI_Cart_shift(info%comm_cart, 0, 1, info%nbr_left, info%nbr_right, ierr)
        end if
        if (ndim >= 2) then
            call MPI_Cart_shift(info%comm_cart, 1, 1, info%nbr_bottom, info%nbr_top, ierr)
        end if
        if (ndim == 3) then
            call MPI_Cart_shift(info%comm_cart, 2, 1, info%nbr_back, info%nbr_front, ierr)
        end if

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

    module subroutine write_decom_to_disk(info,comm)
        class(decomp_info), intent(in) :: info
        type(MPI_Comm), intent(in) :: comm

        character(len=40) :: rank_str
        integer :: rank,size,ierr
        integer :: exitstat,cmdstat

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, size, ierr)

        write(rank_str, '(I0)') rank
        ! Write decomposition to disk
        call execute_command_line("mkdir -p output",.TRUE.,exitstat,cmdstat)
        open(unit=1,file=adjustl('output/domain_rank_' // trim(rank_str)//'.dat'),status='unknown',form='formatted')
        write(1,*) rank
        write(1,*) info%ilow
        write(1,*) info%ihigh
        write(1,*) info%jlow
        write(1,*) info%jhigh
        write(1,*) info%klow
        write(1,*) info%khigh
        close(unit=1)
    end subroutine

end submodule
