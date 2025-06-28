submodule(MOD_MPI_decomposition) SMOD_MPI_halo_exchange
    !! Submodule handling the halo exchange at processor interfaces.
    implicit none
contains

    module subroutine exchange_halos(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
        !! Exchanges halo layers for a multi-variable 2D structured field.
        !!
        !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
        !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.
        class(decomp_info), intent(in) :: info
        !! Object containing the decomposition info
        double precision, allocatable :: dat2D(:,:,:)
        !! Multi-variable 2D array: dat2D(m_var, i, j)
        integer, intent(in) :: m_var        !! Number of variables per grid point
        integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
        integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: comm         !! MPI communicator

        ! local variables
        integer :: i,j, layer, ierr, rank
        integer :: recv_type, send_type
        integer :: sizes(3), subsizes(3), starts(3)
        integer :: nreq
        integer, allocatable :: reqs(:)
        integer :: stat(MPI_STATUS_SIZE)
        integer :: maxreq

        ! Buffer arrays for send and receive — dimensioned for m_var * num_j_points
        double precision, allocatable :: sendbuf_left(:), recvbuf_right(:)
        double precision, allocatable :: sendbuf_right(:), recvbuf_left(:)
        integer :: nj, count,var

        ! local index range of rank
        integer :: ilo,ihi,jlo,jhi

        integer :: offset_i_right_recv, offset_i_left_recv, offset_i_send_left, offset_i_send_right

        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        call info%get_local_block_bounds(ilo,ihi,jlo,jhi)

        !---------------------------------------------------------------
        ! Create derived datatypes for non-contiguous face transfer
        !---------------------------------------------------------------
        ! Overall size of the full array (with ghosts)
        sizes     = [m_var, (1+abs(ihi-ilo)) + 2*num_ghost, (1+abs(jhi-jlo)) + 2*num_ghost]

        ! Size of the subarray to send/recv: (m_var, 1, m_eta)
        subsizes  = [m_var, 1, (1+abs(jhi-jlo))]

        ! These will be set per call — but to define the type we assume j=0 (placeholder)
        starts    = [0, 0, 0]  ! MPI uses 0-based indexing, unlike Fortran

        ! Maximum of 8 requests per ghost layer
        maxreq = 8 * num_ghost
        allocate(reqs(maxreq))
        nreq = 0

        !print *, left,right,top,bottom,rank
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !stop

        ! === LEFT / RIGHT exchange ===
        do layer = 1, num_ghost
            nj = jhi - jlo + 1
            count = m_var * nj

            allocate(sendbuf_left(count))
            allocate(recvbuf_right(count))
            allocate(sendbuf_right(count))
            allocate(recvbuf_left(count))

            ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
            do var = 1, m_var
                do j = jlo, jhi
                    sendbuf_left((var-1)*nj + (j - jlo + 1)) = dat2D(var, ilo + layer - 1, j)
                end do
            end do

            call MPI_Sendrecv(sendbuf_left, count, MPI_DOUBLE_PRECISION, left, 100 + layer, &
                recvbuf_right, count, MPI_DOUBLE_PRECISION, right, 100 + layer, &
                comm, MPI_STATUS_IGNORE, ierr)

            if (right /= MPI_PROC_NULL) then
                do var = 1, m_var
                    do j = jlo, jhi
                        dat2D(var, ihi + layer, j) = recvbuf_right((var-1)*nj + (j - jlo + 1))
                    end do
                end do
            end if

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            do var = 1, m_var
                do j = jlo, jhi
                    sendbuf_right((var-1)*nj + (j - jlo + 1)) = dat2D(var, ihi - layer + 1, j)
                end do
            end do

            call MPI_Sendrecv(sendbuf_right, count, MPI_DOUBLE_PRECISION, right, 200 + layer, &
                recvbuf_left, count, MPI_DOUBLE_PRECISION, left, 200 + layer, &
                comm, MPI_STATUS_IGNORE, ierr)

            if (left /= MPI_PROC_NULL) then
                do var = 1, m_var
                    do j = jlo, jhi
                        dat2D(var, ilo - layer, j) = recvbuf_left((var-1)*nj + (j - jlo + 1))
                    end do
                end do
            end if

            deallocate(sendbuf_left)
            deallocate(recvbuf_right)
            deallocate(sendbuf_right)
            deallocate(recvbuf_left)
        end do

        ! === BOTTOM / TOP exchange ===
        do layer = 1, num_ghost
            ! Receive top ghost
            call MPI_Irecv(dat2D(:, ilo:ihi, jhi+layer), m_var*(1+abs(ihi-ilo)), MPI_DOUBLE_PRECISION, top, 300+layer, comm, reqs(nreq+1), ierr)
            nreq = nreq + 1

            ! Receive bottom ghost
            call MPI_Irecv(dat2D(:, ilo:ihi, jlo-layer), m_var*(1+abs(ihi-ilo)), MPI_DOUBLE_PRECISION, bottom, 400+layer, comm, reqs(nreq+1), ierr)
            nreq = nreq + 1

            ! Send bottom boundary
            call MPI_Isend(dat2D(:, ilo:ihi, jlo+layer-1), m_var*(1+abs(ihi-ilo)), MPI_DOUBLE_PRECISION, bottom, 300+layer, comm, reqs(nreq+1), ierr)
            nreq = nreq + 1

            ! Send top boundary
            call MPI_Isend(dat2D(:, ilo:ihi, jhi+1-layer), m_var*(1+abs(ihi-ilo)), MPI_DOUBLE_PRECISION, top, 400+layer, comm, reqs(nreq+1), ierr)
            nreq = nreq + 1
        end do

        ! Wait for all communications to complete
        call MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE, ierr)

        deallocate(reqs)

    end subroutine exchange_halos
end submodule
