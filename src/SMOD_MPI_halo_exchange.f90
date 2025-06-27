submodule(MOD_MPI_decomposition) SMOD_MPI_halo_exchange
    !! Submodule handling the halo exchange at processor interfaces.
    implicit none
contains

    module subroutine exchange_halos(dat2D, m_var, m_xi, m_eta, num_ghost, left, right, bottom, top, comm)
        !! Exchanges halo layers for a multi-variable 2D structured field.
        !!
        !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
        !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.

        double precision, allocatable :: dat2D(:,:,:)
        !! Multi-variable 2D array: dat2D(m_var, i, j)
        integer, intent(in) :: m_var        !! Number of variables per grid point
        integer, intent(in) :: m_xi         !! Number of *interior* i-cells (excluding ghosts)
        integer, intent(in) :: m_eta        !! Number of *interior* j-cells (excluding ghosts)
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
        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: nj, count,var

        ! local index range of rank
        integer :: ilo,ihi,jlo,jhi

        integer :: offset_i_right_recv, offset_i_left_recv, offset_i_send_left, offset_i_send_right


        logical :: buffered_eta_exchange = .false.
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        call get_local_block_bounds(ilo,ihi,jlo,jhi)

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

            if(buffered_eta_exchange)then
                nj = jhi - jlo + 1
                count = m_var * nj

                allocate(sendbuf(count))
                allocate(recvbuf(count))

                ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
                do var = 1, m_var
                    do j = jlo, jhi
                        sendbuf((var-1)*nj + (j - jlo + 1)) = dat2D(var, ilo + layer - 1, j)
                    end do
                end do

                call MPI_Sendrecv(sendbuf, count, MPI_DOUBLE_PRECISION, left, 100 + layer, &
                    recvbuf, count, MPI_DOUBLE_PRECISION, right, 100 + layer, &
                    comm, MPI_STATUS_IGNORE, ierr)

                if (right /= MPI_PROC_NULL) then
                    do var = 1, m_var
                        do j = jlo, jhi
                            dat2D(var, ihi + layer, j) = recvbuf((var-1)*nj + (j - jlo + 1))
                        end do
                    end do
                end if

                ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
                do var = 1, m_var
                    do j = jlo, jhi
                        sendbuf((var-1)*nj + (j - jlo + 1)) = dat2D(var, ihi - layer + 1, j)
                    end do
                end do

                call MPI_Sendrecv(sendbuf, count, MPI_DOUBLE_PRECISION, right, 200 + layer, &
                    recvbuf, count, MPI_DOUBLE_PRECISION, left, 200 + layer, &
                    comm, MPI_STATUS_IGNORE, ierr)

                if (left /= MPI_PROC_NULL) then
                    do var = 1, m_var
                        do j = jlo, jhi
                            dat2D(var, ilo - layer, j) = recvbuf((var-1)*nj + (j - jlo + 1))
                        end do
                    end do
                end if

                deallocate(sendbuf)
                deallocate(recvbuf)
            else
                !EXPERIMENTAL NOT YET WORKING
                !=== WARNING: Dangerous code segment! Manually prescribes array bounds for MPI array exchange.
                !Huge potential to fuck up and send/recieve out of bounds!! ====!

                !MORE EFFICIENT VERSION: workaround to use a contiguous segment of data for MPI send and recieve
                !achieved through a custom mpi type, cutting data into subarrays
                !short explanation (UPDATE!!!!!):
                !dat2D = pointer to start of data array in memory
                !1 = counter of how many element of type recv_type or send_type are communicated (here one, because that data is "wrapped" into the type)
                !right/left = neighbour rank to right and left
                !100+layer = unique tag
                !rest not important for here

                error stop "Currently experimental version of j-ghost cell exchange. Not working"

                ! Calculate zero-based offsets relative to dat2D's first element (1-num_ghost in Fortran corresponds to 0 in MPI)
                offset_i_right_recv = (ihi + layer) - (1 - num_ghost)
                offset_i_left_recv  = (ilo - layer) - (1 - num_ghost)
                offset_i_send_left  = (ilo + layer - 1) - (1 - num_ghost)
                offset_i_send_right = (ihi - layer + 1) - (1 - num_ghost)

                !ilo_g = lbound(dat2D, 2)
                !jlo_g = lbound(dat2D, 3)

                ! --- Setup receive datatype for column at right ghost (i = ihi + layer) ---
                starts = [0, offset_i_right_recv, jlo - (1 - num_ghost)]
                call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, ierr)
                call MPI_Type_commit(recv_type, ierr)

                call MPI_Irecv(dat2D, 1, recv_type, right, 100 + layer, comm, reqs(nreq + 1), ierr)
                nreq = nreq + 1

                call MPI_Type_free(recv_type, ierr)

                ! --- Setup receive datatype for column at left ghost (i = ilo - layer) ---
                starts = [0, offset_i_left_recv, jlo - (1 - num_ghost)]
                call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, ierr)
                call MPI_Type_commit(recv_type, ierr)

                call MPI_Irecv(dat2D, 1, recv_type, left, 200 + layer, comm, reqs(nreq + 1), ierr)
                nreq = nreq + 1

                call MPI_Type_free(recv_type, ierr)

                ! --- Setup send datatype for column at left boundary (i = ilo + layer - 1) ---
                starts = [0, offset_i_send_left, jlo - (1 - num_ghost)]
                call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_type, ierr)
                call MPI_Type_commit(send_type, ierr)

                call MPI_Isend(dat2D, 1, send_type, left, 100 + layer, comm, reqs(nreq + 1), ierr)
                nreq = nreq + 1

                call MPI_Type_free(send_type, ierr)

                ! --- Setup send datatype for column at right boundary (i = ihi - layer + 1) ---
                starts = [0, offset_i_send_right, jlo - (1 - num_ghost)]
                call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_type, ierr)
                call MPI_Type_commit(send_type, ierr)

                call MPI_Isend(dat2D, 1, send_type, right, 200 + layer, comm, reqs(nreq + 1), ierr)
                nreq = nreq + 1

                call MPI_Type_free(send_type, ierr)
            end if
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
