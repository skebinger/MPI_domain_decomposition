submodule(MOD_MPI_decomposition) SMOD_MPI_halo_exchange
    !! Submodule handling the halo exchange at processor interfaces.
    implicit none
contains

    module subroutine exchange_halos_2D(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
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

        ! Buffer arrays for send and receive — dimensioned for m_var * num_j_points
        double precision, allocatable :: sendbuf_left(:), recvbuf_right(:)
        double precision, allocatable :: sendbuf_right(:), recvbuf_left(:)
        integer :: nj, count,var

        ! local index range of rank
        integer :: ilo,ihi,jlo,jhi

        integer :: offset_i_right_recv, offset_i_left_recv, offset_i_send_left, offset_i_send_right

        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        call info%get_local_block_bounds(ilo,ihi,jlo,jhi)

        !print *, left,right,top,bottom,rank
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! === LEFT / RIGHT exchange ===
        do layer = 1, num_ghost

            ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
            call send_left_right(info, dat2D, m_var, ilo+layer-1, ihi+layer, left, right, 100+layer, comm)

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            call send_left_right(info, dat2D, m_var, ihi-layer+1, ilo-layer, right, left, 200+layer, comm)

        end do

        ! === BOTTOM / TOP exchange ===
        do layer = 1, num_ghost

            ! --- SEND BOTTOM boundary: j = jlo + layer - 1 ---
            call send_top_bottom(info, dat2D, m_var, jlo+layer-1, jhi+layer, bottom, top, 300+layer, comm)

            ! --- SEND TOP boundary: j = jhi - layer + 1 ---
            call send_top_bottom(info, dat2D, m_var, jhi-layer+1, jlo-layer, top, bottom, 400+layer, comm)
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine exchange_halos_2D


    module subroutine send_left_right(info, dat2D, m_var, isend, irecv, rank_origin, rank_target, tag, comm)
        !! Sends the left or right ghost cells of a rank to its neighbour.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var
        integer, intent(in) :: isend,irecv
        integer, intent(in) :: rank_origin,rank_target
        integer, intent(in) :: tag
        integer, intent(in) :: comm

        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: rank, ierr

        ! Check data input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf,recvbuf,info%jlow,info%jhigh,m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        call pack_2d_jslice(dat2D, sendbuf, isend, m_var, info%jlow, info%jhigh, comm)

        call send_1D_array(sendbuf, recvbuf, rank_origin, rank_target, tag, MPI_DOUBLE_PRECISION, comm)

        ! Only unpack data back into 'dat2D' if 'rank_target' is a valid rank target (== internal boundary)
        if(rank_target /= MPI_PROC_NULL) call unpack_2d_jslice(dat2D, recvbuf, irecv, m_var, info%jlow, info%jhigh, comm)

        deallocate(sendbuf,recvbuf)
    end subroutine

    module subroutine send_top_bottom(info, dat2D, m_var, jsend, jrecv, rank_origin, rank_target, tag, comm)
        !! Sends the top or bottom ghost cells of a rank to its neighbour.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var
        integer, intent(in) :: jsend,jrecv
        integer, intent(in) :: rank_origin,rank_target
        integer, intent(in) :: tag
        integer, intent(in) :: comm

        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: rank, ierr

        ! Check data input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf, recvbuf, info%ilow, info%ihigh, m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        call pack_2d_islice(dat2D, sendbuf, jsend, m_var, info%ilow, info%ihigh, comm)

        call send_1D_array(sendbuf, recvbuf, rank_origin, rank_target, tag, MPI_DOUBLE_PRECISION, comm)

        ! Only unpack data back into 'dat2D' if 'rank_target' is a valid rank target (== internal boundary)
        if(rank_target /= MPI_PROC_NULL) call unpack_2d_islice(dat2D, recvbuf, jrecv, m_var, info%ilow, info%ihigh, comm)

        deallocate(sendbuf,recvbuf)
    end subroutine
end submodule
