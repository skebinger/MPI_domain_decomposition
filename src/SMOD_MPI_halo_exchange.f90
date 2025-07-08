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
        type(MPI_Comm), intent(in) :: comm  !! MPI communicator

        ! local variables
        integer :: layer

        ! local index range of rank
        integer :: ilo,ihi,jlo,jhi

        !call MPI_Comm_rank(comm,rank,ierr)

        call info%get_local_block_bounds(ilo,ihi,jlo,jhi)

        ! === LEFT / RIGHT exchange ===
        do layer = 1, num_ghost

            ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
            call send_j_ghost(info, dat2D, m_var, ilo+layer-1, ihi+layer, left, right, 100+layer, comm)

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            call send_j_ghost(info, dat2D, m_var, ihi-layer+1, ilo-layer, right, left, 200+layer, comm)

        end do

        ! === BOTTOM / TOP exchange ===
        do layer = 1, num_ghost

            ! --- SEND BOTTOM boundary: j = jlo + layer - 1 ---
            call send_i_ghost(info, dat2D, m_var, jlo+layer-1, jhi+layer, bottom, top, 300+layer, comm)

            ! --- SEND TOP boundary: j = jhi - layer + 1 ---
            call send_i_ghost(info, dat2D, m_var, jhi-layer+1, jlo-layer, top, bottom, 400+layer, comm)
        end do

    end subroutine exchange_halos_2D


    module subroutine send_j_ghost(info, dat2D, m_var, isend, irecv, rank_target, rank_origin, tag, comm)
        !! Sends the left or right ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: isend !! i-index of the slice of which data is sent
        integer, intent(in) :: irecv !! i-index of the slice to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator

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
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_2d_jslice(dat2D, sendbuf, isend, m_var, info%jlow, info%jhigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_2d_jslice(dat2D, recvbuf, irecv, m_var, info%jlow, info%jhigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_i_ghost(info, dat2D, m_var, jsend, jrecv, rank_target, rank_origin, tag, comm)
        !! Sends the top or bottom ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: jsend !! j-index of the slice of which data is sent
        integer, intent(in) :: jrecv !! j-index of the slice to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator

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
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_2d_islice(dat2D, sendbuf, jsend, m_var, info%ilow, info%ihigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_2d_islice(dat2D, recvbuf, jrecv, m_var, info%ilow, info%ihigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine
end submodule
