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

    module subroutine exchange_halos_1d(info, dat1D, m_var, num_ghost, left, right, comm)
        !! Exchanges halo layers for a multi-variable 1D structured field.
        !!
        !! Performs halo exchange on a field `dat2D(m_var, i)` with ghost cells,
        !! communicating ghost layers with direct MPI neighbors: left and right.
        class(decomp_info), intent(in) :: info
        !! Object containing the decomposition info
        double precision, allocatable :: dat1D(:,:)
        !! Multi-variable 1D array: dat1D(m_var, i)
        integer, intent(in) :: m_var        !! Number of variables per grid point
        integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
        integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
        type(MPI_Comm), intent(in) :: comm  !! MPI communicator

        ! local variables
        integer :: layer

        ! local index range of rank
        integer :: ilo,ihi

        !call MPI_Comm_rank(comm,rank,ierr)

        call info%get_local_block_bounds(ilo,ihi)

        ! === LEFT / RIGHT exchange ===
        do layer = 1, num_ghost

            ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
            call send_i_index_ghost(info, dat1D, m_var, ilo+layer-1, ihi+layer, left, right, 100+layer, comm)

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            call send_i_index_ghost(info, dat1D, m_var, ihi-layer+1, ilo-layer, right, left, 200+layer, comm)

        end do

    end subroutine

    module subroutine exchange_halos_2d(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
        !! Exchanges halo layers for a multi-variable 2D structured field.
        !!
        !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
        !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.
        class(decomp_info), intent(in) :: info
        !! Object containing the decomposition info
        double precision, allocatable :: dat2D(:,:,:)
        !! Multi-variable 2D array: `dat2D(m_var, i, j)`
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
            call send_i_index_ghost(info, dat2D, m_var, ilo+layer-1, ihi+layer, left, right, 100+layer, comm)

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            call send_i_index_ghost(info, dat2D, m_var, ihi-layer+1, ilo-layer, right, left, 200+layer, comm)

        end do

        ! === BOTTOM / TOP exchange ===
        do layer = 1, num_ghost

            ! --- SEND BOTTOM boundary: j = jlo + layer - 1 ---
            call send_j_index_ghost(info, dat2D, m_var, jlo+layer-1, jhi+layer, bottom, top, 300+layer, comm)

            ! --- SEND TOP boundary: j = jhi - layer + 1 ---
            call send_j_index_ghost(info, dat2D, m_var, jhi-layer+1, jlo-layer, top, bottom, 400+layer, comm)
        end do

    end subroutine

    module subroutine exchange_halos_3d(info, dat3D, m_var, num_ghost, left, right, bottom, top, back, front, comm)
        !! Exchanges halo layers for a multi-variable 3D structured field.
        !!
        !! Performs halo exchange on a field `dat3D(m_var, i, j, k)` with ghost cells,
        !! communicating ghost layers with direct MPI neighbors: left, right, bottom, top, back and front.
        class(decomp_info), intent(in) :: info
        !! Object containing the decomposition info
        double precision, allocatable :: dat3D(:,:,:,:)
        !! Multi-variable 2D array: `dat3D(m_var, i, j,k)`
        integer, intent(in) :: m_var        !! Number of variables per grid point
        integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
        integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: back         !! MPI rank of back neighbor (or MPI_PROC_NULL)
        integer, intent(in) :: front        !! MPI rank of front neighbor (or MPI_PROC_NULL)
        type(MPI_Comm), intent(in) :: comm  !! MPI communicator

        ! local variables
        integer :: layer

        ! local index range of rank
        integer :: ilo,ihi,jlo,jhi,klo,khi

        !call MPI_Comm_rank(comm,rank,ierr)

        call info%get_local_block_bounds(ilo,ihi,jlo,jhi,klo,khi)

        ! === LEFT / RIGHT exchange ===
        do layer = 1, num_ghost

            ! --- SEND LEFT boundary: i = ilo + layer - 1 ---
            call send_i_index_ghost(info, dat3D, m_var, ilo+layer-1, ihi+layer, left, right, 100+layer, comm)

            ! --- SEND RIGHT boundary: i = ihi - layer + 1 ---
            call send_i_index_ghost(info, dat3D, m_var, ihi-layer+1, ilo-layer, right, left, 200+layer, comm)

        end do

        ! === BOTTOM / TOP exchange ===
        do layer = 1, num_ghost

            ! --- SEND BOTTOM boundary: j = jlo + layer - 1 ---
            call send_j_index_ghost(info, dat3D, m_var, jlo+layer-1, jhi+layer, bottom, top, 300+layer, comm)

            ! --- SEND TOP boundary: j = jhi - layer + 1 ---
            call send_j_index_ghost(info, dat3D, m_var, jhi-layer+1, jlo-layer, top, bottom, 400+layer, comm)
        end do

        ! === BACK / FRONT exchange ===
        do layer = 1, num_ghost

            ! --- SEND BOTTOM boundary: k = klo + layer - 1 ---
            call send_k_index_ghost(info, dat3D, m_var, klo+layer-1, khi+layer, back, front, 500+layer, comm)

            ! --- SEND TOP boundary: k = khi - layer + 1 ---
            call send_k_index_ghost(info, dat3D, m_var, khi-layer+1, klo-layer, front, back, 600+layer, comm)
        end do

    end subroutine

end submodule
