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

submodule(MOD_MPI_decomposition) SMOD_get
    implicit none
contains

    module subroutine get_local_block_bounds(info, ilo, ihi, jlo, jhi)
        !! Returns the local 2D index bounds for this rank.
        !!
        !! These bounds define the owned rectangular region in the global mesh for
        !! all directional sweeps. Ghost cells must be added externally by the solver.
        class(decomp_info), intent(in) :: info      !! Object containing the decomposition info
        integer, intent(out) :: ilo, ihi, jlo, jhi  !! local indexing range
        ilo = info%ilow
        ihi = info%ihigh
        jlo = info%jlow
        jhi = info%jhigh
    end subroutine get_local_block_bounds

    module function get_cartesian_comm(info) result(comm)
        !! Returns the stored Cartesian communicator
        class(decomp_info), intent(in) :: info !! Object containing the decomposition info
        type(MPI_Comm) :: comm
        comm = info%comm_cart
    end function get_cartesian_comm

    module subroutine get_neighbouring_ranks(info, left, right, bottom, top)
        !! Returns the ranks of neighboring processes
        class(decomp_info), intent(in) :: info              !! Object containing the decomposition info
        integer, intent(out) :: left, right, bottom, top    !! local indexing bounds
        left    = info%nbr_left
        right   = info%nbr_right
        bottom  = info%nbr_bottom
        top     = info%nbr_top
    end subroutine get_neighbouring_ranks
end submodule
