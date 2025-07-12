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

    module subroutine get_local_block_bounds(info, ilo, ihi, jlo, jhi, klo, khi)
        !! Returns the local 2D index bounds for this rank.
        !!
        !! These bounds define the owned rectangular region in the global mesh for
        !! all directional sweeps. Ghost cells must be added externally by the solver.
        class(decomp_info), intent(in) :: info  !! Object containing the decomposition info
        integer, intent(out) :: ilo             !! local indexing range
        integer, intent(out) :: ihi             !! local indexing range
        integer, intent(out), optional :: jlo   !! local indexing range
        integer, intent(out), optional :: jhi   !! local indexing range
        integer, intent(out), optional :: klo   !! local indexing range
        integer, intent(out), optional :: khi   !! local indexing range
        ilo = info%ilow
        ihi = info%ihigh
        if(present(jlo)) jlo = info%jlow
        if(present(jhi)) jhi = info%jhigh
        if(present(klo)) klo = info%klow
        if(present(khi)) khi = info%khigh
    end subroutine get_local_block_bounds

    module function get_cartesian_comm(info) result(comm)
        !! Returns the stored Cartesian communicator
        class(decomp_info), intent(in) :: info  !! Object containing the decomposition info
        type(MPI_Comm) :: comm                  !! MPI cartesian communicator

        comm = info%comm_cart
    end function get_cartesian_comm

    module subroutine get_neighbouring_ranks(info, left, right, bottom, top, back, front)
        !! Returns the ranks of neighboring processes
        class(decomp_info), intent(in) :: info      !! Object containing the decomposition info
        integer, intent(out) :: left                !! Neighbouring ranks
        integer, intent(out) :: right               !! Neighbouring ranks
        integer, intent(out), optional :: bottom    !! Neighbouring ranks
        integer, intent(out), optional :: top       !! Neighbouring ranks
        integer, intent(out), optional :: back      !! Neighbouring ranks
        integer, intent(out), optional :: front     !! Neighbouring ranks

        left                        = info%nbr_left
        right                       = info%nbr_right
        if(present(bottom)) bottom  = info%nbr_bottom
        if(present(top))    top     = info%nbr_top
        if(present(back))   back    = info%nbr_back
        if(present(front))  front   = info%nbr_front
    end subroutine get_neighbouring_ranks
end submodule
