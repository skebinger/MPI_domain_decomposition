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

submodule(MOD_MPI_decomposition) SMOD_helper

    implicit none

contains

    module subroutine allocate_buffers_1d(sendbuf, recvbuf, m_var)
        !! Allocates the send and recieve buffer arrays.
        !!
        !! Takes the two bounds `idxlo` and `idxhi` as lower and upper bounds of one slice of data.
        !! Then `m_var` is the number of slices (= number of field variables).
        double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
        integer, intent(in) :: m_var                                !! Number of field variables

        integer :: count

        count = m_var

        allocate(sendbuf(count))
        allocate(recvbuf(count))
    end subroutine

    module subroutine allocate_buffers_2d(sendbuf, recvbuf, idxlo, idxhi, m_var)
        !! Allocates the send and recieve buffer arrays.
        !!
        !! Takes the two bounds `idxlo` and `idxhi` as lower and upper bounds of one slice of data.
        !! Then `m_var` is the number of slices (= number of field variables).
        double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
        integer, intent(in) :: idxlo                                !! Bounds of array slice
        integer, intent(in) :: idxhi                                !! Bounds of array slice
        integer, intent(in) :: m_var                                !! Number of field variables

        integer :: count, num_idx

        num_idx = idxhi-idxlo+1
        count = m_var*num_idx

        allocate(sendbuf(count))
        allocate(recvbuf(count))
    end subroutine

    module subroutine allocate_buffers_3d(sendbuf, recvbuf, idxlo_i, idxhi_i, idxlo_j, idxhi_j, m_var)
        !! Allocates the send and recieve buffer arrays.
        !!
        !! Takes the two bounds `idxlo_i` and `idxhi_i` as lower and upper bounds of first slicing direction of data.
        !! Takes the two bounds `idxlo_j` and `idxhi_j` as lower and upper bounds of the second direction of the data.
        !! Then `m_var` is the number of slices (= number of field variables).
        double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
        integer, intent(in) :: idxlo_i                              !! Bounds of array slice
        integer, intent(in) :: idxhi_i                              !! Bounds of array slice
        integer, intent(in) :: idxlo_j                              !! Bounds of array slice
        integer, intent(in) :: idxhi_j                              !! Bounds of array slice
        integer, intent(in) :: m_var                                !! Number of field variables

        integer :: count, num_idx_i, num_idx_j

        num_idx_i = idxhi_i-idxlo_i+1
        num_idx_j = idxhi_j-idxlo_j+1
        count = m_var*num_idx_i*num_idx_j

        allocate(sendbuf(count))
        allocate(recvbuf(count))
    end subroutine

    module subroutine deallocate_buffers(sendbuf, recvbuf)
        !! DEallocates the send and recieve buffer arrays.
        double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer

        deallocate(sendbuf)
        deallocate(recvbuf)
    end subroutine

end submodule
