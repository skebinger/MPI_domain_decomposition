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

    module subroutine allocate_buffers(sendbuf, recvbuf, idxlo, idxhi, m_var)
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

    module subroutine deallocate_buffers(sendbuf, recvbuf)
        !! DEallocates the send and recieve buffer arrays.
        double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer

        deallocate(sendbuf)
        deallocate(recvbuf)
    end subroutine

    module subroutine pack_2d_jslice(dat2D, buf, isend, m_var, jlo, jhi, comm)
        !! Compiles a 1D array into which all the data from `dat2D(1:mvar,isend,jlo:jhi)` is packed.
        !!
        !! `dat2D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat2D(:,:,:)   !! Multi-variable 2D array: `dat2D(m_var, i, j)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: isend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat2D`
        integer, intent(in) :: jlo                                  !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                  !! Rank-specific index boundaries of the j-slice
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var, j, num_j

        ! Check input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_j = jhi - jlo + 1

        do var = 1, m_var
            do j = jlo, jhi
                buf((var-1)*num_j + (j - jlo + 1)) = dat2D(var, isend, j)
            end do
        end do
    end subroutine

    module subroutine pack_2d_islice(dat2D, buf, jsend, m_var, ilo, ihi, comm)
        !! Compiles a 1D array into which all the data from `dat2D(1:mvar,ilo:ihi,jsend)` is packed.
        !!
        !! `dat2D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat2D(:,:,:)   !! Multi-variable 2D array: `dat2D(m_var, i, j)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: jsend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat2D`
        integer, intent(in) :: ilo                                  !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                  !! Rank-specific index boundaries of the i-slice
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, num_i

        ! Check input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_i = ihi - ilo + 1

        do var = 1, m_var
            do i = ilo, ihi
                buf((var-1)*num_i + (i - ilo + 1)) = dat2D(var, i, jsend)
            end do
        end do
    end subroutine

    module subroutine unpack_2d_jslice(dat2D, recvbuf, irecv, m_var, jlo, jhi, comm)
        !! Reverses the flattening of `dat2D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat2D`.
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)    !! Multi-variable 2D array: `dat2D(m_var, i, j)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: irecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat2D`
        integer, intent(in) :: jlo                                      !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                      !! Rank-specific index boundaries of the j-slice
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var, j, num_j

        ! Check input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_j = jhi-jlo+1

        do var = 1, m_var
            do j = jlo, jhi
                dat2D(var, irecv, j) = recvbuf((var-1)*num_j + (j - jlo + 1))
            end do
        end do
    end subroutine

    module subroutine unpack_2d_islice(dat2D, recvbuf, jrecv, m_var, ilo, ihi, comm)
        !! Reverses the flattening of `dat2D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat2D`.
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)    !! Multi-variable 2D array: `dat2D(m_var, i, j)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: jrecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat2D`
        integer, intent(in) :: ilo                                      !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                      !! Rank-specific index boundaries of the i-slice
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, num_i

        ! Check input for allocation status
        if(.NOT.allocated(dat2D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_i = ihi-ilo+1

        do var = 1, m_var
            do i = ilo, ihi
                dat2D(var, i, jrecv) = recvbuf((var-1)*num_i + (i - ilo + 1))
            end do
        end do
    end subroutine

end submodule
