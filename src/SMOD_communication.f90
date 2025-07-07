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

submodule(MOD_MPI_decomposition) SMOD_communication
    implicit none
contains

    module subroutine send_1D_array(sendbuf,recvbuf,origin,target,tag,MPI_SEND_TYPE,comm)
        !! Sends a 1D array of data ('sendbuf') from rankg 'origin' to rank 'target'.
        !! The recieved array 'recvbuf' is returned.
        !!
        !! The data to be sent has the type 'MPI_SEND_TYPE', e.g. 'MPI_DOUBLE_PRECISION' and
        !! is sent through the MPI communicator 'comm' with the communication 'tag'.
        double precision, allocatable ,intent(in) :: sendbuf(:)
        double precision, allocatable, intent(inout) :: recvbuf(:)
        integer, intent(in) :: origin, target
        integer, intent(in) :: tag
        type(MPI_Datatype), intent(in) :: MPI_SEND_TYPE
        type(MPI_Comm), intent(in) :: comm

        integer :: l_bound, u_bound, count
        integer :: rank, ierr

        ! Check input for allocation status
        if(.NOT.allocated(sendbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: sendbuf not allocated in send_1D_array"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: recvbuf not allocated in send_1D_array"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        l_bound = lbound(sendbuf,1)
        u_bound = ubound(sendbuf,1)
        count = u_bound-l_bound+1

        call MPI_Sendrecv(sendbuf, count, MPI_SEND_TYPE, origin, tag, &
            recvbuf, count, MPI_SEND_TYPE, target, tag, &
            comm, MPI_STATUS_IGNORE, ierr)
    end subroutine


end submodule
