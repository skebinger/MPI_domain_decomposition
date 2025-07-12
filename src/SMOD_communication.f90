! =============*FORTRAN*============ !
!   _ __ ___  _ __ (_) __| | ___| |  !
!  | `_ ` _ \| `_ \| |/ _` |/ __| |  !
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

    module subroutine send_recv_1D_array_dp(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)
        !! Sends a 1D double precision array of data (`sendbuf`) from the calling rank to rank `rank_target`
        !! and recieves from rank `rank_origin`. The recieved array `recvbuf` is returned.
        !!
        !! Either one of the communications (send/recieve) is only executed if the neighbouring rank
        !! for this communication is a valid rank (=internal boundary).
        double precision, allocatable ,intent(in) :: sendbuf(:)     !! The buffer to be sent
        double precision, allocatable, intent(inout) :: recvbuf(:)  !! The expected recieve buffer
        integer, intent(in) :: rank_target, rank_origin             !! The neighbouring ranks for the communication
        integer, intent(in) :: tag                                  !! The communication tag
        type(MPI_Comm), intent(in) :: comm                          !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_DOUBLE_PRECISION
        integer :: l_bound, u_bound, count
        integer :: rank, ierr

        ! Check input for allocation status
        if(.NOT.allocated(sendbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: sendbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: recvbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if

        l_bound = lbound(sendbuf,1)
        u_bound = ubound(sendbuf,1)
        count = u_bound-l_bound+1

        call MPI_Sendrecv(sendbuf, count, MPI_SEND_TYPE, rank_target, tag, &
            recvbuf, count, MPI_SEND_TYPE, rank_origin, tag, &
            comm, MPI_STATUS_IGNORE, ierr)
    end subroutine

    module subroutine send_recv_scalar_dp(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)
        !! Sends a 1D double precision array of data (`sendbuf`) from the calling rank to rank `rank_target`
        !! and recieves from rank `rank_origin`. The recieved array `recvbuf` is returned.
        !!
        !! Either one of the communications (send/recieve) is only executed if the neighbouring rank
        !! for this communication is a valid rank (=internal boundary).
        double precision, allocatable ,intent(in) :: sendbuf        !! The buffer to be sent
        double precision, allocatable, intent(inout) :: recvbuf     !! The expected recieve buffer
        integer, intent(in) :: rank_target, rank_origin             !! The neighbouring ranks for the communication
        integer, intent(in) :: tag                                  !! The communication tag
        type(MPI_Comm), intent(in) :: comm                          !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_DOUBLE_PRECISION
        integer :: count
        integer :: rank, ierr

        ! Check input for allocation status
        if(.NOT.allocated(sendbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: sendbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: recvbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if

        count = 1

        call MPI_Sendrecv(sendbuf, count, MPI_SEND_TYPE, rank_target, tag, &
            recvbuf, count, MPI_SEND_TYPE, rank_origin, tag, &
            comm, MPI_STATUS_IGNORE, ierr)
    end subroutine

    module subroutine send_recv_1D_array_int(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)
        !! Sends a 1D integer array of data (`sendbuf`) from the calling rank to rank `rank_target`
        !! and recieves from rank `rank_origin`. The recieved array `recvbuf` is returned.
        !!
        !! Either one of the communications (send/recieve) is only executed if the neighbouring rank
        !! for this communication is a valid rank (=internal boundary).
        integer, allocatable ,intent(in) :: sendbuf(:)              !! The buffer to be sent
        integer, allocatable, intent(inout) :: recvbuf(:)           !! The expected recieve buffer
        integer, intent(in) :: rank_target, rank_origin             !! The neighbouring ranks for the communication
        integer, intent(in) :: tag                                  !! The communication tag
        type(MPI_Comm), intent(in) :: comm                          !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_INTEGER
        integer :: l_bound, u_bound, count
        integer :: rank, ierr

        ! Check input for allocation status
        if(.NOT.allocated(sendbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: sendbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: recvbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if

        l_bound = lbound(sendbuf,1)
        u_bound = ubound(sendbuf,1)
        count = u_bound-l_bound+1

        call MPI_Sendrecv(sendbuf, count, MPI_SEND_TYPE, rank_target, tag, &
            recvbuf, count, MPI_SEND_TYPE, rank_origin, tag, &
            comm, MPI_STATUS_IGNORE, ierr)
    end subroutine

    module subroutine send_recv_scalar_int(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)
        !! Sends a 1D integer array of data (`sendbuf`) from the calling rank to rank `rank_target`
        !! and recieves from rank `rank_origin`. The recieved array `recvbuf` is returned.
        !!
        !! Either one of the communications (send/recieve) is only executed if the neighbouring rank
        !! for this communication is a valid rank (=internal boundary).
        integer, allocatable ,intent(in) :: sendbuf                 !! The buffer to be sent
        integer, allocatable, intent(inout) :: recvbuf              !! The expected recieve buffer
        integer, intent(in) :: rank_target, rank_origin             !! The neighbouring ranks for the communication
        integer, intent(in) :: tag                                  !! The communication tag
        type(MPI_Comm), intent(in) :: comm                          !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_INTEGER
        integer :: count
        integer :: rank, ierr

        ! Check input for allocation status
        if(.NOT.allocated(sendbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: sendbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: recvbuf not allocated in send_1D_array"
            call MPI_Abort(comm, 1, ierr)
        end if

        count = 1

        call MPI_Sendrecv(sendbuf, count, MPI_SEND_TYPE, rank_target, tag, &
            recvbuf, count, MPI_SEND_TYPE, rank_origin, tag, &
            comm, MPI_STATUS_IGNORE, ierr)
    end subroutine

    module subroutine find_min_max_scalar_dp(val_local, val_global, operator, comm)
        !! Identifies the min or max of an input double precision scalar variable
        !! across all ranks in the communicator `comm` and returns it.
        !!
        !! The `operator` can take either `MPI_MIN` or `MPI_MAX` as valid input.
        implicit none
        double precision, intent(in)  :: val_local      !! Local value on each rank
        double precision, intent(out) :: val_global     !! Global max (same on all ranks)
        type(MPI_Op), intent(in) :: operator            !! Operator for the communication
        type(MPI_Comm), intent(in)    :: comm           !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_DOUBLE_PRECISION
        integer :: ierr

        ! Perform Allreduce with chosen operation
        call MPI_Allreduce(val_local, val_global, 1, MPI_SEND_TYPE, operator, comm, ierr)
    end subroutine

    module subroutine find_min_max_scalar_int(val_local, val_global, operator, comm)
        !! Identifies the min or max of an input integer scalar variable
        !! across all ranks in the communicator `comm` and returns it.
        !!
        !! The `operator` can take either `MPI_MIN` or `MPI_MAX` as valid input.
        implicit none
        integer, intent(in)  :: val_local               !! Local value on each rank
        integer, intent(out) :: val_global              !! Global max (same on all ranks)
        type(MPI_Op), intent(in) :: operator            !! Operator for the communication
        type(MPI_Comm), intent(in)    :: comm           !! MPI communicator

        type(MPI_Datatype) :: MPI_SEND_TYPE = MPI_INTEGER
        integer :: ierr

        ! Perform Allreduce with chosen operation
        call MPI_Allreduce(val_local, val_global, 1, MPI_SEND_TYPE, operator, comm, ierr)
    end subroutine


end submodule
