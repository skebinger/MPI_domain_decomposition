submodule(MOD_MPI_decomposition) SMOD_send_recieve
    implicit none

contains

    module subroutine send_1d_i_index_ghost(info, dat1D, m_var, isend, irecv, rank_target, rank_origin, tag, comm)
        !! Sends the left or right ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat1D(:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: isend !! i-index of the cell of which data is sent
        integer, intent(in) :: irecv !! i-index of the cell to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator

        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: rank, ierr

        ! Check data input for allocation status
        if(.NOT.allocated(dat1D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf,recvbuf,m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_i_index(dat1D, sendbuf, isend, m_var, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_i_index(dat1D, recvbuf, irecv, m_var, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_2d_i_index_ghost(info, dat2D, m_var, isend, irecv, rank_target, rank_origin, tag, comm)
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
        if(rank_target /= MPI_PROC_NULL) call pack_i_index(dat2D, sendbuf, isend, m_var, info%jlow, info%jhigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_i_index(dat2D, recvbuf, irecv, m_var, info%jlow, info%jhigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_3d_i_index_ghost(info, dat3D, m_var, isend, irecv, rank_target, rank_origin, tag, comm)
        !! Sends the left or right ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: isend !! i-index of the face of which data is sent
        integer, intent(in) :: irecv !! i-index of the face to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator

        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: rank, ierr

        ! Check data input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf, recvbuf, info%jlow, info%jhigh, info%klow, info%khigh, m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_i_index(dat3D, sendbuf, isend, m_var, info%jlow, info%jhigh, info%klow, info%khigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_i_index(dat3D, recvbuf, irecv, m_var, info%jlow, info%jhigh, info%klow, info%khigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_2d_j_index_ghost(info, dat2D, m_var, jsend, jrecv, rank_target, rank_origin, tag, comm)
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
        if(rank_target /= MPI_PROC_NULL) call pack_j_index(dat2D, sendbuf, jsend, m_var, info%ilow, info%ihigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_j_index(dat2D, recvbuf, jrecv, m_var, info%ilow, info%ihigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_3d_j_index_ghost(info, dat3D, m_var, jsend, jrecv, rank_target, rank_origin, tag, comm)
        !! Sends the top or bottom ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)
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
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf, recvbuf, info%ilow, info%ihigh, info%klow, info%khigh, m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_j_index(dat3D, sendbuf, jsend, m_var, info%ilow, info%ihigh, info%klow, info%khigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_j_index(dat3D, recvbuf, jrecv, m_var, info%ilow, info%ihigh, info%klow, info%khigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

    module subroutine send_3d_k_index_ghost(info, dat3D, m_var, ksend, krecv, rank_target, rank_origin, tag, comm)
        !! Sends the top or bottom ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! `target` is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: ksend !! k-index of the slice of which data is sent
        integer, intent(in) :: krecv !! k-index of the slice to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator

        double precision, allocatable :: sendbuf(:), recvbuf(:)
        integer :: rank, ierr

        ! Check data input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Create the send and recieve buffer
        call allocate_buffers(sendbuf, recvbuf, info%ilow, info%ihigh, info%jlow, info%jhigh, m_var)
        ! Pack the ghost slice data to be sent into a 1D contiguous array
        ! Only if the target rank is avalid target (=internal boudnary) do the work and pack it
        if(rank_target /= MPI_PROC_NULL) call pack_k_index(dat3D, sendbuf, ksend, m_var, info%ilow, info%ihigh, info%jlow, info%jhigh, comm)

        call send_rec_1D_array(sendbuf, recvbuf, rank_target, rank_origin, tag, comm)

        ! only execute if target is a valid rank (=internal boundary)
        if(rank_origin /= MPI_PROC_NULL) call unpack_k_index(dat3D, recvbuf, krecv, m_var, info%ilow, info%ihigh, info%jlow, info%jhigh, comm)

        call deallocate_buffers(sendbuf,recvbuf)
    end subroutine

end submodule
