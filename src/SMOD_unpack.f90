submodule(MOD_MPI_decomposition) SMOD_unpack
    implicit none

contains

    module subroutine unpack_1d_i_index(dat1D, recvbuf, irecv, m_var, comm)
        !! Reverses the flattening of `dat1D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat1D`.
        double precision, allocatable, intent(inout) :: dat1D(:,:)      !! Multi-variable 1D array: `dat1D(m_var, i)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: irecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat1D`
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var

        ! Check input for allocation status
        if(.NOT.allocated(dat1D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat1D not allocated in unpack_1d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat1D not allocated in unpack_1d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        do var = 1, m_var
            dat1D(var, irecv) = recvbuf(var)
        end do
    end subroutine

    module subroutine unpack_2d_i_index(dat2D, recvbuf, irecv, m_var, jlo, jhi, comm)
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
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_j = jhi-jlo+1

        do j = jlo, jhi
            do var = 1, m_var
                dat2D(var, irecv, j) = recvbuf((var-1)*num_j + (j - jlo + 1))
            end do
        end do
    end subroutine

    module subroutine unpack_3d_i_index(dat3D, recvbuf, irecv, m_var, jlo, jhi, klo, khi, comm)
        !! Reverses the flattening of `dat3D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat3D`.
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)  !! Multi-variable 3D array: `dat3D(m_var, i, j, k)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: irecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat3D`
        integer, intent(in) :: jlo                                      !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                      !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: klo                                      !! Rank-specific index boundaries of the k-slice
        integer, intent(in) :: khi                                      !! Rank-specific index boundaries of the k-slice
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var, j, k, index

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        index = 1
        do k = klo, khi
            do j = jlo, jhi
                do var = 1, m_var
                    dat3D(var, irecv, j, k) = recvbuf(index)
                    index = index + 1
                end do
            end do
        end do
    end subroutine

    module subroutine unpack_2d_j_index(dat2D, recvbuf, jrecv, m_var, ilo, ihi, comm)
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
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_2d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_i = ihi-ilo+1

        do i = ilo, ihi
            do var = 1, m_var
                dat2D(var, i, jrecv) = recvbuf((var-1)*num_i + (i - ilo + 1))
            end do
        end do
    end subroutine

    module subroutine unpack_3d_j_index(dat3D, recvbuf, jrecv, m_var, ilo, ihi, klo, khi, comm)
        !! Reverses the flattening of `dat3D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat3D`.
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)  !! Multi-variable 3D array: `dat3D(m_var, i, j, k)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: jrecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat3D`
        integer, intent(in) :: ilo                                      !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                      !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: klo                                      !! Rank-specific index boundaries of the k-slice
        integer, intent(in) :: khi                                      !! Rank-specific index boundaries of the k-slice
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, k, index

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        index = 1
        do k = klo, khi
            do i = ilo, ihi
                do var = 1, m_var
                    dat3D(var, i, jrecv, k) = recvbuf(index)
                    index = index + 1
                end do
            end do
        end do
    end subroutine

    module subroutine unpack_3d_jkindex(dat3D, recvbuf, krecv, m_var, ilo, ihi, jlo, jhi, comm)
        !! Reverses the flattening of `dat3D` into a 1D array and unpacks the 1D array
        !! `recbuf` back into `dat3D`.
        double precision, allocatable, intent(inout) :: dat3D(:,:,:,:)  !! Multi-variable 3D array: `dat3D(m_var, i, j, k)`
        double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
        integer, intent(in) :: krecv                                    !! Slice index in recieving rank
        integer, intent(in) :: m_var                                    !! Number of field variables in `dat3D`
        integer, intent(in) :: ilo                                      !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                      !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: jlo                                      !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                      !! Rank-specific index boundaries of the j-slice
        type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, j, index

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.NOT.allocated(recvbuf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat2D not allocated in unpack_3d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        index = 1
        do j = jlo, jhi
            do i = ilo, ihi
                do var = 1, m_var
                    dat3D(var, i, j, krecv) = recvbuf(index)
                    index = index + 1
                end do
            end do
        end do
    end subroutine

end submodule
