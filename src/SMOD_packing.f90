submodule(MOD_MPI_decomposition) SMOD_packing
    implicit none

contains

    module subroutine pack_1d_i_index(dat1D, buf, isend, m_var, comm)
        !! Compiles a 1D array into which all the data from `dat1D(1:mvar,isend)` is packed.
        !!
        !! `dat1D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat1D(:,:)     !! Multi-variable 1D array: `dat1D(m_var, i)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: isend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat1D`
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var

        ! Check input for allocation status
        if(.NOT.allocated(dat1D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat1D not allocated in pack_1d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_1d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        do var = 1, m_var
            buf(var) = dat1D(var, isend)
        end do
    end subroutine

    module subroutine pack_2d_i_index(dat2D, buf, isend, m_var, jlo, jhi, comm)
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
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_i_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_2dpack_2d_i_index_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_j = jhi - jlo + 1

        do j = jlo, jhi
            do var = 1, m_var
                buf((var-1)*num_j + (j - jlo + 1)) = dat2D(var, isend, j)
            end do
        end do
    end subroutine

    module subroutine pack_3d_i_index(dat3D, buf, isend, m_var, jlo, jhi, klo, khi, comm)
        !! Compiles a 1D array into which all the data from `dat3D(1:mvar,isend,jlo:jhi,klo:khi)` is packed.
        !!
        !! `dat3D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat3D(:,:,:,:) !! Multi-variable 3D array: `dat3D(m_var, i, j)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: isend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat3D`
        integer, intent(in) :: jlo                                  !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                  !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: klo                                  !! Rank-specific index boundaries of the k-slice
        integer, intent(in) :: khi                                  !! Rank-specific index boundaries of the k-slice
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var, j, k, idx

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in pack_3d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_3d_jslice"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        idx = 1
        do k = klo, khi
            do j = jlo, jhi
                do var = 1, m_var
                    buf(idx) = dat3D(var, isend, j, k)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine

    module subroutine pack_2d_j_index(dat2D, buf, jsend, m_var, ilo, ihi, comm)
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
            print *, "RANK ", rank, " reports: dat2D not allocated in pack_2d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_2d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        num_i = ihi - ilo + 1

        do i = ilo, ihi
            do var = 1, m_var
                buf((var-1)*num_i + (i - ilo + 1)) = dat2D(var, i, jsend)
            end do
        end do
    end subroutine

    module subroutine pack_3d_j_index(dat3D, buf, jsend, m_var, ilo, ihi, klo, khi, comm)
        !! Compiles a 1D array into which all the data from `dat3D(1:mvar,ilo:ihi,jsend,klo:khi)` is packed.
        !!
        !! `dat3D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat3D(:,:,:,:) !! Multi-variable 3D array: `dat3D(m_var, i, j, k)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: jsend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat3D`
        integer, intent(in) :: ilo                                  !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                  !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: klo                                  !! Rank-specific index boundaries of the k-slice
        integer, intent(in) :: khi                                  !! Rank-specific index boundaries of the k-slice
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, k, index

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in pack_3d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_3d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        index = 1
        do k = klo, khi
            do i = ilo, ihi
                do var = 1, m_var
                    buf(index) = dat3D(var, i, jsend, k)
                    index = index + 1
                end do
            end do
        end do
    end subroutine

    module subroutine pack_3d_k_index(dat3D, buf, ksend, m_var, ilo, ihi, jlo, jhi, comm)
        !! Compiles a 1D array into which all the data from `dat3D(1:mvar,ilo:ihi,jlo:jhi,ksend)` is packed.
        !!
        !! `dat3D` is flattened into a 1D contiguous array.
        double precision, allocatable, intent(in) :: dat3D(:,:,:,:)  !! Multi-variable 3D array: `dat3D(m_var, i, j, k)`
        double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
        integer, intent(in) :: ksend                                !! Slice index of the sending rank
        integer, intent(in) :: m_var                                !! Number of field variables in `dat3D`
        integer, intent(in) :: ilo                                  !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: ihi                                  !! Rank-specific index boundaries of the i-slice
        integer, intent(in) :: jlo                                  !! Rank-specific index boundaries of the j-slice
        integer, intent(in) :: jhi                                  !! Rank-specific index boundaries of the j-slice
        type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator

        integer :: rank, ierr
        integer :: var, i, j, index

        ! Check input for allocation status
        if(.NOT.allocated(dat3D)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: dat3D not allocated in pack_3d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(.NOT.allocated(buf)) then
            call MPI_Comm_rank(comm,rank,ierr)
            print *, "RANK ", rank, " reports: buffer not allocated in pack_3d_j_index"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        index = 1
        do j = jlo, jhi
            do i = ilo, ihi
                do var = 1, m_var
                    buf(index) = dat3D(var, i, j, ksend)
                    index = index + 1
                end do
            end do
        end do
    end subroutine

end submodule
