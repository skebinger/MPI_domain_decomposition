program test_decomposition_driver
    use mpi_f08
    use MOD_MPI_decomposition
    implicit none

    type(decomp_info) :: info
    integer :: rank
    type(MPI_Comm) :: comm
    integer :: ierr

    ! Initialize MPI
    call MPI_Init(ierr)
    comm = MPI_COMM_WORLD
    call MPI_Comm_rank(comm, rank, ierr)

    ! ---- CASE 1: 1D decomposition ----
    if (rank == 0) print *, "=== 1D Decomposition ==="
    call initialize_decomposition(info, m_xi=16, m_eta=0,  m_tau=0,  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    ! ---- CASE 2: 2D decomposition ----
    if (rank == 0) print *, "=== 2D Decomposition ==="
    call initialize_decomposition(info, m_xi=16, m_eta=12, m_tau=0,  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    ! ---- CASE 3: 3D decomposition ----
    if (rank == 0) print *, "=== 3D Decomposition ==="
    call initialize_decomposition(info, m_xi=16, m_eta=12, m_tau=8,  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    ! Finalize MPI
    call MPI_Finalize(ierr)
end program test_decomposition_driver
