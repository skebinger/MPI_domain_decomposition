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

    if (rank == 0) print *, "=== 1D Decomposition but with userdefined processor grid ==="
    call initialize_decomposition(info, m_xi=16, m_eta=0,  m_tau=0, udf_dims=[8,0,0],  comm=comm)
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

    if (rank == 0) print *, "=== 2D Decomposition but with userdefined processor grid ==="
    call initialize_decomposition(info, m_xi=16, m_eta=12, m_tau=0, udf_dims=[2,4,0],  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    call initialize_decomposition(info, m_xi=16, m_eta=12, m_tau=0, udf_dims=[1,8,0],  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    ! ---- CASE 3: 3D decomposition ----
    if (rank == 0) print *, "=== 3D Decomposition ==="
    call initialize_decomposition(info, m_xi=16, m_eta=79, m_tau=8,  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    if (rank == 0) print *, "=== 3D Decomposition but with userdefined processor grid ==="
    call initialize_decomposition(info, m_xi=24, m_eta=12, m_tau=9, udf_dims=[1,4,2],  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)

    call MPI_Barrier(comm)

    call initialize_decomposition(info, m_xi=13, m_eta=12, m_tau=8, udf_dims=[1,2,4],  comm=comm)
    call setup_cartesian_topology(info, comm)
    call print_decomposition_summary(info, comm)
    call print_cartesian_rank_layout(info)


    if (rank == 0) print *, "=== OUTPUT TEST COMPLETED ==="

    ! Finalize MPI
    call MPI_Finalize(ierr)
end program test_decomposition_driver
