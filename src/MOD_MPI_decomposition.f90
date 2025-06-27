
module MOD_MPI_decomposition
    !! (so far 2D) block domain decomposition for MPI-based structured grid solvers.
    !!
    !! This module handles a static 2D block decomposition using MPI. Each MPI rank
    !! owns a rectangular subdomain (i,j) of the global mesh. Ghost cells are not
    !! included in the decomposition ranges and must be handled externally.
    !! Each rank is assigned a portion of the global index range. These ensure
    !! continuous field indizes across ranks without the need to calculate a global
    !! index.
    !!
    !! Setup might be updated some day to also include sliced layout. Currently
    !! only chunked decomposition is supported.
    !!
    !! All sweep directions (xi and eta) operate over the same local block.

    use mpi
    implicit none

    public :: initialize_decomposition
    public :: get_local_block_bounds
    public :: print_decomposition_summary
    public :: setup_cartesian_topology
    public :: get_cartesian_comm
    public :: get_neighbouring_ranks
    public :: exchange_halos

    integer, private :: ilow, ihigh                                 !! Local i-direction bounds
    integer, private :: jlow, jhigh                                 !! Local j-direction bounds

    integer, private :: nbr_top, nbr_bottom, nbr_left, nbr_right    !! neighbouring ranks
    integer, private :: comm_cart                                   !! cartesian communicator

    integer, private :: dims(2)                                     !! Processor grid in i and j directions


    interface
        module subroutine exchange_halos(dat2D, m_var, m_xi, m_eta, num_ghost, left, right, bottom, top, comm)
            !! Exchanges halo layers for a multi-variable 2D structured field.
            !!
            !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
            !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.

            double precision, allocatable :: dat2D(:,:,:)
            !! Multi-variable 2D array: dat2D(m_var, i, j)
            integer, intent(in) :: m_var        !! Number of variables per grid point
            integer, intent(in) :: m_xi         !! Number of *interior* i-cells (excluding ghosts)
            integer, intent(in) :: m_eta        !! Number of *interior* j-cells (excluding ghosts)
            integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
            integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: comm         !! MPI communicator
        end subroutine


        module subroutine initialize_decomposition(m_xi, m_eta, comm)
            !! Initializes the domain decomposition for a 2D grid using MPI.
            !!
            !! The decomposition is performed only by the master rank (rank 0), which computes
            !! a 2D processor grid decomposition of the global domain into chunks.
            !! Each chunk corresponds to a rectangular subdomain of the global mesh.
            !!
            !! Rank 0 then sends the computed local index ranges (bounds) to all other ranks
            !! using point-to-point MPI communication. Each rank receives its assigned local
            !! ranges for both i (xi) and j (eta) directions.
            !!
            !! This approach avoids all ranks calling MPI_Dims_create, thus preventing
            !! potential segmentation faults and ensuring consistent decomposition.
            !!
            !! The subroutine uses a helper `compute_range` to divide each dimension evenly
            !! among the processors assigned along that dimension, distributing any remainder
            !! cells to the lowest ranks.
            !!
            !! Input:
            !! - m_xi: Total number of cells in the i-direction (xi-axis).
            !! - m_eta: Total number of cells in the j-direction (eta-axis).
            !! - comm: MPI communicator over which ranks are defined.
            !!
            !! Output (module private variables set):
            !! - ilow, ihigh: Local i-direction bounds (inclusive).
            !! - jlow, jhigh: Local j-direction bounds (inclusive).

            integer, intent(in) :: m_xi !! Total number of cells in the i-direction (xi-axis)
            integer, intent(in) :: m_eta !! Total number of cells in the j-direction (eta-axis)
            integer, intent(in) :: comm !! MPI communicator (normally 'MPI_COMM_WORLD')
        end subroutine

        module subroutine setup_cartesian_topology(comm_in)
            !! Create a 2D Cartesian communicator and get neighbor ranks.
            !!
            !! The number of processes in each dimension is taken from previous
            !! calculation executed in initialize_decomposition. Sets the global
            !! module variable 'comm_cart'.
            integer, intent(in) :: comm_in !! Input communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_decomposition_summary(comm)
            !! Prints the domain decomposition summary for all MPI ranks.
            !!
            !! Gathers the local index ranges from all ranks and prints a formatted
            !! table on rank 0 showing the (i,j) bounds for each rank's block.
            integer, intent(in) :: comm !! MPI communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_cartesian_rank_layout()
            !! Prints a visual representation of the rank layout.
        end subroutine
    end interface

contains

    subroutine get_local_block_bounds(ilo, ihi, jlo, jhi)
        !! Returns the local 2D index bounds for this rank.
        !!
        !! These bounds define the owned rectangular region in the global mesh for
        !! all directional sweeps. Ghost cells must be added externally by the solver.
        integer :: ilo, ihi, jlo, jhi
        ilo = ilow
        ihi = ihigh
        jlo = jlow
        jhi = jhigh
    end subroutine get_local_block_bounds

    function get_cartesian_comm() result(comm)
        !! Returns the stored Cartesian communicator
        integer :: comm
        comm = comm_cart
    end function get_cartesian_comm

    subroutine get_neighbouring_ranks(left, right, bottom, top)
        !! Returns the ranks of neighboring processes
        integer, intent(out) :: left, right, bottom, top
        left = nbr_left
        right = nbr_right
        bottom = nbr_bottom
        top = nbr_top
    end subroutine get_neighbouring_ranks


end module MOD_MPI_decomposition
