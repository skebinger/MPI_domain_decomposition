
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
    !! The decomposition info is stored in a type 'decomp_info'.
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

    type :: decomp_info
        integer :: ilow, ihigh                                 !! Local i-direction bounds
        integer :: jlow, jhigh                                 !! Local j-direction bounds

        integer :: nbr_top, nbr_bottom, nbr_left, nbr_right    !! neighbouring ranks
        integer :: comm_cart                                   !! cartesian communicator

        integer :: dims(2)                                     !! Processor grid in i and j directions
    contains
        procedure :: initialize_decomposition
        procedure :: exchange_halos
        procedure :: setup_cartesian_topology
        procedure :: print_decomposition_summary
        procedure :: print_cartesian_rank_layout
        procedure :: get_local_block_bounds
        procedure :: get_cartesian_comm
        procedure :: get_neighbouring_ranks
    end type


    interface
        module subroutine initialize_decomposition(info, m_xi, m_eta, comm)
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
            class(decomp_info), intent(out) :: info !! Object containing the decomposition info
            integer, intent(in) :: m_xi             !! Total number of cells in the i-direction (xi-axis)
            integer, intent(in) :: m_eta            !! Total number of cells in the j-direction (eta-axis)
            integer, intent(in) :: comm             !! MPI communicator (normally 'MPI_COMM_WORLD')
        end subroutine

        module subroutine exchange_halos(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
            !! Exchanges halo layers for a multi-variable 2D structured field.
            !!
            !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
            !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.
            class(decomp_info), intent(in) :: info
            !! Object containing the decomposition info
            double precision, allocatable :: dat2D(:,:,:)
            !! Multi-variable 2D array: dat2D(m_var, i, j)
            integer, intent(in) :: m_var        !! Number of variables per grid point
            integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
            integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: comm         !! MPI communicator
        end subroutine

        module subroutine setup_cartesian_topology(info,comm_in)
            !! Create a 2D Cartesian communicator and get neighbor ranks.
            !!
            !! The number of processes in each dimension is taken from previous
            !! calculation executed in initialize_decomposition. Sets the global
            !! module variable 'comm_cart'.
            class(decomp_info) :: info      !! Object containing the decomposition info
            integer, intent(in) :: comm_in  !! Input communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_decomposition_summary(info,comm)
            !! Prints the domain decomposition summary for all MPI ranks.
            !!
            !! Gathers the local index ranges from all ranks and prints a formatted
            !! table on rank 0 showing the (i,j) bounds for each rank's block.
            class(decomp_info) :: info  !! Object containing the decomposition info
            integer, intent(in) :: comm !! MPI communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_cartesian_rank_layout(info)
            !! Prints a visual representation of the rank layout in the form of a chessboard.
            use mpi
            class(decomp_info), intent(in) :: info
        end subroutine

        module subroutine get_local_block_bounds(info,ilo, ihi, jlo, jhi)
            !! Returns the local 2D index bounds for this rank.
            !!
            !! These bounds define the owned rectangular region in the global mesh for
            !! all directional sweeps. Ghost cells must be added externally by the solver.
            class(decomp_info), intent(in) :: info      !! Object containing the decomposition info
            integer, intent(out) :: ilo, ihi, jlo, jhi  !! local indexing range
        end subroutine
        module function get_cartesian_comm(info) result(comm)
            !! Returns the stored Cartesian communicator
            class(decomp_info), intent(in) :: info
            integer :: comm
        end function
        module subroutine get_neighbouring_ranks(info,left, right, bottom, top)
            !! Returns the ranks of neighboring processes
            class(decomp_info), intent(in) :: info              !! Object containing the decomposition info
            integer, intent(out) :: left, right, bottom, top    !! local indexing bounds
        end subroutine
    end interface

contains

end module MOD_MPI_decomposition
