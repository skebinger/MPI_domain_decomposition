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
!  I built this library as a means for a private CFD project in the 
!  hope that it might also be useful for others looking for an 
!  accessible entry into MPI programming in CFD. Feel free to change
!  and improve it as you please.
!  
!  This library is distributed under the BSD 3-Clause License.
! ================================================================= !

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

    use mpi_f08
    implicit none

    type :: decomp_info
        !! Contains the domain decomposition information, such as
        !! indexing bounds and rank neighbours.
        integer :: ilow, ihigh                                 !! Local i-direction bounds
        integer :: jlow, jhigh                                 !! Local j-direction bounds

        integer :: nbr_top, nbr_bottom, nbr_left, nbr_right    !! neighbouring ranks
        type(MPI_Comm) :: comm_cart                            !! cartesian communicator

        integer :: dims(2)                                     !! Processor grid in i and j directions
    contains
        procedure :: initialize_decomposition
        procedure :: exchange_halos_2D
        procedure :: setup_cartesian_topology
        procedure :: print_decomposition_summary
        procedure :: print_cartesian_rank_layout
        procedure :: get_local_block_bounds
        procedure :: get_cartesian_comm
        procedure :: get_neighbouring_ranks
    end type


    !=======================================================================================================!
    ! setup and initialization
    !=======================================================================================================!
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
            type(MPI_Comm), intent(in) :: comm      !! MPI communicator (normally 'MPI_COMM_WORLD')
        end subroutine

        module subroutine setup_cartesian_topology(info, comm_in)
            !! Create a 2D Cartesian communicator and get neighbor ranks.
            !!
            !! The number of processes in each dimension is taken from previous
            !! calculation executed in initialize_decomposition. Sets the global
            !! module variable 'comm_cart'.
            class(decomp_info) :: info      !! Object containing the decomposition info
            type(MPI_Comm), intent(in) :: comm_in  !! Input communicator (usually MPI_COMM_WORLD)
        end subroutine
    end interface

    !=======================================================================================================!
    ! info-output subroutines
    !=======================================================================================================!
    interface
        module subroutine print_decomposition_summary(info, comm)
            !! Prints the domain decomposition summary for all MPI ranks.
            !!
            !! Gathers the local index ranges from all ranks and prints a formatted
            !! table on rank 0 showing the (i,j) bounds for each rank's block.
            class(decomp_info) :: info  !! Object containing the decomposition info
            type(MPI_Comm), intent(in) :: comm !! MPI communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_cartesian_rank_layout(info)
            !! Prints a visual representation of the rank layout in the form of a chessboard.
            class(decomp_info), intent(in) :: info
        end subroutine
    end interface

    !=======================================================================================================!
    ! type-bound getters
    !=======================================================================================================!
    interface
        module subroutine get_local_block_bounds(info, ilo, ihi, jlo, jhi)
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
            type(MPI_Comm) :: comm
        end function
        module subroutine get_neighbouring_ranks(info,left, right, bottom, top)
            !! Returns the ranks of neighboring processes
            class(decomp_info), intent(in) :: info              !! Object containing the decomposition info
            integer, intent(out) :: left, right, bottom, top    !! local indexing bounds
        end subroutine
    end interface

    !=======================================================================================================!
    ! halo exchange
    !=======================================================================================================!
    interface
        module subroutine exchange_halos_2D(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
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
            type(MPI_Comm), intent(in) :: comm  !! MPI communicator
        end subroutine

        module subroutine send_j_ghost(info, dat2D, m_var, isend, irecv, rank_target, rank_origin, tag, comm)
        !! Sends the left or right ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! 'target' is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: isend !! i-index of the slice of which data is sent
        integer, intent(in) :: irecv !! i-index of the slice to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator
        end subroutine

        module subroutine send_i_ghost(info, dat2D, m_var, jsend, jrecv, rank_target, rank_origin, tag, comm)
        !! Sends the top or bottom ghost cells of a rank to its neighbour.
        !! The origin is defined as the neighbouring rank sending data to the calling rank,
        !! 'target' is the rank, that recieves data from this (calling) rank.
        type(decomp_info) :: info
        double precision, allocatable, intent(inout) :: dat2D(:,:,:)
        integer, intent(in) :: m_var !! Number of variables per grid point
        integer, intent(in) :: jsend !! j-index of the slice of which data is sent
        integer, intent(in) :: jrecv !! j-index of the slice to which data is recieved
        integer, intent(in) :: rank_target !! Origin rank of an MPI Transmission
        integer, intent(in) :: rank_origin !! Target rank of an MPI Transmission
        integer, intent(in) :: tag !! Communication tag
        type(MPI_Comm), intent(in) :: comm !! MPI communicator
        end subroutine
    end interface

    !=======================================================================================================!
    ! rank communication subroutines
    !=======================================================================================================!
    interface
        module subroutine send_1D_array(sendbuf,recvbuf,origin,target,tag,MPI_SEND_TYPE,comm)
            double precision, allocatable ,intent(in) :: sendbuf(:)
            double precision, allocatable, intent(inout) :: recvbuf(:)
            integer, intent(in) :: origin, target
            integer, intent(in) :: tag
            type(MPI_Datatype), intent(in) :: MPI_SEND_TYPE
            type(MPI_Comm), intent(in) :: comm
        end subroutine
    end interface

    !=======================================================================================================!
    ! helper subroutines
    !=======================================================================================================!
    interface pack_islice
        module procedure :: pack_2d_islice
        !module procedure :: pack_3d_islice
    end interface
    interface pack_jslice
        module procedure :: pack_2d_jslice
        !module procedure :: pack_3d_jslice
    end interface

    interface
        module subroutine allocate_buffers(sendbuf, recvbuf, idxlo, idxhi, m_var)
            !! Allocates the send and recieve buffer arrays.
            !!
            !! Takes the two bounds 'idxlo' and 'idxhi' as lower and upper bounds of one slice of data.
            !! Then 'm_var' is the number of slices (= number of field variables).
            double precision, allocatable, intent(inout) :: sendbuf(:), recvbuf(:) !! buffer arrays
            integer, intent(in) :: idxlo, idxhi !! Bounds of array slice
            integer, intent(in) :: m_var !! Number of field variables
        end subroutine

        module subroutine pack_2d_jslice(dat2D, buf, isend, m_var, jlo, jhi, comm)
            !! Compiles a 1D array into which all the data from 'dat2D(1:mvar,isend,jlo:jhi)' is packed.
            !!
            !! 'dat2D' is flattened into a 1D contiguous array.
            double precision, allocatable, intent(in) :: dat2D(:,:,:) !! Multi-variable 2D array: dat2D(m_var, i, j)
            double precision, allocatable, intent(inout) :: buf(:) !! Buffer (flattened) 1D array for the mpi send transmission
            integer, intent(in) :: isend !! Slice index of the sending rank
            integer, intent(in) :: m_var !! Number of field variables in 'dat2D'
            integer, intent(in) :: jlo, jhi !! Rank-specific index boundaries of the j-slice
            type(MPI_Comm), intent(in) :: comm !! MPI Communicator
        end subroutine

        module subroutine pack_2d_islice(dat2D, buf, jsend, m_var, ilo, ihi, comm)
            !! Compiles a 1D array into which all the data from 'dat2D(1:mvar,ilo:ihi,jsend)' is packed.
            !!
            !! 'dat2D' is flattened into a 1D contiguous array.
            double precision, allocatable, intent(in) :: dat2D(:,:,:) !! Multi-variable 2D array: dat2D(m_var, i, j)
            double precision, allocatable, intent(inout) :: buf(:) !! Buffer (flattened) 1D array for the mpi send transmission
            integer, intent(in) :: jsend !! Slice index of the sending rank
            integer, intent(in) :: m_var !! Number of field variables in 'dat2D'
            integer, intent(in) :: ilo, ihi !! Rank-specific index boundaries of the i-slice
            type(MPI_Comm), intent(in) :: comm !! MPI Communicator
        end subroutine

        module subroutine unpack_2d_jslice(dat2D, recvbuf, irecv, m_var, jlo, jhi, comm)
            !! Reverses the flattening of 'dat2D' into a 1D array and unpacks the 1D array
            !! 'recbuf' back into 'dat2D'.
            double precision, allocatable, intent(inout) :: dat2D(:,:,:) !! Multi-variable 2D array: dat2D(m_var, i, j)
            double precision, allocatable, intent(in) :: recvbuf(:) !! Flattened slice data
            integer, intent(in) :: irecv !! Slice index in recieving rank
            integer, intent(in) :: m_var !! Number of field variables in 'dat2D'
            integer, intent(in) :: jlo, jhi !! Rank-specific index boundaries of the j-slice
            type(MPI_Comm), intent(in) :: comm !! MPI Communicator
        end subroutine

        module subroutine unpack_2d_islice(dat2D, recvbuf, jrecv, m_var, ilo, ihi, comm)
            !! Reverses the flattening of 'dat2D' into a 1D array and unpacks the 1D array
            !! 'recbuf' back into 'dat2D'.
            double precision, allocatable, intent(inout) :: dat2D(:,:,:) !! Multi-variable 2D array: dat2D(m_var, i, j)
            double precision, allocatable, intent(in) :: recvbuf(:) !! Flattened slice data
            integer, intent(in) :: jrecv !! Slice index in recieving rank
            integer, intent(in) :: m_var !! Number of field variables in 'dat2D'
            integer, intent(in) :: ilo, ihi !! Rank-specific index boundaries of the i-slice
            type(MPI_Comm), intent(in) :: comm !! MPI Communicator
        end subroutine
    end interface

contains

end module MOD_MPI_decomposition
