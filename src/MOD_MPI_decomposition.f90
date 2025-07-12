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
    !! Block domain decomposition for MPI-based structured grid solvers in 1, 2 and 3 
    !! spatial dimensions.
    !!
    !! This library handles a static block decomposition using MPI. Each MPI rank
    !! owns a rectangular subdomain `(i)`, `(i,j)` or `(i,j,k)` of the global mesh. Ghost cells
    !! are not included in the decomposition index ranges and must be accounted for externally.
    !! Each rank is assigned a portion of the global index range. These ensure
    !! continuous field indizes across ranks without the need to calculate a global
    !! index.
    !!
    !! The decomposition info is stored in a type `decomp_info`.
    !!
    !! Setup might be updated some day to also include user defined sliced layout. Currently
    !! only (automatic) chunked decomposition is supported.
    !!
    !! All sweep directions (xi, eta and tau) operate over the same local block.
    !!
    !! Features to be added: load balancing? more detailed tests for allocation status etc.

    use mpi_f08
    implicit none

    type :: decomp_info
        !! Contains the domain decomposition information, such as
        !! indexing bounds and rank neighbours.
        integer :: ilow, ihigh                                  !! Local i-direction bounds
        integer :: jlow, jhigh                                  !! Local j-direction bounds
        integer :: klow, khigh                                  !! Local k-direction bounds

        integer :: nbr_top                                      !! neighbouring ranks
        integer :: nbr_bottom                                   !! neighbouring ranks
        integer :: nbr_left                                     !! neighbouring ranks
        integer :: nbr_right                                    !! neighbouring ranks
        integer :: nbr_front                                    !! neighbouring ranks
        integer :: nbr_back                                     !! neighbouring ranks
        type(MPI_Comm) :: comm_cart                             !! cartesian communicator

        integer :: dims(3)                                      !! Processor grid in i, j and k directions
    contains
        !initialization:
        procedure :: initialize_decomposition
        procedure :: setup_cartesian_topology
        !console output:
        procedure :: print_decomposition_summary
        procedure :: print_cartesian_rank_layout
        !output to disk:
        procedure :: write_decom_to_disk
        !getters:
        procedure :: get_local_block_bounds
        procedure :: get_cartesian_comm
        procedure :: get_neighbouring_ranks
    end type


    !=======================================================================================================!
    ! setup and initialization
    !=======================================================================================================!
    interface
        module subroutine initialize_decomposition(info, m_xi, m_eta, m_tau, comm)
            !! Initializes the domain decomposition for a N-dimensional grid using MPI.
            !!
            !! The decomposition supports 1D, 2D, and 3D domain partitioning depending on the
            !! non-zero inputs among m_xi, m_eta, and m_tau. A zero value for a dimension
            !! disables decomposition in that direction.
            !!
            !! The decomposition is performed only by the master rank (rank 0), which computes
            !! an N-D processor grid decomposition of the global domain into chunks.
            !! Each chunk corresponds to a rectangular subdomain of the global mesh.
            !!
            !! Rank 0 then sends the computed local index ranges (bounds) to all other ranks
            !! using point-to-point MPI communication. Each rank receives its assigned local
            !! ranges for all active directions.
            !!
            !! This approach avoids all ranks calling MPI_Dims_create, thus preventing
            !! potential segmentation faults and ensuring consistent decomposition.
            !!
            !! The subroutine uses a helper `compute_range` to divide each dimension evenly
            !! among the processors assigned along that dimension, distributing any remainder
            !! cells to the lowest ranks.
            !!
            !! Input:
            !! - m_xi: Total number of cells in the i-direction (xi-axis). 0 = inactive.
            !! - m_eta: Total number of cells in the j-direction (eta-axis). 0 = inactive.
            !! - m_tau: Total number of cells in the k-direction (tau-axis). 0 = inactive.
            !! - comm: MPI communicator over which ranks are defined.
            !!
            !! Output (module private variables set):
            !! - ilow, ihigh: Local i-direction bounds (inclusive).
            !! - jlow, jhigh: Local j-direction bounds (inclusive).
            !! - klow, khigh: Local k-direction bounds (inclusive), if active.
            class(decomp_info), intent(out) :: info        !! Object containing the decomposition info
            integer, intent(in) :: m_xi                    !! Total number of cells in the i-direction (xi-axis)
            integer, intent(in) :: m_eta                   !! Total number of cells in the j-direction (eta-axis)
            integer, intent(in) :: m_tau                   !! Total number of cells in the k-direction (tau-axis)
            type(MPI_Comm), intent(in) :: comm             !! MPI communicator (normally 'MPI_COMM_WORLD')
        end subroutine

        module subroutine setup_cartesian_topology(info, comm_in)
            !! Create a 1D, 2D, or 3D Cartesian communicator and get neighbor ranks.
            !!
            !! The number of processes in each dimension is taken from previous
            !! calculation executed in initialize_decomposition. Sets the global
            !! module variable 'comm_cart' and neighbor ranks.
            !!
            !! Neighbor ranks include:
            !! - nbr_left/nbr_right (i-direction)
            !! - nbr_bottom/nbr_top (j-direction)
            !! - nbr_back/nbr_front (k-direction), if active

            class(decomp_info) :: info               !! Object containing the decomposition info
            type(MPI_Comm), intent(in) :: comm_in    !! Input communicator (usually MPI_COMM_WORLD)
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
            !! table on rank 0 showing the (i, j, k) bounds for each rank's block.
            !! Only active dimensions (m_xi, m_eta, m_tau) are shown.
            !!
            !! Output is dimension-aware and adjusts formatting accordingly.

            class(decomp_info) :: info          !! Object containing the decomposition info
            type(MPI_Comm), intent(in) :: comm  !! MPI communicator (usually MPI_COMM_WORLD)
        end subroutine

        module subroutine print_cartesian_rank_layout(info)
            !! Prints a visual representation of the rank layout in the form of a chessboard.
            class(decomp_info), intent(in) :: info !! Object containing the decomposition info
        end subroutine

        module subroutine write_decom_to_disk(info,comm)
            class(decomp_info), intent(in) :: info
            type(MPI_Comm), intent(in) :: comm
        end subroutine
    end interface

    !=======================================================================================================!
    ! type-bound getters
    !=======================================================================================================!
    interface
        module subroutine get_local_block_bounds(info, ilo, ihi, jlo, jhi, klo, khi)
            !! Returns the local 2D index bounds for this rank.
            !!
            !! These bounds define the owned rectangular region in the global mesh for
            !! all directional sweeps. Ghost cells must be added externally by the solver.
            class(decomp_info), intent(in) :: info  !! Object containing the decomposition info
            integer, intent(out) :: ilo             !! local indexing range
            integer, intent(out) :: ihi             !! local indexing range
            integer, intent(out), optional :: jlo   !! local indexing range
            integer, intent(out), optional :: jhi   !! local indexing range
            integer, intent(out), optional :: klo   !! local indexing range
            integer, intent(out), optional :: khi   !! local indexing range
        end subroutine

        module function get_cartesian_comm(info) result(comm)
            !! Returns the stored Cartesian communicator
            class(decomp_info), intent(in) :: info  !! Object containing the decomposition info
            type(MPI_Comm) :: comm                  !! MPI cartesian communicator
        end function

        module subroutine get_neighbouring_ranks(info, left, right, bottom, top, back, front)
            !! Returns the ranks of neighboring processes
            class(decomp_info), intent(in) :: info      !! Object containing the decomposition info
            integer, intent(out) :: left                !! Neighbouring ranks
            integer, intent(out) :: right               !! Neighbouring ranks
            integer, intent(out), optional :: bottom    !! Neighbouring ranks
            integer, intent(out), optional :: top       !! Neighbouring ranks
            integer, intent(out), optional :: back      !! Neighbouring ranks
            integer, intent(out), optional :: front     !! Neighbouring ranks
        end subroutine
    end interface

    !=======================================================================================================!
    ! halo exchange
    !=======================================================================================================!
    interface exchange_halos
        module procedure :: exchange_halos_1d
        module procedure :: exchange_halos_2d
        module procedure :: exchange_halos_3d
    end interface

    interface
        module subroutine exchange_halos_1d(info, dat1D, m_var, num_ghost, left, right, comm)
            !! Exchanges halo layers for a multi-variable 1D structured field.
            !!
            !! Performs halo exchange on a field `dat2D(m_var, i)` with ghost cells,
            !! communicating ghost layers with direct MPI neighbors: left and right.
            class(decomp_info), intent(in) :: info
            !! Object containing the decomposition info
            double precision, allocatable :: dat1D(:,:)
            !! Multi-variable 1D array: dat1D(m_var, i)
            integer, intent(in) :: m_var        !! Number of variables per grid point
            integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
            integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
            type(MPI_Comm), intent(in) :: comm  !! MPI communicator
        end subroutine

        module subroutine exchange_halos_2d(info, dat2D, m_var, num_ghost, left, right, bottom, top, comm)
            !! Exchanges halo layers for a multi-variable 2D structured field.
            !!
            !! Performs halo exchange on a field `dat2D(m_var, i, j)` with ghost cells,
            !! communicating ghost layers with direct MPI neighbors: left, right, bottom, and top.
            class(decomp_info), intent(in) :: info
            !! Object containing the decomposition info
            double precision, allocatable :: dat2D(:,:,:)
            !! Multi-variable 2D array: `dat2D(m_var, i, j)`
            integer, intent(in) :: m_var        !! Number of variables per grid point
            integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
            integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
            type(MPI_Comm), intent(in) :: comm  !! MPI communicator
        end subroutine

        module subroutine exchange_halos_3d(info, dat3D, m_var, num_ghost, left, right, bottom, top, back, front, comm)
            !! Exchanges halo layers for a multi-variable 3D structured field.
            !!
            !! Performs halo exchange on a field `dat3D(m_var, i, j, k)` with ghost cells,
            !! communicating ghost layers with direct MPI neighbors: left, right, bottom, top, back and front.
            class(decomp_info), intent(in) :: info
            !! Object containing the decomposition info
            double precision, allocatable :: dat3D(:,:,:,:)
            !! Multi-variable 2D array: `dat3D(m_var, i, j,k)`
            integer, intent(in) :: m_var        !! Number of variables per grid point
            integer, intent(in) :: num_ghost    !! Number of ghost cells on each side
            integer, intent(in) :: left         !! MPI rank of left neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: right        !! MPI rank of right neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: bottom       !! MPI rank of bottom neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: top          !! MPI rank of top neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: back         !! MPI rank of back neighbor (or MPI_PROC_NULL)
            integer, intent(in) :: front        !! MPI rank of front neighbor (or MPI_PROC_NULL)
            type(MPI_Comm), intent(in) :: comm  !! MPI communicator
        end subroutine
    end interface
    !=======================================================================================================!
    ! sending the ghost indizes
    !=======================================================================================================!
    interface send_i_index_ghost
        module procedure :: send_1d_i_index_ghost
        module procedure :: send_2d_i_index_ghost
        module procedure :: send_3d_i_index_ghost
    end interface

    interface send_j_index_ghost
        module procedure :: send_2d_j_index_ghost
        module procedure :: send_3d_j_index_ghost
    end interface

    interface send_k_index_ghost
        module procedure :: send_3d_k_index_ghost
    end interface

    interface
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
        end subroutine
    end interface

    !=======================================================================================================!
    ! rank communication subroutines
    !=======================================================================================================!
    interface send_rec_1D_array
        !! Generic interface, so that the subroutine call be called with either
        !! integer or double precision arrays.
        module procedure send_recv_1D_array_dp  !! Double precision version of in-/output arrays
        module procedure send_recv_1D_array_int !! Integer version of in-/output arrays
    end interface

    interface find_min_max_scalar
        !! Generic interface, so that the subroutine call be called with either
        !! integer or double precision arrays.
        module procedure find_min_max_scalar_dp !! Double precision version of in-/output scalar
        module procedure find_min_max_scalar_int!! Integer version of in-/output scalar
    end interface

    interface
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
        end subroutine
    end interface

    !=======================================================================================================!
    ! helper subroutines
    !=======================================================================================================!
    interface allocate_buffers
        module procedure :: allocate_buffers_1d
        module procedure :: allocate_buffers_2d
        module procedure :: allocate_buffers_3d
    end interface

    interface pack_i_index
        module procedure :: pack_1d_i_index
        module procedure :: pack_2d_i_index
        module procedure :: pack_3d_i_index
    end interface

    interface pack_j_index
        module procedure :: pack_2d_j_index
        module procedure :: pack_3d_j_index
    end interface

    interface unpack_i_index
        module procedure :: unpack_1d_i_index
        module procedure :: unpack_2d_i_index
        module procedure :: unpack_3d_i_index
    end interface

    interface unpack_j_index
        module procedure :: unpack_2d_j_index
        module procedure :: unpack_3d_j_index
    end interface

    interface pack_k_index
        module procedure :: pack_3d_k_index
    end interface

    interface unpack_k_index
        module procedure :: unpack_3d_jkindex
    end interface

    interface
        module subroutine allocate_buffers_1d(sendbuf, recvbuf, m_var)
            !! Allocates the send and recieve buffer arrays.
            !!
            !! Takes the two bounds `idxlo` and `idxhi` as lower and upper bounds of one slice of data.
            !! Then `m_var` is the number of slices (= number of field variables).
            double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
            double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
            integer, intent(in) :: m_var                                !! Number of field variables
        end subroutine

        module subroutine allocate_buffers_2d(sendbuf, recvbuf, idxlo, idxhi, m_var)
            !! Allocates the send and recieve buffer arrays.
            !!
            !! Takes the two bounds `idxlo` and `idxhi` as lower and upper bounds of one slice of data.
            !! Then `m_var` is the number of slices (= number of field variables).
            double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
            double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
            integer, intent(in) :: idxlo                                !! Bounds of array slice
            integer, intent(in) :: idxhi                                !! Bounds of array slice
            integer, intent(in) :: m_var                                !! Number of field variables
        end subroutine

        module subroutine allocate_buffers_3d(sendbuf, recvbuf, idxlo_i, idxhi_i, idxlo_j, idxhi_j, m_var)
            !! Allocates the send and recieve buffer arrays.
            !!
            !! Takes the two bounds `idxlo_i` and `idxhi_i` as lower and upper bounds of first slicing direction of data.
            !! Takes the two bounds `idxlo_j` and `idxhi_j` as lower and upper bounds of the second direction of the data.
            !! Then `m_var` is the number of slices (= number of field variables).
            double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
            double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
            integer, intent(in) :: idxlo_i                              !! Bounds of array slice
            integer, intent(in) :: idxhi_i                              !! Bounds of array slice
            integer, intent(in) :: idxlo_j                              !! Bounds of array slice
            integer, intent(in) :: idxhi_j                              !! Bounds of array slice
            integer, intent(in) :: m_var                                !! Number of field variables
        end subroutine

        module subroutine deallocate_buffers(sendbuf, recvbuf)
            !! DEallocates the send and recieve buffer arrays.
            double precision, allocatable, intent(inout) :: sendbuf(:)  !! Send buffer
            double precision, allocatable, intent(inout) :: recvbuf(:)  !! Recieve buffer
        end subroutine

        module subroutine pack_1d_i_index(dat1D, buf, isend, m_var, comm)
            !! Compiles a 1D array into which all the data from `dat2D(1:mvar,isend,jlo:jhi)` is packed.
            !!
            !! `dat1D` is flattened into a 1D contiguous array.
            double precision, allocatable, intent(in) :: dat1D(:,:)     !! Multi-variable 1D array: `dat1D(m_var, i)`
            double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
            integer, intent(in) :: isend                                !! Slice index of the sending rank
            integer, intent(in) :: m_var                                !! Number of field variables in `dat1D`
            type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator
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
        end subroutine

        module subroutine pack_3d_j_index(dat3D, buf, jsend, m_var, ilo, ihi, klo, khi, comm)
            !! Compiles a 1D array into which all the data from `dat2D(1:mvar,ilo:ihi,jsend,klo:khi)` is packed.
            !!
            !! `dat3D` is flattened into a 1D contiguous array.
            double precision, allocatable, intent(in) :: dat3D(:,:,:,:)   !! Multi-variable 2D array: `dat3D(m_var, i, j, k)`
            double precision, allocatable, intent(inout) :: buf(:)      !! Buffer (flattened) 1D array for the mpi send transmission
            integer, intent(in) :: jsend                                !! Slice index of the sending rank
            integer, intent(in) :: m_var                                !! Number of field variables in `dat3D`
            integer, intent(in) :: ilo                                  !! Rank-specific index boundaries of the i-slice
            integer, intent(in) :: ihi                                  !! Rank-specific index boundaries of the i-slice
            integer, intent(in) :: klo                                  !! Rank-specific index boundaries of the k-slice
            integer, intent(in) :: khi                                  !! Rank-specific index boundaries of the k-slice
            type(MPI_Comm), intent(in) :: comm                          !! MPI Communicator
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
        end subroutine

        module subroutine unpack_1d_i_index(dat1D, recvbuf, irecv, m_var, comm)
            !! Reverses the flattening of `dat1D` into a 1D array and unpacks the 1D array
            !! `recbuf` back into `dat1D`.
            double precision, allocatable, intent(inout) :: dat1D(:,:)      !! Multi-variable 1D array: `dat1D(m_var, i)`
            double precision, allocatable, intent(in) :: recvbuf(:)         !! Flattened slice data
            integer, intent(in) :: irecv                                    !! Slice index in recieving rank
            integer, intent(in) :: m_var                                    !! Number of field variables in `dat1D`
            type(MPI_Comm), intent(in) :: comm                              !! MPI Communicator
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
        end subroutine
    end interface

contains

end module MOD_MPI_decomposition
