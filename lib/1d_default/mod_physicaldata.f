module mod_physicaldata
   implicit none
   save

   type walloc
      !> ID of a grid block
      integer :: igrid=-1
      !> location of w0-array, 0: cell center, ^D : cell interface in dimension ^D
      integer :: iw0=0
      !> Is e is internal energy or total energy
      logical :: e_is_internal=.false.
      !> If it face a physical boundary
      logical :: is_physical_boundary(2*1)
      !> Variables, normally cell center conservative values
      double precision, dimension(:,:), allocatable :: w
      !> Variables, normally center, temporary state for multi-step scheme
      double precision, dimension(:,:), allocatable :: w1
      !> Variables, normally center, temporary state for multi-step scheme
      double precision, dimension(:,:), allocatable :: w2
      !> Variables, normally center, temporary state for multi-step scheme
      double precision, dimension(:,:), allocatable :: w3
      !> Variables, normally center, temporary state for multi-step scheme
      double precision, dimension(:,:), allocatable :: w4
      !> Variables, normally center, pointer of reference state
      double precision, dimension(:,:), pointer :: wa => Null()
      !> Variables, normally center, pointer of updated state
      double precision, dimension(:,:), pointer :: wb => Null()
      !> Variables, normally center, initial state at the beginning of iteration
      double precision, dimension(:,:), allocatable :: wold
      !> Variables, normally center, for visualization data
      double precision, dimension(:,:), allocatable :: wio
      !> Variables, normally center, one level coarser representation
      double precision, dimension(:,:), allocatable :: wcoarse
      !> Time-independent magnetic field at cell center and cell interface
      double precision, dimension(:,:,:), allocatable :: B0
      !> Time-independent electric current density at cell center
      double precision, dimension(:,:), allocatable :: J0
      !> Cell-center positions
      double precision, dimension(:,:), allocatable :: x
      !> Cell-center positions, one level coarser 
      double precision, dimension(:,:), allocatable :: xcoarse
      !> Cell sizes
      double precision, dimension(:,:), allocatable :: dx
      !> Cell sizes, one level coarser
      double precision, dimension(:,:), allocatable :: dxcoarse
      !> Volumes of a cell
      double precision, dimension(:), allocatable :: dvolume
      !> Volumes of a cell, one level coarser representative
      double precision, dimension(:), allocatable :: dvolumecoarse
      !> Pointers below are used to switch between fine/coarse 
      !> Volumes of a cell, pointer
      double precision, dimension(:), pointer :: dvolumep => Null()
      !> Cell-center positions, pointer
      double precision, dimension(:,:), pointer :: xCCp => Null()
      !> Cell sizes, pointer
      double precision, dimension(:,:), pointer :: dxp => Null()
      !> Areas of cell-center surfaces
      double precision, dimension(:,:), allocatable :: surface
      !> Areas of cell-face surfaces
      double precision, dimension(:,:), allocatable :: surfaceC
   end type walloc



   type walloc_sub
      !> ID of a grid block
      integer :: igrid=-1
      !> Is w in primitive state or not
      logical :: e_is_internal=.false.
      !> Variables, normally center
      double precision, dimension(:), allocatable :: w
      !> Variables for the cornerpositions on the slice 
      double precision, dimension(:), allocatable :: wC
      !> Variables, normally center, one level coarser representative
      double precision, dimension(:), allocatable :: wcoarse
      !> Cell-center positions
      double precision, dimension(:), allocatable :: x
      !> Corner positions on the slice
      double precision, dimension(:), allocatable :: xC
      !> Cell-center positions, one level coarser representative
      double precision, dimension(:), allocatable :: xcoarse
   end type walloc_sub

   type grid_field
      !> Variables new state
      double precision, dimension(:,:), allocatable :: w
      !> Variables old state
      double precision, dimension(:,:), allocatable :: wold
   end type grid_field
   !> Block pointer for using current block
   type(walloc), pointer :: block

   !> array of physical blocks
   type(walloc), dimension(:), allocatable, target :: pw

   !> array of physical blocks in reduced dimension
   type(walloc_sub), dimension(:), allocatable, target :: pw_sub


   double precision, dimension(:), allocatable :: collapsedData


   !> array of physical blocks of meshed fields for particles
   type(grid_field), dimension(:), allocatable, target :: gridvars

end module mod_physicaldata
