module oe_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit

implicit none

private
public :: block_object
public :: cc_par_object
public :: box_object
public :: bc_int_type, bc_string
public :: create_blocks_list
public :: implode_blocks
public :: load_file_cc_par
public :: load_file_grd
public :: replay_splits_on_patches
public :: save_file_cc_par
public :: load_file_icc
public :: popout_blocks_list
public :: save_file_grd
public :: save_file_icc
public :: save_proc_input
public :: update_blocks

! BC parameters
! natural BC
integer(I4P), parameter, public :: BC_NATURAL_WALL                     = -1   !< Viscous wall.
integer(I4P), parameter, public :: BC_NATURAL_SIMMETRY                 = -2   !< Simmetry.
integer(I4P), parameter, public :: BC_NATURAL_INFLOW                   = -3   !< Inflow.
integer(I4P), parameter, public :: BC_NATURAL_INOUTFLOW                = -4   !< In/outflow.
integer(I4P), parameter, public :: BC_NATURAL_ASSIGNED_INFLOW          = -5   !< Assigned inflow.
integer(I4P), parameter, public :: BC_NATURAL_ASSIGNED_PRESSURE        = -6   !< Assigned pressure.
integer(I4P), parameter, public :: BC_NATURAL_ASSIGNED_NORMAL_VELOCITY = -7   !< Assigned normal velocity.
integer(I4P), parameter, public :: BC_NATURAL_ASSIGNED_RIEMANN         = -8   !< Assigned Riemann invariant.
integer(I4P), parameter, public :: BC_NATURAL_EXTRAPOLATED             = -9   !< Extrapolated.
integer(I4P), parameter, public :: BC_NATURAL_MOVING_WALL              = -10  !< Moving wall.
integer(I4P), parameter, public :: BC_NATURAL_INACTIVE_WALL            = -11  !< Inactive wall.
integer(I4P), parameter, public :: BC_NATURAL_EXTRAPOLATED_ALT         = -19  !< Extrapolated (alternative).
integer(I4P), parameter, public :: BC_NATURAL_RCC_RESERVED_DATA        =  19  !< RCC data reserved for BC natural.
! non BC, active cell
integer(I4P), parameter, public :: BC_ACTIVE_CELL                      = 0    !< Non BC, active cell.
! chimera BC, face face-center
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF                  = 20   !< Chimera face.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_I0               = 21   !< Chimera face, centered at i0 face-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_IN               = 22   !< Chimera face, centered at in face-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_J0               = 23   !< Chimera face, centered at j0 face-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_JN               = 24   !< Chimera face, centered at jn face-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_K0               = 25   !< Chimera face, centered at k0 face-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XF_KN               = 26   !< Chimera face, centered at kn face-center.
! chimera BC, cell
integer(I4P), parameter, public :: BC_CHIMERA_CELL                     = 27   !< Chimera cell inside domain.
integer(I4P), parameter, public :: BC_CHIMERA_CELL_INT_WALL            = 28   !< Chimera cell internal wall.
! chimera BC, face cell-center
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC                  = 40   !< Chimera face.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_I0               = 41   !< Chimera face, centered at i0 cell-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_IN               = 42   !< Chimera face, centered at in cell-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_J0               = 43   !< Chimera face, centered at j0 cell-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_JN               = 44   !< Chimera face, centered at jn cell-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_K0               = 45   !< Chimera face, centered at k0 cell-center.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_XC_KN               = 46   !< Chimera face, centered at kn cell-center.
! chimera BC, adjacent
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ                 = 60   !< Adjacent.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_I0              = 61   !< Adjacent along face i0.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_IN              = 62   !< Adjacent along face in.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_J0              = 63   !< Adjacent along face j0.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_JN              = 64   !< Adjacent along face jn.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_K0              = 65   !< Adjacent along face k0.
integer(I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_KN              = 66   !< Adjacent along face kn.
! chimera BC, edge
integer(I4P), parameter, public :: BC_CHIMERA_EDGE                     = 80   !< Edge.

integer(I4P), parameter :: MAX_DONORS=8 !< Maximum number of donors for chimera-like BC.

type :: box_object
   !< Box class.
   integer(I4P) :: btype          !< Box type.
   integer(I4P) :: bblock         !< Box block.
   integer(I4P) :: bgroup         !< Box group.
   real(R8P)    :: nodes(1:3,1:8) !< Box nodes coordinates.
endtype box_object

type :: cc_par_object
   !< File cc.par handler class.
   character(99) :: file_name_input_grd          !< File name of input grd.
   character(99) :: base_name_output             !< File name of output grd.
   logical       :: save_ghost_cells=.true.      !< Save grid with ghost cells.
   logical       :: increase_overlap=.false.     !< Increase overlap between grids.
   logical       :: extend_internal_wall=.false. !< Extend internal wall inside chimera cells.
   integer(I4P)  :: mgl(1:2)                     !< Multi grid levels.
   real(R8P)     :: boundary_layer_thickness     !< Boundary layer thickness.
   real(R8P)     :: numberical_beach             !< Numerical beach extension.
   integer(I4P)  :: blocks_number=0              !< Blocks number.
   integer(I4P)  :: patches_number=0             !< Patches number.
   integer(I4P)  :: boxes_number=0               !< Boxes number.
   integer(I4P)  :: edges_number=0               !< Edges number.
   integer(I4P)  :: circuits_number=0            !< Circuits number.
endtype cc_par_object

type :: patch_object
  !< Patch class.
  logical      :: is_connection=.false. !< Flag to inquire is the patch defines a connection.
  integer(I4P) :: patch_index=0         !< Patch index, local numeration.
  integer(I4P) :: block_index=0         !< Block to which patch belongs, local numeration.
  integer(I4P) :: face_index=0          !< Face index in the overset convention: Imin=>1, Imax=2,Jmin=3, Jmax=4, Kmin=5, Kmax=6.
  integer(I4P) :: boundary_condition=0  !< Boundary condition (or IJK orientation) in the overset convention.
  integer(I4P) :: connect_family=0      !< Index of connected patch or family index of the patch.
  integer(I4P) :: ijk_extents(1:6)=&
                  [0,0,0,0,0,0]         !< IJK extents of the patch in the overset convention: Imin, Imax, Jmin, Jmax, Kmin, Kmax.
  character(99):: comment               !< Patch end-line comment (contain the global numeration index of patch).
endtype patch_object

type :: block_object
   !< Block class.
   integer(I4P)                    :: Ni=0              !< Number of cells in i direction.
   integer(I4P)                    :: Nj=0              !< Number of cells in j direction.
   integer(I4P)                    :: Nk=0              !< Number of cells in k direction.
   integer(I4P)                    :: gc=2              !< Number of ghost cells.
   integer(I4P)                    :: w=0               !< Block weight (work load).
   real(R8P),          allocatable :: nodes(:,:,:,:)    !< Nodes coordinates.
   integer(I4P),       allocatable :: icc(:,:,:)        !< Cell centered icc values.
   integer(I4P),       allocatable :: tcc(:,:,:,:)      !< BC type and index on chimera values [1:2,1-gc:ni,1-gc:nj,1-gc:nk].
   real(R4P),          allocatable :: chimera(:)        !< Chimera values (donors number, bijk-weight for each donor) [1:nchimera].
   type(patch_object), allocatable :: patches(:)        !< Patches on block faces (6 patches or more if faces have split patches).
   integer(I4P)                    :: ab=0              !< Absolute block index.
   integer(I4P)                    :: group=0           !< Index of gruop.
   integer(I4P)                    :: body=0            !< Index of body.
   integer(I4P)                    :: level=0           !< Level for overset algorithms.
   integer(I4P)                    :: priority=0        !< Priority for overset algorithms.
   character(99)                   :: comment           !< Block comment in cc.par.
   integer(I4P)                    :: proc=0            !< Processor assigned to.
   logical                         :: is_loaded=.false. !< Flag for checking if the block is loaded.
   ! splitting data
   integer(I4P), allocatable :: parents(:)        !< List of parents blocks.
   integer(I4P)              :: split_level=0     !< Split level, 0 for original block.
   integer(I4P)              :: split_dir=0       !< Split direction, if = 0 no split.
   contains
      ! public methods
      procedure, pass(self) :: destroy            !< Destroy dynamic memory.
      procedure, pass(self) :: alloc              !< Allocate dynamic memory.
      procedure, pass(self) :: load_dimensions    !< Load block dimensions from file.
      procedure, pass(self) :: load_icc           !< Load block icc from file.
      procedure, pass(self) :: load_nodes         !< Load block nodes from file.
      procedure, pass(self) :: parse_rcc          !< Parse global rcc and store in local tcc/chimera arrays.
      procedure, pass(self) :: sanitize_chimera   !< Sanitize chimera data after a split.
      procedure, pass(self) :: save_block_file    !< Save block data into its own file.
      procedure, pass(self) :: split              !< Split block.
      procedure, pass(self) :: weight             !< Return block weight (work load).
      generic :: assignment(=) => assign_block    !< Assignment operator overloading.
      ! private methods
      procedure, pass(lhs), private :: assign_block !< Operator `=`.
endtype block_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   self%w  = 0
   if (allocated(self%nodes)) deallocate(self%nodes)
   if (allocated(self%icc)) deallocate(self%icc)
   if (allocated(self%tcc)) deallocate(self%tcc)
   if (allocated(self%chimera)) deallocate(self%chimera)
   if (allocated(self%patches)) deallocate(self%patches)
   self%ab        = 0
   self%group     = 0
   self%body      = 0
   self%proc      = 0
   self%level     = 0
   self%priority  = 0
   self%comment   = ''
   self%is_loaded = .false.
   if (allocated(self%parents)) deallocate(self%parents)
   self%split_level = 0
   self%split_dir   = 0
   endsubroutine destroy

   elemental subroutine alloc(self)
   !< Allocate dynamic memory.
   class(block_object), intent(inout) :: self !< Block data.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc)
   allocate(self%nodes(1:3,0-gc:Ni+gc,0-gc:Nj+gc,0-gc:Nk+gc))
   allocate(self%icc(1-gc:Ni+gc,1-gc:Nj+gc,1-gc:Nk+gc))
   allocate(self%tcc(1:2,1-gc:Ni+gc,1-gc:Nj+gc,1-gc:Nk+gc))
   endassociate
   endsubroutine alloc

   subroutine load_dimensions(self, file_unit, ab, allocate_data, gc)
   !< Load block dimensions from file.
   !<
   !< @note The file must be already open and the current record-index must be at the proper block dimensions record.
   !< If ab index is passed it is supposed that this the first time the dimensions is loaded from grd file, thus the block
   !< destroyed and allocated ex novo.
   class(block_object), intent(inout)        :: self          !< Block data.
   integer(I4P),        intent(in)           :: file_unit     !< Logical unit of grd file.
   integer(I4P),        intent(in), optional :: ab            !< Absolute block index.
   logical,             intent(in), optional :: allocate_data !< Sentinel to trigger data allocation.
   integer(I4P),        intent(in), optional :: gc            !< Number of ghost cells.

   if (present(ab)) call self%destroy
   if (present(gc)) self%gc = gc
   read(file_unit, end=10, err=10) self%Ni,self%Nj,self%Nk,self%gc
   10 continue
   if (present(ab)) then
      self%ab = ab
      self%w = self%weight()
   endif
   if (present(allocate_data)) then
      if (allocate_data) call self%alloc
   endif
   endsubroutine load_dimensions

   subroutine load_icc(self, file_unit)
   !< Load block icc from file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block icc record.
   class(block_object), intent(inout) :: self      !< Block data.
   integer(I4P),        intent(in)    :: file_unit !< Logical unit of icc file.
   integer(I4P)                       :: i,j,k     !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,icc=>self%icc)
   read(file_unit)(((icc(i,j,k),i=1-gc,Ni+gc),j=1-gc,Nj+gc),k=1-gc,Nk+gc)
   endassociate
   endsubroutine load_icc

   subroutine load_nodes(self, file_unit)
   !< Load block nodes from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block nodes record.
   class(block_object), intent(inout) :: self      !< Block data.
   integer(I4P),        intent(in)    :: file_unit !< Logical unit of grd file.
   integer(I4P)                       :: i,j,k     !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes)
      read(file_unit)(((nodes(1,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
      read(file_unit)(((nodes(2,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
      read(file_unit)(((nodes(3,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   endassociate
   endsubroutine load_nodes

   subroutine parse_rcc(self, rcc)
   !< Parse global rcc and store in local tcc/chimera arrays.
   class(block_object), intent(inout) :: self                 !< Block data.
   real(R4P),           intent(in)    :: rcc(1:)              !< rcc unstructured array.
   integer(I4P)                       :: nchimera             !< Number of chimera data.
   integer(I4P)                       :: ndonors              !< Number of donors.
   real(R4P)                          :: donors(5,MAX_DONORS) !< Chimera-like BC donors data.
   integer(I4P)                       :: i,j,k,n,p            !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,icc=>self%icc)
   if (allocated(self%chimera)) deallocate(self%chimera)
   self%tcc(1,:,:,:) = BC_ACTIVE_CELL ! set all cell type to active cell
   self%tcc(2,:,:,:) = 0_I4P          ! set all cell pointer-to-chimera to 0
   ! count chimera rcc elements of block and set blocks%tcc
   nchimera = 0
   do k=1-gc, Nk+gc
   do j=1-gc, Nj+gc
   do i=1-gc, Ni+gc
      p = icc(i,j,k)
      if     (p>BC_NATURAL_RCC_RESERVED_DATA) then ! chimera-like BC
         nchimera = nchimera + 1
         call get_donors(p=p, rcc=rcc, ndonors=ndonors)
         self%tcc(1,i,j,k) = nint(rcc(p),I4P)
         self%tcc(2,i,j,k) = nchimera
         nchimera = nchimera + ndonors*5 ! b,i,j,k,weight for each donor
      elseif (p>0.and.p<=BC_NATURAL_RCC_RESERVED_DATA) then ! natural BC
         self%tcc(1,i,j,k) = nint(rcc(p),I4P)
      endif
   enddo
   enddo
   enddo
   if (nchimera>0) then
      ! set blocks%chimera
      allocate(self%chimera(1:nchimera))
      self%chimera = 0._R4P
      nchimera = 0
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         p = icc(i,j,k)
         if     (p>BC_NATURAL_RCC_RESERVED_DATA) then  ! chimera-like BC
            nchimera = nchimera + 1
            call get_donors(p=p, rcc=rcc, ndonors=ndonors, donors=donors)
            self%chimera(nchimera) = real(ndonors,R4P)
            do n=1, ndonors
               self%chimera(nchimera+1+5*(n-1)) = donors(1,n) ! b
               self%chimera(nchimera+2+5*(n-1)) = donors(2,n) ! i
               self%chimera(nchimera+3+5*(n-1)) = donors(3,n) ! j
               self%chimera(nchimera+4+5*(n-1)) = donors(4,n) ! k
               self%chimera(nchimera+5+5*(n-1)) = donors(5,n) ! weight
            enddo
            nchimera = nchimera + ndonors*5 ! b,i,j,k,weight for each donor
         endif
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine parse_rcc

   pure subroutine sanitize_chimera(self, sb)
   !< Sanitize chimera data after a split.
   class(block_object), intent(inout) :: self        !< Block data.
   type(block_object),  intent(in)    :: sb(2)       !< Split blocks.
   integer(I4P)                       :: i,j,k,p,n   !< Counter.
   integer(I4P)                       :: db,di,dj,dk !< Donor indexes.
   integer(I4P)                       :: ndonors     !< Chimera donors number.
   integer(I4P)                       :: o           !< Offset of rcc.

   if (.not.allocated(self%tcc)) return
   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes,icc=>self%icc,tcc=>self%tcc,chimera=>self%chimera)
   do k=1-gc, Nk+gc
   do j=1-gc, Nj+gc
   do i=1-gc, Ni+gc
      select case(tcc(1,i,j,k))
      case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
         p = tcc(2,i,j,k)
         ndonors = nint(chimera(p),I4P)
         do n=1, ndonors
            o = p + 5*(n-1)
            db = nint(chimera(o+1),I4P) ! donor block
            if (db==sb(1)%ab) then
               ! self block has a chimera reference with a split block that must be sanitized, if fall in sb(2)
               di = nint(chimera(o+2),I4P) ; dj = nint(chimera(o+3),I4P) ; dk = nint(chimera(o+4),I4P)
               ! check if reference fall in sb(2) domain
               select case(sb(1)%split_dir)
               case(1)
                  if (di>sb(1)%Ni) then
                     chimera(o+1) = real(db+1,R4P)        ! point to sb(2)%ab = sb(1)%ab+1 = db + 1
                     chimera(o+2) = real(di-sb(1)%Ni,R4P) ! point to sb(2)%i, other 2 indexes are the same
                  endif
               case(2)
                  if (dj>sb(1)%Nj) then
                     chimera(o+1) = real(db+1,R4P)        ! point to sb(2)%ab = sb(1)%ab+1 = db + 1
                     chimera(o+3) = real(dj-sb(1)%Nj,R4P) ! point to sb(2)%j, other 2 indexes are the same
                  endif
               case(3)
                  if (dk>sb(1)%Nk) then
                     chimera(o+1) = real(db+1,R4P)        ! point to sb(2)%ab = sb(1)%ab+1 = db + 1
                     chimera(o+4) = real(dk-sb(1)%Nk,R4P) ! point to sb(2)%k, other 2 indexes are the same
                  endif
               endselect
            elseif (db>sb(1)%ab) then
               ! self block has a chimera reference with a block subsequent to split block, ab must be shifted
               chimera(o+1) = real(db+1,R4P)
               ! other chimera indexes are the same
            else
               ! for chimera references to block previous to split block sanitize is not necessary
            endif
         enddo
      endselect
   enddo
   enddo
   enddo
   endassociate
   endsubroutine sanitize_chimera

   subroutine save_block_file(self, basename, tec)
   !< Save block data into its own file.
   class(block_object), intent(in)           :: self      !< Block data.
   character(*),        intent(in), optional :: basename  !< Basename.
   logical,             intent(in), optional :: tec       !< Save (also) in tecplot (ASCII) format (for debugging).
   character(:), allocatable                 :: bn        !< Basename, local var.
   logical                                   :: tec_      !< Save (also) in tecplot (ASCII) format (for debugging), local var.
   character(len=6)                          :: bstr      !< Block number stringified.
   integer(I4P)                              :: file_unit !< Block unit file.
   integer(I4P)                              :: i,j,k,p,n !< Counter.
   integer(I4P)                              :: o         !< Offset of rcc.

   bn = '' ; if (present(basename)) bn = trim(adjustl(basename))
   tec_ = .false. ; if (present(tec)) tec_ = tec
   write(bstr, '(I6.5)') self%ab
   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes,icc=>self%icc,tcc=>self%tcc,chimera=>self%chimera)
   ! file grid
   open(newunit=file_unit, file=bn//'block-'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write', status='replace')
   write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
   write(file_unit)(((nodes(1,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   write(file_unit)(((nodes(2,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   write(file_unit)(((nodes(3,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   close(file_unit)
   ! file rcc
   open(newunit=file_unit, file='block-rcc'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write', status='replace')
   write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
   do k=1-gc, Nk+gc
   do j=1-gc, Nj+gc
   do i=1-gc, Ni+gc
      select case(tcc(1,i,j,k))
      case(BC_NATURAL_EXTRAPOLATED_ALT:BC_NATURAL_WALL)
         write(file_unit) i,j,k,bc_string(tcc(1,i,j,k))
      case(BC_ACTIVE_CELL)
         ! save nothing for active cell
      case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN,BC_CHIMERA_EDGE)
         p = tcc(2,i,j,k)
         write(file_unit) i,j,k,bc_string(tcc(1,i,j,k))
         write(file_unit) nint(chimera(p),I4P) ! donors number
         do n=1, nint(chimera(p),I4P) ! b,i,j,k,weight for each donor
            o = p + 5*(n-1)
            write(file_unit) nint(chimera(o+1)),nint(chimera(o+2)),nint(chimera(o+3)),nint(chimera(o+4)),chimera(o+5)
         enddo
      case default
         print *, 'error: unknown tcc "',tcc(1,i,j,k),'", b,i,j,k=',self%ab,i,j,k
         stop
      endselect
   enddo
   enddo
   enddo

   if (tec_) then
      ! save block in tecplot format
      open(newunit=file_unit, file='block-'//trim(adjustl(bstr))//'.dat', action='write', status='replace')
      write(file_unit,*) 'TITLE = "Block '//trim(adjustl(bstr))//'"'
      write(file_unit,*) 'VARIABLES = "X", "Y", "Z", "tcc" "b" "i" "j" "k"'
      write(file_unit,*) 'ZONE I=', ni+1+2*gc, ', J=', nj+1+2*gc, ', K=', nk+1+2*gc, &
                         ', DATAPACKING=BLOCK, VARLOCATION=([4,5,6,7,8]=CELLCENTERED)'
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(1,i,j,k) ; enddo ; enddo ; enddo
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(2,i,j,k) ; enddo ; enddo ; enddo
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(3,i,j,k) ; enddo ; enddo ; enddo

      ! tcc
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         write(file_unit,'(E23.15)') real(tcc(1,i,j,k),R8P)
      enddo ; enddo ; enddo
      ! b
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         select case(tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            p = tcc(2,i,j,k)
            write(file_unit,'(E23.15)') real(chimera(p+1),R8P)
         case default
            write(file_unit,'(E23.15)') -33._R8P
         endselect
      enddo ; enddo ; enddo
      ! i
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         select case(tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            p = tcc(2,i,j,k)
            write(file_unit,'(E23.15)') real(chimera(p+2),R8P)
         case default
            write(file_unit,'(E23.15)') -33._R8P
         endselect
      enddo ; enddo ; enddo
      ! j
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         select case(tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            p = tcc(2,i,j,k)
            write(file_unit,'(E23.15)') real(chimera(p+3),R8P)
         case default
            write(file_unit,'(E23.15)') -33._R8P
         endselect
      enddo ; enddo ; enddo
      ! k
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         select case(tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            p = tcc(2,i,j,k)
            write(file_unit,'(E23.15)') real(chimera(p+4),R8P)
         case default
            write(file_unit,'(E23.15)') -33._R8P
         endselect
      enddo ; enddo ; enddo
      close(file_unit)
   endif
   endassociate
   endsubroutine save_block_file

   pure subroutine split(self, mgl, is_split_done, sb, split_data, use_cc_par)
   !< Split block. Split current block (if possible), along a the largest direction, in half: the first
   !< Block substitute current block, the other is added to blocks list.
   class(block_object), intent(in)           :: self          !< Block data.
   integer(I4P),        intent(in)           :: mgl           !< Number of levels of multi-grid to be preserved.
   logical,             intent(out)          :: is_split_done !< Sentinel to check is split has been done.
   type(block_object),  intent(out)          :: sb(2)         !< Split blocks.
   logical,             intent(in), optional :: split_data    !< Sentinel to split also data.
   logical,             intent(in), optional :: use_cc_par    !< Sentinel to enable cc.par use.
   logical                                   :: split_data_   !< Sentinel to split also data, local var.
   logical                                   :: use_cc_par_   !< Sentinel to enable cc.par use, local var.
   integer(I4P)                              :: nadj          !< Number of new adjacent-BC cells.
   integer(I4P)                              :: delta(3)      !< Directions deltas.

   split_data_ = .false. ; if (present(split_data)) split_data_ = split_data
   use_cc_par_ = .false. ; if (present(use_cc_par)) use_cc_par_ = use_cc_par
   is_split_done = .false.
   call sb(1)%destroy
   call sb(2)%destroy
   call find_split(Ni=self%Ni,Nj=self%Nj,Nk=self%Nk,gc=self%gc,mgl=mgl,nadj=nadj,delta=delta,sb=sb)
   if (sb(1)%split_dir>0) then
      sb(1)%w = sb(1)%weight()
      sb(2)%w = sb(2)%weight()
      ! assign group, body, proc tags
      sb%group = self%group
      sb%body  = self%body
      sb%proc  = self%proc
      ! assign absolute block index
      sb(1)%ab = self%ab     ! first  child substitute parent
      sb(2)%ab = self%ab + 1 ! second child substitute subsequent block of parent, shift is necessary in global blocks list
      ! update split history
      sb(1)%split_level = self%split_level + 1
      sb(2)%split_level = self%split_level + 1
      if (allocated(self%parents)) then
         sb(1)%parents = [self%parents,self%ab]
         sb(2)%parents = [self%parents,self%ab]
      else
         sb(1)%parents = [             self%ab]
         sb(2)%parents = [             self%ab]
      endif
      is_split_done = .true.
      if (split_data_) then
         ! split data
         call sb(1)%alloc
         call sb(2)%alloc
         call split_nodes(gc=self%gc,nodes=self%nodes,delta=delta,sb=sb)
         if (use_cc_par_) then
            call split_patches(patches=self%patches,delta=delta,sb=sb)
         else
            call split_tcc(tcc=self%tcc,delta=delta,gc=self%gc,sb=sb)
            call split_chimera(chimera=self%chimera,delta=delta,nadj=nadj,sb_n=1,ab_ob=sb(2)%ab,bs=sb(1))
            call split_chimera(chimera=self%chimera,delta=delta,nadj=nadj,sb_n=2,ab_ob=sb(1)%ab,bs=sb(2))
         endif
      endif
   endif
   contains
      pure subroutine find_split(Ni,Nj,Nk,gc,mgl,nadj,delta,sb)
      !< Find split related data, direction, cells numbers, etc.
      integer(I4P),       intent(in)    :: Ni,Nj,Nk,gc !< Block dimensions.
      integer(I4P),       intent(in)    :: mgl         !< Number of levels of multi-grid to be preserved.
      integer(I4P),       intent(inout) :: nadj        !< Number of new adjacent-BC cells.
      integer(I4P),       intent(inout) :: delta(3)    !< Directions deltas.
      type(block_object), intent(inout) :: sb(2)       !< Split blocks.
      integer(I4P)                      :: maxd,mind   !< Temporary variables.
      integer(I4P)                      :: dir(1:3)    !< Ordered directions list.
      integer(I4P)                      :: N           !< Current direction cells number.
      integer(I4P)                      :: Ns(2)       !< Cells numbers in the two split blocks.
      logical                           :: dir_found   !< Sentil to check direction found.
      integer(I4P)                      :: d           !< Counter.

      sb%split_dir = 0
      ! search for direction with highest number of cells being compatible with MG level
      maxd = maxloc([Ni,Nj,Nk],dim=1) ; mind = minloc([Ni,Nj,Nk],dim=1) ; dir = [maxd,6-maxd-mind,mind]
      dir_found = .false.
      direction_loop: do d=1, 3
         select case(dir(d))
         case(1) ! i direction
            N = Ni
         case(2) ! j direction
            N = Nj
         case(3) ! k direction
            N = Nk
         endselect
         if (N/2**(mgl)<2) cycle direction_loop ! not enough MG levels
         dir_found = .true.
         sb%split_dir = dir(d)
         Ns(1) = (((N)/2**mgl)/2)*(2**mgl) ; Ns(2) = N-Ns(1)
         exit direction_loop
      enddo direction_loop

      if (dir_found) then
         ! assign dimensions and deltas
         sb%Ni = Ni
         sb%Nj = Nj
         sb%Nk = Nk
         sb%gc = gc
         delta = 0
         select case(sb(1)%split_dir)
         case(1) ! i direction
            sb(1)%Ni = Ns(1) ; sb(2)%Ni = Ns(2)
            nadj = gc*(Nj+2*gc)*(Nk+2*gc)
            delta(1) = 1
         case(2) ! j direction
            sb(1)%Nj = Ns(1) ; sb(2)%Nj = Ns(2)
            nadj = (Ni+2*gc)*gc*(Nk+2*gc)
            delta(2) = 1
         case(3) ! k direction
            sb(1)%Nk = Ns(1) ; sb(2)%Nk = Ns(2)
            nadj = (Ni+2*gc)*(Nj+2*gc)*gc
            delta(3) = 1
         endselect
      endif
      endsubroutine find_split

      pure subroutine split_chimera(chimera,delta,nadj,sb_n,ab_ob,bs)
      real(R4P),          intent(in)    :: chimera(:)    !< Parent chimera data.
      integer(I4P),       intent(in)    :: delta(3)      !< Deltas.
      integer(I4P),       intent(in)    :: nadj          !< Number of new adjacent-BC cells.
      integer(I4P),       intent(in)    :: sb_n          !< Current split block number, 1 or 2.
      integer(I4P),       intent(in)    :: ab_ob         !< Absolute block index of other split block.
      type(block_object), intent(inout) :: bs            !< Split block.
      integer(I4P)                      :: nchimera      !< Number of chimera data for split block.
      integer(I4P)                      :: n,o,p,ndonors !< Counter.
      integer(I4P)                      :: i,j,k         !< Counter.
      ! local parameters
      integer(I4P), parameter :: BC_ADJ0(3)=[BC_CHIMERA_FACE_ADJ_I0,&
                                             BC_CHIMERA_FACE_ADJ_J0,&
                                             BC_CHIMERA_FACE_ADJ_K0]          !< List of BC for new adj of sb(2).
      integer(I4P), parameter :: BC_ADJN(3)=[BC_CHIMERA_FACE_ADJ_IN,&
                                             BC_CHIMERA_FACE_ADJ_JN,&
                                             BC_CHIMERA_FACE_ADJ_KN]          !< List of BC for new adj of sb(1).
      integer(I4P), parameter :: BC_ADJ(3,2)=reshape([BC_ADJN,BC_ADJ0],[3,2]) !< List of BC for new adj.

      associate(Ni=>bs%Ni,Nj=>bs%Nj,Nk=>bs%Nk,gc=>bs%gc)
      ! count chimera data of parent
      nchimera = 0
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         select case(bs%tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN,BC_CHIMERA_EDGE)
            p = bs%tcc(2,i,j,k)
            ndonors = nint(chimera(p))
            nchimera = nchimera + 1 + ndonors*5 ! b,i,j,k,weight for each donor
         endselect
      enddo
      enddo
      enddo
      allocate(bs%chimera(nchimera+nadj*6)) ! chimera data from self + new adjacent chimera data (1 donor)
      ! assign old chimera data
      nchimera = 0
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         select case(bs%tcc(1,i,j,k))
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN,BC_CHIMERA_EDGE)
            nchimera = nchimera + 1
            p = bs%tcc(2,i,j,k)               ! point to parent chimera array
            bs%tcc(2,i,j,k) = nchimera        ! point to new split block chimera array
            bs%chimera(nchimera) = chimera(p) ! assigno donors number
            do n=1, nint(chimera(p),I4P) ! b,i,j,k,weight for each donor
               o = p + 5*(n-1)
               ! Check if donor references the parent block (intra-block reference)
               if (nint(chimera(o+1),R4P)==bs%ab+(1-sb_n)) then
                  ! Intra-block reference: donor is in parent block being split
                  ! Determine which sub-block the donor falls into
                  select case(bs%split_dir)
                  case(1) ! Split in i direction
                     if (nint(chimera(o+2),R4P)>sb(1)%Ni) then
                        ! Donor falls in sb(2)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n)+1,R4P)       ! Point to sb(2)%ab
                        bs%chimera(nchimera+2) = chimera(o+2) - real(sb(1)%Ni,R4P) ! Adjust i index
                        bs%chimera(nchimera+3) = chimera(o+3)                      ! j unchanged
                        bs%chimera(nchimera+4) = chimera(o+4)                      ! k unchanged
                     else
                        ! Donor falls in sb(1)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n),R4P)         ! Point to sb(1)%ab
                        bs%chimera(nchimera+2) = chimera(o+2)                      ! i unchanged
                        bs%chimera(nchimera+3) = chimera(o+3)                      ! j unchanged
                        bs%chimera(nchimera+4) = chimera(o+4)                      ! k unchanged
                     endif
                  case(2) ! Split in j direction
                     if (nint(chimera(o+3),R4P)>sb(1)%Nj) then
                        ! Donor falls in sb(2)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n)+1,R4P)       ! Point to sb(2)%ab
                        bs%chimera(nchimera+2) = chimera(o+2)                      ! i unchanged
                        bs%chimera(nchimera+3) = chimera(o+3) - real(sb(1)%Nj,R4P) ! Adjust j index
                        bs%chimera(nchimera+4) = chimera(o+4)                      ! k unchanged
                     else
                        ! Donor falls in sb(1)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n),R4P)         ! Point to sb(1)%ab
                        bs%chimera(nchimera+2) = chimera(o+2)                      ! i unchanged
                        bs%chimera(nchimera+3) = chimera(o+3)                      ! j unchanged
                        bs%chimera(nchimera+4) = chimera(o+4)                      ! k unchanged
                     endif
                  case(3) ! Split in k direction
                     if (nint(chimera(o+4),R4P)>sb(1)%Nk) then
                        ! Donor falls in sb(2)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n)+1,R4P)       ! Point to sb(2)%ab
                        bs%chimera(nchimera+2) = chimera(o+2)                      ! i unchanged
                        bs%chimera(nchimera+3) = chimera(o+3)                      ! j unchanged
                        bs%chimera(nchimera+4) = chimera(o+4) - real(sb(1)%Nk,R4P) ! Adjust k index
                     else
                        ! Donor falls in sb(1)
                        bs%chimera(nchimera+1) = real(bs%ab+(1-sb_n),R4P)         ! Point to sb(1)%ab
                        bs%chimera(nchimera+2) = chimera(o+2)                      ! i unchanged
                        bs%chimera(nchimera+3) = chimera(o+3)                      ! j unchanged
                        bs%chimera(nchimera+4) = chimera(o+4)                      ! k unchanged
                     endif
                  endselect
               elseif (nint(chimera(o+1),R4P)>bs%ab+(1-sb_n)) then
                  ! Reference to a block subsequent to the split one, ab index must be shifted
                  bs%chimera(nchimera+1) = chimera(o+1) + 1
                  bs%chimera(nchimera+2) = chimera(o+2)
                  bs%chimera(nchimera+3) = chimera(o+3)
                  bs%chimera(nchimera+4) = chimera(o+4)
               else
                  ! Reference to a block before the split one, no change needed
                  bs%chimera(nchimera+1) = chimera(o+1)
                  bs%chimera(nchimera+2) = chimera(o+2)
                  bs%chimera(nchimera+3) = chimera(o+3)
                  bs%chimera(nchimera+4) = chimera(o+4)
               endif
               bs%chimera(nchimera+5) = chimera(o+5)  ! Weight always unchanged
               nchimera = nchimera + 5  ! increment inside loop for each donor
            enddo
         endselect
      enddo
      enddo
      enddo
      ! assign new chimera-adjacent data
      do k=1-gc, Nk+gc - (Nk+gc)*delta(3)
      do j=1-gc, Nj+gc - (Nj+gc)*delta(2)
      do i=1-gc, Ni+gc - (Ni+gc)*delta(1)
         nchimera = nchimera + 1
         bs%tcc(1,i+(Ni+gc)*delta(1)*(2-sb_n),j+(Nj+gc)*delta(2)*(2-sb_n),k+(Nk+gc)*delta(3)*(2-sb_n)) = BC_ADJ(bs%split_dir,sb_n)
         bs%tcc(2,i+(Ni+gc)*delta(1)*(2-sb_n),j+(Nj+gc)*delta(2)*(2-sb_n),k+(Nk+gc)*delta(3)*(2-sb_n)) = nchimera
         bs%chimera(nchimera  ) = 1._R4P
         bs%chimera(nchimera+1) = real(ab_ob        ,R4P)
         bs%chimera(nchimera+2) = real(i+gc*delta(1),R4P)
         bs%chimera(nchimera+3) = real(j+gc*delta(2),R4P)
         bs%chimera(nchimera+4) = real(k+gc*delta(3),R4P)
         bs%chimera(nchimera+5) = 1._R4P
         nchimera = nchimera + 5
      enddo
      enddo
      enddo
      ! ! count chimera data of parent
      ! nchimera = 0
      ! do k=1-gc, Nk+gc
      ! do j=1-gc, Nj+gc
      ! do i=1-gc, Ni+gc
      !    select case(bs%tcc(1,i,j,k))
      !    case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN,BC_CHIMERA_EDGE)
      !       p = bs%tcc(2,i,j,k)
      !       ndonors = nint(chimera(p))
      !       nchimera = nchimera + 1 + ndonors*5 ! b,i,j,k,weight for each donor
      !    endselect
      ! enddo
      ! enddo
      ! enddo
      ! allocate(bs%chimera(nchimera+nadj*6)) ! chimera data from self + new adjacent chimera data (1 donor)
      ! ! assign old chimera data
      ! nchimera = 0
      ! do k=1-gc, Nk+gc
      ! do j=1-gc, Nj+gc
      ! do i=1-gc, Ni+gc
      !    select case(bs%tcc(1,i,j,k))
      !    case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN,BC_CHIMERA_EDGE)
      !       nchimera = nchimera + 1
      !       p = bs%tcc(2,i,j,k)               ! point to parent chimera array
      !       bs%tcc(2,i,j,k) = nchimera        ! point to new split block chimera array
      !       bs%chimera(nchimera) = chimera(p) ! assigno donors number
      !       do n=1, nint(chimera(p),I4P) ! b,i,j,k,weight for each donor
      !          o = p + 5*(n-1)
      !          if (nint(chimera(o+1),R4P)>bs%ab+(1-sb_n)) then
      !             ! reference to a block subsequent to the split one, ab index must be shifted
      !             bs%chimera(nchimera+1) = chimera(o+1) + 1
      !          else
      !             bs%chimera(nchimera+1) = chimera(o+1)
      !          endif
      !          bs%chimera(nchimera+2) = chimera(o+2)
      !          bs%chimera(nchimera+3) = chimera(o+3)
      !          bs%chimera(nchimera+4) = chimera(o+4)
      !          bs%chimera(nchimera+5) = chimera(o+5)
      !          nchimera = nchimera + 5  ! increment inside loop for each donor
      !       enddo
      !    endselect
      ! enddo
      ! enddo
      ! enddo
      ! ! assign new chimera-adjacent data
      ! do k=1-gc, Nk+gc - (Nk+gc)*delta(3)
      ! do j=1-gc, Nj+gc - (Nj+gc)*delta(2)
      ! do i=1-gc, Ni+gc - (Ni+gc)*delta(1)
      !    nchimera = nchimera + 1
      !    bs%tcc(1,i+(Ni+gc)*delta(1)*(2-sb_n),j+(Nj+gc)*delta(2)*(2-sb_n),k+(Nk+gc)*delta(3)*(2-sb_n)) = BC_ADJ(bs%split_dir,sb_n)
      !    bs%tcc(2,i+(Ni+gc)*delta(1)*(2-sb_n),j+(Nj+gc)*delta(2)*(2-sb_n),k+(Nk+gc)*delta(3)*(2-sb_n)) = nchimera
      !    bs%chimera(nchimera  ) = 1._R4P
      !    bs%chimera(nchimera+1) = real(ab_ob        ,R4P)
      !    bs%chimera(nchimera+2) = real(i+gc*delta(1),R4P)
      !    bs%chimera(nchimera+3) = real(j+gc*delta(2),R4P)
      !    bs%chimera(nchimera+4) = real(k+gc*delta(3),R4P)
      !    bs%chimera(nchimera+5) = 1._R4P
      !    nchimera = nchimera + 5
      ! enddo
      ! enddo
      ! enddo
      endassociate
      endsubroutine split_chimera

      pure subroutine split_patches(patches,delta,sb)
      !< Split patches.
      type(patch_object), intent(in)    :: patches(1:) !< Patches on block faces (6 patches or more if faces have split patches).
      integer(I4P),       intent(in)    :: delta(3)    !< Deltas.
      type(block_object), intent(inout) :: sb(2)       !< Split blocks.
      integer(I4P)                      :: p1, p2, p   !< Counter.

  ! logical      :: is_connection=.false. !< Flag to inquire is the patch defines a connection.
  ! integer(I4P) :: patch_index=0         !< Patch index, local numeration.
  ! integer(I4P) :: block_index=0         !< Block to which patch belongs, local numeration.
  ! integer(I4P) :: face_index=0          !< Face index in the overset convention: Imin=>1, Imax=2,Jmin=3, Jmax=4, Kmin=5, Kmax=6.
  ! integer(I4P) :: boundary_condition=0  !< Boundary condition (or IJK orientation) in the overset convention.
  ! integer(I4P) :: connect_family=0      !< Index of connected patch or family index of the patch.
  ! integer(I4P) :: ijk_extents(1:6)=&
  !                 [0,0,0,0,0,0]         !< IJK extents of the patch in the overset convention: Imin, Imax, Jmin, Jmax, Kmin, Kmax.
      sb(1)%patches = patches
      sb(2)%patches = patches
      p1 = 0
      p2 = 0
      if     (delta(1)==1_I4P) then
         do p=1, size(patches,dim=1)
            if     (patches(p)%face_index==1_I4P) then
               ! patch is YZ min
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(1)%patches(p1+1) = patches(p)
               sb(1)%patches(p1+1)%face_index = 2_I4P
               sb(1)%patches(p1+1)%boundary_condition = -1 ! to be implemented
               sb(1)%patches(p1+1)%connect_family = sb(2)%ab
               sb(1)%patches(p1+1)%ijk_extents(1) = sb(1)%ni
               sb(1)%patches(p1+1)%ijk_extents(2) = sb(1)%ni
            elseif (patches(p)%face_index==2_I4P) then
               ! patch is YZ max
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(2)%patches(p1-1) = patches(p)
               sb(2)%patches(p1-1)%face_index = 1_I4P
               sb(2)%patches(p1-1)%boundary_condition = -1 ! to be implemented
               sb(2)%patches(p1-1)%connect_family = sb(1)%ab
               sb(2)%patches(p1-1)%ijk_extents(1) = 0_I4P
               sb(2)%patches(p1-1)%ijk_extents(2) = 0_I4P
            else
               ! patch is XY or XZ
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               sb(1)%patches(p1)%ijk_extents(2) = sb(1)%ni
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               sb(2)%patches(p2)%ijk_extents(2) = sb(2)%ni
            endif
         enddo
      elseif (delta(2)==1_I4P) then
         do p=1, size(patches,dim=1)
            if     (patches(p)%face_index==3_I4P) then
               ! patch is XZ min
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(1)%patches(p1+1) = patches(p)
               sb(1)%patches(p1+1)%face_index = 4_I4P
               sb(1)%patches(p1+1)%boundary_condition = -1 ! to be implemented
               sb(1)%patches(p1+1)%connect_family = sb(2)%ab
               sb(1)%patches(p1+1)%ijk_extents(1) = sb(1)%nj
               sb(1)%patches(p1+1)%ijk_extents(2) = sb(1)%nj
            elseif (patches(p)%face_index==4_I4P) then
               ! patch is XZ max
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(2)%patches(p1-1) = patches(p)
               sb(2)%patches(p1-1)%face_index = 3_I4P
               sb(2)%patches(p1-1)%boundary_condition = -1 ! to be implemented
               sb(2)%patches(p1-1)%connect_family = sb(1)%ab
               sb(2)%patches(p1-1)%ijk_extents(1) = 0_I4P
               sb(2)%patches(p1-1)%ijk_extents(2) = 0_I4P
            else
               ! patch is XY or YZ
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               sb(1)%patches(p1)%ijk_extents(4) = sb(1)%nj
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               sb(2)%patches(p2)%ijk_extents(4) = sb(2)%nj
            endif
         enddo
      elseif (delta(3)==1_I4P) then
         do p=1, size(patches,dim=1)
            if     (patches(p)%face_index==5_I4P) then
               ! patch is XY min
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(1)%patches(p1+1) = patches(p)
               sb(1)%patches(p1+1)%face_index = 6_I4P
               sb(1)%patches(p1+1)%boundary_condition = -1 ! to be implemented
               sb(1)%patches(p1+1)%connect_family = sb(2)%ab
               sb(1)%patches(p1+1)%ijk_extents(1) = sb(1)%nk
               sb(1)%patches(p1+1)%ijk_extents(2) = sb(1)%nk
            elseif (patches(p)%face_index==6_I4P) then
               ! patch is XY max
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               ! assuming only 6 patch ever, the following is not general enough
               sb(2)%patches(p1-1) = patches(p)
               sb(2)%patches(p1-1)%face_index = 5_I4P
               sb(2)%patches(p1-1)%boundary_condition = -1 ! to be implemented
               sb(2)%patches(p1-1)%connect_family = sb(1)%ab
               sb(2)%patches(p1-1)%ijk_extents(1) = 0_I4P
               sb(2)%patches(p1-1)%ijk_extents(2) = 0_I4P
            else
               ! patch is XZ or YZ
               p1 = p1 + 1
               sb(1)%patches(p1) = patches(p)
               sb(1)%patches(p1)%ijk_extents(6) = sb(1)%nj
               p2 = p2 + 1
               sb(2)%patches(p2) = patches(p)
               sb(2)%patches(p2)%ijk_extents(6) = sb(2)%nj
            endif
         enddo
      endif
      endsubroutine split_patches

      pure subroutine split_nodes(gc,nodes,delta,sb)
      !< Split nodes.
      integer(I4P),       intent(in)    :: gc                          !< Number of ghost cells.
      real(R8P),          intent(in)    :: nodes(1:,0-gc:,0-gc:,0-gc:) !< Nodes coordinates.
      integer(I4P),       intent(in)    :: delta(3)                    !< Deltas.
      type(block_object), intent(inout) :: sb(2)                       !< Split blocks.
      integer(I4P)                      :: i,j,k                       !< Counter.

      sb(1)%nodes = 0._R8P
      do k=0-gc, sb(1)%Nk+gc
      do j=0-gc, sb(1)%Nj+gc
      do i=0-gc, sb(1)%Ni+gc
         sb(1)%nodes(:,i,j,k) = nodes(:,i,j,k)
      enddo
      enddo
      enddo
      sb(2)%nodes = 0._R8P
      do k=0-gc, sb(2)%Nk+gc
      do j=0-gc, sb(2)%Nj+gc
      do i=0-gc, sb(2)%Ni+gc
         sb(2)%nodes(:,i,j,k) = nodes(:,i+sb(1)%Ni*delta(1),j+sb(1)%Nj*delta(2),k+sb(1)%Nk*delta(3))
      enddo
      enddo
      enddo
      endsubroutine split_nodes

      pure subroutine split_tcc(tcc,delta,gc,sb)
      !< Split tcc.
      integer(I4P),       intent(in)    :: delta(3)                  !< Deltas.
      integer(I4P),       intent(in)    :: gc                        !< Number of ghost cells.
      integer(I4P),       intent(in)    :: tcc(1:,1-gc:,1-gc:,1-gc:) !< Parent tcc data.
      type(block_object), intent(inout) :: sb(2)                     !< Split blocks.
      integer(I4P)                      :: i,j,k                     !< Counter.

      sb(1)%tcc = 0
      do k=1-gc, sb(1)%Nk+gc
      do j=1-gc, sb(1)%Nj+gc
      do i=1-gc, sb(1)%Ni+gc
         sb(1)%tcc(:,i,j,k) = tcc(:,i,j,k)
      enddo
      enddo
      enddo
      sb(2)%tcc = 0
      do k=1-gc, sb(2)%Nk+gc
      do j=1-gc, sb(2)%Nj+gc
      do i=1-gc, sb(2)%Ni+gc
         sb(2)%tcc(:,i,j,k) = tcc(:,i+sb(1)%Ni*delta(1),j+sb(1)%Nj*delta(2),k+sb(1)%Nk*delta(3))
      enddo
      enddo
      enddo
      endsubroutine split_tcc
   endsubroutine split

   elemental function weight(self)
   !< Return block weight (work load).
   class(block_object), intent(in) :: self !< Block data.
   integer(I4P)                    :: weight !< Block weight (work load).

   weight = self%Ni*self%Nj*self%Nk
   endfunction weight

   ! private methods
   pure subroutine assign_block(lhs, rhs)
   !< Operator `=`
   class(block_object), intent(inout) :: lhs !< Left hand side.
   type(block_object),  intent(in)    :: rhs !< Right hand side.

   call lhs%destroy
                               lhs%Ni          = rhs%Ni
                               lhs%Nj          = rhs%Nj
                               lhs%Nk          = rhs%Nk
                               lhs%gc          = rhs%gc
                               lhs%w           = rhs%w
   if (allocated(rhs%nodes  )) lhs%nodes       = rhs%nodes
   if (allocated(rhs%icc    )) lhs%icc         = rhs%icc
   if (allocated(rhs%tcc    )) lhs%tcc         = rhs%tcc
   if (allocated(rhs%chimera)) lhs%chimera     = rhs%chimera
   if (allocated(rhs%patches)) lhs%patches     = rhs%patches
                               lhs%ab          = rhs%ab
                               lhs%group       = rhs%group
                               lhs%body        = rhs%body
                               lhs%level       = rhs%level
                               lhs%priority    = rhs%priority
                               lhs%comment     = rhs%comment
                               lhs%proc        = rhs%proc
                               lhs%is_loaded   = rhs%is_loaded
   if (allocated(rhs%parents)) lhs%parents     = rhs%parents
                               lhs%split_level = rhs%split_level
                               lhs%split_dir   = rhs%split_dir
   endsubroutine assign_block

   ! non TBP
   pure function bc_int_type(bc_string)
   !< Return BC integer-type-parameter given a string tag (UPPER or lower case).
   character(*), intent(in) :: bc_string   !< String tag (UPPER or lower case).
   integer(I4P)             :: bc_int_type !< BC integer-type-parameter.

   select case(trim(adjustl(bc_string)))
   case('BC_NATURAL_WALL'                    ,'bc_natural_wall'                    );bc_int_type=BC_NATURAL_WALL
   case('BC_NATURAL_SIMMETRY'                ,'bc_natural_simmetry'                );bc_int_type=BC_NATURAL_SIMMETRY
   case('BC_NATURAL_INFLOW'                  ,'bc_natural_inflow'                  );bc_int_type=BC_NATURAL_INFLOW
   case('BC_NATURAL_INOUTFLOW'               ,'bc_natural_inoutflow'               );bc_int_type=BC_NATURAL_INOUTFLOW
   case('BC_NATURAL_ASSIGNED_INFLOW'         ,'bc_natural_assigned_inflow'         );bc_int_type=BC_NATURAL_ASSIGNED_INFLOW
   case('BC_NATURAL_ASSIGNED_PRESSURE'       ,'bc_natural_assigned_pressure'       );bc_int_type=BC_NATURAL_ASSIGNED_PRESSURE
   case('BC_NATURAL_ASSIGNED_NORMAL_VELOCITY','bc_natural_assigned_normal_velocity');bc_int_type=BC_NATURAL_ASSIGNED_NORMAL_VELOCITY
   case('BC_NATURAL_ASSIGNED_RIEMANN'        ,'bc_natural_assigned_riemann'        );bc_int_type=BC_NATURAL_ASSIGNED_RIEMANN
   case('BC_NATURAL_EXTRAPOLATED'            ,'bc_natural_extrapolated'            );bc_int_type=BC_NATURAL_EXTRAPOLATED
   case('BC_NATURAL_MOVING_WALL'             ,'bc_natural_moving_wall'             );bc_int_type=BC_NATURAL_MOVING_WALL
   case('BC_NATURAL_INACTIVE_WALL'           ,'bc_natural_inactive_wall'           );bc_int_type=BC_NATURAL_INACTIVE_WALL
   case('BC_NATURAL_EXTRAPOLATED_ALT'        ,'bc_natural_extrapolated_alt'        );bc_int_type=BC_NATURAL_EXTRAPOLATED_ALT
   case('BC_ACTIVE_CELL'                     ,'bc_active_cell'                     );bc_int_type=BC_ACTIVE_CELL
   case('BC_CHIMERA_FACE_XF'                 ,'bc_chimera_face_xf'                 );bc_int_type=BC_CHIMERA_FACE_XF
   case('BC_CHIMERA_FACE_XF_I0'              ,'bc_chimera_face_xf_i0'              );bc_int_type=BC_CHIMERA_FACE_XF_I0
   case('BC_CHIMERA_FACE_XF_IN'              ,'bc_chimera_face_xf_in'              );bc_int_type=BC_CHIMERA_FACE_XF_IN
   case('BC_CHIMERA_FACE_XF_J0'              ,'bc_chimera_face_xf_j0'              );bc_int_type=BC_CHIMERA_FACE_XF_J0
   case('BC_CHIMERA_FACE_XF_JN'              ,'bc_chimera_face_xf_jn'              );bc_int_type=BC_CHIMERA_FACE_XF_JN
   case('BC_CHIMERA_FACE_XF_K0'              ,'bc_chimera_face_xf_k0'              );bc_int_type=BC_CHIMERA_FACE_XF_K0
   case('BC_CHIMERA_FACE_XF_KN'              ,'bc_chimera_face_xf_kn'              );bc_int_type=BC_CHIMERA_FACE_XF_KN
   case('BC_CHIMERA_CELL'                    ,'bc_chimera_cell'                    );bc_int_type=BC_CHIMERA_CELL
   case('BC_CHIMERA_CELL_INT_WALL'           ,'bc_chimera_cell_int_wall'           );bc_int_type=BC_CHIMERA_CELL_INT_WALL
   case('BC_CHIMERA_FACE_XC'                 ,'bc_chimera_face_xc'                 );bc_int_type=BC_CHIMERA_FACE_XC
   case('BC_CHIMERA_FACE_XC_I0'              ,'bc_chimera_face_xc_i0'              );bc_int_type=BC_CHIMERA_FACE_XC_I0
   case('BC_CHIMERA_FACE_XC_IN'              ,'bc_chimera_face_xc_in'              );bc_int_type=BC_CHIMERA_FACE_XC_IN
   case('BC_CHIMERA_FACE_XC_J0'              ,'bc_chimera_face_xc_j0'              );bc_int_type=BC_CHIMERA_FACE_XC_J0
   case('BC_CHIMERA_FACE_XC_JN'              ,'bc_chimera_face_xc_jn'              );bc_int_type=BC_CHIMERA_FACE_XC_JN
   case('BC_CHIMERA_FACE_XC_K0'              ,'bc_chimera_face_xc_k0'              );bc_int_type=BC_CHIMERA_FACE_XC_K0
   case('BC_CHIMERA_FACE_XC_KN'              ,'bc_chimera_face_xc_kn'              );bc_int_type=BC_CHIMERA_FACE_XC_KN
   case('BC_CHIMERA_FACE_ADJ'                ,'bc_chimera_face_adj'                );bc_int_type=BC_CHIMERA_FACE_ADJ
   case('BC_CHIMERA_FACE_ADJ_I0'             ,'bc_chimera_face_adj_i0'             );bc_int_type=BC_CHIMERA_FACE_ADJ_I0
   case('BC_CHIMERA_FACE_ADJ_IN'             ,'bc_chimera_face_adj_in'             );bc_int_type=BC_CHIMERA_FACE_ADJ_IN
   case('BC_CHIMERA_FACE_ADJ_J0'             ,'bc_chimera_face_adj_j0'             );bc_int_type=BC_CHIMERA_FACE_ADJ_J0
   case('BC_CHIMERA_FACE_ADJ_JN'             ,'bc_chimera_face_adj_jn'             );bc_int_type=BC_CHIMERA_FACE_ADJ_JN
   case('BC_CHIMERA_FACE_ADJ_K0'             ,'bc_chimera_face_adj_k0'             );bc_int_type=BC_CHIMERA_FACE_ADJ_K0
   case('BC_CHIMERA_FACE_ADJ_KN'             ,'bc_chimera_face_adj_kn'             );bc_int_type=BC_CHIMERA_FACE_ADJ_KN
   case('BC_CHIMERA_EDGE'                    ,'bc_chimera_edge'                    );bc_int_type=BC_CHIMERA_EDGE
   endselect
   endfunction bc_int_type

   pure function bc_string(bc_int_type)
   !< Return string tag (UPPER case) given a BC integer-type-parameter.
   integer(I4P), intent(in)  :: bc_int_type !< BC integer-type-parameter.
   character(:), allocatable :: bc_string   !< String tag (UPPER case).

   select case(bc_int_type)
   case(BC_NATURAL_WALL                    ) ; bc_string = 'BC_NATURAL_WALL'
   case(BC_NATURAL_SIMMETRY                ) ; bc_string = 'BC_NATURAL_SIMMETRY'
   case(BC_NATURAL_INFLOW                  ) ; bc_string = 'BC_NATURAL_INFLOW'
   case(BC_NATURAL_INOUTFLOW               ) ; bc_string = 'BC_NATURAL_INOUTFLOW'
   case(BC_NATURAL_ASSIGNED_INFLOW         ) ; bc_string = 'BC_NATURAL_ASSIGNED_INFLOW'
   case(BC_NATURAL_ASSIGNED_PRESSURE       ) ; bc_string = 'BC_NATURAL_ASSIGNED_PRESSURE'
   case(BC_NATURAL_ASSIGNED_NORMAL_VELOCITY) ; bc_string = 'BC_NATURAL_ASSIGNED_NORMAL_VELOCITY'
   case(BC_NATURAL_ASSIGNED_RIEMANN        ) ; bc_string = 'BC_NATURAL_ASSIGNED_RIEMANN'
   case(BC_NATURAL_EXTRAPOLATED            ) ; bc_string = 'BC_NATURAL_EXTRAPOLATED'
   case(BC_NATURAL_MOVING_WALL             ) ; bc_string = 'BC_NATURAL_MOVING_WALL'
   case(BC_NATURAL_INACTIVE_WALL           ) ; bc_string = 'BC_NATURAL_INACTIVE_WALL'
   case(BC_NATURAL_EXTRAPOLATED_ALT        ) ; bc_string = 'BC_NATURAL_EXTRAPOLATED_ALT'
   case(BC_ACTIVE_CELL                     ) ; bc_string = 'BC_ACTIVE_CELL'
   case(BC_CHIMERA_FACE_XF                 ) ; bc_string = 'BC_CHIMERA_FACE_XF'
   case(BC_CHIMERA_FACE_XF_I0              ) ; bc_string = 'BC_CHIMERA_FACE_XF_I0'
   case(BC_CHIMERA_FACE_XF_IN              ) ; bc_string = 'BC_CHIMERA_FACE_XF_IN'
   case(BC_CHIMERA_FACE_XF_J0              ) ; bc_string = 'BC_CHIMERA_FACE_XF_J0'
   case(BC_CHIMERA_FACE_XF_JN              ) ; bc_string = 'BC_CHIMERA_FACE_XF_JN'
   case(BC_CHIMERA_FACE_XF_K0              ) ; bc_string = 'BC_CHIMERA_FACE_XF_K0'
   case(BC_CHIMERA_FACE_XF_KN              ) ; bc_string = 'BC_CHIMERA_FACE_XF_KN'
   case(BC_CHIMERA_CELL                    ) ; bc_string = 'BC_CHIMERA_CELL'
   case(BC_CHIMERA_CELL_INT_WALL           ) ; bc_string = 'BC_CHIMERA_CELL_INT_WALL'
   case(BC_CHIMERA_FACE_XC                 ) ; bc_string = 'BC_CHIMERA_FACE_XC'
   case(BC_CHIMERA_FACE_XC_I0              ) ; bc_string = 'BC_CHIMERA_FACE_XC_I0'
   case(BC_CHIMERA_FACE_XC_IN              ) ; bc_string = 'BC_CHIMERA_FACE_XC_IN'
   case(BC_CHIMERA_FACE_XC_J0              ) ; bc_string = 'BC_CHIMERA_FACE_XC_J0'
   case(BC_CHIMERA_FACE_XC_JN              ) ; bc_string = 'BC_CHIMERA_FACE_XC_JN'
   case(BC_CHIMERA_FACE_XC_K0              ) ; bc_string = 'BC_CHIMERA_FACE_XC_K0'
   case(BC_CHIMERA_FACE_XC_KN              ) ; bc_string = 'BC_CHIMERA_FACE_XC_KN'
   case(BC_CHIMERA_FACE_ADJ                ) ; bc_string = 'BC_CHIMERA_FACE_ADJ'
   case(BC_CHIMERA_FACE_ADJ_I0             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_I0'
   case(BC_CHIMERA_FACE_ADJ_IN             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_IN'
   case(BC_CHIMERA_FACE_ADJ_J0             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_J0'
   case(BC_CHIMERA_FACE_ADJ_JN             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_JN'
   case(BC_CHIMERA_FACE_ADJ_K0             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_K0'
   case(BC_CHIMERA_FACE_ADJ_KN             ) ; bc_string = 'BC_CHIMERA_FACE_ADJ_KN'
   case(BC_CHIMERA_EDGE                    ) ; bc_string = 'BC_CHIMERA_EDGE'
   endselect
   endfunction bc_string

   recursive subroutine blocks_quick_sort(blocks,blocks_list,first,last)
   !< Order blocks list in decreasing-workload-order by quick sort algorithm.
   type(block_object), intent(in)    :: blocks(1:)      !< Blocks data.
   integer(I4P),       intent(inout) :: blocks_list(1:) !< Blocks ordered list.
   integer(I4P),       intent(in)    :: first           !< First sort index.
   integer(I4P),       intent(in)    :: last            !< Last sort index.
   integer(I4P)                      :: i,j             !< Counter.
   integer(I4P)                      :: pivot           !< Pivot.
   integer(I4P)                      :: temp            !< Temporary buffer.

   i = first
   j = last
   pivot = blocks_list((first + last)/2)

   do
      do while (blocks(blocks_list(i))%w > blocks(pivot)%w)
         i = i + 1
      enddo
      do while (blocks(blocks_list(j))%w < blocks(pivot)%w)
         j = j - 1
      enddo
      if (i <= j) then
         temp = blocks_list(i)
         blocks_list(i) = blocks_list(j)
         blocks_list(j) = temp
         i = i + 1
         j = j - 1
      endif
      if (i > j) exit
   enddo
   if (first < j) call blocks_quick_sort(blocks=blocks,blocks_list=blocks_list,first=first,last=j   )
   if (i < last ) call blocks_quick_sort(blocks=blocks,blocks_list=blocks_list,first=i    ,last=last)
   endsubroutine blocks_quick_sort

   subroutine create_blocks_list(blocks, blocks_list)
   !< Create blocks list, workload-decreasing-ordered.
   type(block_object), intent(in)                 :: blocks(1:)     !< Blocks data.
   integer(I4P),       intent(inout), allocatable :: blocks_list(:) !< Blocks list.
   integer(I4P)                                   :: blocks_number  !< Blocks number.
   integer(I4P)                                   :: b              !< Counter.

   blocks_number = size(blocks,dim=1)
   if (blocks_number>0) then
      blocks_list = [(b,b=1,blocks_number)]
      call blocks_quick_sort(blocks=blocks,blocks_list=blocks_list,first=1,last=blocks_number)
   endif
   endsubroutine create_blocks_list

   pure subroutine get_donors(p, rcc, ndonors, donors)
   !< Return donors data given pointer to rcc.
   integer(I4P), intent(in)            :: p             !< Pointer to rcc.
   real(R4P),    intent(in)            :: rcc(1:)       !< rcc unstructured array.
   integer(I4P), intent(out)           :: ndonors       !< Number of donors.
   real(R4P),    intent(out), optional :: donors(1:,1:) !< Donors data.
   integer(I4P)                        :: n             !< Counter.

   if (p<=BC_NATURAL_RCC_RESERVED_DATA) then
      ndonors = 0
      if (present(donors)) donors = 0._R4P
      return
   endif
   ndonors = nint(rcc(p+1))
   if (present(donors)) then
      do n = 1, ndonors
         donors(1,n) = rcc(p + 1 + 5*(n-1) + 1)
         donors(2,n) = rcc(p + 1 + 5*(n-1) + 2)
         donors(3,n) = rcc(p + 1 + 5*(n-1) + 3)
         donors(4,n) = rcc(p + 1 + 5*(n-1) + 4)
         donors(5,n) = rcc(p + 1 + 5*(n-1) + 5)
      enddo
   endif
   endsubroutine get_donors

   subroutine implode_blocks(blocks, rcc)
   !< Implode previously exploded blocks in order to use legacy Xnavis/Xall inputs, transitional debugging mode.
   type(block_object), intent(inout)              :: blocks(1:)    !< Blocks data.
   real(R4P),          intent(inout), allocatable :: rcc(:)        !< rcc unstructured array.
   integer(I4P)                                   :: nchimera      !< Number of chimera data.
   integer(I4P)                                   :: ndonors       !< Number of donors.
   integer(I4P)                                   :: Ni,Nj,Nk,gc   !< Block dimensions.
   integer(I4P)                                   :: b,i,j,k,n,o,p !< Counter.

   if (allocated(rcc)) deallocate(rcc)
   ! count rcc elements and set blocks%icc
   nchimera = BC_NATURAL_RCC_RESERVED_DATA ! start from reserved natural BC data
   do b=1, size(blocks,dim=1)
      Ni = blocks(b)%Ni ; Nj = blocks(b)%Nj ; Nk = blocks(b)%Nk ; gc = blocks(b)%gc
      blocks(b)%icc = 0_I4P
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         if     (blocks(b)%tcc(1,i,j,k)<0) then
            blocks(b)%icc(i,j,k) = -blocks(b)%tcc(1,i,j,k) ! point to BC natural reserved elements
         elseif (blocks(b)%tcc(1,i,j,k)>0) then
            nchimera = nchimera + 1          ! BC type
            blocks(b)%icc(i,j,k) = nchimera  ! update icc pointer
            nchimera = nchimera + 1          ! BC chimera donors
            p = blocks(b)%tcc(2,i,j,k)
            ndonors = nint(blocks(b)%chimera(p))
            nchimera = nchimera + ndonors*5  ! b,i,j,k,weight for each donor
         endif
      enddo
      enddo
      enddo
   enddo
   ! at least the natural BC reserved data is present
   allocate(rcc(1:nchimera))
   do p=1, BC_NATURAL_RCC_RESERVED_DATA
      rcc(p) = -real(p,R4P)
   enddo
   ! set other elements of rcc
   do b=1, size(blocks,dim=1)
      Ni = blocks(b)%Ni ; Nj = blocks(b)%Nj ; Nk = blocks(b)%Nk ; gc = blocks(b)%gc
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         if  (blocks(b)%tcc(1,i,j,k)>0) then
            p = blocks(b)%icc(i,j,k)
            o = blocks(b)%tcc(2,i,j,k)
            ndonors = nint(blocks(b)%chimera(o))
               rcc(p            ) = real(blocks(b)%tcc(1,i,j,k),R4P) ! BC type
               rcc(p+1          ) = real(ndonors,R4P)                ! BC chimera donors
            do n=1, ndonors                                          ! for each donor
               rcc(p+1+5*(n-1)+1) = blocks(b)%chimera(o+5*(n-1)+1)   ! b
               rcc(p+1+5*(n-1)+2) = blocks(b)%chimera(o+5*(n-1)+2)   ! i
               rcc(p+1+5*(n-1)+3) = blocks(b)%chimera(o+5*(n-1)+3)   ! j
               rcc(p+1+5*(n-1)+4) = blocks(b)%chimera(o+5*(n-1)+4)   ! k
               rcc(p+1+5*(n-1)+5) = blocks(b)%chimera(o+5*(n-1)+5)   ! weight
            enddo
         endif
      enddo
      enddo
      enddo
   enddo
   endsubroutine implode_blocks

   subroutine load_file_cc_par(file_name,cc_par,blocks,boxes)
   !< Load file cc.par.
   character(*),        intent(in)                 :: file_name !< File name
   type(cc_par_object), intent(inout)              :: cc_par    !< File cc.par handler.
   type(block_object),  intent(inout)              :: blocks(:) !< Blocks data.
   type(box_object),    intent(inout), allocatable :: boxes(:)  !< Boxes.
   integer(I4P)                                    :: file_unit !< File unit.

   open(newunit=file_unit, file=trim(adjustl(file_name)))
   call load_header
   print '(A)', '  cc.par: header loaded'
   call load_blocks
   print '(A,I0)', '  cc.par: blocks loaded, n=', cc_par%blocks_number
   call load_patches
   print '(A,I0)', '  cc.par: patches loaded, n=', cc_par%patches_number
   call load_edges
   print '(A)', '  cc.par: edges loaded'
   call load_boxes
   print '(A,I0)', '  cc.par: boxes loaded, n=', cc_par%boxes_number
   call load_circuits
   print '(A)', '  cc.par: circuits loaded'
   close(file_unit)
   contains
      subroutine load_header()
      !< Load header.

      read(file_unit,*) cc_par%file_name_input_grd
      read(file_unit,*) cc_par%base_name_output
      read(file_unit,*) cc_par%save_ghost_cells
      read(file_unit,*) cc_par%increase_overlap
      read(file_unit,*) cc_par%extend_internal_wall
      read(file_unit,*)
      read(file_unit,*) cc_par%mgl(1), cc_par%mgl(2)
      read(file_unit,*)
      read(file_unit,*) cc_par%boundary_layer_thickness
      read(file_unit,*) cc_par%numberical_beach
      read(file_unit,*)
      endsubroutine load_header

      subroutine load_blocks()
      !< Load blocks section.
      integer(I4P) :: b !< Counter.

      read(file_unit, *) cc_par%blocks_number
      if (cc_par%blocks_number /= size(blocks,dim=1)) then
         write(stderr, '(A)')'error: grd and cc.par have different number of blocks'
         write(stderr, '(A,I9)')'       grd blocks number:   ', size(blocks,dim=1)
         write(stderr, '(A,I9)')'       cc.par blocks number:', cc_par%blocks_number
         stop
      endif
      read(file_unit, *)
      do b=1, cc_par%blocks_number
         read(file_unit, *) blocks(b)%level, blocks(b)%group, blocks(b)%priority, blocks(b)%comment
      enddo
      read(file_unit, *)
      endsubroutine load_blocks

      subroutine load_patches()
      !< Load patches section.
      integer(I4P)                    :: p                !< Counter.
      integer(I4P)                    :: b                !< Counter.
      integer(I4P)                    :: bp               !< Counter for per-block patch index.
      integer(I4P),       allocatable :: patches_count(:) !< Number of patches per block.
      type(patch_object), allocatable :: patches(:)       !< Patches temporary buffer.

      read(file_unit, *) cc_par%patches_number
      read(file_unit, *)
      if (cc_par%patches_number>0) then
         allocate(patches(1:cc_par%patches_number))
         do p=1, cc_par%patches_number
            read(file_unit, *) patches(p)%block_index,        &
                               patches(p)%face_index,         &
                               patches(p)%boundary_condition, &
                               patches(p)%connect_family,     &
                               patches(p)%ijk_extents(1:6),   &
                               patches(p)%comment
            patches(p)%patch_index  = p
            patches(p)%is_connection = patches(p)%connect_family > 0
         enddo
         read(file_unit, *)
         ! distribute patches to blocks
         allocate(patches_count(1:cc_par%blocks_number), source=0)
         do p=1, cc_par%patches_number
            b = patches(p)%block_index
            patches_count(b) = patches_count(b) + 1
         enddo
         do b=1, cc_par%blocks_number
            allocate(blocks(b)%patches(1:patches_count(b)))
         enddo
         patches_count = 0
         do p=1, cc_par%patches_number
            b = patches(p)%block_index
            patches_count(b) = patches_count(b) + 1
            blocks(b)%patches(patches_count(b)) = patches(p)
         enddo
         deallocate(patches_count)
         deallocate(patches)
      endif
      endsubroutine load_patches

      subroutine load_edges()
      !< Load edges section.
      !< Not yet supported, just a placeholder.

      read(file_unit, *) cc_par%edges_number ! keep this always 0
      read(file_unit, *)
      endsubroutine load_edges

      subroutine load_boxes()
      !< Load boxes section.
      integer(I4P) :: b,n   !< Counter.
      integer(I4P) :: ios    !< I/O status.

      if (allocated(boxes)) deallocate(boxes)
      read(file_unit, *, iostat=ios) cc_par%boxes_number
      if (ios /= 0) then
         cc_par%boxes_number = 0
         return
      endif
      read(file_unit, *, iostat=ios)
      if (cc_par%boxes_number>0) then
         allocate(boxes(1:cc_par%boxes_number))
         do b=1, cc_par%boxes_number
            read(file_unit, *, iostat=ios) boxes(b)%btype, boxes(b)%bblock, boxes(b)%bgroup
            if (ios /= 0) then
               write(stderr, '(A,I0)') 'warning: error reading box ', b
               cc_par%boxes_number = b - 1
               return
            endif
            do n=1, 8 ! box xyz nodes loop
               read(file_unit, *) boxes(b)%nodes(1,n), boxes(b)%nodes(2,n), boxes(b)%nodes(3,n)
            enddo
         enddo
      endif
      endsubroutine load_boxes

      subroutine load_circuits()
      !< Load circuits section.
      !< Not yet supported, just a placeholder.
      integer(I4P) :: ios !< I/O status.

      read(file_unit, *, iostat=ios) cc_par%circuits_number ! keep this always 0
      if (ios /= 0) then
         cc_par%circuits_number = 0
         return
      endif
      read(file_unit, *, iostat=ios)
      endsubroutine load_circuits
   endsubroutine load_file_cc_par

   subroutine load_file_grd(file_name,blocks,blocks_number,gc)
   !< Load file grd.
   character(*),       intent(in)                 :: file_name     !< File name
   type(block_object), intent(inout), allocatable :: blocks(:)     !< Blocks data.
   integer(I4P),       intent(out)                :: blocks_number !< Blocks number.
   integer(I4P),       intent(in), optional       :: gc            !< Number of ghost cells.
   integer(I4P)                                   :: file_unit     !< File unit.
   integer(I4P)                                   :: b             !< Counter.

   if (allocated(blocks)) deallocate(blocks)
   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
   read(file_unit, end=10, err=10) blocks_number
   allocate(blocks(1:blocks_number))
   do b=1, blocks_number
      call blocks(b)%load_dimensions(file_unit=file_unit, ab=b, allocate_data=.true., gc=gc)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_nodes(file_unit=file_unit)
   enddo
   10 close(file_unit)
   endsubroutine load_file_grd

   subroutine load_file_icc(file_name,blocks,blocks_number,rcc)
   !< Load file icc.
   character(*),       intent(in)                 :: file_name          !< File name
   type(block_object), intent(inout)              :: blocks(1:)         !< Blocks data.
   integer(I4P),       intent(in)                 :: blocks_number      !< Blocks number.
   real(R4P),          intent(inout), allocatable :: rcc(:)             !< rcc unstructured array.
   integer(I4P)                                   :: file_unit          !< File unit.
   integer(I4P)                                   :: unstruct_dimension !< Dimension of unstructured array rcc.
   integer(I4P)                                   :: b, i               !< Counter.

   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
   read(file_unit) b
   if (b/=blocks_number) then
      write(stderr, "(A)")'error: grd and icc have different number of blocks'
      stop
   endif
   do b=1, blocks_number
      call blocks(b)%load_dimensions(file_unit=file_unit)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_icc(file_unit=file_unit)
   enddo
   read(file_unit) unstruct_dimension
   allocate(rcc(1:unstruct_dimension))
   read(file_unit) (rcc(i),i=1,unstruct_dimension)
   close(file_unit)
   endsubroutine load_file_icc

   pure subroutine popout_blocks_list(blocks_list)
   !< Pop-out first element of blocks list.
   integer(I4P), intent(inout), allocatable :: blocks_list(:)  !< Blocks list.
   integer(I4P), allocatable                :: blocks_list_(:) !< Blocks list, local var.
   integer(I4P)                             :: b               !< Counter.

   if (size(blocks_list,dim=1)>1) then
      ! pop out the block assigned
      allocate(blocks_list_(1:size(blocks_list,dim=1)-1))
      do b=2, size(blocks_list,dim=1)
         blocks_list_(b-1) = blocks_list(b)
      enddo
      call move_alloc(from=blocks_list_, to=blocks_list)
   else
      ! no other blocks in list, destroy the list
      deallocate(blocks_list)
   endif
   if (allocated(blocks_list_)) deallocate(blocks_list_)
   endsubroutine popout_blocks_list

   subroutine save_file_grd(file_name, blocks)
   !< Save file grd.
   character(*),       intent(in) :: file_name     !< File name
   type(block_object), intent(in) :: blocks(1:)    !< Blocks data.
   integer(I4P)                   :: blocks_number !< Blocks number.
   integer(I4P)                   :: file_unit     !< File unit.
   integer(I4P)                   :: Ni,Nj,Nk,gc   !< Block dimensions.
   integer(I4P)                   :: b,i,j,k       !< Counter.

   blocks_number = size(blocks,dim=1)
   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='write', status='replace')
   write(file_unit) blocks_number
   do b=1, blocks_number
      write(file_unit) blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk,blocks(b)%gc
   enddo
   do b=1, blocks_number
      Ni = blocks(b)%Ni ; Nj = blocks(b)%Nj ; Nk = blocks(b)%Nk ; gc = blocks(b)%gc
      write(file_unit)(((blocks(b)%nodes(1,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
      write(file_unit)(((blocks(b)%nodes(2,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
      write(file_unit)(((blocks(b)%nodes(3,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   enddo
   close(file_unit)
   endsubroutine save_file_grd

   subroutine save_file_icc(file_name, blocks, rcc)
   !< Save file icc.
   character(*),       intent(in) :: file_name          !< File name
   type(block_object), intent(in) :: blocks(1:)         !< Blocks data.
   real(R4P),          intent(in) :: rcc(1:)            !< rcc unstructured array.
   integer(I4P)                   :: file_unit          !< File unit.
   integer(I4P)                   :: blocks_number      !< Blocks number.
   integer(I4P)                   :: unstruct_dimension !< Dimension of unstructured array rcc.
   integer(I4P)                   :: Ni,Nj,Nk,gc        !< Block dimensions.
   integer(I4P)                   :: b,i,j,k            !< Counter.

   blocks_number = size(blocks,dim=1)
   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='write', status='replace')
   write(file_unit) blocks_number
   do b=1, blocks_number
      write(file_unit) blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk,blocks(b)%gc
   enddo
   do b=1, blocks_number
      Ni = blocks(b)%Ni ; Nj = blocks(b)%Nj ; Nk = blocks(b)%Nk ; gc = blocks(b)%gc
      write(file_unit)(((blocks(b)%icc(i,j,k),i=1-gc,Ni+gc),j=1-gc,Nj+gc),k=1-gc,Nk+gc)
   enddo
   unstruct_dimension = size(rcc,dim=1)
   write(file_unit) unstruct_dimension
   write(file_unit) (rcc(i),i=1,unstruct_dimension)
   close(file_unit)
   endsubroutine save_file_icc

   subroutine save_proc_input(blocks, file_name)
   !< Save proc.input file.
   class(block_object), intent(in)           :: blocks(1:)    !< Blocks.
   character(*),        intent(in), optional :: file_name     !< File name.
   integer(I4P)                              :: file_unit     !< Logical file unit.
   integer(I4P)                              :: blocks_number !< Blocks number.
   integer(I4P)                              :: b             !< Counter.

   blocks_number = size(blocks,dim=1)
   if (present(file_name)) then
      open(newunit=file_unit, file=trim(adjustl(file_name)), action='write', status='replace')
   else
      open(newunit=file_unit, file='proc.input',             action='write', status='replace')
   endif
   write(file_unit, *) 'generated by overset-exploded'
   write(file_unit, *) ! record skipped
   write(file_unit, *) ! record skipped
   write(file_unit, *) blocks_number, ' number of blocks'
   write(file_unit, *) ! record skipped
   write(file_unit, *) ! record skipped
   write(file_unit, *) ' block - group - body - processor [history of splits, parents list]'
   do b=1, blocks_number
      if (allocated(blocks(b)%parents)) then
         write(file_unit,*) b, blocks(b)%group, blocks(b)%body, blocks(b)%proc, '!', blocks(b)%parents
      else
         write(file_unit,*) b, blocks(b)%group, blocks(b)%body, blocks(b)%proc, '! original block, no splits'
      endif
   enddo
   close(file_unit)
   endsubroutine save_proc_input

   subroutine update_blocks(blocks, sb, blocks_number, use_cc_par)
   !< Update blocks data after a block split.
   type(block_object), intent(inout), allocatable :: blocks(:)     !< Blocks data.
   type(block_object), intent(in)                 :: sb(1:)        !< Split blocks.
   integer(I4P),       intent(out)                :: blocks_number !< Blocks number.
   logical,            intent(in)                 :: use_cc_par    !< Use cc.par instead of icc.
   type(block_object), allocatable                :: blocks_(:)    !< New blocks data.
   integer(I4P)                                   :: b             !< Counter.

   blocks_number = size(blocks,dim=1)
   ! sanitize old chimera data
   do b=1, blocks_number
      if (b==sb(1)%ab) cycle ! split block does not need to be santized, it is replaced by sb
      if (.not.use_cc_par) call blocks(b)%sanitize_chimera(sb=sb)
   enddo
   ! ab shift
   do b=sb(2)%ab, blocks_number
      blocks(b)%ab = blocks(b)%ab + 1
   enddo
   ! create new blocks data
   allocate(blocks_(1:blocks_number+1))
   do b=1, sb(1)%ab - 1
      blocks_(b) = blocks(b)
      call blocks(b)%destroy  ! free memory immediately
   enddo
   ! destroy the block that was split before replacing it
   call blocks(sb(1)%ab)%destroy
   blocks_(sb(1)%ab) = sb(1)
   blocks_(sb(2)%ab) = sb(2)
   do b=sb(1)%ab + 1, blocks_number
      blocks_(b+1) = blocks(b)
      call blocks(b)%destroy  ! free memory immediately
   enddo
   call move_alloc(from=blocks_, to=blocks)
   blocks_number = blocks_number + 1
   endsubroutine update_blocks

   subroutine replay_splits_on_patches(blocks, splits, splits_dir, splits_nijk, boxes)
   !< Replay split history on patch data.
   !< For each split, distributes patches to children, splits partner patches,
   !< creates new internal connection patches at split interfaces, and shifts indices.
   type(block_object), intent(inout), allocatable :: blocks(:)     !< Blocks data.
   integer(I4P),       intent(in)                 :: splits(:)     !< Splits history (block indices at time of split).
   integer(I4P),       intent(in)                 :: splits_dir(:) !< Splits direction history.
   integer(I4P),       intent(in)                 :: splits_nijk(:)!< Splits first-child size history.
   type(box_object),   intent(inout), allocatable :: boxes(:)      !< Boxes data.
   integer(I4P)                                   :: s             !< Split counter.
   integer(I4P)                                   :: B             !< Block index being split.
   integer(I4P)                                   :: d             !< Split direction.
   integer(I4P)                                   :: Ns1           !< First child size in split direction.
   integer(I4P)                                   :: Nd            !< Parent size in split direction.
   integer(I4P)                                   :: nb            !< Number of blocks.
   integer(I4P)                                   :: np_total      !< Total patches count.
   type(patch_object), allocatable                :: all_patches(:)!< Flat array of all patches.
   type(patch_object), allocatable                :: new_patches(:)!< Temporary for building new patches.
   type(patch_object)                             :: p1, p2        !< Temporary patches for splitting.
   type(block_object), allocatable                :: blocks_(:)    !< Temporary blocks for insertion.
   type(box_object),   allocatable                :: boxes_(:)     !< Temporary boxes.
   integer(I4P)                                   :: b_idx, p_idx  !< Counters.
   integer(I4P)                                   :: ip, jp        !< Patch loop counters.
   integer(I4P)                                   :: new_count     !< New patches count.
   integer(I4P)                                   :: partner_block !< Partner block index.
   integer(I4P)                                   :: partner_patch !< Partner patch index in flat array.
   integer(I4P)                                   :: orient_code   !< 3-digit orientation code.
   integer(I4P)                                   :: o1, o2, o3    !< Orientation digits.
   integer(I4P)                                   :: partner_axis  !< Partner axis corresponding to split direction.
   integer(I4P)                                   :: partner_sign  !< Sign of partner axis mapping (+1 or -1).
   integer(I4P)                                   :: split_pos_partner !< Split position on partner extent.
   integer(I4P)                                   :: ext_lo, ext_hi    !< Extent low/high for partner in split axis.
   integer(I4P)                                   :: face_lo, face_hi  !< Face indices for split direction.
   integer(I4P)                                   :: Nj_face, Nk_face  !< Face dimensions for new internal patches.
   integer(I4P)                                   :: dim_perp1, dim_perp2 !< Perpendicular dimension sizes.
   integer(I4P)                                   :: pc                   !< Patches count per block.
   integer(I4P), allocatable                      :: patches_count(:)     !< Patches count per block.
   integer(I4P), allocatable                      :: pos_to_global(:)     !< Mapping: new_patches position -> global index.

   nb = size(blocks, dim=1)

   do s = 1, size(splits, dim=1)
      B   = splits(s)
      d   = splits_dir(s)
      Ns1 = splits_nijk(s)

      ! Get parent size in split direction
      select case(d)
      case(1) ; Nd = blocks(B)%Ni
      case(2) ; Nd = blocks(B)%Nj
      case(3) ; Nd = blocks(B)%Nk
      endselect

      ! face indices: face_lo = min face (2*d-1), face_hi = max face (2*d)
      face_lo = 2*d - 1
      face_hi = 2*d

      ! Collect all patches into a flat array
      np_total = 0
      do b_idx = 1, nb
         if (allocated(blocks(b_idx)%patches)) np_total = np_total + size(blocks(b_idx)%patches)
      enddo

      allocate(all_patches(np_total))
      ip = 0
      do b_idx = 1, nb
         if (allocated(blocks(b_idx)%patches)) then
            do jp = 1, size(blocks(b_idx)%patches)
               ip = ip + 1
               all_patches(ip) = blocks(b_idx)%patches(jp)
            enddo
         endif
      enddo

      ! Process patches: distribute patches on block B to children
      ! child 1 stays at index B, child 2 goes to index B+1
      ! Allocate new_patches with generous size (original + potential splits + 2 internal)
      allocate(new_patches(np_total * 2 + 2))
      new_count = 0

      ! Two-pass approach: process block-B patches first (which may mark partners as sentinels),
      ! then copy non-B patches. This avoids order-dependent sentinel issues.
      ! Pass 1: process patches on block B
      do ip = 1, np_total
         if (all_patches(ip)%block_index == B) then
            call distribute_patch_to_children(all_patches(ip), d, Ns1, Nd, B, face_lo, face_hi, &
                                              new_patches, new_count, &
                                              all_patches, np_total)
         endif
      enddo
      ! Pass 2: copy non-B patches, skipping sentinels (block_index set to -1 by split_partner_patch)
      do ip = 1, np_total
         if (all_patches(ip)%block_index /= B .and. all_patches(ip)%block_index > 0) then
            new_count = new_count + 1
            new_patches(new_count) = all_patches(ip)
            if (new_patches(new_count)%block_index > B) then
               new_patches(new_count)%block_index = new_patches(new_count)%block_index + 1
            endif
         endif
      enddo

      ! Add two internal connection patches at split interface
      ! Child 1 (block B) gets a face_hi patch connecting to child 2 (B+1)
      new_count = new_count + 1
      new_patches(new_count)%is_connection = .true.
      new_patches(new_count)%block_index = B
      new_patches(new_count)%face_index = face_hi
      new_patches(new_count)%boundary_condition = 135  ! identity orientation
      new_patches(new_count)%connect_family = -(new_count + 1) ! temporary marker, will be fixed
      select case(d)
      case(1)
         new_patches(new_count)%ijk_extents = [Ns1, Ns1, 0, blocks(B)%Nj, 0, blocks(B)%Nk]
      case(2)
         new_patches(new_count)%ijk_extents = [0, blocks(B)%Ni, Ns1, Ns1, 0, blocks(B)%Nk]
      case(3)
         new_patches(new_count)%ijk_extents = [0, blocks(B)%Ni, 0, blocks(B)%Nj, Ns1, Ns1]
      endselect
      new_patches(new_count)%comment = '! split-internal'
      p_idx = new_count  ! remember index for cross-linking

      ! Child 2 (block B+1) gets a face_lo patch connecting to child 1 (B)
      new_count = new_count + 1
      new_patches(new_count)%is_connection = .true.
      new_patches(new_count)%block_index = B + 1
      new_patches(new_count)%face_index = face_lo
      new_patches(new_count)%boundary_condition = 135  ! identity orientation
      new_patches(new_count)%connect_family = p_idx    ! links to the child-1 patch
      select case(d)
      case(1)
         new_patches(new_count)%ijk_extents = [0, 0, 0, blocks(B)%Nj, 0, blocks(B)%Nk]
      case(2)
         new_patches(new_count)%ijk_extents = [0, blocks(B)%Ni, 0, 0, 0, blocks(B)%Nk]
      case(3)
         new_patches(new_count)%ijk_extents = [0, blocks(B)%Ni, 0, blocks(B)%Nj, 0, 0]
      endselect
      new_patches(new_count)%comment = '! split-internal'

      ! Fix cross-link for child-1 internal patch
      new_patches(p_idx)%connect_family = new_count

      ! Renumber connect_family references and assign sequential patch_index.
      ! NOTE: do NOT reassign patch_index before this call - renumber_connect_families
      ! needs the old patch_index values to build the old-to-new mapping.
      ! Pre-linked patches (from split_partner_patch) have patch_index=0 and are skipped.
      call renumber_connect_families(new_patches, new_count)

      ! Insert child 2 block in blocks array
      nb = nb + 1
      allocate(blocks_(nb))
      do b_idx = 1, B
         blocks_(b_idx) = blocks(b_idx)
      enddo
      ! Child 2: copy metadata from parent, adjust dimensions
      blocks_(B+1) = blocks(B)
      blocks_(B+1)%ab = B + 1
      select case(d)
      case(1)
         blocks_(B)%Ni   = Ns1
         blocks_(B+1)%Ni = Nd - Ns1
      case(2)
         blocks_(B)%Nj   = Ns1
         blocks_(B+1)%Nj = Nd - Ns1
      case(3)
         blocks_(B)%Nk   = Ns1
         blocks_(B+1)%Nk = Nd - Ns1
      endselect
      ! Shift ab for all subsequent blocks
      do b_idx = B + 2, nb
         blocks_(b_idx) = blocks(b_idx - 1)
         blocks_(b_idx)%ab = b_idx
      enddo
      deallocate(blocks)
      call move_alloc(from=blocks_, to=blocks)

      ! Shift box indices
      if (allocated(boxes)) then
         do b_idx = 1, size(boxes)
            if (boxes(b_idx)%bblock > B) then
               boxes(b_idx)%bblock = boxes(b_idx)%bblock + 1
            endif
         enddo
      endif

      ! Redistribute patches to blocks
      allocate(patches_count(nb))
      patches_count = 0
      do ip = 1, new_count
         b_idx = new_patches(ip)%block_index
         if (b_idx >= 1 .and. b_idx <= nb) then
            patches_count(b_idx) = patches_count(b_idx) + 1
         endif
      enddo
      do b_idx = 1, nb
         if (allocated(blocks(b_idx)%patches)) deallocate(blocks(b_idx)%patches)
         if (patches_count(b_idx) > 0) allocate(blocks(b_idx)%patches(patches_count(b_idx)))
      enddo
      patches_count = 0
      do ip = 1, new_count
         b_idx = new_patches(ip)%block_index
         if (b_idx >= 1 .and. b_idx <= nb) then
            patches_count(b_idx) = patches_count(b_idx) + 1
            blocks(b_idx)%patches(patches_count(b_idx)) = new_patches(ip)
         endif
      enddo
      deallocate(patches_count)

      ! Remap connect_family from new_patches positions to global block-by-block indices.
      ! The patches were redistributed to blocks, but their connect_family still references
      ! positions in new_patches. The global ordering (block by block) differs from new_patches order.
      allocate(pos_to_global(new_count), source=0)
      pc = 0
      do b_idx = 1, nb
         if (.not.allocated(blocks(b_idx)%patches)) cycle
         do jp = 1, size(blocks(b_idx)%patches)
            pc = pc + 1
            if (blocks(b_idx)%patches(jp)%patch_index >= 1 .and. &
                blocks(b_idx)%patches(jp)%patch_index <= new_count) then
               pos_to_global(blocks(b_idx)%patches(jp)%patch_index) = pc
            endif
         enddo
      enddo
      ! Update connect_family using the mapping, and assign final sequential patch_index
      pc = 0
      do b_idx = 1, nb
         if (.not.allocated(blocks(b_idx)%patches)) cycle
         do jp = 1, size(blocks(b_idx)%patches)
            pc = pc + 1
            if (blocks(b_idx)%patches(jp)%is_connection .and. &
                blocks(b_idx)%patches(jp)%connect_family >= 1 .and. &
                blocks(b_idx)%patches(jp)%connect_family <= new_count) then
               blocks(b_idx)%patches(jp)%connect_family = &
                  pos_to_global(blocks(b_idx)%patches(jp)%connect_family)
            endif
            blocks(b_idx)%patches(jp)%patch_index = pc
         enddo
      enddo
      deallocate(pos_to_global)

      deallocate(all_patches)
      deallocate(new_patches)
   enddo
   contains
      subroutine distribute_patch_to_children(patch, d, Ns1, Nd, B, face_lo, face_hi, &
                                              new_patches, new_count, all_patches, np_total)
      !< Distribute a patch on block B to child 1 (B) and child 2 (B+1).
      type(patch_object), intent(in)    :: patch           !< Patch to distribute.
      integer(I4P),       intent(in)    :: d               !< Split direction (1=i, 2=j, 3=k).
      integer(I4P),       intent(in)    :: Ns1             !< First child size in split direction.
      integer(I4P),       intent(in)    :: Nd              !< Parent size in split direction.
      integer(I4P),       intent(in)    :: B               !< Block index being split.
      integer(I4P),       intent(in)    :: face_lo         !< Min face index (2*d-1).
      integer(I4P),       intent(in)    :: face_hi         !< Max face index (2*d).
      type(patch_object), intent(inout) :: new_patches(:)  !< New patches array.
      integer(I4P),       intent(inout) :: new_count       !< Current count of new patches.
      type(patch_object), intent(inout) :: all_patches(:)  !< All patches (for partner splitting).
      integer(I4P),       intent(in)    :: np_total        !< Total patches count.
      integer(I4P)                      :: ext_lo, ext_hi  !< Extent in split direction.
      integer(I4P)                      :: d_lo, d_hi      !< Direction extent indices (1-based pair).

      ! Determine extent indices in ijk_extents for the split direction
      ! ijk_extents = [imin, imax, jmin, jmax, kmin, kmax]
      d_lo = 2*d - 1  ! index into ijk_extents for min of split direction
      d_hi = 2*d      ! index into ijk_extents for max of split direction

      ext_lo = patch%ijk_extents(d_lo)
      ext_hi = patch%ijk_extents(d_hi)

      if (patch%face_index == face_lo) then
         ! Min face perpendicular to d -> goes to child 1, extents unchanged
         new_count = new_count + 1
         new_patches(new_count) = patch
         new_patches(new_count)%block_index = B
      elseif (patch%face_index == face_hi) then
         ! Max face perpendicular to d -> goes to child 2, shift extent
         new_count = new_count + 1
         new_patches(new_count) = patch
         new_patches(new_count)%block_index = B + 1
         new_patches(new_count)%ijk_extents(d_lo) = ext_lo - Ns1
         new_patches(new_count)%ijk_extents(d_hi) = ext_hi - Ns1
         ! If connection, update partner to point to new patch index
         if (patch%is_connection .and. patch%connect_family > 0) then
            call update_partner_block_ref(all_patches, np_total, patch%connect_family, B, B+1)
         endif
      else
         ! Parallel face: check extent range in split direction
         if (ext_hi <= Ns1) then
            ! Fully within child 1
            new_count = new_count + 1
            new_patches(new_count) = patch
            new_patches(new_count)%block_index = B
         elseif (ext_lo >= Ns1) then
            ! Fully within child 2
            new_count = new_count + 1
            new_patches(new_count) = patch
            new_patches(new_count)%block_index = B + 1
            new_patches(new_count)%ijk_extents(d_lo) = ext_lo - Ns1
            new_patches(new_count)%ijk_extents(d_hi) = ext_hi - Ns1
            if (patch%is_connection .and. patch%connect_family > 0) then
               call update_partner_block_ref(all_patches, np_total, patch%connect_family, B, B+1)
            endif
         else
            ! Spans Ns1 -> split into two patches
            ! Part 1: child 1 (ext_lo to Ns1)
            new_count = new_count + 1
            new_patches(new_count) = patch
            new_patches(new_count)%block_index = B
            new_patches(new_count)%ijk_extents(d_hi) = Ns1

            ! Part 2: child 2 (Ns1 to ext_hi, shifted)
            new_count = new_count + 1
            new_patches(new_count) = patch
            new_patches(new_count)%block_index = B + 1
            new_patches(new_count)%ijk_extents(d_lo) = 0
            new_patches(new_count)%ijk_extents(d_hi) = ext_hi - Ns1

            ! If this is a connection, also split the partner patch
            if (patch%is_connection .and. patch%connect_family > 0) then
               call split_partner_patch(all_patches, np_total, patch%connect_family, &
                                        patch%boundary_condition, d, Ns1, ext_lo, ext_hi, Nd, &
                                        new_patches, new_count, B)
            endif
         endif
      endif
      endsubroutine distribute_patch_to_children

      subroutine update_partner_block_ref(all_patches, np_total, partner_idx, old_block, new_block)
      !< Update partner patch to reference new block index when a connection patch moves to child 2.
      type(patch_object), intent(inout) :: all_patches(:) !< All patches.
      integer(I4P),       intent(in)    :: np_total       !< Total patches count.
      integer(I4P),       intent(in)    :: partner_idx    !< Partner patch index.
      integer(I4P),       intent(in)    :: old_block      !< Old block index.
      integer(I4P),       intent(in)    :: new_block      !< New block index.
      ! No-op: partner's connect_family references patch index, not block
      ! The renumbering step will handle this
      endsubroutine update_partner_block_ref

      subroutine split_partner_patch(all_patches, np_total, partner_idx, orient, &
                                     d, Ns1, ext_lo, ext_hi, Nd, &
                                     new_patches, new_count, B)
      !< Split partner patch when a connection patch spanning the split position is split.
      type(patch_object), intent(inout) :: all_patches(:) !< All patches.
      integer(I4P),       intent(in)    :: np_total       !< Total patches count.
      integer(I4P),       intent(in)    :: partner_idx    !< Partner patch index in all_patches.
      integer(I4P),       intent(in)    :: orient         !< 3-digit orientation code of the connection.
      integer(I4P),       intent(in)    :: d              !< Split direction on block B (1=i, 2=j, 3=k).
      integer(I4P),       intent(in)    :: Ns1            !< First child size in split direction.
      integer(I4P),       intent(in)    :: ext_lo         !< Original low extent on block B.
      integer(I4P),       intent(in)    :: ext_hi         !< Original high extent on block B.
      integer(I4P),       intent(in)    :: Nd             !< Parent size in split direction.
      type(patch_object), intent(inout) :: new_patches(:) !< New patches array.
      integer(I4P),       intent(inout) :: new_count      !< Current count of new patches.
      integer(I4P),       intent(in)    :: B              !< Block being split.
      integer(I4P)                      :: o(3)           !< Orientation digits.
      integer(I4P)                      :: p_axis         !< Partner axis (1=i,2=j,3=k).
      integer(I4P)                      :: p_sign         !< Partner axis sign (+1 or -1).
      integer(I4P)                      :: p_d_lo, p_d_hi !< Partner extent indices.
      integer(I4P)                      :: p_ext_lo       !< Partner extent low.
      integer(I4P)                      :: p_ext_hi       !< Partner extent high.
      integer(I4P)                      :: split_pos      !< Split position on partner.
      integer(I4P)                      :: child1_pidx    !< Index of child-1 patch in new_patches.
      integer(I4P)                      :: child2_pidx    !< Index of child-2 patch in new_patches.

      if (partner_idx < 1 .or. partner_idx > np_total) return

      ! Decode orientation: each digit tells which axis of block B maps to which axis of partner
      o(1) = orient / 100         ! hundreds digit: partner's i-axis relation
      o(2) = mod(orient/10, 10)   ! tens digit: partner's j-axis relation
      o(3) = mod(orient, 10)      ! ones digit: partner's k-axis relation

      ! Find which partner axis corresponds to direction d on block B
      ! o(axis) encodes: 1=+I, 2=+J, 3=+K, 4=-I, 5=-J, 6=-K
      ! We need: which axis index in partner corresponds to d
      call find_partner_axis(o, d, p_axis, p_sign)

      if (p_axis == 0) return ! should not happen

      p_d_lo = 2*p_axis - 1
      p_d_hi = 2*p_axis

      p_ext_lo = all_patches(partner_idx)%ijk_extents(p_d_lo)
      p_ext_hi = all_patches(partner_idx)%ijk_extents(p_d_hi)

      ! Compute split position on partner
      if (p_sign > 0) then
         ! Same direction: split at position (Ns1 - ext_lo) from partner's low extent
         split_pos = p_ext_lo + (Ns1 - ext_lo)
      else
         ! Reversed: split at position from partner's high extent
         split_pos = p_ext_hi - (Ns1 - ext_lo)
      endif

      ! Mark original partner as consumed (will be replaced by two new patches)
      ! We do this by not copying the original partner in the main loop -
      ! Instead, we add two split copies here

      ! child1_pidx and child2_pidx are the indices of the patches in new_patches
      ! that correspond to child 1 and child 2 of block B
      child1_pidx = new_count - 1  ! the child-1 patch we just added
      child2_pidx = new_count      ! the child-2 patch we just added

      ! Partner part 1: connects to child 1 of B
      new_count = new_count + 1
      new_patches(new_count) = all_patches(partner_idx)
      if (new_patches(new_count)%block_index > B) then
         new_patches(new_count)%block_index = new_patches(new_count)%block_index + 1
      endif
      if (p_sign > 0) then
         new_patches(new_count)%ijk_extents(p_d_hi) = split_pos
      else
         new_patches(new_count)%ijk_extents(p_d_lo) = split_pos
      endif

      ! Partner part 2: connects to child 2 of B (at B+1)
      new_count = new_count + 1
      new_patches(new_count) = all_patches(partner_idx)
      if (new_patches(new_count)%block_index > B) then
         new_patches(new_count)%block_index = new_patches(new_count)%block_index + 1
      endif
      if (p_sign > 0) then
         new_patches(new_count)%ijk_extents(p_d_lo) = split_pos
      else
         new_patches(new_count)%ijk_extents(p_d_hi) = split_pos
      endif

      ! Set up direct cross-links between the 4 split patches.
      ! child1 <-> partner_part1, child2 <-> partner_part2
      new_patches(child1_pidx)%connect_family  = new_count - 1  ! partner part 1
      new_patches(child2_pidx)%connect_family  = new_count      ! partner part 2
      new_patches(new_count - 1)%connect_family = child1_pidx   ! child 1
      new_patches(new_count)%connect_family     = child2_pidx   ! child 2
      ! Mark all 4 as pre-linked (patch_index=0) so renumber_connect_families skips them
      new_patches(child1_pidx)%patch_index  = 0
      new_patches(child2_pidx)%patch_index  = 0
      new_patches(new_count - 1)%patch_index = 0
      new_patches(new_count)%patch_index     = 0

      ! Mark original partner patch so it won't be copied in the main loop
      all_patches(partner_idx)%block_index = -1  ! sentinel to skip
      endsubroutine split_partner_patch

      subroutine find_partner_axis(o, d, p_axis, p_sign)
      !< Find which partner axis corresponds to direction d on the source block.
      !< o(axis) encodes: 1=+I, 2=+J, 3=+K, 4=-I, 5=-J, 6=-K
      !< This means: axis 'axis' of the partner maps to direction o(axis) on the source block.
      !< We want: which partner axis maps to direction d on the source block.
      integer(I4P), intent(in)  :: o(3)    !< Orientation digits.
      integer(I4P), intent(in)  :: d       !< Direction on source block (1=i, 2=j, 3=k).
      integer(I4P), intent(out) :: p_axis  !< Partner axis (1=i, 2=j, 3=k).
      integer(I4P), intent(out) :: p_sign  !< Sign (+1 or -1).
      integer(I4P)              :: a       !< Counter.
      integer(I4P)              :: mapped_dir !< Direction that this partner axis maps to.
      integer(I4P)              :: mapped_sign !< Sign of the mapping.

      p_axis = 0
      p_sign = 1
      do a = 1, 3
         if (o(a) <= 3) then
            mapped_dir  = o(a)
            mapped_sign = 1
         else
            mapped_dir  = o(a) - 3
            mapped_sign = -1
         endif
         if (mapped_dir == d) then
            p_axis = a
            p_sign = mapped_sign
            return
         endif
      enddo
      endsubroutine find_partner_axis

      subroutine renumber_connect_families(patches, n)
      !< Renumber connect_family references after patch splitting.
      !< Patches with patch_index==0 are "pre-linked" (their connect_family already points
      !< to the correct position in the patches array) and are skipped during remapping.
      !< Other patches have connect_family pointing to old patch indices that need remapping.
      type(patch_object), intent(inout) :: patches(:) !< Patches array.
      integer(I4P),       intent(in)    :: n          !< Number of patches.
      integer(I4P)                      :: ip         !< Counter.
      integer(I4P)                      :: old_idx    !< Old patch index.
      integer(I4P)                      :: cf         !< Connect family.
      integer(I4P)                      :: max_old    !< Maximum old patch index.
      integer(I4P), allocatable         :: old_to_new(:) !< Mapping: old patch index -> new position.
      integer(I4P), allocatable         :: saved_old(:)  !< Saved old patch indices.

      ! Save old patch indices before reassignment
      allocate(saved_old(n))
      do ip = 1, n
         saved_old(ip) = patches(ip)%patch_index
      enddo

      ! Find maximum old patch index (only from non-pre-linked patches)
      max_old = 0
      do ip = 1, n
         if (saved_old(ip) > max_old) max_old = saved_old(ip)
      enddo
      if (max_old == 0) max_old = n

      ! Build old-to-new mapping: each non-pre-linked patch maps its old index to its new position.
      ! With pre-linking, each old index should map to exactly one new position.
      allocate(old_to_new(max_old), source=0)
      do ip = 1, n
         old_idx = saved_old(ip)
         if (old_idx >= 1 .and. old_idx <= max_old) then
            old_to_new(old_idx) = ip
         endif
      enddo

      ! Assign sequential patch_index
      do ip = 1, n
         patches(ip)%patch_index = ip
      enddo

      ! Update connect_family references
      do ip = 1, n
         if (.not. patches(ip)%is_connection) cycle
         if (saved_old(ip) == 0) cycle  ! pre-linked patch, connect_family already correct
         cf = patches(ip)%connect_family
         if (cf < 1 .or. cf > max_old) cycle
         if (old_to_new(cf) > 0) then
            patches(ip)%connect_family = old_to_new(cf)
         endif
      enddo

      deallocate(saved_old)
      deallocate(old_to_new)
      endsubroutine renumber_connect_families
   endsubroutine replay_splits_on_patches

   subroutine save_file_cc_par(file_name, cc_par, blocks, boxes)
   !< Save file cc.par.
   character(*),        intent(in)              :: file_name      !< File name.
   type(cc_par_object), intent(in)              :: cc_par         !< File cc.par handler.
   type(block_object),  intent(in)              :: blocks(:)      !< Blocks data.
   type(box_object),    intent(in), allocatable :: boxes(:)       !< Boxes data.
   integer(I4P)                                 :: file_unit      !< File unit.
   integer(I4P)                                 :: blocks_number  !< Blocks number.
   integer(I4P)                                 :: patches_number !< Patches number.

   blocks_number = size(blocks, dim=1)
   open(newunit=file_unit, file=trim(adjustl(file_name)), action='write', status='replace')
   call save_header
   call save_blocks
   call save_patches
   call save_edges
   call save_boxes
   call save_circuits
   close(file_unit)
   contains
      subroutine save_header()
      !< Save header section.

      write(file_unit, '(A)') "'split-balanced-cc.grd'"
      write(file_unit, '(A)') "'"//trim(cc_par%base_name_output)//"'"
      write(file_unit, '(A)') merge('.true. ', '.false.', cc_par%save_ghost_cells)
      write(file_unit, '(A)') merge('.true. ', '.false.', cc_par%increase_overlap)
      write(file_unit, '(A)') merge('.true. ', '.false.', cc_par%extend_internal_wall)
      write(file_unit, '(A)') ''
      write(file_unit, '(I3,1X,I3)') cc_par%mgl(1), cc_par%mgl(2)
      write(file_unit, '(A)') ''
      write(file_unit, '(ES12.4)') cc_par%boundary_layer_thickness
      write(file_unit, '(ES12.4)') cc_par%numberical_beach
      write(file_unit, '(A)') ''
      endsubroutine save_header

      subroutine save_blocks()
      !< Save blocks section.
      integer(I4P) :: b !< Counter.

      write(file_unit, '(I0,A)') blocks_number, char(9)//'! blocks number'
      write(file_unit, '(A)') ''
      do b = 1, blocks_number
         write(file_unit, '(I9,1X,I9,1X,I9,1X,A)') blocks(b)%level, blocks(b)%group, blocks(b)%priority, &
            '! '//trim(adjustl(blocks(b)%comment))
      enddo
      write(file_unit, '(A)') ''
      endsubroutine save_blocks

      subroutine save_patches()
      !< Save patches section.
      integer(I4P) :: b, p     !< Counters.
      integer(I4P) :: global_p !< Global patch counter.

      patches_number = 0
      do b = 1, blocks_number
         if (allocated(blocks(b)%patches)) patches_number = patches_number + size(blocks(b)%patches)
      enddo

      write(file_unit, '(I0,A)') patches_number, char(9)//'! patches number'
      write(file_unit, '(A)') ''
      global_p = 0
      do b = 1, blocks_number
         if (.not.allocated(blocks(b)%patches)) cycle
         do p = 1, size(blocks(b)%patches)
            global_p = global_p + 1
            write(file_unit, '(10(I9,1X),A,I9)')        &
               blocks(b)%patches(p)%block_index,        &
               blocks(b)%patches(p)%face_index,         &
               blocks(b)%patches(p)%boundary_condition, &
               blocks(b)%patches(p)%connect_family,     &
               blocks(b)%patches(p)%ijk_extents(1:6),   &
               '! ', global_p
         enddo
      enddo
      write(file_unit, '(A)') ''
      endsubroutine save_patches

      subroutine save_edges()
      !< Save edges section (not yet supported).

      write(file_unit, '(I0,A)') 0, '! edges number'
      write(file_unit, '(A)') ''
      endsubroutine save_edges

      subroutine save_boxes()
      !< Save boxes section.
      integer(I4P) :: b, n         !< Counter.
      integer(I4P) :: boxes_number !< Boxes number.

      if (allocated(boxes)) then
         boxes_number = size(boxes)
      else
         boxes_number = 0
      endif

      write(file_unit, '(I0,A)') boxes_number, ' ! boxes number'
      write(file_unit, '(A)') ''
      if (boxes_number > 0) then
         do b = 1, boxes_number
            write(file_unit, '(I9,1X,I9,1X,I9,1X,A)') boxes(b)%btype, boxes(b)%bblock, boxes(b)%bgroup, &
               ' ! type, block, group associated'
            do n = 1, 8
               write(file_unit, '(3(ES23.12,1X))') boxes(b)%nodes(1,n), boxes(b)%nodes(2,n), boxes(b)%nodes(3,n)
            enddo
            write(file_unit, '(A)') ''
         enddo
      endif
      endsubroutine save_boxes

      subroutine save_circuits()
      !< Save circuits section (not yet supported).
      ! Only write if there are no boxes trailing (file ends after boxes)
      ! Actually, the original format doesn't have circuits after boxes in this file
      ! But to be safe and compatible, omit circuits if the original didn't have them
      endsubroutine save_circuits
   endsubroutine save_file_cc_par
endmodule oe_block_object

module oe_process_object
!< Overset-Exploded, definition of class process_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32

implicit none

private
public :: process_object

type :: process_object
   !< Class process object.
   integer(I4P)              :: nb=0        !< Number of assigned blocks.
   integer(I4P), allocatable :: blocks(:)   !< List of assigned blocks.
   integer(I4P)              :: w=0         !< Process workload.
   integer(I4P)              :: unbalance=0 !< Unbalance respect ideal workload.
   contains
      ! public methods
      procedure, pass(self) :: assign_block !< Assign block to process.
      procedure, pass(self) :: initialize   !< Initialize process data.
endtype process_object
contains
   ! public methods
   pure subroutine assign_block(self, ab, wb, ideal_workload)
   !< Assign block to process.
   class(process_object), intent(inout) :: self           !< Process.
   integer(I4P),          intent(in)    :: ab             !< Absolute block index.
   integer(I4P),          intent(in)    :: wb             !< Block weight.
   integer(I4P),          intent(in)    :: ideal_workload !< Ideal process workload.

   self%nb        = self%nb + 1
   self%blocks    = [self%blocks,ab]
   self%w         = self%w + wb
   self%unbalance = nint(((ideal_workload*1._R8P-self%w*1._R8P)/ideal_workload*1._R8P)*100._R8P)
   endsubroutine assign_block

   elemental subroutine initialize(self)
   !< Initialize process data.
   class(process_object), intent(inout) :: self !< Process.

   self%nb = 0_I4P
   self%blocks = [0_I4P]
   self%w = 0_I4P
   self%unbalance = 0_I4P
   endsubroutine initialize
endmodule oe_process_object

module oe_library
!< Overset-Exploded, module library.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object
use oe_process_object

implicit none

private
public :: balance_workload
public :: str, strz

! interface for auxiliary procedures
interface str
  !< Convert number (real and integer) to string (number to string type casting).
  procedure str_I4P, str_a_I4P
endinterface

interface strz
  !< Convert integer, to string, prefixing with the right number of zeros (integer to string type casting with zero padding).
  procedure strz_I4P
endinterface
contains
   ! public methods
   subroutine balance_workload(file_name,use_cc_par,mgl,max_unbalance,procs_number,save_bsplit_par,&
                              processes,splits,splits_dir,splits_nijk)
   !< Balance workload distributing blocks (eventually splitted) over processes.
   !< The balanced per-process blocks assignment is returned as well the eventually blocks-splits.
   !< From file grd (or icc) only the blocks dimensions are loaded, other data are not allocated.
   character(*),         intent(in)                 :: file_name          !< File name
   logical,              intent(in)                 :: use_cc_par         !< Sentinel to enable use of cc.par.
   integer(I4P),         intent(in)                 :: mgl                !< Multigrid level to be preserved.
   integer(I4P),         intent(in)                 :: max_unbalance      !< Maximum processes unbalancing in percent.
   integer(I4P),         intent(in)                 :: procs_number       !< Number of processes for load balancing.
   logical,              intent(in)                 :: save_bsplit_par    !< Sentinel to save bsplit.par.
   type(process_object), intent(inout)              :: processes(0:)      !< Processes data.
   integer(I4P),         intent(inout), allocatable :: splits(:)          !< Splits history.
   integer(I4P),         intent(inout), allocatable :: splits_dir(:)      !< Splits direction history.
   integer(I4P),         intent(inout), allocatable :: splits_nijk(:)     !< Splits nijk history.
   type(block_object), allocatable                  :: blocks(:)          !< Blocks data.
   integer(I4P)                                     :: blocks_number      !< Blocks number.
   integer(I4P)                                     :: file_unit          !< File unit.
   integer(I4P)                                     :: total_blocks_weight!< Total blocks weight.
   integer(I4P)                                     :: ideal_proc_workload!< Ideal process weight for load balancing.
   logical                                          :: is_split_done      !< Sentinel to check is split has been done.
   type(block_object)                               :: sb(2)              !< Split blocks.
   integer(I4P), allocatable                        :: blocks_list(:)     !< Blocks (unassigned) list (decreasing-workload) ordered.
   integer(I4P)                                     :: gc                 !< Number of ghost cells in input file.
   integer(I4P)                                     :: b, bb, p           !< Counter.

   gc = 2 ; if (use_cc_par) gc = 0

   splits      = [0_I4P]
   splits_dir  = [0_I4P]
   splits_nijk = [0_I4P]

   ! load original blocks dimensions
   if (allocated(blocks)) deallocate(blocks)
   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
   read(file_unit, end=10, err=10) blocks_number
   allocate(blocks(1:blocks_number))
   do b=1, blocks_number
      call blocks(b)%load_dimensions(file_unit=file_unit, ab=b, gc=gc)
   enddo
   10 close(file_unit)

   ! load balancing
   print '(A)', 'work load balancing stats'
   total_blocks_weight = 0
   do b=1, blocks_number
      print '(A)', '    block "'//trim(strz(b,9))//'" weight: '//trim(str(blocks(b)%w,.true.))
      total_blocks_weight = total_blocks_weight + blocks(b)%w
   enddo
   ideal_proc_workload = total_blocks_weight / procs_number
   print '(A)', 'total work load "'//trim(str(total_blocks_weight))//'"'
   print '(A)', 'ideal work load for np "'//trim(strz(procs_number,6))//'" processes: '//trim(str(ideal_proc_workload,.true.))

   call create_blocks_list(blocks=blocks, blocks_list=blocks_list)
   print '(A)', 'blocks list in decreasing-workload-order'
   do b=1, blocks_number
      bb = blocks_list(b)
      print '(A)', '  block "'//trim(strz(bb,9))//'" weight: '//trim(str(blocks(bb)%w,.true.))//&
         ' Ni,Nj,Nk: '//trim(str([blocks(bb)%Ni,blocks(bb)%Nj,blocks(bb)%Nk]))
   enddo

   ! assign blocks to processes
   assign_blocks_loop : do while(allocated(blocks_list))
      p = minloc(processes(0:)%w,dim=1)-1 ! process with minimum workload
      b = blocks_list(1)                  ! first blocks in unassigned list, the current biggest block
      if (processes(p)%w+blocks(b)%w<=ideal_proc_workload*(100._R8P+max_unbalance)/100._R8P) then
         blocks(b)%proc = p
         call processes(p)%assign_block(ab=b, wb=blocks(b)%w, ideal_workload=ideal_proc_workload)
         call popout_blocks_list(blocks_list=blocks_list)
      else
         print '(A)', 'block "'//trim(strz(blocks(b)%ab,9))//'" must be split to be insert into process '//trim(strz(p,6))
         call blocks(b)%split(mgl=mgl, is_split_done=is_split_done, sb=sb)
         if (is_split_done) then
            print '(A)', '   block "'//trim(strz(b,9))//'" split'
            print '(A)', '      first split block  (ni,nj,nk) '//trim(str([sb(1)%Ni,sb(1)%Nj,sb(1)%Nk]))
            print '(A)', '      second split block (ni,nj,nk) '//trim(str([sb(2)%Ni,sb(2)%Nj,sb(2)%Nk]))
            print '(A)', '      first block parents list      '//trim(str(sb(1)%parents,.true.))
            print '(A)', '      second block parents list     '//trim(str(sb(2)%parents,.true.))
            print '(A)', '      update blocks data'
            ! recreate unassigned blocks list and reset processes data, thus the blocks assignment restart
            call update_blocks(blocks=blocks, sb=sb, blocks_number=blocks_number, use_cc_par=use_cc_par)
            call create_blocks_list(blocks=blocks, blocks_list=blocks_list)
            call processes%initialize
            ! update splits history
            splits     = [splits, sb(1)%ab]
            splits_dir = [splits_dir, sb(1)%split_dir]
            select case(sb(1)%split_dir)
            case(1)
               splits_nijk = [splits_nijk, sb(1)%ni]
            case(2)
               splits_nijk = [splits_nijk, sb(1)%nj]
            case(3)
               splits_nijk = [splits_nijk, sb(1)%nk]
            endselect
         else
            print '(A)', 'block "'//trim(strz(blocks(b)%ab,9))//'" split failed, assigned anyway to process '//trim(strz(p,6))
            blocks(b)%proc = p
            call processes(p)%assign_block(ab=b, wb=blocks(b)%w, ideal_workload=ideal_proc_workload)
            call popout_blocks_list(blocks_list=blocks_list)
         endif
      endif
   enddo assign_blocks_loop
   if (size(splits,dim=1)>1) then ! trim out first block set to 0 for convenience
      splits      = splits(     2:)
      splits_dir  = splits_dir( 2:)
      splits_nijk = splits_nijk(2:)
   else ! no splits necessary, destroy splits history
      deallocate(splits)
      deallocate(splits_dir)
      deallocate(splits_nijk)
   endif

   print '(A)', 'processes workload'
   do p=0, procs_number-1
      print '(A)', '  proc '//trim(strz(p,6))//&
                   ' unbalancing '//trim(str(processes(p)%unbalance))//&
                   '% assigned blocks '//trim(str(processes(p)%blocks(2:),.true.))
   enddo
   if (save_bsplit_par.and.allocated(splits)) call save_file_bsplit_par
   contains
      subroutine save_file_bsplit_par
      !< Save file bsplit.par.
      integer(I4P) :: file_unit !< File unit.
      integer(I4P) :: s         !< Counter.

      open(newunit=file_unit, file='bsplit.par', action='write', status='replace')
      write(file_unit,'(A)') "'none' ! basename file solution"
      write(file_unit,'(A)') "'cc'   ! basename grd file"
      write(file_unit,'(A)') "'cc'   ! basename icc file"
      write(file_unit,'(A)') trim(str(mgl))//" ! multigrid level"
      write(file_unit,'(A)') "0 ! variables time n"
      write(file_unit,'(A)') "0 ! variables time n-1"
      write(file_unit,'(A)') "5 ! debug level"
      do s=1, size(splits,dim=1)
         write(file_unit,'(A)') trim(str([splits(s),splits_dir(s),splits_nijk(s)]))
      enddo
      close(file_unit)
      endsubroutine save_file_bsplit_par
   endsubroutine balance_workload

   ! string procedures
   elemental function str_I4P(n, no_sign) result(str)
   !< Converting integer to string.
   integer(I4P), intent(in)           :: n       !< Integer to be converted.
   logical,      intent(in), optional :: no_sign !< Flag for leaving out the sign.
   character(11)                      :: str     !< Returned string containing input number plus padding zeros.

   write(str, '(I11)') n             ! Casting of n to string.
   str = adjustl(trim(str))          ! Removing white spaces.
   if (n>=0_I4P) str='+'//trim(str)  ! Prefixing plus if n>0.
   if (present(no_sign)) str=str(2:) ! Leaving out the sign.
   endfunction str_I4P

   elemental function strz_I4P(n, nz_pad) result(str)
   !< Convert integer to string, prefixing with the right number of zeros.
   integer(I4P), intent(in)           :: n      !< Integer to be converted.
   integer(I4P), intent(in), optional :: nz_pad !< Number of zeros padding.
   character(11)                      :: str    !< Returned string containing input number plus padding zeros.

   write(str,'(I11.10)') n                      ! Casting of n to string.
   str=str(2:)                                  ! Leaving out the sign.
   if (present(nz_pad)) str=str(11-nz_pad:11-1) ! Leaving out the extra zeros padding
   endfunction strz_I4P

   pure function str_a_I4P(n, no_sign, separator, delimiters) result(str)
   !< Convert integer array to string.
   integer(I4P), intent(in)           :: n(:)            !< Integer array to be converted.
   logical,      intent(in), optional :: no_sign         !< Flag for leaving out the sign.
   character(1), intent(in), optional :: separator       !< Eventual separator of array values.
   character(*), intent(in), optional :: delimiters(1:2) !< Eventual delimiters of array values.
   character(len=:), allocatable      :: str             !< Returned string containing input number.
   character(11)                      :: strn            !< String containing of element of input array number.
   character(len=1)                   :: sep             !< Array values separator
   integer                            :: i               !< Counter.

   str = ''
   sep = ','
   if(present(separator)) sep = separator
   if (present(no_sign)) then
     do i=1,size(n)
       strn = str_I4P(no_sign=no_sign, n=n(i))
       str = str//sep//trim(strn)
     enddo
   else
     do i=1,size(n)
       strn = str_I4P(n=n(i))
       str = str//sep//trim(strn)
     enddo
   endif
   str = trim(str(2:))
   if (present(delimiters)) str = delimiters(1)//str//delimiters(2)
   endfunction str_a_I4P
endmodule oe_library

program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object
use oe_process_object
use oe_library

implicit none

character(len=99)                 :: file_name_grd        !< Grid file name.
character(len=99)                 :: file_name_icc        !< Icc file name.
character(len=99)                 :: file_name_input      !< Input file name: grd, icc or infocc.out.
character(len=99)                 :: file_name_proc_input !< Name of proc.input file.
logical                           :: save_block_tecplot   !< Save blocks also in tecplot (ASCII) format.
logical                           :: save_imploded        !< Save imploded blocks after explosion.
logical                           :: save_exploded        !< Save exploded blocks.
logical                           :: save_bsplit_par      !< Save bsplit.par.
logical                           :: use_cc_par           !< Use cc.par input instead of overset output.
character(len=99)                 :: exploded_basename    !< Exploded files basebame.
integer(I4P)                      :: mgl                  !< Multigrid level to be preserved.
integer(I4P)                      :: blocks_number        !< Number of blocks contained into the files.
type(block_object), allocatable   :: blocks(:)            !< Blocks data.
type(cc_par_object)               :: cc_par               !< File cc.par handler.
type(box_object), allocatable     :: boxes(:)             !< Boxes.
real(R4P), allocatable            :: rcc(:)               !< rcc unstructured array.
integer(I4P)                      :: procs_number         !< Number of processes for load balancing.
type(process_object), allocatable :: processes(:)         !< Processes data.
integer(I4P), allocatable         :: splits(:)            !< Splits history.
integer(I4P), allocatable         :: splits_dir(:)        !< Splits direction history.
integer(I4P), allocatable         :: splits_nijk(:)       !< Splits nijk history.
logical                           :: is_split_done        !< Sentinel to check is split has been done.
type(block_object)                :: sb(2)                !< Split blocks.
integer(I4P)                      :: max_unbalance        !< Maximum processes unbalancing in percent.
integer(I4P)                      :: i,b,bb,p             !< Counter.
character(99)                     :: fname_cc_par         !< Overset input cc.par file name.

call parse_command_line(fgrd=file_name_grd,ficc=file_name_icc,fpci=file_name_proc_input,                   &
                        stec=save_block_tecplot,simp=save_imploded,sexp=save_exploded,sbsp=save_bsplit_par,&
                        uccp=use_cc_par,fccp=fname_cc_par,ebn=exploded_basename,np=procs_number,mu=max_unbalance,mgl=mgl)

if (use_cc_par) then
   ! use grid file (from geogrd without ghost cells) or infocc.out from infocc
   if     (is_file_found(file_name_grd).and.is_file_found(fname_cc_par)) then
      file_name_input = trim(adjustl(file_name_grd))
   else
      write(stderr, '(A)')'error: file "'//trim(adjustl(file_name_grd))//'" or '//&
                                      '"'//trim(adjustl(fname_cc_par))//'" not found!'
      stop
   endif
else
   ! use grid or icc file from overset
   if     (is_file_found(file_name_grd)) then
      file_name_input = trim(adjustl(file_name_grd))
   elseif (is_file_found(file_name_icc)) then
      file_name_input = trim(adjustl(file_name_icc))
   else
      write(stderr, '(A)')'error: file "'//trim(adjustl(file_name_grd))//'" or '//&
                                      '"'//trim(adjustl(file_name_icc))//'" not found!'
      stop
   endif
endif

allocate(processes(0:procs_number-1))
call processes%initialize

print '(A)', 'perform work load balancing using file '//trim(adjustl(file_name_input))
call balance_workload(file_name=file_name_input,use_cc_par=use_cc_par,mgl=mgl,                              &
                      max_unbalance=max_unbalance,procs_number=procs_number,save_bsplit_par=save_bsplit_par,&
                      processes=processes,splits=splits,splits_dir=splits_dir,splits_nijk=splits_nijk)

if (use_cc_par) then
   print '(A)', 'load blocks grid from file '//trim(adjustl(file_name_input))
   call load_file_grd(file_name=file_name_input,blocks=blocks,blocks_number=blocks_number,gc=0)
   print '(A)', 'load parameter file '//trim(adjustl(fname_cc_par))
   call load_file_cc_par(file_name=fname_cc_par,cc_par=cc_par,blocks=blocks,boxes=boxes)
else
   print '(A)', 'load blocks grid from file '//trim(adjustl(file_name_input))
   call load_file_grd(file_name=file_name_input,blocks=blocks,blocks_number=blocks_number)
   print '(A)', 'load icc file '//trim(adjustl(file_name_icc))
   call load_file_icc(file_name=file_name_icc,blocks=blocks,blocks_number=blocks_number,rcc=rcc)
   print '(A)', 'parse global rcc and create block-local-rcc'
   do b=1, blocks_number
      call blocks(b)%parse_rcc(rcc=rcc)
      print '(A)', 'block '//trim(str(b,.true.))//' BC chimera cells number: '//trim(str(size(blocks(b)%chimera,dim=1)))
   enddo
endif

if (use_cc_par) then
   ! if (allocated(splits)) then
   !    print '(A)', 'replay splits on patches'
   !    flush(6)
   !    call replay_splits_on_patches(blocks=blocks, splits=splits, splits_dir=splits_dir, &
   !                                  splits_nijk=splits_nijk, boxes=boxes)
   !    blocks_number = size(blocks, dim=1)
   !    print '(A)', 'new blocks number after split replay: '//trim(str(blocks_number,.true.))
   ! endif
   ! print '(A)', 'save split cc.par file'
   ! call save_file_cc_par('split-balanced-cc.par', cc_par, blocks, boxes)
endif

if (allocated(splits)) then
   print '(A)', 'split blocks'
   do b=1, size(splits,dim=1)
      print '(A)', '  block '//trim(str(splits(b),.true.))
      call blocks(splits(b))%split(mgl=mgl, is_split_done=is_split_done, sb=sb, split_data=.true., use_cc_par=use_cc_par)
      ! debuuuuggg
      print*, 'cazzo b',splits(b),splits_dir(b)
      print*, 'cazzo                       p,          b,            face,     bc,         conn,         ijk'
      do p=1, size(blocks(splits(b))%patches,dim=1)
         print*, 'cazzo block        ',blocks(splits(b))%patches(p)%patch_index,&
                                       blocks(splits(b))%patches(p)%block_index,&
                                       blocks(splits(b))%patches(p)%face_index,&
                                       blocks(splits(b))%patches(p)%boundary_condition,&
                                       blocks(splits(b))%patches(p)%connect_family,&
                                       blocks(splits(b))%patches(p)%ijk_extents
         print*, 'cazzo split block 1',sb(1)%patches(p)%patch_index,&
                                       sb(1)%patches(p)%block_index,&
                                       sb(1)%patches(p)%face_index,&
                                       sb(1)%patches(p)%boundary_condition,&
                                       sb(1)%patches(p)%connect_family,&
                                       sb(1)%patches(p)%ijk_extents
         print*, 'cazzo split block 2',sb(2)%patches(p)%patch_index,&
                                       sb(2)%patches(p)%block_index,&
                                       sb(2)%patches(p)%face_index,&
                                       sb(2)%patches(p)%boundary_condition,&
                                       sb(2)%patches(p)%connect_family,&
                                       sb(2)%patches(p)%ijk_extents
      enddo
  ! integer(I4P) :: patch_index=0         !< Patch index, local numeration.
  ! integer(I4P) :: block_index=0         !< Block to which patch belongs, local numeration.
  ! integer(I4P) :: face_index=0          !< Face index in the overset convention: Imin=>1, Imax=2,Jmin=3, Jmax=4, Kmin=5, Kmax=6.
  ! integer(I4P) :: boundary_condition=0  !< Boundary condition (or IJK orientation) in the overset convention.
  ! integer(I4P) :: connect_family=0      !< Index of connected patch or family index of the patch.
      stop
      if (is_split_done) then
         call update_blocks(blocks=blocks, sb=sb, blocks_number=blocks_number, use_cc_par=use_cc_par)
         print '(A)', '     new blocks number '//trim(str(blocks_number,.true.))
         ! free sb memory immediately after update_blocks copies the data
         call sb(1)%destroy
         call sb(2)%destroy
      else
         write(stderr,'(A)')'error: unable to split block "'//trim(str(splits(b)))//'"'
         stop
      endif
   enddo
endif
! update blocks to processes assignment
do p=0, procs_number - 1
   do b=1, size(processes(p)%blocks,dim=1)
      bb = processes(p)%blocks(b)
      if (bb>0) blocks(bb)%proc = p
   enddo
enddo
call save_proc_input(blocks=blocks, file_name=file_name_proc_input)

if (save_exploded) then
   print '(A)', 'save exploded blocks'
   do b=1, blocks_number
      call blocks(b)%save_block_file(basename=exploded_basename, tec=save_block_tecplot)
   enddo
endif

if (save_imploded) then
   print '(A)', 'implode exploded blocks'
   if (.not.use_cc_par) then
      call implode_blocks(blocks=blocks, rcc=rcc)
      print '(A)', 'total chimera elements after implosion'//trim(str(size(rcc,dim=1)))
      do b=1, blocks_number
         print '(A)', 'block '//trim(str(b,.true.))//' BC chimera cells number: '//trim(str(size(blocks(b)%chimera,dim=1)))
      enddo
      call save_file_icc(file_name='split-balanced-cc', blocks=blocks, rcc=rcc)
   endif
   call save_file_grd(file_name='split-balanced-cc.grd', blocks=blocks)
endif
contains
   subroutine parse_command_line(fgrd,ficc,fpci,stec,simp,sexp,sbsp,uccp,fccp,ebn,np,mu,mgl)
   !< Parse command line inputs.
   character(*), intent(out) :: fgrd      !< Grid file name.
   character(*), intent(out) :: ficc      !< Icc file name.
   character(*), intent(out) :: fpci      !< Name of proc.input file.
   logical,      intent(out) :: stec      !< Save blocks also in tecplot (ASCII) format.
   logical,      intent(out) :: simp      !< Save imploded blocks after explosion.
   logical,      intent(out) :: sexp      !< Save exploded blocks.
   logical,      intent(out) :: sbsp      !< Save bsplit.par.
   logical,      intent(out) :: uccp      !< Use cc.par overset input instead of overset output.
   character(*), intent(out) :: fccp      !< File name of cc.par.
   character(*), intent(out) :: ebn       !< Exploded files basename.
   integer(I4P), intent(out) :: np        !< Number of processes.
   integer(I4P), intent(out) :: mu        !< Maximum processes unbalancing in percent.
   integer(I4P), intent(out) :: mgl       !< Multigrid level to be preserved.
   integer(I4P)              :: na        !< Number of command line arguments.
   character(len=99)         :: ca_buffer !< Command argument buffer.
   integer(I4P)              :: a         !< Counter.

   ! defaults
   fgrd = 'cc.01.grd'
   ficc = 'cc.01'
   fpci = 'proc.input'
   stec = .false.
   simp = .false.
   sexp = .false.
   sbsp = .false.
   uccp = .false.
   ebn  = 'exploded-'
   np   = 1_I4P
   mu   = 1_I4P
   mgl  = 4_I4P

   na = command_argument_count()
   a = 1
   ca_loop : do
      call get_command_argument(a, ca_buffer)
      select case(trim(adjustl(ca_buffer)))
      case('-grd')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         fgrd = trim(adjustl(ca_buffer))
      case('-icc')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         ficc = trim(adjustl(ca_buffer))
      case('-proc-input')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         fpci = trim(adjustl(ca_buffer))
      case('-np')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         read(ca_buffer,*) np
      case('-max-unbalance')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         read(ca_buffer,*) mu
      case('-mgl')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         read(ca_buffer,*) mgl
      case('-tec')
         stec = .true.
      case('-save-imploded')
         simp = .true.
      case('-save-exploded')
         sexp = .true.
      case('-save-bsplit_par')
         sbsp = .true.
      case('-exploded-basename')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         ebn = trim(adjustl(ca_buffer))
      case('-use-cc-par')
         uccp = .true.
         a = a + 1
         call get_command_argument(a, ca_buffer)
         fccp = trim(adjustl(ca_buffer))
      case('-h','--help')
         call print_help
         stop
      case default
         write(stderr, "(A)")'error: command line argument "'//trim(adjustl(ca_buffer))//'" unkwnown!'
         call print_help
         stop
      endselect
      a = a + 1
      if (a>na) exit ca_loop
   enddo ca_loop
   endsubroutine parse_command_line

   subroutine print_help
   !< Print help message.
   write(*, '(A)')'overset-exploded: overset post-processor, automatic blocks-splitting, load-balancing, blocks-explosion'
   write(*, '(A)')'usage:'
   write(*, '(A)')'   overset-exploded [args]'
   write(*, '(A)')'args list:'
   write(*, '(A)')'   -grd file_name_grd               => GRD file name, default "cc.01.grd"'
   write(*, '(A)')'   -icc file_name_icc               => ICC file name, default "cc.01"'
   write(*, '(A)')'   -proc-input file_name_proc_input => proc.input file name, default "proc.input"'
   write(*, '(A)')'   -np processes_number             => number of processes for load balancing, default 1'
   write(*, '(A)')'   -max-unbalance mu                => maximum processes unbalancing in percent, default 1%'
   write(*, '(A)')'   -mgl mgl                         => multigrid level to be preserved, default 4'
   write(*, '(A)')'   -tec                             => enable tecplot output for debug, default .false.'
   write(*, '(A)')'   -save-imploded                   => save imploded blocks after explosion, default .false.'
   write(*, '(A)')'   -save-exploded                   => save exploded blocks, default .false.'
   write(*, '(A)')'   -exploded-basename               => exploded files basename, default "exploded-"'
   write(*, '(A)')'   -use-cc-par file_name_cc_par     => use cc.par overset input insteod of overset output, default .false.'
   write(*, '(A)')'   -h, --help                       => print this help message'
   write(*, '(A)')'examples:'
   write(*, '(A)')'   overset-exploded -np 32'
   write(*, '(A)')'   overset-exploded -grd cc.02.grd -icc cc.02 -np 16'
   write(*, '(A)')'   overset-exploded -grd geo.grd -use-cc-par -np 16'
   write(*, '(A)')'   overset-exploded -np 16 -max-unbalance 4'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16 -save-imploded'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16 -save-exploded'
   write(*, '(A)')'   overset-exploded -grd cc.03.grd -icc cc.03 -np 2 -tec -max-unbalance 3'
   endsubroutine print_help

   function is_file_found(file_name) result(is_found)
   !< Inquire is the file path is valid and the file is found.
   character(*), intent(in) :: file_name !< File name.
   logical                  :: is_found  !< Inquiring result.

   inquire(file=trim(adjustl(file_name)), exist=is_found)
   endfunction is_file_found
endprogram overset_exploded
