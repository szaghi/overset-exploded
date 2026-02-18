module overset_exploded_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use            :: overset_exploded_patch_object

implicit none

private
public :: block_object
public :: bc_int_type, bc_string

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
   integer(I4P)                    :: np=0              !< Patches number.
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
      procedure, pass(self) :: sanitize_patches   !< Sanitize patches data after a split.
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
   self%np        = 0
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

   pure subroutine sanitize_patches(self, sb)
   !< Sanitize patches data after a split.
   class(block_object), intent(inout) :: self              !< Block data.
   type(block_object),  intent(in)    :: sb(2)             !< Split blocks.
   integer(I4P)                       :: p,np1             !< Counter.
   integer(I4P)                       :: face_min,face_max !< Min and max face index in split direction.

   if (.not.allocated(self%patches)) return
   np1 = size(sb(1)%patches,dim=1)
   do p=1, size(self%patches,dim=1)
      ! shift connection family
      ! if (self%patches(p)%is_connection) then
      !    if     (self%patches(p)%connect_family < sb(1)%patches(1  )%patch_index) then
      !       ! do nothing, patch index is not moved
      !    elseif (self%patches(p)%connect_family > sb(1)%patches(np1)%patch_index) then
      !       ! shift connection patch
      !       self%patches(p)%connect_family = self%patches(p)%connect_family + 6
      !    else
      !       ! there is a connection to a patch that could have been split
      !       face_min = 1 + (sb(1)%split_dir-1) * 2
      !       face_max = 2 + (sb(1)%split_dir-1) * 2
      !       if      (self%patches(p)%face_index==face_min) then
      !          ! do nothing, patch index is not moved
      !       elseif  (self%patches(p)%face_index==face_max) then
      !          ! shift connection patch
      !          self%patches(p)%connect_family = self%patches(p)%connect_family + 6
      !       else
      !          ! the connection is to a patch that has been split, we need to ammend connection extents and add new sub-patch
      !       endif
      !    endif
      ! endif
      ! shift patch index
      ! if (self%patches(p)%patch_index>sb(1)%patches(np1)%patch_index) then
      !    self%patches(p)%patch_index = self%patches(p)%patch_index + 6
      ! endif
      ! shift patch block index
      if (self%patches(p)%block_index>sb(1)%ab) then
         self%patches(p)%block_index = self%patches(p)%block_index + 1
      endif
   enddo
   endsubroutine sanitize_patches

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

   subroutine split(self, mgl, is_split_done, sb, split_data, use_cc_par)
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
      endassociate
      endsubroutine split_chimera

      subroutine split_patches(patches,delta,sb)
      !< Split patches.
      !< Distributes parent patches to the two sub-blocks and creates new adjacent connection patches at the split interface.
      !< Perpendicular face patches are clipped/split based on their extent overlap with each sub-block.
      !< New adjacent patches store -(sibling_ab * 10 + face_index) in connect_family. The sanitize_patches
      !< subroutine uses this to find the partner patch on the sibling block and establish cross-references.
      type(patch_object), intent(in)    :: patches(1:) !< Patches on block faces (6 patches or more if faces have split patches).
      integer(I4P),       intent(in)    :: delta(3)    !< Deltas.
      type(block_object), intent(inout) :: sb(2)       !< Split blocks.
      integer(I4P)                      :: p1, p2, p   !< Counter.
      integer(I4P)                      :: face_min    !< Face index of min face in split direction (1, 3, or 5).
      integer(I4P)                      :: face_max    !< Face index of max face in split direction (2, 4, or 6).
      integer(I4P)                      :: ext_lo      !< Index of low extent in split direction (1, 3, or 5).
      integer(I4P)                      :: ext_hi      !< Index of high extent in split direction (2, 4, or 6).
      integer(I4P)                      :: ns1         !< Number of cells in split direction for sb(1).
      integer(I4P)                      :: np          !< Number of parent patches.
      type(patch_object), allocatable   :: patches1(:) !< Patches buffer for sb(1).
      type(patch_object), allocatable   :: patches2(:) !< Patches buffer for sb(2).

      np = size(patches, dim=1)

      ! determine split direction parameters
      if     (delta(1)==1_I4P) then
         face_min = 1 ; face_max = 2 ; ext_lo = 1 ; ext_hi = 2 ; ns1 = sb(1)%Ni
      elseif (delta(2)==1_I4P) then
         face_min = 3 ; face_max = 4 ; ext_lo = 3 ; ext_hi = 4 ; ns1 = sb(1)%Nj
      else
         face_min = 5 ; face_max = 6 ; ext_lo = 5 ; ext_hi = 6 ; ns1 = sb(1)%Nk
      endif

!     ! first pass: count patches for each sub-block (1 extra each for new adjacent patch at split interface)
!     p1 = 1
!     p2 = 1
!     do p=1, np
!        if     (patches(p)%face_index==face_min) then
!           p1 = p1 + 1 ! min face belongs entirely to sb(1)
!        elseif (patches(p)%face_index==face_max) then
!           p2 = p2 + 1 ! max face belongs entirely to sb(2)
!        else
!           ! perpendicular face: check extent overlap with each sub-block
!           if (patches(p)%ijk_extents(ext_lo) < ns1) p1 = p1 + 1
!           if (patches(p)%ijk_extents(ext_hi) > ns1) p2 = p2 + 1
!        endif
!     enddo
!     print*,'cazzo np,p1,p2',np,p1,p2
!
!     allocate(patches1(1:p1))
!     allocate(patches2(1:p2))
!
!     ! second pass: distribute patches to sub-blocks
!     p1 = 0
!     p2 = 0
!     do p=1, np
!        if     (patches(p)%face_index==face_min) then
!           ! min face patch belongs entirely to sb(1), sb(2) has the new connection in this patch
!           p1 = p1 + 1
!           patches1(p1)             = patches(p)
!           patches1(p1)%block_index = sb(1)%ab
!           ! Update extents to sb(1) dimensions (min face: extent in split dir stays 0)
!           patches1(p1)%ijk_extents(1) = min(patches1(p1)%ijk_extents(1), sb(1)%Ni)
!           patches1(p1)%ijk_extents(2) = min(patches1(p1)%ijk_extents(2), sb(1)%Ni)
!           patches1(p1)%ijk_extents(3) = min(patches1(p1)%ijk_extents(3), sb(1)%Nj)
!           patches1(p1)%ijk_extents(4) = min(patches1(p1)%ijk_extents(4), sb(1)%Nj)
!           patches1(p1)%ijk_extents(5) = min(patches1(p1)%ijk_extents(5), sb(1)%Nk)
!           patches1(p1)%ijk_extents(6) = min(patches1(p1)%ijk_extents(6), sb(1)%Nk)
!
!           p2 = p2 + 1
!           patches2(p2)%is_connection       = .true.
!           patches2(p2)%patch_index         = patches(p)%patch_index
!           patches2(p2)%block_index         = sb(2)%ab
!           patches2(p2)%face_index          = face_min
!           patches2(p2)%boundary_condition  = -123
!           patches2(p2)%connect_family      = patches(p)%patch_index + 1
!           patches2(p2)%ijk_extents         = [0, sb(2)%Ni, 0, sb(2)%Nj, 0, sb(2)%Nk]
!           patches2(p2)%ijk_extents(ext_lo) = 0
!           patches2(p2)%ijk_extents(ext_hi) = 0
!        elseif (patches(p)%face_index==face_max) then
!           ! max face patch belongs entirely to sb(2), sb(1) has the new connection in this patch
!           p2 = p2 + 1
!           patches2(p2)             = patches(p)
!           patches2(p2)%block_index = sb(2)%ab
!           ! Update extents to sb(2) dimensions
!           ! First shift the split-direction extents, then cap all to sb(2) dimensions
!           patches2(p2)%ijk_extents(ext_lo) = patches2(p2)%ijk_extents(ext_lo) - ns1
!           patches2(p2)%ijk_extents(ext_hi) = patches2(p2)%ijk_extents(ext_hi) - ns1
!           patches2(p2)%ijk_extents(1) = min(patches2(p2)%ijk_extents(1), sb(2)%Ni)
!           patches2(p2)%ijk_extents(2) = min(patches2(p2)%ijk_extents(2), sb(2)%Ni)
!           patches2(p2)%ijk_extents(3) = min(patches2(p2)%ijk_extents(3), sb(2)%Nj)
!           patches2(p2)%ijk_extents(4) = min(patches2(p2)%ijk_extents(4), sb(2)%Nj)
!           patches2(p2)%ijk_extents(5) = min(patches2(p2)%ijk_extents(5), sb(2)%Nk)
!           patches2(p2)%ijk_extents(6) = min(patches2(p2)%ijk_extents(6), sb(2)%Nk)
!
!           p1 = p1 + 1
!           patches1(p1)%is_connection       = .true.
!           patches1(p1)%patch_index         = patches(p)%patch_index
!           patches1(p1)%block_index         = sb(1)%ab
!           patches1(p1)%face_index          = face_max
!           patches1(p1)%boundary_condition  = -123
!           patches1(p1)%connect_family      = patches(p)%patch_index + 5
!           patches1(p1)%ijk_extents         = [0, sb(1)%Ni, 0, sb(1)%Nj, 0, sb(1)%Nk]
!           patches1(p1)%ijk_extents(ext_lo) = ns1
!           patches1(p1)%ijk_extents(ext_hi) = ns1
        ! first pass: count patches for each sub-block
        p1 = 0
        p2 = 0
        do p=1, np
           if     (patches(p)%face_index==face_min) then
              p1 = p1 + 1 ! min face belongs entirely to sb(1)
           elseif (patches(p)%face_index==face_max) then
              p2 = p2 + 1 ! max face belongs entirely to sb(2)
           else
              ! perpendicular face: check extent overlap with each sub-block
              if (patches(p)%ijk_extents(ext_lo) < ns1) p1 = p1 + 1
              if (patches(p)%ijk_extents(ext_hi) > ns1) p2 = p2 + 1
           endif
        enddo
        p1 = p1 + 1 ! one adjacency at face_max of sb(1)
        p2 = p2 + 1 ! one adjacency at face_min of sb(2)
        print*,'cazzo np,p1,p2',np,p1,p2

        allocate(patches1(1:p1))
        allocate(patches2(1:p2))

        ! second pass: distribute patches to sub-blocks
        p1 = 0
        p2 = 0
        do p=1, np
           if     (patches(p)%face_index==face_min) then
              ! min face patch belongs entirely to sb(1)
              p1 = p1 + 1
              patches1(p1)             = patches(p)
              patches1(p1)%block_index = sb(1)%ab
              ! Update extents to sb(1) dimensions (min face: extent in split dir stays 0)
              patches1(p1)%ijk_extents(1) = min(patches1(p1)%ijk_extents(1), sb(1)%Ni)
              patches1(p1)%ijk_extents(2) = min(patches1(p1)%ijk_extents(2), sb(1)%Ni)
              patches1(p1)%ijk_extents(3) = min(patches1(p1)%ijk_extents(3), sb(1)%Nj)
              patches1(p1)%ijk_extents(4) = min(patches1(p1)%ijk_extents(4), sb(1)%Nj)
              patches1(p1)%ijk_extents(5) = min(patches1(p1)%ijk_extents(5), sb(1)%Nk)
              patches1(p1)%ijk_extents(6) = min(patches1(p1)%ijk_extents(6), sb(1)%Nk)

           elseif (patches(p)%face_index==face_max) then
              ! max face patch belongs entirely to sb(2)
              p2 = p2 + 1
              patches2(p2)             = patches(p)
              patches2(p2)%block_index = sb(2)%ab
              ! First shift the split-direction extents, then cap all to sb(2) dimensions
              patches2(p2)%ijk_extents(ext_lo) = patches2(p2)%ijk_extents(ext_lo) - ns1
              patches2(p2)%ijk_extents(ext_hi) = patches2(p2)%ijk_extents(ext_hi) - ns1
              patches2(p2)%ijk_extents(1) = min(patches2(p2)%ijk_extents(1), sb(2)%Ni)
              patches2(p2)%ijk_extents(2) = min(patches2(p2)%ijk_extents(2), sb(2)%Ni)
              patches2(p2)%ijk_extents(3) = min(patches2(p2)%ijk_extents(3), sb(2)%Nj)
              patches2(p2)%ijk_extents(4) = min(patches2(p2)%ijk_extents(4), sb(2)%Nj)
              patches2(p2)%ijk_extents(5) = min(patches2(p2)%ijk_extents(5), sb(2)%Nk)
              patches2(p2)%ijk_extents(6) = min(patches2(p2)%ijk_extents(6), sb(2)%Nk)
         else
            ! perpendicular face: clip extents to each sub-block range
            if (patches(p)%ijk_extents(ext_lo) < ns1) then
               ! overlaps with sb(1): extent range [ext_lo, min(ext_hi, ns1)]
               p1 = p1 + 1
               patches1(p1)             = patches(p)
               ! patches1(p1)%block_index = -sb(1)%ab
               if (patches(p)%ijk_extents(ext_hi) > ns1) then
                  patches1(p1)%block_index = -sb(1)%ab ! spans boundary -> needs remote split in connected patch
               else
                  patches1(p1)%block_index =  sb(1)%ab ! fits within sb(1) -> no remote split
               endif
               if (patches1(p1)%ijk_extents(ext_hi) > ns1) patches1(p1)%ijk_extents(ext_hi) = ns1
               ! Cap extents to sb(1) dimensions
               patches1(p1)%ijk_extents(1) = min(patches1(p1)%ijk_extents(1), sb(1)%Ni)
               patches1(p1)%ijk_extents(2) = min(patches1(p1)%ijk_extents(2), sb(1)%Ni)
               patches1(p1)%ijk_extents(3) = min(patches1(p1)%ijk_extents(3), sb(1)%Nj)
               patches1(p1)%ijk_extents(4) = min(patches1(p1)%ijk_extents(4), sb(1)%Nj)
               patches1(p1)%ijk_extents(5) = min(patches1(p1)%ijk_extents(5), sb(1)%Nk)
               patches1(p1)%ijk_extents(6) = min(patches1(p1)%ijk_extents(6), sb(1)%Nk)
            endif
            if (patches(p)%ijk_extents(ext_hi) > ns1) then
               ! overlaps with sb(2): shift extents to sb(2) local coordinates
               p2 = p2 + 1
               patches2(p2)             = patches(p)
               ! patches2(p2)%block_index = -sb(2)%ab
               if (patches(p)%ijk_extents(ext_lo) < ns1) then
                  patches2(p2)%block_index = -sb(2)%ab ! spans boundary -> needs remote split in connected patch
               else
                  patches2(p2)%block_index =  sb(2)%ab ! fits within sb(2) -> no remote split
               endif
               patches2(p2)%ijk_extents(ext_lo) = max(patches(p)%ijk_extents(ext_lo), ns1) - ns1
               patches2(p2)%ijk_extents(ext_hi) = patches(p)%ijk_extents(ext_hi) - ns1
               ! Cap extents to sb(2) dimensions
               patches2(p2)%ijk_extents(1) = min(patches2(p2)%ijk_extents(1), sb(2)%Ni)
               patches2(p2)%ijk_extents(2) = min(patches2(p2)%ijk_extents(2), sb(2)%Ni)
               patches2(p2)%ijk_extents(3) = min(patches2(p2)%ijk_extents(3), sb(2)%Nj)
               patches2(p2)%ijk_extents(4) = min(patches2(p2)%ijk_extents(4), sb(2)%Nj)
               patches2(p2)%ijk_extents(5) = min(patches2(p2)%ijk_extents(5), sb(2)%Nk)
               patches2(p2)%ijk_extents(6) = min(patches2(p2)%ijk_extents(6), sb(2)%Nk)
            endif
         endif
      enddo

      ! create one adjacency at face_max of sb(1) connecting to sb(2)
      p1 = p1 + 1
      patches1(p1)%is_connection       = .true.
      patches1(p1)%patch_index         = 0
      patches1(p1)%block_index         = sb(1)%ab
      patches1(p1)%face_index          = face_max
      patches1(p1)%boundary_condition  = -123
      patches1(p1)%connect_family      = 0
      patches1(p1)%ijk_extents         = [0, sb(1)%Ni, 0, sb(1)%Nj, 0, sb(1)%Nk]
      patches1(p1)%ijk_extents(ext_lo) = ns1
      patches1(p1)%ijk_extents(ext_hi) = ns1

      ! create one adjacency at face_min of sb(2) connecting to sb(1)
      p2 = p2 + 1
      patches2(p2)%is_connection       = .true.
      patches2(p2)%patch_index         = 0
      patches2(p2)%block_index         = sb(2)%ab
      patches2(p2)%face_index          = face_min
      patches2(p2)%boundary_condition  = -123
      patches2(p2)%connect_family      = 0
      patches2(p2)%ijk_extents         = [0, sb(2)%Ni, 0, sb(2)%Nj, 0, sb(2)%Nk]
      patches2(p2)%ijk_extents(ext_lo) = 0
      patches2(p2)%ijk_extents(ext_hi) = 0

      ! assign patch arrays to sub-blocks
      sb(1)%patches = patches1
      sb(2)%patches = patches2
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
                               lhs%np          = rhs%np
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
endmodule overset_exploded_block_object
