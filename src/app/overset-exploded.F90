module oe_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit

implicit none

private
public :: block_object
public :: bc_int_type, bc_string
public :: create_blocks_list
public :: implode_blocks
public :: load_file_grd
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

type :: block_object
   !<  Block class.
   integer(I4P)              :: Ni=0              !< Number of cells in i direction.
   integer(I4P)              :: Nj=0              !< Number of cells in j direction.
   integer(I4P)              :: Nk=0              !< Number of cells in k direction.
   integer(I4P)              :: gc=2              !< Number of ghost cells.
   integer(I4P)              :: w=0               !< Block weight (work load).
   real(R8P),    allocatable :: nodes(:,:,:,:)    !< Nodes coordinates.
   integer(I4P), allocatable :: icc(:,:,:)        !< Cell centered icc values.
   integer(I4P), allocatable :: tcc(:,:,:,:)      !< BC type and eventual index on chimera values [1:2,1-gc:ni,1-gc:nj,1-gc:nk].
   real(R4P),    allocatable :: chimera(:)        !< Chimera values (donors number, bijk-weight for each donor) [1:nchimera].
   integer(I4P)              :: ab=0              !< Absolute block index.
   integer(I4P)              :: group=0           !< Index of gruop.
   integer(I4P)              :: body=0            !< Index of body.
   integer(I4P)              :: proc=0            !< Processor assigned to.
   logical                   :: is_loaded=.false. !< Flag for checking if the block is loaded.
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
   self%ab    = 0
   self%group = 0
   self%body  = 0
   self%proc  = 0
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

   subroutine load_dimensions(self, file_unit, ab)
   !< Load block dimensions from file.
   !<
   !< @note The file must be already open and the current record-index must be at the proper block dimensions record.
   !< If ab index is passed it is supposed that this the first time the dimensions is loaded from grd file, thus the block
   !< destroyed and allocated ex novo.
   class(block_object), intent(inout)        :: self      !< Block data.
   integer(I4P),        intent(in)           :: file_unit !< Logical unit of grd file.
   integer(I4P),        intent(in), optional :: ab        !< Absolute block index.

   if (present(ab)) call self%destroy
   read(file_unit, end=10, err=10) self%Ni,self%Nj,self%Nk,self%gc
   10 continue
   if (present(ab)) then
      call self%alloc
      self%ab = ab
      self%w = self%weight()
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
   print*, 'cazzo ni,nj,nk,gc',ni,nj,nk,gc
      read(file_unit)(((icc(i,j,k),i=1-gc,Ni+gc),j=1-gc,Nj+gc),k=1-gc,Nk+gc)
   print*, 'cazzo icc',count(icc>0.and.icc<BC_NATURAL_RCC_RESERVED_DATA)
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
   logical, allocatable               :: on_rcc(:)

   do p=1, 25
      print*, 'cazzo rcc aaaa ',p,rcc(p)
   enddo

   allocate(on_rcc(1:size(rcc,dim=1)))
   on_rcc = .false.
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
      ! set blocks%tcc
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
            on_rcc(p  ) = .true.
            on_rcc(p+1) = .true.
            do n=1, ndonors
               self%chimera(nchimera+1+5*(n-1)) = donors(1,n) ! b
               self%chimera(nchimera+2+5*(n-1)) = donors(2,n) ! i
               self%chimera(nchimera+3+5*(n-1)) = donors(3,n) ! j
               self%chimera(nchimera+4+5*(n-1)) = donors(4,n) ! k
               self%chimera(nchimera+5+5*(n-1)) = donors(5,n) ! weight
               on_rcc(p+1+5*(n-1)+1) = .true.
               on_rcc(p+1+5*(n-1)+2) = .true.
               on_rcc(p+1+5*(n-1)+3) = .true.
               on_rcc(p+1+5*(n-1)+4) = .true.
               on_rcc(p+1+5*(n-1)+5) = .true.
            enddo
            nchimera = nchimera + ndonors*5 ! b,i,j,k,weight for each donor
         elseif (p>0.and.p<=BC_NATURAL_RCC_RESERVED_DATA) then ! natural BC
            on_rcc(p) = .true.
         endif
      enddo
      enddo
      enddo
   endif
   endassociate
   print*, 'cazzo block ',self%ab,'on_rcc elements      ',count(on_rcc)
   print*, 'cazzo block ',self%ab,'on_rcc elements miss ',count(.not.on_rcc)
   endsubroutine parse_rcc

   pure subroutine sanitize_chimera(self, sb)
   !< Sanitize chimera data after a split.
   class(block_object), intent(inout) :: self        !< Block data.
   type(block_object),  intent(in)    :: sb(2)       !< Split blocks.
   integer(I4P)                       :: i,j,k,p,n   !< Counter.
   integer(I4P)                       :: db,di,dj,dk !< Donor indexes.
   integer(I4P)                       :: ndonors     !< Chimera donors number.
   integer(I4P)                       :: o           !< Offset of rcc.

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

   pure subroutine split(self, mgl, is_split_done, sb)
   !< Split block. Split current block (if possible), along a the largest direction, in half: the first
   !< Block substitute current block, the other is added to blocks list.
   class(block_object), intent(in)  :: self          !< Block data.
   integer(I4P),        intent(in)  :: mgl           !< Number of levels of multi-grid to be preserved.
   logical,             intent(out) :: is_split_done !< Sentinel to check is split has been done.
   type(block_object),  intent(out) :: sb(2)         !< Split blocks.
   integer(I4P)                     :: nadj          !< Number of new adjacent-BC cells.
   integer(I4P)                     :: delta(3)      !< Directions deltas.

   is_split_done = .false.
   call sb(1)%destroy
   call sb(2)%destroy
   call find_split(Ni=self%Ni,Nj=self%Nj,Nk=self%Nk,gc=self%gc,mgl=mgl,nadj=nadj,delta=delta,sb=sb)
   if (sb(1)%split_dir>0) then
      ! alloc split blocks
      call sb(1)%alloc
      call sb(2)%alloc
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
      ! split data
      call split_nodes(nodes=self%nodes,delta=delta,gc=self%gc,sb=sb)
      call split_tcc(tcc=self%tcc,delta=delta,gc=self%gc,sb=sb)
      call split_chimera(chimera=self%chimera,delta=delta,nadj=nadj,sb_n=1,ab_ob=sb(2)%ab,bs=sb(1))
      call split_chimera(chimera=self%chimera,delta=delta,nadj=nadj,sb_n=2,ab_ob=sb(1)%ab,bs=sb(2))
      is_split_done = .true.
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
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
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
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            nchimera = nchimera + 1
            p = bs%tcc(2,i,j,k)               ! point to parent chimera array
            bs%tcc(2,i,j,k) = nchimera        ! point to new split block chimera array
            bs%chimera(nchimera) = chimera(p) ! assigno donors number
            do n=1, nint(chimera(p),I4P) ! b,i,j,k,weight for each donor
               o = p + 5*(n-1)
               if (nint(chimera(o+1),R4P)>bs%ab+(1-sb_n)) then
                  ! reference to a block subsequent to the split one, ab index must be shifted
                  bs%chimera(nchimera+1) = chimera(o+1) + 1
               else
                  bs%chimera(nchimera+1) = chimera(o+1)
               endif
               bs%chimera(nchimera+2) = chimera(o+2)
               bs%chimera(nchimera+3) = chimera(o+3)
               bs%chimera(nchimera+4) = chimera(o+4)
               bs%chimera(nchimera+5) = chimera(o+5)
            enddo
            nchimera = nchimera + 5
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

      pure subroutine split_nodes(nodes,delta,gc,sb)
      !< Split nodes.
      integer(I4P),       intent(in)    :: delta(3)                    !< Deltas.
      integer(I4P),       intent(in)    :: gc                          !< Number of ghost cells.
      real(R8P),          intent(in)    :: nodes(1:,0-gc:,0-gc:,0-gc:) !< Nodes coordinates.
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
                               lhs%ab          = rhs%ab
                               lhs%group       = rhs%group
                               lhs%body        = rhs%body
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
      do n = 1, min(ndonors, size(donors,dim=2))
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

   subroutine load_file_grd(file_name,blocks,blocks_number)
   !< Load file grd.
   character(*),       intent(in)                 :: file_name     !< File name
   type(block_object), intent(inout), allocatable :: blocks(:)     !< Blocks data.
   integer(I4P),       intent(out)                :: blocks_number !< Blocks number.
   integer(I4P)                                   :: file_unit     !< File unit.
   integer(I4P)                                   :: b             !< Counter.

   if (allocated(blocks)) deallocate(blocks)
   open(newunit=file_unit, file=trim(adjustl(file_name)), form='unformatted', action='read')
   read(file_unit, end=10, err=10) blocks_number
   allocate(blocks(1:blocks_number))
   do b=1, blocks_number
      call blocks(b)%load_dimensions(ab=b, file_unit=file_unit)
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
         write(file_unit,*) b, blocks(b)%group, blocks(b)%body, blocks(b)%proc, blocks(b)%parents
      else
         write(file_unit,*) b, blocks(b)%group, blocks(b)%body, blocks(b)%proc, ' original block, no splits'
      endif
   enddo
   close(file_unit)
   endsubroutine save_proc_input

   subroutine update_blocks(blocks, sb, blocks_number)
   !< Update blocks data after a block split.
   type(block_object), intent(inout), allocatable :: blocks(:)     !< Blocks data.
   type(block_object), intent(in)                 :: sb(1:)        !< Split blocks.
   integer(I4P),       intent(out)                :: blocks_number !< Blocks number.
   type(block_object), allocatable                :: blocks_(:)    !< New blocks data.
   integer(I4P)                                   :: b             !< Counter.

   blocks_number = size(blocks,dim=1)
   ! sanitize old chimera data
   do b=1, blocks_number
      if (b==sb(1)%ab) cycle ! split block does not need to be santized, it is replaced by sb
      call blocks(b)%sanitize_chimera(sb=sb)
   enddo
   ! ab shift
   do b=sb(2)%ab, blocks_number
      blocks(b)%ab = blocks(b)%ab + 1
   enddo
   ! create new blocks data
   print*,'blocks number ', blocks_number
   allocate(blocks_(1:blocks_number+1))
   do b=1, sb(1)%ab - 1
      blocks_(b) = blocks(b)
      call blocks(b)%destroy  ! free memory immediately
   enddo
   blocks_(sb(1)%ab) = sb(1)
   blocks_(sb(2)%ab) = sb(2)
   do b=sb(1)%ab + 1, blocks_number
      blocks_(b+1) = blocks(b)
      call blocks(b)%destroy  ! free memory immediately
   enddo
   call move_alloc(from=blocks_, to=blocks)
   if (allocated(blocks_)) deallocate(blocks_)
   blocks_number = blocks_number + 1
   endsubroutine update_blocks
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

program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object
use oe_process_object

implicit none

character(len=99)                 :: file_name_grd        !< Grid file name.
character(len=99)                 :: file_name_icc        !< Icc file name.
character(len=99)                 :: file_name_proc_input !< Name of proc.input file.
logical                           :: save_block_tecplot   !< Save blocks also in tecplot (ASCII) format.
logical                           :: save_imploded        !< Save imploded blocks after explosion.
logical                           :: save_exploded        !< Save exploded blocks.
character(len=99)                 :: exploded_basename    !< Exploded files basebame.
integer(I4P)                      :: mgl                  !< Multigrid level to be preserved.
integer(I4P)                      :: blocks_number        !< Number of blocks contained into the files.
type(block_object), allocatable   :: blocks(:)            !< Blocks data.
integer(I4P)                      :: total_blocks_weight  !< Total blocks weight.
real(R4P), allocatable            :: rcc(:)               !< rcc unstructured array.
logical                           :: is_split_done        !< Sentinel to check is split has been done.
type(block_object)                :: sb(2)                !< Split blocks.
integer(I4P), allocatable         :: blocks_list(:)       !< Blocks (unassigned) list (decreasing-workload) ordered.
integer(I4P)                      :: procs_number         !< Number of processes for load balancing.
type(process_object), allocatable :: processes(:)         !< Processes data.
integer(I4P)                      :: ideal_proc_workload  !< Ideal process weight for load balancing.
integer(I4P)                      :: max_unbalance        !< Maximum processes unbalancing in percent.
integer(I4P)                      :: i,b,bb,p             !< Counter.

! interface for auxiliary procedures
interface str
  !< Convert number (real and integer) to string (number to string type casting).
  procedure str_I4P, str_a_I4P
endinterface

interface strz
  !< Convert integer, to string, prefixing with the right number of zeros (integer to string type casting with zero padding).
  procedure strz_I4P
endinterface

call parse_command_line(fgrd=file_name_grd,ficc=file_name_icc,fpci=file_name_proc_input,&
                        stec=save_block_tecplot,simp=save_imploded,sexp=save_exploded,  &
                        ebn=exploded_basename,np=procs_number,mu=max_unbalance,mgl=mgl)

if ((.not.is_file_found(file_name_grd)).or.(.not.is_file_found(file_name_icc))) then
   write(stderr, '(A)')'error: file "'//trim(adjustl(file_name_grd))//'" or '//&
                                   '"'//trim(adjustl(file_name_icc))//'" not found!'
   stop
endif

allocate(processes(0:procs_number-1))
call processes%initialize

print '(A)', 'load grd file '//trim(adjustl(file_name_grd))
call load_file_grd(file_name=file_name_grd,blocks=blocks,blocks_number=blocks_number)

print '(A)', 'load icc file '//trim(adjustl(file_name_icc))
call load_file_icc(file_name=file_name_icc,blocks=blocks,blocks_number=blocks_number,rcc=rcc)
print '(A)', 'finish load input files'

! parse global rcc
print '(A)', 'parse global rcc and create block-local-rcc'
print '(A)', 'total chimera elements '//trim(str(size(rcc,dim=1)))
do b=1, blocks_number
   call blocks(b)%parse_rcc(rcc=rcc)
   print '(A)', 'block '//trim(str(b,.true.))//' BC chimera cells number: '//trim(str(size(blocks(b)%chimera,dim=1)))
enddo
print '(A)', 'finish parse global rcc'

! load balancing
print '(A)', 'load balancing stats'
total_blocks_weight = 0
do b=1, blocks_number
   print '(A)', '    block "'//trim(strz(b,9))//'" weight: '//trim(str(blocks(b)%w,.true.))
   total_blocks_weight = total_blocks_weight + blocks(b)%w
enddo
ideal_proc_workload = total_blocks_weight / procs_number
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
         call update_blocks(blocks=blocks, sb=sb, blocks_number=blocks_number)
         call create_blocks_list(blocks=blocks, blocks_list=blocks_list)
         call processes%initialize
      else
         print '(A)', 'block "'//trim(strz(blocks(b)%ab,9))//'" split failed, assigned anyway to process '//trim(strz(p,6))
         blocks(b)%proc = p
         call processes(p)%assign_block(ab=b, wb=blocks(b)%w, ideal_workload=ideal_proc_workload)
         call popout_blocks_list(blocks_list=blocks_list)
      endif
   endif
enddo assign_blocks_loop

print '(A)', 'processes workload'
do p=0, procs_number-1
print '(A)', '  proc '//trim(strz(p,6))//&
             ' unbalancing '//trim(str(processes(p)%unbalance))//&
             '% assigned blocks '//trim(str(processes(p)%blocks(2:),.true.))
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
   call implode_blocks(blocks=blocks, rcc=rcc)
   print '(A)', 'total chimera elements after implosion'//trim(str(size(rcc,dim=1)))
   do b=1, blocks_number
      print '(A)', 'block '//trim(str(b,.true.))//' BC chimera cells number: '//trim(str(size(blocks(b)%chimera,dim=1)))
   enddo
   print '(A)', 'save imploded blocks in legacy overset format (split and load-balanced)'
   ! call save_file_grd(file_name='split-balanced-'//trim(adjustl(file_name_grd)), blocks=blocks)
   ! call save_file_icc(file_name='split-balanced-'//trim(adjustl(file_name_icc)), blocks=blocks, rcc=rcc)
endif
contains
   subroutine parse_command_line(fgrd,ficc,fpci,stec,simp,sexp,ebn,np,mu,mgl)
   !< Parse command line inputs.
   character(*), intent(out) :: fgrd      !< Grid file name.
   character(*), intent(out) :: ficc      !< Icc file name.
   character(*), intent(out) :: fpci      !< Name of proc.input file.
   logical,      intent(out) :: stec      !< Save blocks also in tecplot (ASCII) format.
   logical,      intent(out) :: simp      !< Save imploded blocks after explosion.
   logical,      intent(out) :: sexp      !< Save exploded blocks.
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
      case('-exploded-basename')
         a = a + 1
         call get_command_argument(a, ca_buffer)
         ebn = trim(adjustl(ca_buffer))
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
   write(*, '(A)')'   -h, --help                       => print this help message'
   write(*, '(A)')'examples:'
   write(*, '(A)')'   overset-exploded -np 32'
   write(*, '(A)')'   overset-exploded -grd cc.02.grd -icc cc.02 -np 16'
   write(*, '(A)')'   overset-exploded -np 16 -max-unbalance 4'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16 -save-imploded'
   write(*, '(A)')'   overset-exploded -np 16 -proc-input proc.input-pes16 -save-exploded'
   write(*, '(A)')'   overset-exploded -grd cc.03.grd -icc cc.03 -np 2 -tec -max-unbalance 3'
   endsubroutine print_help

   ! files procedures
   function is_file_found(file_name) result(is_found)
   !< Inquire is the file path is valid and the file is found.
   character(*), intent(in) :: file_name !< File name.
   logical                  :: is_found  !< Inquiring result.

   inquire(file=trim(adjustl(file_name)), exist=is_found)
   endfunction is_file_found

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
endprogram overset_exploded
