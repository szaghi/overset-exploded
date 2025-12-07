module oe_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32

implicit none

private
public :: block_object
public :: bc_int_type, bc_string

! BC parameters
! natural BC
integer(kind=I4P), parameter, public :: BC_NATURAL_WALL                     = -1   !< Viscous wall.
integer(kind=I4P), parameter, public :: BC_NATURAL_SIMMETRY                 = -2   !< Simmetry.
integer(kind=I4P), parameter, public :: BC_NATURAL_INFLOW                   = -3   !< Inflow.
integer(kind=I4P), parameter, public :: BC_NATURAL_INOUTFLOW                = -4   !< In/outflow.
integer(kind=I4P), parameter, public :: BC_NATURAL_ASSIGNED_INFLOW          = -5   !< Assigned inflow.
integer(kind=I4P), parameter, public :: BC_NATURAL_ASSIGNED_PRESSURE        = -6   !< Assigned pressure.
integer(kind=I4P), parameter, public :: BC_NATURAL_ASSIGNED_NORMAL_VELOCITY = -7   !< Assigned normal velocity.
integer(kind=I4P), parameter, public :: BC_NATURAL_ASSIGNED_RIEMANN         = -8   !< Assigned Riemann invariant.
integer(kind=I4P), parameter, public :: BC_NATURAL_EXTRAPOLATED             = -9   !< Extrapolated.
integer(kind=I4P), parameter, public :: BC_NATURAL_MOVING_WALL              = -10  !< Moving wall.
integer(kind=I4P), parameter, public :: BC_NATURAL_INACTIVE_WALL            = -11  !< Inactive wall.
integer(kind=I4P), parameter, public :: BC_NATURAL_EXTRAPOLATED_ALT         = -19  !< Extrapolated (alternative).
! non BC, active cell
integer(kind=I4P), parameter, public :: BC_ACTIVE_CELL                      = 0    !< Non BC, active cell.
! chimera BC, face face-center
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF                  = 20   !< Chimera face.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_I0               = 21   !< Chimera face, centered at i0 face-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_IN               = 22   !< Chimera face, centered at in face-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_J0               = 23   !< Chimera face, centered at j0 face-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_JN               = 24   !< Chimera face, centered at jn face-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_K0               = 25   !< Chimera face, centered at k0 face-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XF_KN               = 26   !< Chimera face, centered at kn face-center.
! chimera BC, cell
integer(kind=I4P), parameter, public :: BC_CHIMERA_CELL                     = 27   !< Chimera cell inside domain.
integer(kind=I4P), parameter, public :: BC_CHIMERA_CELL_INT_WALL            = 28   !< Chimera cell internal wall.
! chimera BC, face cell-center
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC                  = 40   !< Chimera face.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_I0               = 41   !< Chimera face, centered at i0 cell-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_IN               = 42   !< Chimera face, centered at in cell-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_J0               = 43   !< Chimera face, centered at j0 cell-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_JN               = 44   !< Chimera face, centered at jn cell-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_K0               = 45   !< Chimera face, centered at k0 cell-center.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_XC_KN               = 46   !< Chimera face, centered at kn cell-center.
! chimera BC, adjacent
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ                 = 60   !< Adjacent.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_I0              = 61   !< Adjacent along face i0.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_IN              = 62   !< Adjacent along face in.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_J0              = 63   !< Adjacent along face j0.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_JN              = 64   !< Adjacent along face jn.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_K0              = 65   !< Adjacent along face k0.
integer(kind=I4P), parameter, public :: BC_CHIMERA_FACE_ADJ_KN              = 66   !< Adjacent along face kn.
! edge BC
integer(kind=I4P), parameter, public :: BC_EDGE                             = 80   !< Edge.

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
endtype block_object

contains
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

   pure subroutine parse_rcc(self, rcc)
   !< Parse global rcc and store in local tcc/chimera arrays.
   class(block_object), intent(inout) :: self      !< Block data.
   real(R4P),           intent(in)    :: rcc(1:)   !< rcc unstructured array.
   integer(I4P)                       :: nchimera  !< Number of chimera data.
   integer(I4P)                       :: ndonors   !< Number of donors.
   integer(I4P)                       :: i,j,k,n,p !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,icc=>self%icc)
   if (allocated(self%chimera)) deallocate(self%chimera)
   self%tcc(1,:,:,:) = BC_ACTIVE_CELL ! set all cell type to active cell
   nchimera = 0
   do k=1-gc, Nk+gc
   do j=1-gc, Nj+gc
   do i=1-gc, Ni+gc
      p = icc(i,j,k)
      if (p>0) then
         self%tcc(1,i,j,k) = int(rcc(p),I4P)
         select case(int(rcc(p),I4P))
         case(BC_NATURAL_EXTRAPOLATED_ALT:BC_NATURAL_WALL)
            ! no additional data is necessary
         case(BC_ACTIVE_CELL)
            ! no additional data is necessary
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            ! for chimera BC it is necessary to store other data
            nchimera = nchimera + 1
            self%tcc(2,i,j,k) = nchimera
            ndonors = nint(rcc(p+1),I4P)    ! donors number
            nchimera = nchimera + ndonors*5 ! b,i,j,k,weight for each donor
         case(BC_EDGE)
            ! no additional data is necessary
         case default
            ! print *, 'error: unknown tcc "',self%tcc(1,i,j,k),'", ab,i,j,k=',self%ab,i,j,k
         endselect
      endif
   enddo
   enddo
   enddo
   if (nchimera>0) then
      allocate(self%chimera(1:nchimera))
      nchimera = 0
      do k=1-gc, Nk+gc
      do j=1-gc, Nj+gc
      do i=1-gc, Ni+gc
         p = icc(i,j,k)
         if (p>0) then
            select case(int(rcc(p),I4P))
            case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
               nchimera = nchimera + 1
               self%chimera(nchimera) = rcc(p+1) ! donors number
               do n=1, nint(rcc(p+1),I4P)
                  self%chimera(nchimera+1+5*(n-1)) = rcc(p+2+5*(n-1)) ! b
                  self%chimera(nchimera+2+5*(n-1)) = rcc(p+3+5*(n-1)) ! i
                  self%chimera(nchimera+3+5*(n-1)) = rcc(p+4+5*(n-1)) ! j
                  self%chimera(nchimera+4+5*(n-1)) = rcc(p+5+5*(n-1)) ! k
                  self%chimera(nchimera+5+5*(n-1)) = rcc(p+6+5*(n-1)) ! weight
               enddo
               nchimera = nchimera + nint(rcc(p+1),I4P)*5 ! b,i,j,k,weight for each donor
            endselect
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
               ! self block has a chimera reference with a splitted block that must be sanitized, if fall in sb(2)
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

   subroutine save_block_file(self, b, rcc, tec)
   !< Save block data into its own file.
   class(block_object), intent(in)           :: self      !< Block data.
   integer(I4P),        intent(in)           :: b         !< Current block number, global numeration.
   real(R4P),           intent(in)           :: rcc(1:)   !< rcc unstructured array.
   logical,             intent(in), optional :: tec       !< Save (also) in tecplot (ASCII) format (for debugging).
   logical                                   :: tec_      !< Save (also) in tecplot (ASCII) format (for debugging), local var.
   character(len=6)                          :: bstr      !< Block number stringified.
   integer(I4P)                              :: file_unit !< Block unit file.
   integer(I4P)                              :: i,j,k,p,n !< Counter.
   integer(I4P)                              :: o         !< Offset of rcc.

   tec_ = .false. ; if (present(tec)) tec_ = tec
   write(bstr, '(I6.5)') b
   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes,icc=>self%icc,tcc=>self%tcc,chimera=>self%chimera)
   ! file grid
   open(newunit=file_unit, file='block-'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write', status='replace')
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
      case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
         p = tcc(2,i,j,k)
         write(file_unit) i,j,k,bc_string(tcc(1,i,j,k))
         write(file_unit) nint(chimera(p),I4P) ! donors number
         do n=1, nint(chimera(p),I4P) ! b,i,j,k,weight for each donor
            o = p + 5*(n-1)
            write(file_unit) nint(chimera(o+1)),nint(chimera(o+2)),nint(chimera(o+3)),nint(chimera(o+4)),chimera(o+5)
         enddo
      case(BC_EDGE)
         write(file_unit) i,j,k,bc_string(tcc(1,i,j,k))
      case default
         print *, 'error: unknown tcc "',tcc(1,i,j,k),'", b,i,j,k=',b,i,j,k
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
                  ! reference to a block subsequent to the splitted one, ab index must be shifted
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
         sb(2)%tcc(:,i,j,k) = tcc(:,i+(1+sb(1)%Ni)*delta(1),j+(1+sb(1)%Nj)*delta(2),k+(1+sb(1)%Nk)*delta(3))
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
   case('BC_EDGE'                            ,'bc_edge'                            );bc_int_type=BC_EDGE
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
   case(BC_EDGE                            ) ; bc_string = 'BC_EDGE'
   endselect
   endfunction bc_string
endmodule oe_block_object

module oe_process_object
!< Overset-Exploded, definition of class process_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32

implicit none

private
public :: process_object

type :: process_object
   !< Class process object.
   integer(I4P)              :: id = 0      !< Process ID.
   integer(I4P), allocatable :: blocks(:)   !< List of assigned blocks.
   integer(I4P)              :: w = 0       !< Process workload.
   integer(I4P)              :: unbalance=0 !< Unbalance respect ideal workload.
endtype process_object
endmodule oe_process_object

program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object
use oe_process_object

implicit none

character(len=99)                 :: ca_buffer            !< Command argument buffer.
character(len=99)                 :: file_name_grd        !< Grid file name.
character(len=99)                 :: file_name_icc        !< Icc file name.
integer(I4P)                      :: file_unit_grd        !< Grid unit file.
integer(I4P)                      :: file_unit_icc        !< Icc unit file.
integer(I4P)                      :: na                   !< Number of command line arguments.
logical                           :: save_block_tecplot   !< Save blocks also in tecplot (ASCII) format.
integer(I4P)                      :: blocks_number=0      !< Number of blocks contained into the files.
type(block_object), allocatable   :: blocks(:)            !< Blocks data.
integer(I4P)                      :: unstruct_dimension   !< Dimension of unstructured array of rcc.
integer(I4P)                      :: total_blocks_weight  !< Total blocks weight.
integer(I4P)                      :: ideal_block_weight   !< Ideal block weight for load balancing.
real(R4P), allocatable            :: rcc(:)               !< rcc unstructured array.
logical                           :: is_split_done        !< Sentinel to check is split has been done.
type(block_object)                :: sb(2)                !< Split blocks.
integer(I4P), allocatable         :: blocks_olist(:)      !< Blocks ordered (decreasing-workload) list.
integer(I4P), allocatable         :: blocks_ulist(:)      !< Blocks unassigned list.
integer(I4P)                      :: procs_number=1       !< Number of processes for load balancing.
type(process_object), allocatable :: processes(:)         !< Processes data.
integer(I4P), allocatable         :: processes_olist(:)   !< Processes ordered (increasing-workload) list.
integer(I4P)                      :: ideal_proc_weight    !< Ideal process weight for load balancing.
integer(I4P)                      :: proc_unbalance       !< Current process unbalancing in percent.
integer(I4P)                      :: max_unbalance=4      !< Maximum unbalancing in percent.
integer(I4P)                      :: i,b,bb,bbb,p         !< Counter.

na = command_argument_count()
if (na<3) then
   write(stderr, "(A)")'error: you must pass file_name_grd, file_name_icc and procs number as command line arguments'
   write(stderr, "(A)")'example: overset-exploded cc.01.grd cc.01 16'
   write(stderr, "(A)")'if you want to save blocks in tecplot format pass also "save_tecplot"'
   write(stderr, "(A)")'example: overset-exploded cc.01.grd cc.01 16 save_tecplot'
   stop
else
   save_block_tecplot = .false.
   call get_command_argument(1, file_name_grd)
   call get_command_argument(2, file_name_icc)
   call get_command_argument(3, ca_buffer)
   read(ca_buffer,*) procs_number
   if (na==4) then
      call get_command_argument(4, ca_buffer)
      save_block_tecplot = (trim(adjustl(ca_buffer))=='save_tecplot')
   endif
   allocate(processes(0:procs_number-1))
   ! assign processes ID
   do p=0, procs_number-1
      processes(p)%id = p
      processes(p)%blocks = [0_I4P]
   enddo
   allocate(processes_olist(0:procs_number-1))
   processes_olist = [(p,b=0,procs_number-1)]
endif

if (is_file_found(file_name_grd).and.is_file_found(file_name_icc)) then
   ! grd file
   print *, 'load grd file ',trim(adjustl(file_name_grd))
   open(newunit=file_unit_grd, file=trim(adjustl(file_name_grd)), form='unformatted', action='read')
   read(file_unit_grd, end=10, err=10) blocks_number
   allocate(blocks(1:blocks_number))
   do b=1, blocks_number
      call blocks(b)%load_dimensions(ab=b, file_unit=file_unit_grd)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_nodes(file_unit=file_unit_grd)
   enddo
   10 close(file_unit_grd)

   ! icc file
   print *, 'load icc file ',trim(adjustl(file_name_icc))
   open(newunit=file_unit_icc, file=trim(adjustl(file_name_icc)), form='unformatted', action='read')
   read(file_unit_icc) b
   if (b/=blocks_number) then
      write(stderr, "(A)")'error: grd and icc have different number of blocks'
      stop
   endif
   do b=1, blocks_number
      call blocks(b)%load_dimensions(file_unit=file_unit_icc)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_icc(file_unit=file_unit_icc)
   enddo
   read(file_unit_icc) unstruct_dimension
   allocate(rcc(1:unstruct_dimension))
   read(file_unit_icc) (rcc(i),i=1,unstruct_dimension)
   close(file_unit_icc)
   print *, 'finish load input files'

   ! parse global rcc
   print *, 'parse global rcc'
   do b=1, blocks_number
      print *, '  block ',b
      call blocks(b)%parse_rcc(rcc=rcc)
   enddo
   print *, 'finish parse global rcc'
   ! stop

   ! load balancing
   print *, 'load balancing stats'
   total_blocks_weight = 0
   do b=1, blocks_number
      print *, '    block ',b, ' weight ', blocks(b)%w
      total_blocks_weight = total_blocks_weight + blocks(b)%w
   enddo
   ideal_block_weight = total_blocks_weight / blocks_number
   ideal_proc_weight = total_blocks_weight / procs_number
   print *, '  ideal work load for np "',procs_number,'" processes: ', ideal_proc_weight

   print *, 'order blocks in decreasing-workload-order'
   blocks_olist = [(b,b=1,blocks_number)] ; call blocks_quick_sort(bs=blocks,bl=blocks_olist,first=1,last=blocks_number)
   do b=1, blocks_number
      print *, '    block ',blocks(blocks_olist(b))%ab, ' weight ', blocks(blocks_olist(b))%w, &
               ' Ni,Nj,Nk ', blocks(blocks_olist(b))%Ni, blocks(blocks_olist(b))%Nj, blocks(blocks_olist(b))%Nk
   enddo

   ! assign blocks to processes
   blocks_ulist = blocks_olist
   do while(allocated(blocks_ulist))
      p = minloc(processes(0:)%w,dim=1)-1 ! process with minimum workload
      proc_unbalance = unbalance(ideal_proc_weight,processes(p)%w+blocks(blocks_ulist(1))%w)
      if (processes(p)%w+blocks(blocks_ulist(1))%w<=ideal_proc_weight) then
         ! processes(p)%blocks    = [processes(p)%blocks,blocks_ulist(1)]
         ! processes(p)%w         = processes(p)%w + blocks(blocks_ulist(1))%w
         ! processes(p)%unbalance = unbalance(ideal_proc_weight,processes(p)%w)
         ! if (size(blocks_ulist,dim=1)>1) then
         !    blocks_ulist = blocks_ulist(2:)
         ! else
         !    deallocate(blocks_ulist)
         ! endif
      else
         print *, 'block ',blocks(blocks_ulist(1))%ab,' must be split to be insert into process ',&
               p,' process workload ',processes(p)%w,' new block workload ',blocks(blocks_ulist(1))%w,' unbalancing ',proc_unbalance
      endif
         processes(p)%blocks    = [processes(p)%blocks,blocks_ulist(1)]
         processes(p)%w         = processes(p)%w + blocks(blocks_ulist(1))%w
         processes(p)%unbalance = unbalance(ideal_proc_weight,processes(p)%w)
         if (size(blocks_ulist,dim=1)>1) then
            blocks_ulist = blocks_ulist(2:)
         else
            deallocate(blocks_ulist)
         endif
   enddo
   print *, 'processes workload'
   do p=0, procs_number-1
      print *, '  process ',p,' workload ',processes(p)%w,' unbalancing ',processes(p)%unbalance,'%',&
               ' assigned blocks',processes(p)%blocks(2:)
   enddo

   stop

   b = 1
   split_loop : do
      if (blocks(b)%w>ideal_block_weight) then
         print *, '    block ',b,' nijk ',blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk,' weight ',blocks(b)%w,'>',ideal_block_weight
         print *, '    try to split block ',b
         call blocks(b)%split(mgl=2, is_split_done=is_split_done, sb=sb)
         if (is_split_done) then
            print *, '    block ',b,' splitted'
            print *, '      first split block  (ab,ni,nj,nk) ',sb(1)%ab,sb(1)%Ni,sb(1)%Nj,sb(1)%Nk
            print *, '      second split block (ab,ni,nj,nk) ',sb(2)%ab,sb(2)%Ni,sb(2)%Nj,sb(2)%Nk
            print *, '      first block parents list         ',sb(1)%parents
            print *, '      second block parents list        ',sb(2)%parents
            print *, '      update blocks list'
            ! sanitize old chimera data
            do bb=1, blocks_number
               if (bb==b) cycle ! splitted block does not need to be santized, it is replaced by sb
               call blocks(bb)%sanitize_chimera(sb=sb)
            enddo
            ! ab shift
            do bb=sb(2)%ab, blocks_number
               blocks(bb)%ab = blocks(bb)%ab + 1
            enddo
            ! update blocks list
            if     (sb(1)%ab==1) then ! first block in list has been splitted
               blocks = [sb,blocks(2:blocks_number)]
            elseif (sb(1)%ab==blocks_number) then ! last  block in list has been splitted
               blocks = [blocks(1:blocks_number-1),sb]
            else
               blocks = [blocks(1:sb(1)%ab-1),sb,blocks(sb(1)%ab+1:blocks_number)]
            endif
            blocks_number = blocks_number + 1
            if (sb(1)%w>ideal_block_weight) cycle split_loop
         endif
      endif
      if (b==blocks_number) exit split_loop
      b = b + 1
   enddo split_loop
   ! stop

   print *, 'save exploded blocks'
   do b=1, blocks_number
      call blocks(b)%save_block_file(b=b, rcc=rcc, tec=save_block_tecplot)
   enddo
else
   write(stderr, "(A)")'error: file "'//trim(adjustl(file_name_grd))//'" or '//&
                                   '"'//trim(adjustl(file_name_icc))//'" not found!'
endif

contains
   recursive subroutine blocks_quick_sort(bs,bl,first,last)
   !< Order blocks list in decreasing-workload-order by quick sort algorithm.
   type(block_object), intent(in)    :: bs(1:) !< Blocks data.
   integer(I4P),       intent(inout) :: bl(1:) !< Blocks ordered list.
   integer(I4P),       intent(in)    :: first  !< First sort index.
   integer(I4P),       intent(in)    :: last   !< Last sort index.
   integer(I4P)                      :: i,j    !< Counter.
   integer(I4P)                      :: pivot  !< Pivot.
   integer(I4P)                      :: temp   !< Temporary buffer.

   i = first
   j = last
   pivot = bl((first + last)/2)

   do
      do while (bs(bl(i))%w > bs(pivot)%w)
          i = i + 1
      enddo
      do while (bs(bl(j))%w < bs(pivot)%w)
          j = j - 1
      enddo
      if (i <= j) then
          temp = bl(i)
          bl(i) = bl(j)
          bl(j) = temp
          i = i + 1
          j = j - 1
      endif
      if (i > j) exit
   enddo
   if (first < j) call blocks_quick_sort(bs=bs,bl=bl,first=first,last=j   )
   if (i < last ) call blocks_quick_sort(bs=bs,bl=bl,first=i    ,last=last)
   endsubroutine blocks_quick_sort

   function is_file_found(file_name) result(is_found)
   !< Inquire is the file path is valid and the file is found.
   character(*), intent(in) :: file_name !< File name.
   logical                  :: is_found  !< Inquiring result.

   inquire(file=trim(adjustl(file_name)), exist=is_found)
   endfunction is_file_found

   pure function unbalance(ideal_wload,wload)
   !< Return the unbalance of workload.
   integer(I4P), intent(in) :: ideal_wload !< Ideal workload.
   integer(I4P), intent(in) :: wload       !< Workload.
   integer(I4P)             :: unbalance   !< Workload unbalancing.

   unbalance = nint(((ideal_wload*1._R8P-wload*1._R8P)/ideal_wload*1._R8P)*100._R8P)
   endfunction unbalance
endprogram overset_exploded
