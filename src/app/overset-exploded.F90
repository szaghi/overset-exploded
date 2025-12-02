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
   integer(I4P)              :: weight=0          !< Block weight.
   real(R8P),    allocatable :: nodes(:,:,:,:)    !< Nodes coordinates.
   integer(I4P), allocatable :: tcc(:,:,:,:)      !< BC type and eventual index on chimera values [1:2,1-gc:ni,1-gc:nj,1-gc:nk].
   real(R4P),    allocatable :: chimera(:)        !< Chimera values (type,donors number, bijk-weight for each donor) [1:nchimera].
   integer(I4P), allocatable :: adj(:,:)          !< New adjacent-BC indexes for split blocks [1:4,nadj].
   integer(I4P)              :: nadj=0            !< Number of new adjacent-BC cells.
   integer(I4P)              :: ab=0              !< Absolute block index.
   integer(I4P)              :: group=0           !< Index of gruop.
   integer(I4P)              :: body=0            !< Index of body.
   integer(I4P)              :: proc=0            !< Processor assigned to.
   ! old icc structure
   integer(I4P), allocatable :: icc(:,:,:)        !< Cell centered icc values.
   logical                   :: is_loaded=.false. !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: alloc           !< Allocate dynamic memory.
      procedure, pass(self) :: load_dimensions !< Load block dimensions from file.
      procedure, pass(self) :: load_icc        !< Load block icc from file.
      procedure, pass(self) :: load_nodes      !< Load block nodes from file.
      procedure, pass(self) :: parse_rcc       !< Parse global rcc and store in local tcc/chimera arrays.
      procedure, pass(self) :: save_block_file !< Save block data into its own file.
      procedure, pass(self) :: split           !< Split block.
endtype block_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_object), intent(inout) :: self !< Block data.

   self%Ni     = 0
   self%Nj     = 0
   self%Nk     = 0
   self%gc     = 2
   self%weight = 0
   if (allocated(self%nodes)) deallocate(self%nodes)
   if (allocated(self%icc)) deallocate(self%icc)
   if (allocated(self%tcc)) deallocate(self%tcc)
   if (allocated(self%chimera)) deallocate(self%chimera)
   if (allocated(self%adj)) deallocate(self%adj)
   self%nadj  = 0
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
   if (self%nadj>0) allocate(self%adj(1:4,self%nadj))
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
      self%weight = self%Ni*self%Nj*self%Nk
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
         case(BC_CHIMERA_FACE_XF:BC_CHIMERA_FACE_ADJ_KN)
            ! for chimera BC it is necessary to store other data
            nchimera = nchimera + 1
            self%tcc(2,i,j,k) = nchimera
            ndonors = nint(rcc(p+1),I4P)    ! donors number
            nchimera = nchimera + ndonors*5 ! b,i,j,k,weight for each donor
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
         print *, 'error: unknown icc "',int(rcc(p),I4P),'", b,i,j,k=',b,i,j,k
         stop
      endselect
   enddo
   enddo
   enddo

   if (tec_) then
      ! save block in tecplot format
      open(newunit=file_unit, file='block-'//trim(adjustl(bstr))//'.dat', action='write', status='replace')
      write(file_unit,*) 'TITLE = "Block '//trim(adjustl(bstr))//'"'
      write(file_unit,*) 'VARIABLES = "X", "Y", "Z", "tcc" "rcc"'
      write(file_unit,*) 'ZONE I=', ni+1+2*gc, ', J=', nj+1+2*gc, ', K=', nk+1+2*gc, &
                         ', DATAPACKING=BLOCK, VARLOCATION=([4,5]=CELLCENTERED)'
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(1,i,j,k) ; enddo ; enddo ; enddo
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(2,i,j,k) ; enddo ; enddo ; enddo
      do k=0-gc, Nk+gc ; do j=0-gc, Nj+gc ; do i=0-gc, Ni+gc ; write(file_unit,'(E23.15)') nodes(3,i,j,k) ; enddo ; enddo ; enddo

      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         write(file_unit,'(E23.15)') real(tcc(1,i,j,k),R8P)
      enddo ; enddo ; enddo
      do k=1-gc, Nk+gc ; do j=1-gc, Nj+gc ; do i=1-gc, Ni+gc
         p = icc(i,j,k)
         if (p>0) then
            write(file_unit,'(E23.15)') real(rcc(p),R8P)
         else
            write(file_unit,'(E23.15)') real(icc(i,j,k),R8P)
         endif
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
   integer(I4P)                     :: dir(1:3)      !< Ordered directions list.
   integer(I4P)                     :: N             !< Current direction cells number.
   integer(I4P)                     :: Ns(2)         !< Cells numbers in the two split blocks.
   logical                          :: dir_found     !< Sentil to check direction found.
   integer(I4P)                     :: delta(3)      !< Directions deltas.
   integer(I4P)                     :: i,j,k,b,d,ijk !< Counter.

   call sb(1)%destroy
   call sb(2)%destroy
   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc)
   ! search for direction with highest number of cells being compatible with MG level
   dir = hdirs(Ni=Ni,Nj=Nj,Nk=Nk)
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
      select case(dir(d))
      case(1) ! i direction
         sb(1)%Ni = Ns(1) ; sb(2)%Ni = Ns(2)
         sb%nadj = gc*(Nj+2*gc)*(Nk+2*gc)
         delta(1) = 1
      case(2) ! j direction
         sb(1)%Nj = Ns(1) ; sb(2)%Nj = Ns(2)
         sb%nadj = (Ni+2*gc)*gc*(Nk+2*gc)
         delta(2) = 1
      case(3) ! k direction
         sb(1)%Nk = Ns(1) ; sb(2)%Nk = Ns(2)
         sb%nadj = (Ni+2*gc)*(Nj+2*gc)*gc
         delta(3) = 1
      endselect
      ! alloc split blocks
      call sb(1)%alloc
      call sb(2)%alloc
      ! assign group, body, proc tags
      sb%group = self%group
      sb%body  = self%body
      sb%proc  = self%proc
      ! assign nodes
      do k=0-gc, sb(1)%Nk+gc
      do j=0-gc, sb(1)%Nj+gc
      do i=0-gc, sb(1)%Ni+gc
         sb(1)%nodes(:,i,j,k) = self%nodes(:,i,j,k)
      enddo
      enddo
      enddo
      do k=0-gc, sb(2)%Nk+gc
      do j=0-gc, sb(2)%Nj+gc
      do i=0-gc, sb(2)%Ni+gc
         sb(2)%nodes(:,i,j,k) = self%nodes(:,i+sb(1)%Ni*delta(1),j+sb(1)%Nj*delta(2),k+sb(1)%Nk*delta(3))
      enddo
      enddo
      enddo
      ! assign icc
      do k=1-gc, sb(1)%Nk+gc
      do j=1-gc, sb(1)%Nj+gc
      do i=1-gc, sb(1)%Ni+gc
         sb(1)%icc(i,j,k) = self%icc(i,j,k)
      enddo
      enddo
      enddo
      do k=1-gc, sb(2)%Nk+gc
      do j=1-gc, sb(2)%Nj+gc
      do i=1-gc, sb(2)%Ni+gc
         sb(2)%icc(i,j,k) = self%icc(i+(1+sb(1)%Ni)*delta(1),j+(1+sb(1)%Nj)*delta(2),k+(1+sb(1)%Nk)*delta(3))
      enddo
      enddo
      enddo
      ! assign absolute block index
      sb(1)%ab = self%ab     ! first  child substitute parent
      sb(2)%ab = self%ab + 1 ! second child substitute subsequent block of parent, shift is necessary in global blocks list
      ! create new adjacent icc along the split direction
      ijk = 0
      do k=1-gc, Nk+gc - (Nk+gc)*delta(3)
      do j=1-gc, Nj+gc - (Nj+gc)*delta(2)
      do i=1-gc, Ni+gc - (Ni+gc)*delta(1)
         ijk = ijk + 1
         sb(1)%icc(i+(sb(1)%Ni+gc)*delta(1),j+(sb(1)%Nj+gc)*delta(2),k+(sb(1)%Nk+gc)*delta(3)) = - ijk
         sb(2)%icc(i                       ,j                       ,k                       ) = - ijk
         sb(1)%adj(1:4,ijk) = [sb(2)%ab,i+gc*delta(1)      ,j+gc*delta(2)      ,k+gc*delta(3)      ] ! adjacent indexes, b,i,j,k
         sb(2)%adj(1:4,ijk) = [sb(1)%ab,i+sb(1)%Ni*delta(1),j+sb(1)%Nj*delta(2),k+sb(1)%Nk*delta(3)] ! adjacent indexes, b,i,j,k
      enddo
      enddo
      enddo
   endif
   endassociate
   contains
      pure function hdirs(Ni,Nj,Nk)
      !< Return the list of directions ordered from the one with highest-number of cells to the lowest one.
      integer(I4P), intent(in) :: Ni,Nj,Nk   !< Number of cells along each direction.
      integer(I4P)             :: hdirs(1:3) !< list of ordered directions.
      integer(I4P)             :: maxd,mind  !< Temporary variables.

      maxd = maxloc([Ni,Nj,Nk],dim=1) ; mind = minloc([Ni,Nj,Nk],dim=1) ; hdirs = [maxd,6-maxd-mind,mind]
      endfunction hdirs
   endsubroutine split

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

program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object

implicit none

character(len=99)               :: ca_buffer            !< Command argument buffer.
character(len=99)               :: file_name_grd        !< Grid file name.
character(len=99)               :: file_name_icc        !< Icc file name.
integer(I4P)                    :: file_unit_grd        !< Grid unit file.
integer(I4P)                    :: file_unit_icc        !< Icc unit file.
integer(I4P)                    :: blocks_number=0      !< Number of blocks contained into the files.
type(block_object), allocatable :: blocks(:)            !< Blocks contained into the files.
logical                         :: is_split_done        !< Sentinel to check is split has been done.
type(block_object)              :: sb(2)                !< Split blocks.
integer(I4P)                    :: unstruct_dimension   !< Dimension of unstructured array of rcc.
integer(I4P)                    :: procs_number         !< Number of processes for load balancing.
integer(I4P)                    :: total_blocks_weight  !< Total blocks weight.
integer(I4P)                    :: ideal_block_weight   !< Ideal block weight for load balancing.
integer(I4P)                    :: nblocks_proc         !< Ideal number of blocks per process.
real(R4P), allocatable          :: rcc(:)               !< rcc unstructured array.
integer(I4P)                    :: na                   !< Number of command line arguments.
logical                         :: save_block_tecplot   !< Save blocks also in tecplot (ASCII) format.
integer(I4P)                    :: i,b                  !< Counter.

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

   ! load balancing
   nblocks_proc = blocks_number / procs_number
   print *, 'load balancing stats'
   print *, '  ideal number of blocks per process', nblocks_proc
   total_blocks_weight = 0
   do b=1, blocks_number
      print *, '    block ',b, ' weight ', blocks(b)%weight
      total_blocks_weight = total_blocks_weight + blocks(b)%weight
   enddo
   ideal_block_weight = total_blocks_weight / blocks_number
   print *, '  ideal block weight for load balancing ', ideal_block_weight
   print *, '  list of blocks to be splitted'
   do b=1, blocks_number
      if (blocks(b)%weight>ideal_block_weight) &
      print *, '    block ',b,' nijk ',blocks(b)%Ni,blocks(b)%Nj,blocks(b)%Nk,' weight ',blocks(b)%weight,'>',ideal_block_weight
   enddo

   ! call blocks(2)%split(mgl=2, is_split_done=is_split_done, sb=sb)
   ! if (is_split_done) then
   !    print *, 'parent block           ',2,blocks(2)%Ni,blocks(2)%Nj,blocks(2)%Nk
   !    print *, 'first split block      ',sb(1)%ab,sb(1)%Ni,sb(1)%Nj,sb(1)%Nk
   !    print *, 'second split block     ',sb(2)%ab,sb(2)%Ni,sb(2)%Nj,sb(2)%Nk
   !    print *, 'first split block adj  ',sb(1)%adj(1:4,1),sb(1)%adj(1:4,2)
   !    print *, 'second split block adj ',sb(2)%adj(1:4,1),sb(2)%adj(1:4,2)
   !    ! ab shift
   !    do b=sb(2)%ab, blocks_number
   !       blocks(b)%ab = blocks(b)%ab + 1
   !    enddo
   !    ! update blocks list
   !    if     (sb(1)%ab==1) then ! first block in list has been splitted
   !       blocks = [sb,blocks(2:blocks_number)]
   !    elseif (sb(1)%ab==blocks_number) then ! last  block in list has been splitted
   !       blocks = [blocks(1:blocks_number-1),sb]
   !    else
   !       blocks = [blocks(1:sb(1)%ab-1),sb,blocks(sb(1)%ab+1:blocks_number)]
   !    endif
   !    blocks_number = blocks_number + 1
   ! endif
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
   function is_file_found(file_name) result(is_found)
   !< Inquire is the file path is valid and the file is found.
   character(*), intent(in) :: file_name !< File name.
   logical                  :: is_found  !< Inquiring result.

   inquire(file=trim(adjustl(file_name)), exist=is_found)
   endfunction is_file_found
endprogram overset_exploded
