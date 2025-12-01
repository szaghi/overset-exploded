module oe_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32

implicit none

private
public :: block_object

type :: block_object
   !<  Block class.
   integer(I4P)              :: Ni=0              !< Number of cells in i direction.
   integer(I4P)              :: Nj=0              !< Number of cells in j direction.
   integer(I4P)              :: Nk=0              !< Number of cells in k direction.
   integer(I4P)              :: gc=2              !< Number of ghost cells.
   real(R8P),    allocatable :: nodes(:,:,:,:)    !< Nodes coordinates.
   integer(I4P), allocatable :: icc(:,:,:)        !< Cell centered icc values.
   integer(I4P), allocatable :: adj(:,:)          !< New adjacent-BC indexes for split blocks [1:4,nadj].
   integer(I4P)              :: nadj=0            !< Number of new adjacent-BC cells.
   integer(I4P)              :: ab=0              !< Absolute block index.
   integer(I4P)              :: group=0           !< Index of gruop.
   integer(I4P)              :: body=0            !< Index of body.
   integer(I4P)              :: proc=0            !< Processor assigned to.
   logical                   :: is_loaded=.false. !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy dynamic memory.
      procedure, pass(self) :: alloc           !< Allocate dynamic memory.
      procedure, pass(self) :: load_dimensions !< Load block dimensions from file.
      procedure, pass(self) :: load_icc        !< Load block icc from file.
      procedure, pass(self) :: load_nodes      !< Load block nodes from file.
      procedure, pass(self) :: save_block_file !< Save block data into its own file.
      procedure, pass(self) :: split           !< Split block.
      procedure, pass(self) :: traslate        !< Traslate block nodes by a given traslation vector.
endtype block_object

contains
   elemental subroutine destroy(self)
   !< Destroy dynamic memory.
   class(block_object), intent(inout) :: self !< Block data.

   self%Ni = 0
   self%Nj = 0
   self%Nk = 0
   self%gc = 2
   if (allocated(self%nodes)) deallocate(self%nodes)
   if (allocated(self%icc)) deallocate(self%icc)
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

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nadj=>self%nadj)
   allocate(self%nodes(1:3,0-gc:Ni+gc,0-gc:Nj+gc,0-gc:Nk+gc))
   allocate(self%icc(1-gc:Ni+gc,1-gc:Nj+gc,1-gc:Nk+gc))
   if (nadj>0) allocate(self%adj(1:4,1:nadj))
   endassociate
   endsubroutine alloc

   subroutine load_dimensions(self, ab, file_unit)
   !< Load block dimensions from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_object), intent(inout) :: self      !< Block data.
   integer(I4P),        intent(in)    :: ab        !< Absolute block index.
   integer(I4P),        intent(in)    :: file_unit !< Logical unit of grd file.

   call self%destroy
   read(file_unit, end=10, err=10) self%Ni, self%Nj, self%Nk, self%gc
   10 call self%alloc
   self%ab = ab
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

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, nodes => self%nodes)
      read(file_unit)(((nodes(1,i, j, k), i=0-gc,Ni+gc), j=0-gc,Nj+gc), k=0-gc,Nk+gc)
      read(file_unit)(((nodes(2,i, j, k), i=0-gc,Ni+gc), j=0-gc,Nj+gc), k=0-gc,Nk+gc)
      read(file_unit)(((nodes(3,i, j, k), i=0-gc,Ni+gc), j=0-gc,Nj+gc), k=0-gc,Nk+gc)
   endassociate
   endsubroutine load_nodes

   subroutine save_block_file(self, b, rcc)
   !< Save block data into its own file.
   class(block_object), intent(in) :: self      !< Block data.
   integer(I4P),        intent(in) :: b         !< Current block number, global numeration.
   real(R4P),           intent(in) :: rcc(1:)   !< rcc unstructured array.
   character(len=6)                :: bstr      !< Block number stringified.
   integer(I4P)                    :: file_unit !< Block unit file.
   integer(I4P)                    :: i,j,k,p,n !< Counter.
   integer(I4P)                    :: offset    !< Offset of rcc.

   write(bstr, '(I6.5)') b
   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes,icc=>self%icc)
   ! file grid
   open(newunit=file_unit, file='block-'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write')
   write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
   write(file_unit)(((nodes(1,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   write(file_unit)(((nodes(2,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   write(file_unit)(((nodes(3,i,j,k),i=0-gc,Ni+gc),j=0-gc,Nj+gc),k=0-gc,Nk+gc)
   close(file_unit)
   ! file rcc
   open(newunit=file_unit, file='block-rcc-'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write')
   write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
   do k=1-gc, Nk+gc
   do j=1-gc, Nj+gc
   do i=1-gc, Ni+gc
      p = icc(i,j,k)
      if (p>0) then
         select case(int(rcc(p),I4P))
         case(-1 )
            write(file_unit) i,j,k,'WALL'
         case(-2 )
            write(file_unit) i,j,k,'SIMMETRY'
         case(-3 )
            write(file_unit) i,j,k,'INFLOW'
         case(-4 )
            write(file_unit) i,j,k,'INOUTFLOW'
         case(-5 )
            write(file_unit) i,j,k,'ASSIGNED-INFLOW'
         case(-6 )
            write(file_unit) i,j,k,'ASSIGNED-PRESSURE'
         case(-7 )
            write(file_unit) i,j,k,'ASSIGNED-NORMAL-VELOCITY'
         case(-8 )
            write(file_unit) i,j,k,'ASSIGNED-RIEMANN'
         case(-9 )
            write(file_unit) i,j,k,'EXTRAPOLATED'
         case(-10)
            write(file_unit) i,j,k,'MOVING-WALL'
         case(-11)
            write(file_unit) i,j,k,'INACTIVE-WALL'
         case(-19)
            write(file_unit) i,j,k,'EXTRAPOLATED-ALT'
         case(0)
            ! save nothing for active cell
         case(20:46)
            write(file_unit) i,j,k,'CHIMERA'
            write(file_unit) nint(rcc(p),I4P),nint(rcc(p+1),I4P) ! type, donors number
            do n=1, nint(rcc(p+1),I4P)
               offset = p + 1 + 5*(n-1)
               write(file_unit) nint(rcc(offset+1)), nint(rcc(offset+2)), nint(rcc(offset+3)), nint(rcc(offset+4)) ! b,i,j,k donor
               write(file_unit) rcc(offset+5) ! donor weight
            enddo
         case(60:66)
            write(file_unit) i,j,k,'ADJACENT'
            write(file_unit) nint(rcc(p+2)), nint(rcc(p+3)), nint(rcc(p+4)), nint(rcc(p+5)) ! b,i,j,k adjacent
         case(80)
         case default
            print *, 'error: unknown icc "',int(rcc(p),I4P),'", b,i,j,k=',b,i,j,k
            stop
         endselect
      elseif (p<0) then
         ! this is a cell of a new block coming from a split, the new adjacent BC are in adj array
         write(file_unit) i,j,k,'ADJACENT'
         write(file_unit) self%adj(1:4,-p) ! b,i,j,k adjacent
      endif
   enddo
   enddo
   enddo
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

   pure subroutine traslate(self, traslation)
   !< Traslate block nodes by a given traslation vector.
   class(block_object), intent(inout) :: self          !< Block data.
   real(R8P),           intent(in)    :: traslation(3) !< Traslation vector.
   integer(I4P)                       :: i,j,k         !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes)
   if (allocated(self%nodes)) then
      do k=1-self%gc,self%Nk+self%gc
      do j=1-self%gc,self%Nj+self%gc
      do i=1-self%gc,self%Ni+self%gc
         self%nodes(1:3,i,j,k) = self%nodes(1:3,i,j,k) + traslation(1:3)
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine traslate
endmodule oe_block_object

program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use oe_block_object

implicit none

character(len=99)               :: file_name_grd        !< Grid file name.
character(len=99)               :: file_name_icc        !< Icc file name.
integer(I4P)                    :: file_unit_grd        !< Grid unit file.
integer(I4P)                    :: file_unit_icc        !< Icc unit file.
integer(I4P)                    :: blocks_number=0      !< Number of blocks contained into the files.
type(block_object), allocatable :: blocks(:)            !< Blocks contained into the files.
logical                         :: is_split_done        !< Sentinel to check is split has been done.
type(block_object)              :: sb(2)                !< Split blocks.
integer(I4P)                    :: unstruct_dimension=0 !< Dimension of unstructured array of rcc.
real(R4P), allocatable          :: rcc(:)               !< rcc unstructured array.
integer(I4P)                    :: na                   !< Number of command line arguments.
integer(I4P)                    :: i,b                  !< Counter.

na = command_argument_count()
if (na==0) then
   write(stderr, "(A)")'error: you must pass file_name_grd and file_name_icc as command line arguments'
   write(stderr, "(A)")'example: overset-exploded cc.01.grd cc.01'
   stop
else
   call get_command_argument(1, file_name_grd)
   call get_command_argument(2, file_name_icc)
endif

if (is_file_found(file_name_grd).and.is_file_found(file_name_icc)) then
   ! grd file
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
   open(newunit=file_unit_icc, file=trim(adjustl(file_name_icc)), form='unformatted', action='read')
   read(file_unit_icc) b
   if (b/=blocks_number) then
      write(stderr, "(A)")'error: grd and icc have different number of blocks'
      stop
   endif
   do b=1, blocks_number
      call blocks(b)%load_dimensions(ab=b, file_unit=file_unit_icc)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_icc(file_unit=file_unit_icc)
   enddo
   read(file_unit_icc) unstruct_dimension
   allocate(rcc(1:unstruct_dimension))
   read(file_unit_icc) (rcc(i),i=1,unstruct_dimension)
   close(file_unit_icc)

   call blocks(2)%split(mgl=2, is_split_done=is_split_done, sb=sb)
   if (is_split_done) then
      print *, 'parent block           ',2,blocks(2)%Ni,blocks(2)%Nj,blocks(2)%Nk
      print *, 'first split block      ',sb(1)%ab,sb(1)%Ni,sb(1)%Nj,sb(1)%Nk
      print *, 'second split block     ',sb(2)%ab,sb(2)%Ni,sb(2)%Nj,sb(2)%Nk
      print *, 'first split block adj  ',sb(1)%adj(1:4,1),sb(1)%adj(1:4,2)
      print *, 'second split block adj ',sb(2)%adj(1:4,1),sb(2)%adj(1:4,2)
      ! ab shift
      do b=sb(2)%ab, blocks_number
         blocks(b)%ab = blocks(b)%ab + 1
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
   endif
   ! stop

   do b=1, blocks_number
      call blocks(b)%save_block_file(b=b, rcc=rcc)
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
