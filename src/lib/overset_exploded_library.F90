module overset_exploded_library
!< Overset-Exploded, module library.

!< When a block is split and replaced by 2 split blocks all the patches must be sanitized.
!<```
!<    +---------------+---------------+---------------+
!<    |      P04      |      P08      |      P12      |
!<    |               |               |               |
!<    |      B01      |      B02      |      B03      |
!<    |               |               |               |
!<    |P01         P02|P05         P06|P09         P10|
!<    |               |               |               |
!<    |               |               |               |
!<    |               |               |               |
!<    |      P03      |      P07      |      P11      |
!<    +---------------+---------------+---------------+
!<    |      P16      |      P20      |      P24      |
!<    |               |               |               |
!<    |      B04      |      B05      |      B06      |
!<    |               |               |               |
!<    |P13         P14|P17   SB2   P18|P21         P22|
!<    |         SP1402|_______________|SP2102         |
!<    |         SP1401|      SB1      |SP2101         |
!<    |               |               |               |
!<    |      P15      |      P19      |      P23      |
!<    +---------------+---------------+---------------+
!<    |      P28      |      P32      |      P36      |
!<    |               |               |               |
!<    |      B07      |      B08      |      B09      |
!<    |               |               |               |
!<    |P25         P26|P29         P30|P33         P34|
!<    |               |               |               |
!<    |               |               |               |
!<    |               |               |               |
!<    |      P27      |      P31      |      P35      |
!< j  +---------------+---------------+---------------+
!< ^
!< |
!< |
!< o---->i
!<
!<```
!< Asssuming that the block B05 is split in j direction (SB1 and SB2 are the split blocks created), the patches P14 and P21 must be
!< split accordingly because they are connections between blocks B04, B05 and B06. Thus the rules are that all patch references to
!< patches between P15-P20 must be shifted of 1 place; all references to P21 must be shifted of 7 places; all references to patches
!< higher than P21 must be shift of 8 places. The ruls is to add 1 place for each split patch and 6 places for each split block.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use            :: overset_exploded_block_object
use            :: overset_exploded_patch_object
use            :: overset_exploded_process_object

implicit none

private
public :: balance_workload
public :: create_blocks_list
public :: implode_blocks
public :: load_file_grd
public :: load_file_icc
public :: popout_blocks_list
public :: save_file_grd
public :: save_file_icc
public :: save_proc_input
public :: str, strz
public :: update_blocks

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
   ! public procedures
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

   print '(A)', 'perform work load balancing using file '//trim(adjustl(file_name))

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
   type(block_object), intent(inout), allocatable :: blocks(:)        !< Blocks data.
   type(block_object), intent(inout)              :: sb(1:)           !< Split blocks.
   integer(I4P),       intent(out)                :: blocks_number    !< Blocks number.
   logical,            intent(in)                 :: use_cc_par       !< Use cc.par instead of icc.
   type(block_object), allocatable                :: blocks_(:)       !< New blocks data.
   type(patch_object)                             :: split_patches(2) !< Split patches couple.
   type(patch_object)                             :: new_patch        !< New patch.
   logical                                        :: found_new_patch  !< Sentinel to check new patch.
   integer(I4P)                                   :: sbm_ijk          !< Split block minimum ijk patch extent.
   integer(I4P)                                   :: new_ijk          !< New patch minimum ijk extent.
   integer(I4P)                                   :: b, p, pp, ppp    !< Counter.

   if (use_cc_par.and.allocated(sb(1)%patches)) then
      call sanitize_patches(blocks=blocks, sb=sb)
   else
      call sanitize_chimera(blocks=blocks, sb=sb)
   endif
   blocks_number = size(blocks,dim=1)
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
   if (use_cc_par.and.allocated(sb(1)%patches)) then
      ! update block_index of patches
      do b=1, blocks_number
         do p=1, size(blocks(b)%patches, dim=1)
            blocks(b)%patches(p)%block_index = blocks(b)%ab
         enddo
      enddo
      ! update patch_index of patches
      pp = 0
      do b=1, blocks_number
         do p=1, size(blocks(b)%patches, dim=1)
            pp = pp + 1
            blocks(b)%patches(p)%patch_index = pp
         enddo
      enddo
      ! fix internal adjacencies cross-references (BC=-123)
      do p=1, size(blocks(sb(1)%ab)%patches, dim=1)
         if (blocks(sb(1)%ab)%patches(p)%boundary_condition == -123) then
            ! this is sb(1)'s face_max adj, find sb(2)'s face_min adj
            do pp=1, size(blocks(sb(2)%ab)%patches, dim=1)
               if (blocks(sb(2)%ab)%patches(pp)%boundary_condition == -123) then
                  blocks(sb(1)%ab)%patches(p)%connect_family  = blocks(sb(2)%ab)%patches(pp)%patch_index
                  blocks(sb(2)%ab)%patches(pp)%connect_family = blocks(sb(1)%ab)%patches(p)%patch_index
                  exit
               endif
            enddo
            exit
         endif
      enddo
      ! change sign to ijk orientation of the new split block patch
      do p=1, size(blocks(sb(1)%ab)%patches, dim=1)
         if (blocks(sb(1)%ab)%patches(p)%boundary_condition == -123) then
            blocks(sb(1)%ab)%patches(p)%boundary_condition = 123
            exit
         endif
      enddo
      do p=1, size(blocks(sb(2)%ab)%patches, dim=1)
         if (blocks(sb(2)%ab)%patches(p)%boundary_condition == -123) then
            blocks(sb(2)%ab)%patches(p)%boundary_condition = 123
            exit
         endif
      enddo
   endif
   endsubroutine update_blocks

   ! private procedures
   subroutine sanitize_chimera(blocks, sb)
   !< Sanitize chimera data after a block split.
   type(block_object), intent(inout) :: blocks(:) !< Blocks data.
   type(block_object), intent(in)    :: sb(1:)    !< Split blocks.
   integer(I4P)                      :: b         !< Counter.

   do b=1, size(blocks,dim=1)
      if (b==sb(1)%ab) cycle ! split block does not need to be santized, it is replaced by sb
      call blocks(b)%sanitize_chimera(sb=sb)
   enddo
   endsubroutine sanitize_chimera

   subroutine sanitize_patches(blocks, sb)
   !< Sanitize all patches after a block split.
   type(block_object), intent(inout) :: blocks(:)         !< Blocks data.
   type(block_object), intent(inout) :: sb(1:)            !< Split blocks.
   integer(I4P)                      :: blocks_number     !< Blocks number.
   integer(I4P), allocatable         :: remote_splits(:)  !< Patch indexes of split patches.
   integer(I4P)                      :: last_parent_patch !< Last parent patch index.
   integer(I4P)                      :: b,p,pp            !< Counter.
   integer(I4P)                      :: temp              !< Temporary buffer for sorting remote_splits.
   logical                           :: swapped           !< Sentinel for sorting remote_splits.

   blocks_number = size(blocks,dim=1)
   remote_splits = [0] ! start with a fake patch index to avoid 2 steps count-allocation/filling
   do p=1, size(sb(1)%patches,dim=1)
      if (sb(1)%patches(p)%block_index<0.and.sb(1)%patches(p)%is_connection) then
         remote_splits = [remote_splits, sb(1)%patches(p)%connect_family]
      endif
   enddo
   if (size(remote_splits,dim=1)>1) then
      remote_splits = remote_splits(2:) ! trim out the start fake patch index
      ! sort remote_splits ascending
      sort_split_patches : do
         swapped = .false.
         do p=1, size(remote_splits,dim=1)-1
            if (remote_splits(p) > remote_splits(p+1)) then
               temp = remote_splits(p)
               remote_splits(p) = remote_splits(p+1)
               remote_splits(p+1) = temp
               swapped = .true.
            endif
         enddo
         if (.not.swapped) exit sort_split_patches
      enddo sort_split_patches
      call create_remote_split_patches(blocks=blocks, sb=sb, remote_splits=remote_splits)
      last_parent_patch = blocks(sb(1)%ab)%patches(size(blocks(sb(1)%ab)%patches,dim=1))%patch_index
      do b=1, blocks_number
         if (b == sb(1)%ab) cycle ! parent block will be replaced by sb
         call shift_patches(patches=blocks(b)%patches, remote_splits=remote_splits, last_parent_patch=last_parent_patch)
      enddo
      call shift_patches(patches=sb(1)%patches, remote_splits=remote_splits, last_parent_patch=last_parent_patch)
      call shift_patches(patches=sb(2)%patches, remote_splits=remote_splits, last_parent_patch=last_parent_patch)

      do p=1, size(sb(2)%patches, dim=1)
         if (sb(2)%patches(p)%block_index < 0 .and. sb(2)%patches(p)%is_connection) then
            sb(2)%patches(p)%connect_family = sb(2)%patches(p)%connect_family + 1
         endif
      enddo

      ! fix copies: set connect_family = original's shifted value + size(sb(1)%patches)
      do b=1, blocks_number
         if (b == sb(1)%ab) cycle
         do p=1, size(blocks(b)%patches, dim=1)
            if (blocks(b)%patches(p)%connect_family < 0) then
               do pp=1, size(blocks(b)%patches, dim=1)
                  if (blocks(b)%patches(pp)%patch_index == -blocks(b)%patches(p)%connect_family &
                      .and. blocks(b)%patches(pp)%connect_family > 0) then
                     blocks(b)%patches(p)%connect_family = blocks(b)%patches(pp)%connect_family &
                                                         + size(sb(1)%patches, dim=1)
                     exit
                  endif
               enddo
            endif
         enddo
      enddo
   endif
   endsubroutine sanitize_patches

   subroutine shift_patches(patches, remote_splits, last_parent_patch)
   !< Shift connect_family of all connection patches of a block.
   type(patch_object), intent(inout) :: patches(:)        !< Block patches.
   integer(I4P),       intent(in)    :: remote_splits(:)  !< Sorted ascending remote split patch indices.
   integer(I4P),       intent(in)    :: last_parent_patch !< Last patch index of parent block.
   integer(I4P)                      :: p, pp, shift      !< Counter and shift.

   do p=1, size(patches,dim=1)
      if (patches(p)%is_connection .and. patches(p)%boundary_condition /= -123) then
         shift = 0
         do pp=size(remote_splits,dim=1), 1, -1
            if (patches(p)%connect_family > remote_splits(pp)) then
               shift = pp
               exit
            endif
         enddo
         if (patches(p)%connect_family > last_parent_patch) shift = shift + 6
         patches(p)%connect_family = patches(p)%connect_family + shift
      endif
   enddo
   endsubroutine shift_patches

   subroutine create_remote_split_patches(blocks, sb, remote_splits)
   !< Create and clip copies of remote split patches in adjacent blocks.
   !< For each remote split, the original patch is clipped to match sb(1) and a new copy is
   !< created and clipped to match sb(2). The copy's connect_family is set negative as sentinel
   !< so that shift_patches will skip it; it must be fixed up afterwards.
   type(block_object), intent(inout) :: blocks(:)        !< Blocks data.
   type(block_object), intent(in)    :: sb(1:)           !< Split blocks.
   integer(I4P),       intent(in)    :: remote_splits(:) !< Sorted ascending remote split patch indices.
   type(patch_object)                :: new_patch        !< Copy of remote split patch.
   integer(I4P)                      :: ns1              !< Cells in split direction for sb(1).
   integer(I4P)                      :: orient_digit     !< Orientation digit for split direction.
   integer(I4P)                      :: dir_lo, dir_hi   !< Extent indices of remote split direction.
   integer(I4P)                      :: mid              !< Split point in remote block coordinates.
   integer(I4P)                      :: b, p, r, pp      !< Counters.

   ! compute ns1
   select case(sb(1)%split_dir)
   case(1); ns1 = sb(1)%Ni
   case(2); ns1 = sb(1)%Nj
   case(3); ns1 = sb(1)%Nk
   end select

   do b=1, size(blocks, dim=1)
      if (b == sb(1)%ab) cycle
      loop_block_patches : do p=1, size(blocks(b)%patches, dim=1)
         do r=1, size(remote_splits, dim=1)
            if (blocks(b)%patches(p)%patch_index == remote_splits(r)) then
               new_patch = blocks(b)%patches(p)
               ! find the sb(1) perpendicular patch connecting to this remote patch
               do pp=1, size(sb(1)%patches, dim=1)
                  if (sb(1)%patches(pp)%connect_family == remote_splits(r)) then
                     select case(sb(1)%split_dir)
                     case(1); orient_digit = sb(1)%patches(pp)%boundary_condition / 100
                     case(2); orient_digit = mod(sb(1)%patches(pp)%boundary_condition / 10, 10)
                     case(3); orient_digit = mod(sb(1)%patches(pp)%boundary_condition, 10)
                     end select
                     exit
                  endif
               enddo
               ! determine which extent pair to clip in the remote block
               select case(orient_digit)
               case(1,2); dir_lo = 1; dir_hi = 2
               case(3,4); dir_lo = 3; dir_hi = 4
               case(5,6); dir_lo = 5; dir_hi = 6
               end select
               ! clip extents
               if (mod(orient_digit, 2) == 1) then
                  ! forward: original keeps low half (sb(1)), copy keeps high half (sb(2))
                  mid = blocks(b)%patches(p)%ijk_extents(dir_lo) + ns1
                  blocks(b)%patches(p)%ijk_extents(dir_hi) = mid
                  new_patch%ijk_extents(dir_lo) = mid
               else
                  ! reversed: original keeps high half (sb(1)), copy keeps low half (sb(2))
                  mid = blocks(b)%patches(p)%ijk_extents(dir_hi) - ns1
                  blocks(b)%patches(p)%ijk_extents(dir_lo) = mid
                  new_patch%ijk_extents(dir_hi) = mid
               endif
               ! negative sentinel: shift_patches will skip this, fixed up afterwards
               new_patch%connect_family = -remote_splits(r)
               blocks(b)%patches = [blocks(b)%patches(1:p), new_patch, blocks(b)%patches(p+1:)]
               exit loop_block_patches
            endif
         enddo
      enddo loop_block_patches
   enddo
   endsubroutine create_remote_split_patches

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
endmodule overset_exploded_library
