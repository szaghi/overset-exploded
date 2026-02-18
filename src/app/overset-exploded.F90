program overset_exploded
!< Overset-Exploded program, convert overset output files into exploded per-block files.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use            :: overset_exploded_box_object
use            :: overset_exploded_block_object
use            :: overset_exploded_cc_par_object
use            :: overset_exploded_process_object
use            :: overset_exploded_library

implicit none

character(len=99)                 :: file_name_grd        !< Grid file name.
character(len=99)                 :: file_name_icc        !< Icc file name.
character(len=99)                 :: file_name_input      !< Input file name: grd, icc or infocc.out.
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

call parse_command_line(fgrd=file_name_grd,ficc=file_name_icc,                                             &
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

call balance_workload(file_name=file_name_input,use_cc_par=use_cc_par,mgl=mgl,                              &
                      max_unbalance=max_unbalance,procs_number=procs_number,save_bsplit_par=save_bsplit_par,&
                      processes=processes,splits=splits,splits_dir=splits_dir,splits_nijk=splits_nijk)

if (use_cc_par) then
   print '(A)', 'load blocks grid from file '//trim(adjustl(file_name_input))
   call load_file_grd(file_name=file_name_input,blocks=blocks,blocks_number=blocks_number,gc=0)
   print '(A)', 'load parameter file '//trim(adjustl(fname_cc_par))
   call cc_par%load_file(file_name=fname_cc_par,blocks=blocks,boxes=boxes)
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

if (allocated(splits)) then
   print '(A)', 'split blocks'
   splits_loop : do b=1, size(splits,dim=1)
      print '(A)', '  block '//trim(str(splits(b),.true.))
      ! print*, 'cazzo patches before split'
      ! do p=1, size(blocks(splits(b)-2)%patches,dim=1)
      !    print*, 'cazzo ',trim(blocks(splits(b)-2)%patches(p)%description())
      ! enddo
      ! do p=1, size(blocks(splits(b)-1)%patches,dim=1)
      !    print*, 'cazzo ',trim(blocks(splits(b)-1)%patches(p)%description())
      ! enddo
      ! do p=1, size(blocks(splits(b))%patches,dim=1)
      !    print*, 'cazzo ',trim(blocks(splits(b))%patches(p)%description())
      ! enddo
      ! do p=1, size(blocks(splits(b)+1)%patches,dim=1)
      !    print*, 'cazzo ',trim(blocks(splits(b)+1)%patches(p)%description())
      ! enddo
      call blocks(splits(b))%split(mgl=mgl, is_split_done=is_split_done, sb=sb, split_data=.true., use_cc_par=use_cc_par)
      if (is_split_done) then
         call update_blocks(blocks=blocks, sb=sb, blocks_number=blocks_number, use_cc_par=use_cc_par)
         ! print*, 'cazzo split dir ', sb(1)%split_dir
         ! print*, 'cazzo patches after split'
         ! do p=1, size(blocks(splits(b)-2)%patches,dim=1)
         !    print*, 'cazzo ',trim(blocks(splits(b)-2)%patches(p)%description())
         ! enddo
         ! do p=1, size(blocks(splits(b)-1)%patches,dim=1)
         !    print*, 'cazzo ',trim(blocks(splits(b)-1)%patches(p)%description())
         ! enddo
         ! do p=1, size(blocks(splits(b))%patches,dim=1)
         !    print*, 'cazzo ',trim(blocks(splits(b))%patches(p)%description())
         ! enddo
         ! do p=1, size(blocks(splits(b)+1)%patches,dim=1)
         !    print*, 'cazzo ',trim(blocks(splits(b)+1)%patches(p)%description())
         ! enddo
         ! do p=1, size(blocks(splits(b)+2)%patches,dim=1)
         !    print*, 'cazzo ',trim(blocks(splits(b)+2)%patches(p)%description())
         ! enddo
         ! stop
         print '(A)', '     new blocks number '//trim(str(blocks_number,.true.))
         ! free sb memory immediately after update_blocks copies the data
         call sb(1)%destroy
         call sb(2)%destroy
      else
         write(stderr,'(A)')'error: unable to split block "'//trim(str(splits(b)))//'"'
         stop
      endif
   enddo splits_loop
endif

if (use_cc_par) call cc_par%save_file(file_name='split-balanced-cc.par',blocks=blocks,boxes=boxes)

! update blocks to processes assignment
do p=0, procs_number - 1
   do b=1, size(processes(p)%blocks,dim=1)
      bb = processes(p)%blocks(b)
      if (bb>0) blocks(bb)%proc = p
   enddo
enddo
call save_proc_input(blocks=blocks, file_name='split-balanced-proc.input')

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
   subroutine parse_command_line(fgrd,ficc,stec,simp,sexp,sbsp,uccp,fccp,ebn,np,mu,mgl)
   !< Parse command line inputs.
   character(*), intent(out) :: fgrd      !< Grid file name.
   character(*), intent(out) :: ficc      !< Icc file name.
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
