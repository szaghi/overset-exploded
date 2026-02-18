module overset_exploded_cc_par_object
!< Overset-Exploded, definition of class cc_par_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit
use            :: overset_exploded_block_object
use            :: overset_exploded_box_object
use            :: overset_exploded_patch_object

implicit none

private
public :: cc_par_object

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
   contains
      ! public methods
      procedure, pass(self) :: load_file !< Load cc.par file.
      procedure, pass(self) :: save_file !< Save cc.par file.
endtype cc_par_object

contains
   ! public methods
   subroutine load_file(self, file_name, blocks, boxes)
   !< Load file cc.par.
   class(cc_par_object), intent(inout)              :: self      !< File cc.par handler.
   character(*),         intent(in)                 :: file_name !< File name
   type(block_object),   intent(inout)              :: blocks(:) !< Blocks data.
   type(box_object),     intent(inout), allocatable :: boxes(:)  !< Boxes.
   integer(I4P)                                     :: file_unit !< File unit.

   open(newunit=file_unit, file=trim(adjustl(file_name)))
   call load_header
   call load_blocks
   call load_patches
   call load_edges
   call load_boxes
   close(file_unit)
   contains
      subroutine load_header()
      !< Load header.

      read(file_unit,*) self%file_name_input_grd
      read(file_unit,*) self%base_name_output
      read(file_unit,*) self%save_ghost_cells
      read(file_unit,*) self%increase_overlap
      read(file_unit,*) self%extend_internal_wall
      read(file_unit,*)
      read(file_unit,*) self%mgl(1), self%mgl(2)
      read(file_unit,*)
      read(file_unit,*) self%boundary_layer_thickness
      read(file_unit,*) self%numberical_beach
      read(file_unit,*)
      endsubroutine load_header

      subroutine load_blocks()
      !< Load blocks section.
      integer(I4P) :: b !< Counter.

      read(file_unit, *) self%blocks_number
      if (self%blocks_number /= size(blocks,dim=1)) then
         write(stderr, '(A)')'error: grd and cc.par have different number of blocks'
         write(stderr, '(A,I9)')'       grd blocks number:   ', size(blocks,dim=1)
         write(stderr, '(A,I9)')'       cc.par blocks number:', self%blocks_number
         stop
      endif
      read(file_unit, *)
      do b=1, self%blocks_number
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

      read(file_unit, *) self%patches_number
      read(file_unit, *)
      if (self%patches_number>0) then
         allocate(patches(1:self%patches_number))
         do p=1, self%patches_number
            call patches(p)%load_from_file(file_unit=file_unit, patch_index=p)
         enddo
         read(file_unit, *)
         ! distribute patches to blocks
         allocate(patches_count(1:self%blocks_number), source=0)
         do p=1, self%patches_number
            b = patches(p)%block_index
            patches_count(b) = patches_count(b) + 1
         enddo
         do b=1, self%blocks_number
            allocate(blocks(b)%patches(1:patches_count(b)))
         enddo
         patches_count = 0
         do p=1, self%patches_number
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

      read(file_unit, *) self%edges_number ! keep this always 0
      read(file_unit, *)
      endsubroutine load_edges

      subroutine load_boxes()
      !< Load boxes section.
      integer(I4P) :: b !< Counter.

      if (allocated(boxes)) deallocate(boxes)
      read(file_unit, *) self%boxes_number
      read(file_unit, *)
      if (self%boxes_number>0) then
         allocate(boxes(1:self%boxes_number))
         do b=1, self%boxes_number
            call boxes(b)%load_from_file(file_unit=file_unit)
            read(file_unit, *)
         enddo
      endif
      endsubroutine load_boxes
   endsubroutine load_file

   subroutine save_file(self, file_name, blocks, boxes)
   !< Save file cc.par.
   class(cc_par_object), intent(in)              :: self           !< File cc.par handler.
   character(*),         intent(in)              :: file_name      !< File name.
   type(block_object),   intent(in)              :: blocks(:)      !< Blocks data.
   type(box_object),     intent(in), allocatable :: boxes(:)       !< Boxes data.
   integer(I4P)                                  :: file_unit      !< File unit.
   integer(I4P)                                  :: blocks_number  !< Blocks number.
   integer(I4P)                                  :: patches_number !< Patches number.

   print '(A)', 'save split cc.par file named '//trim(adjustl(file_name))
   blocks_number = size(blocks, dim=1)
   open(newunit=file_unit, file=trim(adjustl(file_name)), action='write', status='replace')
   call save_header
   call save_blocks
   call save_patches
   call save_edges
   call save_boxes
   close(file_unit)
   contains
      subroutine save_header()
      !< Save header section.

      write(file_unit, '(A)') "'split-balanced-cc.grd'"
      write(file_unit, '(A)') "'"//trim(self%base_name_output)//"'"
      write(file_unit, '(A)') merge('.true. ', '.false.', self%save_ghost_cells)
      write(file_unit, '(A)') merge('.true. ', '.false.', self%increase_overlap)
      write(file_unit, '(A)') merge('.true. ', '.false.', self%extend_internal_wall)
      write(file_unit, '(A)') ''
      write(file_unit, '(I3,1X,I3)') self%mgl(1), self%mgl(2)
      write(file_unit, '(A)') ''
      write(file_unit, '(ES12.4)') self%boundary_layer_thickness
      write(file_unit, '(ES12.4)') self%numberical_beach
      write(file_unit, '(A)') ''
      endsubroutine save_header

      subroutine save_blocks()
      !< Save blocks section.
      integer(I4P) :: b !< Counter.

      write(file_unit, '(I0,A)') blocks_number, ' ! blocks number'
      write(file_unit, '(A)') ''
      do b = 1, blocks_number
         write(file_unit, '(I9,1X,I9,1X,I9,1X,A)') blocks(b)%level, blocks(b)%group, blocks(b)%priority, &
            '! '//trim(adjustl(blocks(b)%comment))
      enddo
      write(file_unit, '(A)') ''
      endsubroutine save_blocks

      subroutine save_patches()
      !< Save patches section.
      integer(I4P) :: b, p !< Counters.

      patches_number = 0
      do b = 1, blocks_number
         if (allocated(blocks(b)%patches)) patches_number = patches_number + size(blocks(b)%patches)
      enddo

      write(file_unit, '(I0,A)') patches_number, ' ! patches number'
      write(file_unit, '(A)') ''
      do b = 1, blocks_number
         if (.not.allocated(blocks(b)%patches)) cycle
         do p = 1, size(blocks(b)%patches)
            call blocks(b)%patches(p)%save_into_file(file_unit=file_unit)
         enddo
      enddo
      write(file_unit, '(A)') ''
      endsubroutine save_patches

      subroutine save_edges()
      !< Save edges section (not yet supported).

      write(file_unit, '(I0,A)') 0, ' ! edges number'
      write(file_unit, '(A)') ''
      endsubroutine save_edges

      subroutine save_boxes()
      !< Save boxes section.
      integer(I4P) :: b            !< Counter.
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
            call boxes(b)%save_into_file(file_unit=file_unit)
            write(file_unit, '(A)') ''
         enddo
      endif
      endsubroutine save_boxes
   endsubroutine save_file
endmodule overset_exploded_cc_par_object
