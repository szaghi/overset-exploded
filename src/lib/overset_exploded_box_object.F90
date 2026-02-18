module overset_exploded_box_object
!< Overset-Exploded, definition of class box_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit

implicit none

private
public :: box_object

type :: box_object
   !< Box class.
   integer(I4P) :: btype          !< Box type.
   integer(I4P) :: bblock         !< Box block.
   integer(I4P) :: bgroup         !< Box group.
   real(R8P)    :: nodes(1:3,1:8) !< Box nodes coordinates.
   contains
     ! public methods
     procedure, pass(self) :: load_from_file !< Load box from (opened) file unit.
     procedure, pass(self) :: save_into_file !< Save box into (opened) file unit.
endtype box_object

contains
   ! public methods
   subroutine load_from_file(self, file_unit)
   !< Load box from (opened) file unit.
   class(box_object), intent(inout) :: self      !< Box data.
   integer(I4P),      intent(in)    :: file_unit !< File unit.
   integer(I4P)                     :: n         !< Counter.

   read(file_unit, *) self%btype, self%bblock, self%bgroup
   do n=1, 8 ! box xyz nodes loop
      read(file_unit, *) self%nodes(1,n), self%nodes(2,n), self%nodes(3,n)
   enddo
   endsubroutine load_from_file

   subroutine save_into_file(self, file_unit)
   !< Save box into (opened) file unit.
   class(box_object), intent(in) :: self      !< Box data.
   integer(I4P),      intent(in) :: file_unit !< File unit.
   integer(I4P)                  :: n         !< Counter.

   write(file_unit, '(I9,1X,I9,1X,I9,1X,A)') self%btype, self%bblock, self%bgroup, ' ! type, block, group associated'
   do n = 1, 8
      write(file_unit, '(3(ES23.12,1X))') self%nodes(1,n), self%nodes(2,n), self%nodes(3,n)
   enddo
   endsubroutine save_into_file
endmodule overset_exploded_box_object
