module overset_exploded_process_object
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
endmodule overset_exploded_process_object
