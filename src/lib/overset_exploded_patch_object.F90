module overset_exploded_patch_object
!< Overset-Exploded, definition of class patch_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32, stderr=>error_unit

implicit none

private
public :: patch_object

type :: patch_object
  !< Patch class.
  logical      :: is_connection=.false.       !< Flag to inquire is the patch defines a connection.
  integer(I4P) :: patch_index=0               !< Patch index, local numeration.
  integer(I4P) :: block_index=0               !< Block to which patch belongs, local numeration.
  integer(I4P) :: face_index=0                !< Face index: Imin=>1, Imax=2,Jmin=3, Jmax=4, Kmin=5, Kmax=6.
  integer(I4P) :: boundary_condition=0        !< Boundary condition (or IJK orientation) in the overset convention.
  integer(I4P) :: connect_family=0            !< Index of connected patch or family index of the patch.
  integer(I4P) :: ijk_extents(6)=[0,0,0,0,0,0]!< IJK extents of the patch: Imin,Imax,Jmin,Jmax,Kmin,Kmax.
  contains
     ! public methods
     procedure, pass(self) :: description    !< Return pretty printed description.
     procedure, pass(self) :: load_from_file !< Load patch from (opened) file unit.
     procedure, pass(self) :: save_into_file !< Save patch into (opened) file unit.
endtype patch_object

contains
   ! public methods
   function description(self) result(dsc)
   !< Return pretty printed description.
   class(patch_object), intent(in) :: self !< Patch data.
   character(:), allocatable       :: dsc  !< Pretty printed description.

   dsc = repeat(' ',150)
   write(dsc, '(10(I9,1X),A3,I9)') self%block_index,        &
                                   self%face_index,         &
                                   self%boundary_condition, &
                                   self%connect_family,     &
                                   self%ijk_extents(1:6),   &
                                   ' ! ',                   &
                                   self%patch_index
   endfunction description

   subroutine load_from_file(self, file_unit, patch_index)
   !< Load patch from (opened) file unit.
   class(patch_object), intent(inout) :: self        !< Patch data.
   integer(I4P),        intent(in)    :: file_unit   !< File unit.
   integer(I4P),        intent(in)    :: patch_index !< Current patch index (global).

   self%patch_index = patch_index
   read(file_unit, *) self%block_index,        &
                      self%face_index,         &
                      self%boundary_condition, &
                      self%connect_family,     &
                      self%ijk_extents(1:6)
   ! a patch is a connection only if connect_family > 0 AND BC >= 100 (connection orientation code)
   ! patches with BC < 100 (wall=1, outflow=40, etc.) are not connections even if connect_family is non-zero
   self%is_connection = ((self%connect_family>0).and.(self%boundary_condition>=100))
   endsubroutine load_from_file

   subroutine save_into_file(self, file_unit)
   !< Save patch into (opened) file unit.
   class(patch_object), intent(in) :: self      !< Patch data.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit, '(10(I9,1X),A,I9)') self%block_index,        &
                                        self%face_index,         &
                                        self%boundary_condition, &
                                        self%connect_family,     &
                                        self%ijk_extents(1:6),   &
                                        ' ! ', self%patch_index
   endsubroutine save_into_file
endmodule overset_exploded_patch_object
