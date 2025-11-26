module oe_block_object
!< Overset-Exploded, definition of class block_object.

use, intrinsic :: iso_fortran_env, only : R4P=>real32, R8P=>real64, I4P=>int32

implicit none

private
public :: block_object

type :: block_object
   !<  Block class.
   integer(I4P)              :: Ni=0                   !< Number of cells in i direction.
   integer(I4P)              :: Nj=0                   !< Number of cells in j direction.
   integer(I4P)              :: Nk=0                   !< Number of cells in k direction.
   integer(I4P)              :: gc(6)=[2,2,2,2,2,2]    !< Number of ghost cells.
   real(R8P),    allocatable :: nodes(:,:,:,:)         !< Nodes coordinates.
   real(R8P),    allocatable :: centers(:,:,:,:)       !< Centers coordinates.
   real(R8P),    allocatable :: extents(:,:)           !< Box extents, [min, max].
   integer(I4P), allocatable :: icc(:,:,:)             !< Cell centered icc values.
   real(R4P),    allocatable :: rcc(:,:,:)             !< Cell centered rcc values.
   integer(I4P), allocatable :: tcc(:,:,:)             !< Cell centered tcc values.
   logical                   :: is_loaded=.false.      !< Flag for checking if the block is loaded.
   contains
      ! public methods
      procedure, pass(self) :: destroy             !< Destroy dynamic memory.
      procedure, pass(self) :: alloc               !< Allocate dynamic memory.
      procedure, pass(self) :: compute_cc          !< Compute cells centered rcc and tcc from unstructured rcc.
      procedure, pass(self) :: get_patches_extents !< Return the patches extents of given patch boundary conditions.
      procedure, pass(self) :: load_dimensions     !< Load block dimensions from file.
      procedure, pass(self) :: load_icc            !< Load block icc from file.
      procedure, pass(self) :: load_nodes          !< Load block nodes from file.
      procedure, pass(self) :: save_block_file     !< Save block data into its own file.
      procedure, pass(self) :: traslate            !< Traslate block nodes by a given traslation vector.
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
   if (allocated(self%centers)) deallocate(self%centers)
   if (allocated(self%extents)) deallocate(self%extents)
   if (allocated(self%icc)) deallocate(self%icc)
   if (allocated(self%rcc)) deallocate(self%rcc)
   if (allocated(self%tcc)) deallocate(self%tcc)
   self%is_loaded = .false.
   endsubroutine destroy

   elemental subroutine alloc(self, is_centers_to_allocate, is_extents_to_allocate)
   !< Allocate dynamic memory.
   class(block_object), intent(inout)        :: self                   !< Block data.
   logical,             intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,             intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc)
   allocate(self%nodes(1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
   if (present(is_centers_to_allocate)) then
      if (is_centers_to_allocate) allocate(self%centers(1:3,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
   endif
   if (present(is_extents_to_allocate)) then
      if (is_extents_to_allocate) then
         allocate(self%extents(1:3,1:2)) ! min-max
      endif
   endif
   allocate(self%icc(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
   allocate(self%rcc(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
   allocate(self%tcc(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
   endassociate
   endsubroutine alloc

   pure subroutine compute_cc(self, rcc)
   !< Compute cells/nodes centered rcc and tcc from unstructured rcc.
   class(block_object), intent(inout) :: self             !< Block data.
   real(R4P),           intent(in)    :: rcc(1:)          !< Unstructured rcc.
   integer(I4P)                       :: i,j,k            !< Counter.

   do k=1, self%Nk
   do j=1, self%Nj
   do i=1, self%Ni
      self%rcc(i,j,k) = 0
      if (self%icc(i,j,k) > 0)       self%rcc(i,j,k) = nint(rcc(self%icc(i,j,k)))
      if (self%rcc(i,j,k) > 28._R4P) self%rcc(i,j,k) = 0._R4P
   enddo
   enddo
   enddo
   do k=0, self%Nk+1
   do j=0, self%Nj+1
   do i=0, self%Ni+1
      self%tcc(i,j,k) = 0
      if (self%icc(i,j,k)>0) self%tcc(i,j,k) = abs(nint(rcc(self%icc(i,j,k))))
   enddo
   enddo
   enddo
   endsubroutine compute_cc

   pure subroutine get_patches_extents(self, patch, patches_extents, offset)
   !< Return the patches extents of given patch boundary conditions.
   class(block_object),       intent(in)           :: self                   !< Block data.
   integer(I4P),              intent(in)           :: patch                  !< Patch bc to be found.
   integer(I4P), allocatable, intent(inout)        :: patches_extents(:,:)   !< Patches extents, [np, 13]. The second index means
                                                                             !+ 0  => patch face (1,2,3,4,5,6);
                                                                             !+ 1  => cell i-min;
                                                                             !+ 2  => cell i-max;
                                                                             !+ 3  => cell j-min;
                                                                             !+ 4  => cell j-max;
                                                                             !+ 5  => cell k-min;
                                                                             !+ 6  => cell k-max;
                                                                             !+ 7  => node i-min;
                                                                             !+ 8  => node i-max;
                                                                             !+ 9  => node j-min;
                                                                             !+ 10 => node j-max;
                                                                             !+ 11 => node k-min;
                                                                             !+ 12 => node k-max;
   integer(I4P),              intent(in), optional :: offset                 !< Offset from patch.
   integer(I4P)                                    :: np                     !< Number of patches found.
   integer(I4P)                                    :: offset_                !< Offset from patch, local variable.
   integer(I4P), parameter                         :: ci1=1,  ci2=2,  cj1=3, &
                                                      cj2=4,  ck1=5,  ck2=6, &
                                                      ni1=7,  ni2=8,  nj1=9, &
                                                      nj2=10, nk1=11, nk2=12 !< Named indexes.

   if (allocated(patches_extents)) deallocate(patches_extents)
   if (.not.allocated(self%tcc)) return
   offset_ = 0 ; if (present(offset)) offset_ = offset
   np = 0
   associate(tcc=>self%tcc, Ni=>self%Ni, Nj=>self%Nj, Nk=>self%Nk)
      if (any(tcc(0    , :    , :   ) == patch).or.any(tcc(0     , :    , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(Ni+1 , :    , :   ) == patch).or.any(tcc(Ni+1  , :    , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:    , 0    , :   ) == patch).or.any(tcc(:     , 0    , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:    , Nj+1 , :   ) == patch).or.any(tcc(:     , Nj+1 , :   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:    , :    , 0   ) == patch).or.any(tcc(:     , :    , 0   ) == patch+10_I4P)) np = np + 1
      if (any(tcc(:    , :    , Nk+1) == patch).or.any(tcc(:     , :    , Nk+1) == patch+10_I4P)) np = np + 1
      if (np > 0) then
         allocate(patches_extents(1:np, 0:12))
         np = 0
         if (any(tcc(0    , :    , :   ) == patch).or.any(tcc(0     , :    , :   ) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 1
            patches_extents(np, ci1) = 1 + offset_                  ; patches_extents(np, ci2) = 1 + offset_
            patches_extents(np, cj1) = 1                            ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1                            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = 0 + offset_                  ; patches_extents(np, ni2) = 0 + offset_
            patches_extents(np, nj1) = patches_extents(np, cj1) - 1 ; patches_extents(np, nj2) = patches_extents(np, cj2)
            patches_extents(np, nk1) = patches_extents(np, ck1) - 1 ; patches_extents(np, nk2) = patches_extents(np, ck2)
         endif
         if (any(tcc(Ni+1 , :    , :   ) == patch).or.any(tcc(Ni+1  , :    , :   ) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 2
            patches_extents(np, ci1) = Ni - offset_                 ; patches_extents(np, ci2) = Ni - offset_
            patches_extents(np, cj1) = 1                            ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1                            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = Ni - offset_                 ; patches_extents(np, ni2) = Ni - offset_
            patches_extents(np, nj1) = patches_extents(np, cj1) - 1 ; patches_extents(np, nj2) = patches_extents(np, cj2)
            patches_extents(np, nk1) = patches_extents(np, ck1) - 1 ; patches_extents(np, nk2) = patches_extents(np, ck2)
         endif
         if (any(tcc(:    , 0    , :   ) == patch).or.any(tcc(:     , 0    , :   ) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 3
            patches_extents(np, ci1) = 1                            ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1 + offset_                  ; patches_extents(np, cj2) = 1 + offset_
            patches_extents(np, ck1) = 1                            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = patches_extents(np, ci1) - 1 ; patches_extents(np, ni2) = patches_extents(np, ci2)
            patches_extents(np, nj1) = 0 + offset_                  ; patches_extents(np, nj2) = 0 + offset_
            patches_extents(np, nk1) = patches_extents(np, ck1) - 1 ; patches_extents(np, nk2) = patches_extents(np, ck2)
         endif
         if (any(tcc(:    , Nj+1 , :   ) == patch).or.any(tcc(:     , Nj+1 , :   ) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 4
            patches_extents(np, ci1) = 1                            ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = Nj - offset_                 ; patches_extents(np, cj2) = Nj - offset_
            patches_extents(np, ck1) = 1                            ; patches_extents(np, ck2) = Nk
            patches_extents(np, ni1) = patches_extents(np, ci1) - 1 ; patches_extents(np, ni2) = patches_extents(np, ci2)
            patches_extents(np, nj1) = Nj - offset_                 ; patches_extents(np, nj2) = Nj - offset_
            patches_extents(np, nk1) = patches_extents(np, ck1) - 1 ; patches_extents(np, nk2) = patches_extents(np, ck2)
         endif
         if (any(tcc(:    , :    , 0   ) == patch).or.any(tcc(:     , :    , 0   ) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 5
            patches_extents(np, ci1) = 1               ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1                                         ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = 1 + offset_                  ; patches_extents(np, ck2) = 1 + offset_
            patches_extents(np, ni1) = patches_extents(np, ci1) - 1 ; patches_extents(np, ni2) = patches_extents(np, ci2)
            patches_extents(np, nj1) = patches_extents(np, cj1) - 1 ; patches_extents(np, nj2) = patches_extents(np, cj2)
            patches_extents(np, nk1) = 0 + offset_                  ; patches_extents(np, nk2) = 0 + offset_
         endif
         if (any(tcc(:    , :    , Nk+1) == patch).or.any(tcc(:     , :    , Nk+1) == patch+10_I4P)) then
            np = np + 1
            patches_extents(np, 0) = 6
            patches_extents(np, ci1) = 1                            ; patches_extents(np, ci2) = Ni
            patches_extents(np, cj1) = 1                            ; patches_extents(np, cj2) = Nj
            patches_extents(np, ck1) = Nk - offset_                 ; patches_extents(np, ck2) = Nk - offset_
            patches_extents(np, ni1) = patches_extents(np, ci1) - 1 ; patches_extents(np, ni2) = patches_extents(np, ci2)
            patches_extents(np, nj1) = patches_extents(np, cj1) - 1 ; patches_extents(np, nj2) = patches_extents(np, cj2)
            patches_extents(np, nk1) = Nk - offset_                 ; patches_extents(np, nk2) = Nk - offset_
         endif
      endif
   endassociate
   endsubroutine get_patches_extents

   subroutine load_dimensions(self, file_unit, is_centers_to_allocate, is_extents_to_allocate)
   !< Load block dimensions from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block dimensions record.
   class(block_object), intent(inout)        :: self                   !< Block data.
   integer(I4P),        intent(in)           :: file_unit              !< Logical unit of grd file.
   logical,             intent(in), optional :: is_centers_to_allocate !< Flag to allocate also centers array.
   logical,             intent(in), optional :: is_extents_to_allocate !< Flag to allocate also extents arrays.

   call self%destroy
   read(file_unit, end=10, err=10) self%Ni, self%Nj, self%Nk, self%gc
   10 call self%alloc(is_centers_to_allocate=is_centers_to_allocate, is_extents_to_allocate=is_extents_to_allocate)
   endsubroutine load_dimensions

   subroutine load_icc(self, file_unit)
   !< Load block icc from file.
   !<
   !< @note The icc file must be already open and the current record-index must be at the proper block icc record.
   class(block_object), intent(inout) :: self      !< Block data.
   integer(I4P),        intent(in)    :: file_unit !< Logical unit of icc file.
   integer(I4P)                       :: i,j,k     !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,icc=>self%icc)
      read(file_unit)(((icc(i,j,k),i=1-gc(1),Ni+gc(2)),j=1-gc(3),Nj+gc(4)),k=1-gc(5),Nk+gc(6))
   endassociate
   endsubroutine load_icc

   subroutine load_nodes(self, file_unit)
   !< Load block nodes from file.
   !<
   !< @note The grd file must be already open and the current record-index must be at the proper block nodes record.
   !<
   !< @note If *centers* and *extents* are allocated they are computed from nodes values.
   class(block_object), intent(inout) :: self      !< Block data.
   integer(I4P),        intent(in)    :: file_unit !< Logical unit of grd file.
   integer(I4P)                       :: i,j,k     !< Counter.

   associate(Ni => self%Ni, Nj => self%Nj, Nk => self%Nk, gc => self%gc, nodes => self%nodes)
      read(file_unit)(((nodes(1,i, j, k), i=0-gc(1),Ni+gc(2)), j=0-gc(3),Nj+gc(4)), k=0-gc(5),Nk+gc(6))
      read(file_unit)(((nodes(2,i, j, k), i=0-gc(1),Ni+gc(2)), j=0-gc(3),Nj+gc(4)), k=0-gc(5),Nk+gc(6))
      read(file_unit)(((nodes(3,i, j, k), i=0-gc(1),Ni+gc(2)), j=0-gc(3),Nj+gc(4)), k=0-gc(5),Nk+gc(6))
   endassociate
   if (allocated(self%centers)) then
      do k=1-self%gc(5),self%Nk+self%gc(6)
      do j=1-self%gc(3),self%Nj+self%gc(4)
      do i=1-self%gc(1),self%Ni+self%gc(2)
         self%centers(1:3,i,j,k) = (self%nodes(1:3,i,   j,   k  ) + &
                                    self%nodes(1:3,i-1, j,   k  ) + &
                                    self%nodes(1:3,i  , j-1, k  ) + &
                                    self%nodes(1:3,i  , j  , k-1) + &
                                    self%nodes(1:3,i-1, j-1, k-1) + &
                                    self%nodes(1:3,i  , j-1, k-1) + &
                                    self%nodes(1:3,i-1, j  , k-1) + &
                                    self%nodes(1:3,i-1, j-1, k  )) * 0.125_R8P
      enddo
      enddo
      enddo
   endif
   if (allocated(self%extents)) self%extents = compute_extents(i_extents=[0,self%Ni], &
                                                               j_extents=[0,self%Nj], &
                                                               k_extents=[0,self%Nk])
   contains
      function compute_extents(i_extents, j_extents, k_extents) result(extents)
      !< Compute (sub)block extents provided indexes extents.
      integer(I4P), intent(in) :: i_extents(2) !< Index "i" extents (min, max) of the (sub)block.
      integer(I4P), intent(in) :: j_extents(2) !< Index "j" extents (min, max) of the (sub)block.
      integer(I4P), intent(in) :: k_extents(2) !< Index "k" extents (min, max) of the (sub)block.
      real(R8P)                :: extents(3,2) !< Block extents.

      extents(1,1) = minval(self%nodes(1,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))
      extents(1,2) = maxval(self%nodes(1,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))

      extents(2,1) = minval(self%nodes(2,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))
      extents(2,2) = maxval(self%nodes(2,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))

      extents(3,1) = minval(self%nodes(3,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))
      extents(3,2) = maxval(self%nodes(3,i_extents(1):i_extents(2), j_extents(1):j_extents(2), k_extents(1):k_extents(2)))
      endfunction compute_extents
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
   write(file_unit)(((nodes(1,i,j,k),i=0-gc(1),Ni+gc(2)),j=0-gc(3),Nj+gc(4)),k=0-gc(5),Nk+gc(6))
   write(file_unit)(((nodes(2,i,j,k),i=0-gc(1),Ni+gc(2)),j=0-gc(3),Nj+gc(4)),k=0-gc(5),Nk+gc(6))
   write(file_unit)(((nodes(3,i,j,k),i=0-gc(1),Ni+gc(2)),j=0-gc(3),Nj+gc(4)),k=0-gc(5),Nk+gc(6))
   close(file_unit)
   ! file rcc
   open(newunit=file_unit, file='block-rcc-'//trim(adjustl(bstr))//'.blk', form='unformatted', action='write')
   write(file_unit) self%Ni, self%Nj, self%Nk, self%gc
   do k=1-gc(1), Nk+gc(2)
   do j=1-gc(3), Nj+gc(4)
   do i=1-gc(5), Ni+gc(6)
      p = icc(i,j,k)
      if (p<0) then ! natural BC
         select case(p)
         case(-1 )
            write(file_unit) 'WALL'
         case(-2 )
            write(file_unit) 'SIMMETRY'
         case(-3 )
            write(file_unit) 'INFLOW'
         case(-4 )
            write(file_unit) 'INOUTFLOW'
         case(-5 )
            write(file_unit) 'ASSIGNED-INFLOW'
         case(-6 )
            write(file_unit) 'ASSIGNED-PRESSURE'
         case(-7 )
            write(file_unit) 'ASSIGNED-NORMAL-VELOCITY'
         case(-8 )
            write(file_unit) 'ASSIGNED-RIEMANN'
         case(-9 )
            write(file_unit) 'EXTRAPOLATED'
         case(-10)
            write(file_unit) 'MOVING-WALL'
         case(-11)
            write(file_unit) 'INACTIVE-WALL'
         case(-19)
            write(file_unit) 'EXTRAPOLATED-ALT'
         endselect
      elseif (p==0) then
         write(file_unit) 'ACTIVE-CELL'
      else ! chimera-like BC (chimera, adjacent...)
         write(file_unit) nint(rcc(p),I4P) ! chimera type
         write(file_unit) nint(rcc(p+1),I4P) ! donors number
         do n=1, nint(rcc(p+1),I4P)
            offset = p + 1 + 5*(n-1)
            write(file_unit) nint(rcc(offset+1)), nint(rcc(offset+2)), nint(rcc(offset+3)), nint(rcc(offset+4)) ! b,i,j,k donor
            write(file_unit) rcc(offset+5) ! donor weight
         enddo
      endif
   enddo
   enddo
   enddo
   endassociate
   endsubroutine save_block_file

   pure subroutine traslate(self, traslation)
   !< Traslate block nodes by a given traslation vector.
   class(block_object), intent(inout) :: self          !< Block data.
   real(R8P),           intent(in)    :: traslation(3) !< Traslation vector.
   integer(I4P)                       :: i,j,k         !< Counter.

   associate(Ni=>self%Ni,Nj=>self%Nj,Nk=>self%Nk,gc=>self%gc,nodes=>self%nodes)
   if (allocated(self%nodes)) then
      do k=1-self%gc(5),self%Nk+self%gc(6)
      do j=1-self%gc(3),self%Nj+self%gc(4)
      do i=1-self%gc(1),self%Ni+self%gc(2)
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
      call blocks(b)%load_dimensions(file_unit=file_unit_grd,      &
                                     is_centers_to_allocate=.true.,&
                                     is_extents_to_allocate=.true.)
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
      call blocks(b)%load_dimensions(file_unit=file_unit_icc)
   enddo
   do b=1, blocks_number
      call blocks(b)%load_icc(file_unit=file_unit_icc)
   enddo
   read(file_unit_icc) unstruct_dimension
   allocate(rcc(1:unstruct_dimension))
   read(file_unit_icc) (rcc(i),i=1,unstruct_dimension)
   close(file_unit_icc)
   do b=1, blocks_number
      call blocks(b)%compute_cc(rcc=rcc)
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
