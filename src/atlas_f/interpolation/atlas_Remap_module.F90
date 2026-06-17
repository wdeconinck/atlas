! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Remap_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_char, c_null_ptr, c_associated, c_size_t, c_int, c_loc

implicit none

public :: atlas_Remap

private

!------------------------------------------------------------------------------
type atlas_Remap

! Purpose :
! -------
!   *atlas_Remap* : Cache and execute interpolations between blocked fields

! Methods :
! -------
!   setup_interpolation : create and cache an interpolation
!   interpolate         : interpolate a field or fieldset
!   clear               : remove cached interpolations

!------------------------------------------------------------------------------
  type(c_ptr) :: c_ptr = c_null_ptr

contains
  procedure, public  :: setup_interpolation
  procedure, public  :: clear
  procedure, private :: interpolate_field
  procedure, private :: interpolate_fieldset
  generic, public    :: interpolate => interpolate_field, interpolate_fieldset
end type atlas_Remap

interface
  type(c_ptr) function atlas__Remap__new() bind(C, name="atlas__Remap__new")
    use, intrinsic :: iso_c_binding, only : c_ptr
  end function atlas__Remap__new

  subroutine atlas__Remap__delete(this) bind(C, name="atlas__Remap__delete")
    use, intrinsic :: iso_c_binding, only : c_ptr
    type(c_ptr), value :: this
  end subroutine atlas__Remap__delete

  subroutine atlas__Remap__setup_interpolation(this, method, method_len, source, target) &
      bind(C, name="atlas__Remap__setup_interpolation")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
    type(c_ptr), value :: this
    type(c_ptr), value :: method
    integer(c_size_t), value :: method_len
    type(c_ptr), value :: source
    type(c_ptr), value :: target
  end subroutine atlas__Remap__setup_interpolation

  subroutine atlas__Remap__interpolate_field(this, method, method_len, source, target, on_device) &
      bind(C, name="atlas__Remap__interpolate_field")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_int
    type(c_ptr), value :: this
    type(c_ptr), value :: method
    integer(c_size_t), value :: method_len
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__Remap__interpolate_field

  subroutine atlas__Remap__interpolate_fieldset(this, method, method_len, source, target, on_device) &
      bind(C, name="atlas__Remap__interpolate_fieldset")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_int
    type(c_ptr), value :: this
    type(c_ptr), value :: method
    integer(c_size_t), value :: method_len
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__Remap__interpolate_fieldset
end interface

!========================================================
contains
!========================================================

subroutine setup_interpolation(this, method, source, target)
  use atlas_FunctionSpace_module, only : atlas_FunctionSpace
  class(atlas_Remap), intent(inout) :: this
  character(len=*,kind=c_char), target, intent(in) :: method
  class(atlas_FunctionSpace), intent(in) :: source
  class(atlas_FunctionSpace), intent(in) :: target

  if (.not. c_associated(this%c_ptr)) then
    this%c_ptr = atlas__Remap__new()
  endif

  call atlas__Remap__setup_interpolation(this%c_ptr, c_loc(method), int(len_trim(method), c_size_t), &
      source%c_ptr(), target%c_ptr())
end subroutine setup_interpolation

subroutine clear(this)
  class(atlas_Remap), intent(inout) :: this

  if (c_associated(this%c_ptr)) then
    call atlas__Remap__delete(this%c_ptr)
    this%c_ptr = c_null_ptr
  endif
end subroutine clear

subroutine interpolate_field(this, method, source, target, on_device)
  use atlas_Field_module, only : atlas_Field
  class(atlas_Remap), intent(inout) :: this
  character(len=*,kind=c_char), target, intent(in) :: method
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  if (.not. c_associated(this%c_ptr)) then
    this%c_ptr = atlas__Remap__new()
  endif

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__Remap__interpolate_field(this%c_ptr, c_loc(method), int(len_trim(method), c_size_t), &
      source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine interpolate_field

subroutine interpolate_fieldset(this, method, source, target, on_device)
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_Remap), intent(inout) :: this
  character(len=*,kind=c_char), target, intent(in) :: method
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  if (.not. c_associated(this%c_ptr)) then
    this%c_ptr = atlas__Remap__new()
  endif

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__Remap__interpolate_fieldset(this%c_ptr, c_loc(method), int(len_trim(method), c_size_t), &
      source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine interpolate_fieldset

end module atlas_Remap_module
