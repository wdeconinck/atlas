! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Trace_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated

implicit none

public :: atlas_Trace


private

!-----------------------------
! atlas_Trace                !
!-----------------------------

! With GNU Fortran 10, assign constructor results to a scalar before broadcasting to an array.
type :: atlas_Trace
  ! Allocatable storage is initialized reliably for Intel OpenMP PRIVATE arrays.
  type(c_ptr), allocatable, private :: ptr
  ! Constructor results transfer their initial reference on assignment, as in FCKIT.
  logical, private :: return_value = .false.
contains
  procedure, public :: init => atlas_Trace__init
  procedure, public :: final => atlas_Trace__final
  procedure, public :: c_ptr => atlas_Trace__c_ptr
  procedure, public :: is_null => atlas_Trace__is_null
! Public methods

  procedure, public :: running
  procedure, public :: start
  procedure, public :: stop
  procedure, public :: pause
  procedure, public :: resume
  procedure, public :: elapsed

  procedure, private :: assignment_operator => atlas_Trace__assign
  generic, public :: assignment(=) => assignment_operator
#if FCKIT_HAVE_FINAL
  final :: atlas_Trace__final_auto
#endif
end type

interface atlas_Trace
  module procedure atlas_Trace__loc
  module procedure atlas_Trace__labels_1
  module procedure atlas_Trace__labels_2
  module procedure atlas_Trace__labels_3
  module procedure atlas_Trace__labels_4
  module procedure atlas_Trace__labels_5
end interface

!========================================================
contains
!========================================================

function atlas_Trace__loc(file,line,title) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  call this%init(file, line, title)
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function

function atlas_Trace__labels_1(file,line,title,label) result(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label
  this%ptr = new_atlas_Trace_labels_view_1(file, len_trim(file, kind=c_size_t), line, &
    & title, len_trim(title, kind=c_size_t), label, len_trim(label, kind=c_size_t))
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function
function atlas_Trace__labels_2(file,line,title,label1,label2) result(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  this%ptr = new_atlas_Trace_labels_view_2(file, len_trim(file, kind=c_size_t), line, &
    & title, len_trim(title, kind=c_size_t), label1, len_trim(label1, kind=c_size_t), &
    & label2, len_trim(label2, kind=c_size_t))
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function
function atlas_Trace__labels_3(file,line,title,label1,label2,label3) result(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  this%ptr = new_atlas_Trace_labels_view_3(file, len_trim(file, kind=c_size_t), line, &
    & title, len_trim(title, kind=c_size_t), label1, len_trim(label1, kind=c_size_t), &
    & label2, len_trim(label2, kind=c_size_t), label3, len_trim(label3, kind=c_size_t))
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function
function atlas_Trace__labels_4(file,line,title,label1,label2,label3,label4) result(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  character(len=*) , intent(in) :: label4
  this%ptr = new_atlas_Trace_labels_view_4(file, len_trim(file, kind=c_size_t), line, &
    & title, len_trim(title, kind=c_size_t), label1, len_trim(label1, kind=c_size_t), &
    & label2, len_trim(label2, kind=c_size_t), label3, len_trim(label3, kind=c_size_t), &
    & label4, len_trim(label4, kind=c_size_t))
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function
function atlas_Trace__labels_5(file,line,title,label1,label2,label3,label4,label5) result(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  character(len=*) , intent(in) :: label4
  character(len=*) , intent(in) :: label5
  this%ptr = new_atlas_Trace_labels_view_5(file, len_trim(file, kind=c_size_t), line, &
    & title, len_trim(title, kind=c_size_t), label1, len_trim(label1, kind=c_size_t), &
    & label2, len_trim(label2, kind=c_size_t), label3, len_trim(label3, kind=c_size_t), &
    & label4, len_trim(label4, kind=c_size_t), label5, len_trim(label5, kind=c_size_t))
  call return_atlas_Trace(this%ptr)
  this%return_value = .true.
end function
!-------------------------------------------------------------------------------

subroutine atlas_Trace__init(this, file, line, title)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_Trace_c_binding
  use fckit_exception_module, only : fckit_exception
  class(atlas_Trace), intent(inout) :: this
  character(len=*), intent(in) :: file, title
  integer, intent(in) :: line
  if (allocated(this%ptr)) then
    call fckit_exception%abort("atlas_Trace is already initialized; call %final() before %init()", file, line)
  endif
  this%ptr = new_atlas_Trace_view(file, len_trim(file, kind=c_size_t), line, title, &
    & len_trim(title, kind=c_size_t))
  this%return_value = .false.
end subroutine

impure elemental subroutine atlas_Trace__final(this)
  use atlas_Trace_c_binding
  class(atlas_Trace), intent(inout) :: this
  if (allocated(this%ptr)) then
    call release_atlas_Trace(this%ptr)
    deallocate(this%ptr)
  endif
  this%return_value = .false.
end subroutine

impure elemental subroutine atlas_Trace__assign(lhs, rhs)
  use atlas_Trace_c_binding
  class(atlas_Trace), intent(inout) :: lhs
#if defined(__NVCOMPILER_MAJOR__) && \
    (__NVCOMPILER_MAJOR__ < 25 || (__NVCOMPILER_MAJOR__ == 25 && __NVCOMPILER_MINOR__ < 9))
  ! NVHPC before 25.9 can bypass host-associated assignment unless RHS is polymorphic.
  ! Required with tested versions 22.11, 23.7, 24.1, 24.5, 24.7, 24.9, 25.1, 25.3, 25.5 and 25.7.
  class(atlas_Trace), intent(in) :: rhs
#else
  type(atlas_Trace), intent(in) :: rhs
#endif
  if (c_associated(lhs%c_ptr(), rhs%c_ptr())) return
  if (allocated(rhs%ptr)) call retain_atlas_Trace(rhs%ptr)
  call lhs%final()
  if (allocated(rhs%ptr)) lhs%ptr = rhs%ptr
end subroutine

function atlas_Trace__c_ptr(this) result(ptr)
  class(atlas_Trace), intent(in) :: this
  type(c_ptr) :: ptr
  ptr = c_null_ptr
  if (allocated(this%ptr)) ptr = this%ptr
end function

logical function atlas_Trace__is_null(this)
  class(atlas_Trace), intent(in) :: this
  atlas_Trace__is_null = .not. allocated(this%ptr)
end function

!-------------------------------------------------------------------------------

function running( this )
  use atlas_Trace_c_binding
  logical :: running
  class(atlas_Trace) :: this
  if( atlas_Trace__running( this%c_ptr() ) == 0 ) then
    running = .False.
  else
    running = .True.
  endif
end function

!-------------------------------------------------------------------------------

function elapsed( this )
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_Trace_c_binding
  real(c_double) :: elapsed
  class(atlas_Trace) :: this
  elapsed = atlas_Trace__elapsed( this%c_ptr() )
end function

!-------------------------------------------------------------------------------

subroutine start( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__start( this%c_ptr() )
end subroutine

!-------------------------------------------------------------------------------

subroutine stop( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__stop( this%c_ptr() )
end subroutine

!-------------------------------------------------------------------------------

subroutine pause( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__pause( this%c_ptr() )
end subroutine

!-------------------------------------------------------------------------------

subroutine resume( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__resume( this%c_ptr() )
end subroutine

!-------------------------------------------------------------------------------

#if FCKIT_HAVE_FINAL
ATLAS_FINAL subroutine atlas_Trace__final_auto(this)
  type(atlas_Trace), intent(inout) :: this
  if (.not. this%return_value) call this%final()
end subroutine
#endif

end module atlas_Trace_module

