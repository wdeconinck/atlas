! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

program test_trace_openmp_private_array
  use atlas_module, only: atlas_library, atlas_trace
  use omp_lib, only: omp_get_thread_num
  implicit none
#if !defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  type(atlas_trace) :: traces(3)
#endif
  integer :: repetition, completed

  call atlas_library%initialise()
  completed=0
  do repetition=1,32
#if defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
    ! Intel Classic 2021.10 does not initialize the PRIVATE array reliably.
!$omp parallel num_threads(2) reduction(+:completed)
    call run_with_local_array(completed)
#else
!$omp parallel num_threads(2) private(traces) reduction(+:completed)
    call run_intervals(traces, completed)
#endif
!$omp end parallel
  enddo
  if (completed /= 32*2*8*2) error stop 'unexpected number of trace intervals'
  call atlas_library%finalise()

contains

#if defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
  recursive subroutine run_with_local_array(completed)
    integer, intent(inout) :: completed
    type(atlas_trace) :: traces(3)

    call run_intervals(traces, completed)
  end subroutine
#endif

  recursive subroutine run_intervals(traces, completed)
    type(atlas_trace), intent(inout) :: traces(3)
    integer, intent(inout) :: completed
    integer :: trace_index, interval

    do trace_index=1,size(traces)
      if (trace_index == 3) cycle
      traces(trace_index)=atlas_trace('test_trace_openmp_private_array.F90',22,'private array interval')
      call traces(trace_index)%pause()
    enddo
    do interval=1,8
      do trace_index=1,2
        call traces(trace_index)%resume()
        if (traces(trace_index)%running() .neqv. (omp_get_thread_num() == 0)) error stop 'unexpected trace state'
        completed=completed+1
        call traces(trace_index)%pause()
      enddo
    enddo
    do trace_index=2,1,-1
      call traces(trace_index)%final()
    enddo
  end subroutine
end program test_trace_openmp_private_array
