! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! Build with ENABLE_SANDBOX=ON. Run a single process, with threads pinned to cores:
! ATLAS_TRACE=0 ATLAS_TRACE_BARRIERS=0 ATLAS_TRACE_REPORT=0 \
! OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close \
! atlas-benchmark-trace-fortran [iterations-per-thread=10000] [max-threads=OMP_NUM_THREADS]
!
! Weak-scaling stress test: every thread repeatedly instruments almost no work.
! All cases include the same volatile increment. Results are medians of three
! samples, excluding thread creation and persistent-trace setup/finalization.
! no-trace makes no trace calls in the timed loop, representing instrumentation
! preprocessed out. net-ns/op subtracts this baseline at the same team size.
! retain/release isolates registry locking/lookups without per-iteration allocation.
! It does not isolate mutex cost from hash-table cost. Only thread zero records
! timings, but all threads exercise Fortran ownership. ATLAS_TRACE=0 mutes logging,
! not timer collection. Use the same build type and binding when comparing builds.

program atlas_benchmark_trace_fortran
  use iso_c_binding, only: c_ptr
  use iso_fortran_env, only: int64, real64
  use omp_lib
  use atlas_module, only: atlas_library, atlas_trace
  use atlas_Trace_c_binding, only: retain_atlas_Trace, release_atlas_Trace
  implicit none
  integer, parameter :: baseline=1, initialize=2, construct=3, copy_handle=4, ownership=5, reuse=6
  character(len=17), parameter :: names(6) = [character(len=17) :: &
    'no-trace', 'init/final', 'constructor/final', 'copy/final', 'retain/release', 'resume/pause']
  integer :: iterations, max_threads, threads, mode, sample
  real(real64) :: samples(3), elapsed, baseline_time, operations, serial_time(6)

  iterations = 10000
  max_threads = omp_get_max_threads()
  if (command_argument_count() > 2) error stop 'usage: atlas-benchmark-trace-fortran [iterations] [max-threads]'
  call read_positive_argument(1, iterations)
  call read_positive_argument(2, max_threads)
  call require_disabled('ATLAS_TRACE')
  call require_disabled('ATLAS_TRACE_BARRIERS')
  call require_disabled('ATLAS_TRACE_REPORT')
  call omp_set_dynamic(.false.)
  call atlas_library%initialise()

  print '(a,i0)', '# iterations per thread: ', iterations
  print '(a)', '# median of 3 samples; ns/op uses total operations across the team'
  print '(a)', '# thread-ns/op = wall time / iterations; ideal throughput scaling has slowdown = 1'
  print '(a)', 'threads,case,seconds,ns/op,net-ns/op,thread-ns/op,scaling-slowdown,vs-no-trace'
  threads = 1
  do
    operations = real(iterations, real64) * real(threads, real64)
    do mode=baseline,reuse
      call measure(mode, min(iterations, 100), threads, elapsed)
      do sample=1,size(samples)
        call measure(mode, iterations, threads, samples(sample))
      enddo
      elapsed = sum(samples) - minval(samples) - maxval(samples)
      if (mode == baseline) baseline_time = elapsed
      if (threads == 1) serial_time(mode) = elapsed
      print '(i0,a,a,6(a,es13.6))', threads, ',', trim(names(mode)), ',', elapsed, &
        ',', elapsed * 1.e9_real64 / operations, ',', (elapsed-baseline_time) * 1.e9_real64 / operations, &
        ',', elapsed * 1.e9_real64 / real(iterations,real64), ',', elapsed / serial_time(mode), &
        ',', elapsed / baseline_time
    enddo
    if (threads == max_threads) exit
    threads = threads + min(threads, max_threads-threads)
  enddo
  call atlas_library%finalise()

contains

  subroutine read_positive_argument(index, value)
    integer, intent(in) :: index
    integer, intent(inout) :: value
    character(len=64) :: argument
    integer :: status
    if (command_argument_count() < index) return
    call get_command_argument(index, argument, status=status)
    if (status /= 0) error stop 'cannot read benchmark argument'
    read(argument, *, iostat=status) value
    if (status /= 0) error stop 'benchmark arguments must be positive integers'
    if (value <= 0) error stop 'benchmark arguments must be positive integers'
  end subroutine

  subroutine require_disabled(name)
    character(len=*), intent(in) :: name
    character(len=16) :: value
    integer :: status
    call get_environment_variable(name, value, status=status)
    if (status == 1) return
    if (status /= 0 .or. trim(value) /= '0') then
      print '(a)', 'Set '//name//'=0 to avoid logging, reporting or barrier overhead'
      error stop 'invalid benchmark environment'
    endif
  end subroutine

  subroutine measure(mode, iterations, threads, seconds)
    integer, intent(in) :: mode, iterations, threads
    real(real64), intent(out) :: seconds
    integer(int64) :: completed
    real(real64) :: started
    completed = 0_int64
!$omp parallel num_threads(threads) reduction(+:completed)
    block
      type(atlas_trace) :: trace, copy
      type(c_ptr) :: ptr
      integer :: i
      integer(int64), volatile :: work
      if (omp_get_num_threads() /= threads) error stop 'OpenMP did not supply the requested team size'
      work = 0_int64
      if (mode >= copy_handle) then
        call trace%init('trace-benchmark', 1, 'persistent trace')
        call trace%pause()
        ptr = trace%c_ptr()
      endif
!$omp barrier
!$omp master
      started = omp_get_wtime()
!$omp end master
!$omp barrier
      select case(mode)
      case(baseline)
        do i=1,iterations
          work = work + 1_int64
        enddo
      case(initialize)
        do i=1,iterations
          call trace%init('trace-benchmark', 2, 'short trace')
          work = work + 1_int64
          call trace%final()
        enddo
      case(construct)
        do i=1,iterations
          trace = atlas_trace('trace-benchmark', 3, 'short trace')
          work = work + 1_int64
          call trace%final()
        enddo
      case(copy_handle)
        do i=1,iterations
          copy = trace
          work = work + 1_int64
          call copy%final()
        enddo
      case(ownership)
        do i=1,iterations
          call retain_atlas_Trace(ptr)
          work = work + 1_int64
          call release_atlas_Trace(ptr)
        enddo
      case(reuse)
        do i=1,iterations
          call trace%resume()
          work = work + 1_int64
          call trace%pause()
        enddo
      end select
!$omp barrier
!$omp master
      seconds = omp_get_wtime() - started
!$omp end master
!$omp barrier
      if (mode >= copy_handle) call trace%final()
      completed = completed + work
    end block
!$omp end parallel
    if (completed /= int(iterations,int64) * int(threads,int64)) error stop 'incorrect benchmark operation count'
    if (seconds <= 0._real64) error stop 'benchmark timer did not advance; increase iterations'
  end subroutine
end program atlas_benchmark_trace_fortran
