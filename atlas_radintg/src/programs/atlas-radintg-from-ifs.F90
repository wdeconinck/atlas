program main

use type_geometry
use yomprad
use ATLAS_REMAP_MODULE
use ATLAS_RELAYOUT_MODULE
use atlas_module
use pluto_module
implicit none

type(GEOMETRY) :: ydgeometry
type(RADIATION_GRID_STRUCT) :: radgrid

type(atlas_functionspace_BlockStructuredColumns) :: ifs_blocked_fs
type(atlas_functionspace_BlockStructuredColumns) :: rad_blocked_fs

type(atlas_Remap) :: remap

type(atlas_FieldSet) :: ifs_fields, rad_fields

real(jprb), pointer :: ifs_field_data(:,:,:)
real(jprb), pointer :: rad_field_data(:,:,:)

real(jprb), allocatable :: ZRGP_IFS(:,:,:)
real(jprb), allocatable :: ZRGP_RAD(:,:,:)

type(atlas_Trace) :: trace

integer :: nlev, nproma_rad, nproma_ifs

! --- Program initialization

call atlas_init()

trace = atlas_Trace("atlas-radintg-from-ifs.F90",__LINE__,"main")

! --- Setup IFS data structures

! We assume ydgeometry and radgrid are already set up in the IFS
call setup_geometry(ydgeometry,"O32")
call print_geometry(ydgeometry, 0)
call setup_radgrid(radgrid, "O16", lmatch_ifs=.true.)
! We should extract following still from IFS data structures
nlev = 20       ! grid dependent, but we can set it to a fixed value for this test
nproma_ifs = 16 ! CPU: 16, GPU: 32
nproma_rad = 32 ! CPU: 64??, GPU: 8192

! --- Setup Atlas function spaces based only on IFS data structures

call setup_ifs_blocked_fs(ifs_blocked_fs, ydgeometry, nproma=nproma_ifs)
call setup_rad_blocked_fs(rad_blocked_fs, radgrid, nproma=nproma_rad)

write(0,*) "Finished setting up function spaces"
write(0,*) "IFS blocked function space:"
write(0,*) "  nblks: ", ifs_blocked_fs%nblks()
write(0,*) "  nproma: ", ifs_blocked_fs%nproma()
write(0,*) "  nblks*nproma: ", ifs_blocked_fs%size()
write(0,*) "RAD blocked function space:"
write(0,*) "  nblks: ", rad_blocked_fs%nblks()
write(0,*) "  nproma: ", rad_blocked_fs%nproma()
write(0,*) "  nblks*nproma: ", rad_blocked_fs%size()

!--- Setup interpolations.

! The setup_interpolation may be expensive. It could be part of suecrad; only setup the needed ones.
CALL remap%setup_interpolation("bilinear", ifs_blocked_fs, rad_blocked_fs)
CALL remap%setup_interpolation("bicubic",  ifs_blocked_fs, rad_blocked_fs)
CALL remap%setup_interpolation("cons",     ifs_blocked_fs, rad_blocked_fs)

!--- Create fieldsets

ifs_fields = atlas_FieldSet("ifs_fields")
rad_fields = atlas_FieldSet("rad_fields")

! Add Atlas-managed fields
call ifs_fields%add(ifs_blocked_fs%create_field(name="ifs_field", KIND=ATLAS_REAL(JPRB), levels=nlev))
call rad_fields%add(rad_blocked_fs%create_field(name="rad_field", KIND=ATLAS_REAL(JPRB), levels=nlev))

! And/OR wrap existing data arrays in Atlas fields
ALLOCATE(ZRGP_IFS(ifs_blocked_fs%nproma(), nlev, ifs_blocked_fs%nblks()))
ALLOCATE(ZRGP_RAD(rad_blocked_fs%nproma(), nlev, rad_blocked_fs%nblks()))
block
  type(atlas_Field) :: zrgp_ifs_field, zrgp_rad_field
  zrgp_ifs_field = atlas_Field(name="zrgp_ifs_field", data=ZRGP_IFS)
  zrgp_rad_field = atlas_Field(name="zrgp_rad_field", data=ZRGP_RAD)
  call zrgp_ifs_field%set_functionspace(ifs_blocked_fs)
  call zrgp_rad_field%set_functionspace(rad_blocked_fs)
  call zrgp_ifs_field%set_levels(nlev)
  call zrgp_rad_field%set_levels(nlev)
  call ifs_fields%add(zrgp_ifs_field)
  call rad_fields%add(zrgp_rad_field)
  call zrgp_ifs_field%final()
  call zrgp_rad_field%final()
end block

!--- Initialize ifs_fields

call ifs_fields%data("ifs_field",ifs_field_data)
call rad_fields%data("rad_field",rad_field_data)

ifs_field_data(:,:,:) = 1.0_jprb

write(0,'(3A,F5.3,A,F5.3)') "ifs_fields checksum = ", ifs_blocked_fs%checksum(ifs_fields), &
  & ",    max = ", MAXVAL(ifs_field_data(:,:,1:ifs_blocked_fs%nblks()-1)), ",    min = ", MINVAL(ifs_field_data(:,:,1:ifs_blocked_fs%nblks()-1))


! --- Now we have ifs_fields and rad_fields fieldsets allocated and filled in.

! 1. Set the default host allocator to a memory pool for intermediate temporary fields
call pluto%scope%push()
call pluto%host%set_default_resource("pluto::host_pool_resource")

! 2. Test the interpolations, and print checksums and min/max values for the results. Note that the checksums, min and max values are not expected to be the same as in the IFS, as the interpolation is not conservative and we have different blocking and halo sizes, but they should be reasonable for the given input (e.g. all positive, and not orders of magnitude different from the input).

if (interpolation_is_copy(ifs_blocked_fs, rad_blocked_fs)) then
  ! In case the grids of the IFS and RAD function spaces are the same, we can also test the copy blocked to blocked without interpolation, which should be exact.
  CALL COPY_BLOCKED_TO_BLOCKED(ifs_fields, rad_fields, on_device=.false.)
  write(0,'(3A,F5.3,A,F5.3)') "After copy     interpolation, rad_fields checksum = ", rad_blocked_fs%checksum(rad_fields), &
    & ",    max = ", MAXVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1)), ",    min = ", MINVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1))
else
  CALL remap%interpolate("bilinear", ifs_fields, rad_fields, on_device=.false.)
  write(0,'(3A,F5.3,A,F5.3)') "After bilinear interpolation, rad_fields checksum = ", rad_blocked_fs%checksum(rad_fields), &
    & ",    max = ", MAXVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1)), ",    min = ", MINVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1))

  CALL remap%interpolate("bicubic",  ifs_fields, rad_fields, on_device=.false.)
  write(0,'(3A,F5.3,A,F5.3)') "After bicubic  interpolation, rad_fields checksum = ", rad_blocked_fs%checksum(rad_fields), &
    & ",    max = ", MAXVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1)), ",    min = ", MINVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1))

  CALL remap%interpolate("cons",  ifs_fields, rad_fields, on_device=.false.)
  write(0,'(3A,F5.3,A,F5.3)') "After cons     interpolation, rad_fields checksum = ", rad_blocked_fs%checksum(rad_fields), &
    & ",    max = ", MAXVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1)), ",    min = ", MINVAL(rad_field_data(:,:,1:rad_blocked_fs%nblks()-1))
endif

! 3. Restore the default host allocator
call pluto%scope%pop()

! --- Finalize program
call ifs_fields%final()
call rad_fields%final()

call remap%clear()
call trace%final()
call atlas_final()

contains

function interpolation_is_copy(ifs_blocked_fs, rad_blocked_fs) result(is_copy)
    logical :: is_copy
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: ifs_blocked_fs
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: rad_blocked_fs
    type(atlas_Grid) :: ifsgrid_atlas, radgrid_atlas
    type(atlas_Grid) :: ifs_grid, rad_grid
    ifs_grid = ifs_blocked_fs%grid()
    rad_grid = rad_blocked_fs%grid()
    is_copy = (ifs_grid%uid() == rad_grid%uid())
    call ifs_grid%final()
    call rad_grid%final()
end function

subroutine setup_geometry(ydgeometry, gridname)
    use atlas_module
    use fckit_mpi_module
    implicit none
    type(GEOMETRY), intent(inout) :: ydgeometry
    character(len=*), intent(in) :: gridname
    type(atlas_StructuredGrid) :: grid
    type(atlas_Partitioner) :: partitioner
    type(atlas_GridDistribution) :: distribution
    integer(ATLAS_KIND_IDX), allocatable :: nb_pts_per_partition(:)
    integer :: jp

    grid = atlas_grid(gridname)
    partitioner = atlas_partitioner("equal_regions")
    distribution = partitioner%partition(grid)
    nb_pts_per_partition = distribution%nb_pts()

    ! TODO make this more straightforward, though it is explicit as is
    allocate(YDGEOMETRY%YRMP%NGLOBALPROC(grid%size()))
    do jp = 1, grid%size()
        YDGEOMETRY%YRMP%NGLOBALPROC(jp) = distribution%partition(jp) + 1
    end do

    ! allocate(YDGEOMETRY%YRGEM%NLOENG(grid%ny())) ! no need with direct assignment on next line
    YDGEOMETRY%YRGEM%NLOENG  = grid%nx_array()
    YDGEOMETRY%YRDIM%NDGLG   = grid%ny()
    YDGEOMETRY%YRGEM%NGPTOTG = grid%size()
    YDGEOMETRY%YRGEM%NGPTOT  = nb_pts_per_partition(fckit_mpi%rank() + 1)
    call distribution%final
    call partitioner%final
    call grid%final
end subroutine

subroutine print_geometry(ydgeometry,unit)
    implicit none
    type(GEOMETRY), intent(in) :: ydgeometry
    integer, intent(in) :: unit
    write(unit,*) "YDGEOMETRY%YRGEM%NGPTOTG     : ", YDGEOMETRY%YRGEM%NGPTOTG
    write(unit,*) "YDGEOMETRY%YRDIM%NDGLG       : ", YDGEOMETRY%YRDIM%NDGLG
    write(unit,*) "YDGEOMETRY%YRGEM%NLOENG      : ", YDGEOMETRY%YRGEM%NLOENG
    write(unit,*) "YDGEOMETRY%YRGEM%NGPTOT      : ", YDGEOMETRY%YRGEM%NGPTOT
    write(unit,*) "YDGEOMETRY%YRMP%NGLOBALPROC  : [ ", YDGEOMETRY%YRMP%NGLOBALPROC(1:MIN(3,YDGEOMETRY%YRGEM%NGPTOTG)), &
                         & " ... ", YDGEOMETRY%YRMP%NGLOBALPROC(MAX(1,YDGEOMETRY%YRGEM%NGPTOTG-2):YDGEOMETRY%YRGEM%NGPTOTG), " ]"
end subroutine

subroutine setup_radgrid(radgrid, gridname, lmatch_ifs)
    implicit none
    type(RADIATION_GRID_STRUCT), intent(inout) :: radgrid
    character(len=*), intent(in) :: gridname
    logical, intent(in) :: lmatch_ifs
    type(atlas_StructuredGrid) :: grid
    type(atlas_Partitioner) :: partitioner
    type(atlas_GridDistribution) :: distribution
    type(atlas_functionspace_StructuredColumns) :: fs
    type(atlas_Field) :: lonlat_deg
    REAL(KIND=JPRD), pointer :: lonlat_deg_array(:,:)
    REAL(KIND=JPRD), PARAMETER :: RPI = 2.0_JPRD * ASIN(1.0_JPRD)
    REAL(KIND=JPRD), PARAMETER :: DEG_TO_RAD = RPI / 180.0_JPRD

    integer :: jp

    grid = atlas_grid(gridname)
    call setup_distribution(distribution, grid, ydgeometry, lmatch_ifs)
    allocate(radgrid%NRGRI(grid%ny()))
    radgrid%NRGRI  = grid%nx_array()
    radgrid%NLOENG = grid%nx_array() ! as NRGRI but coud be extended with extra latitudes over poles
    radgrid%NGPTOTG = grid%size()
    radgrid%NDGLG = grid%ny()

    allocate(radgrid%NGLOBALPROC(grid%size()))
    do jp = 1, grid%size()
        radgrid%NGLOBALPROC(jp) = distribution%partition(jp) + 1
    end do

    fs = atlas_functionspace_StructuredColumns(grid=grid, distribution=distribution)
    radgrid%NGPTOT = fs%size()
    lonlat_deg = fs%lonlat()
    call lonlat_deg%data(lonlat_deg_array)
    ALLOCATE(RADGRID%GELAM(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GELAT(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GESLO(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GECLO(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GEMU (RADGRID%NGPTOT))
    DO jp = 1, RADGRID%NGPTOT
        RADGRID%GELAM(jp) = lonlat_deg_array(1, jp) * DEG_TO_RAD
        RADGRID%GELAT(jp) = lonlat_deg_array(2, jp) * DEG_TO_RAD
        RADGRID%GESLO(jp) = SIN(RADGRID%GELAM(jp))
        RADGRID%GECLO(jp) = COS(RADGRID%GELAM(jp))
        RADGRID%GEMU(jp)  = SIN(RADGRID%GELAT(jp))
    END DO

    call distribution%final
    call partitioner%final
    call grid%final

end subroutine

subroutine setup_distribution(distribution, grid, ydgeometry, lmatch_ifs)
    type(atlas_GridDistribution), intent(inout) :: distribution
    type(atlas_StructuredGrid), intent(in) :: grid
    type(GEOMETRY), intent(in) :: ydgeometry
    logical, intent(in) :: lmatch_ifs

    type(atlas_StructuredGrid) :: grid_ifs
    type(atlas_Partitioner) :: partitioner
    type(atlas_functionspace_StructuredColumns) :: fs_ifs
    type(atlas_GridDistribution) :: distribution_ifs

    if (lmatch_ifs) then
      grid_ifs = ATLAS_REDUCEDGAUSSIANGRID( YDGEOMETRY%YRGEM%NLOENG(1:YDGEOMETRY%YRDIM%NDGLG) )
      distribution_ifs = atlas_griddistribution(YDGEOMETRY%YRMP%NGLOBALPROC, PART0=1)
      fs_ifs = atlas_functionspace_StructuredColumns(grid=grid_ifs, distribution=distribution_ifs)
      partitioner = atlas_MatchingPartitioner(fs_ifs)
      distribution = partitioner%partition(grid)
      call partitioner%final
      call fs_ifs%final
      call distribution_ifs%final
      call grid_ifs%final
    else
      partitioner = atlas_partitioner("equal_regions")
      distribution = partitioner%partition(grid)
      call partitioner%final
    end if
end subroutine

subroutine setup_ifs_blocked_fs(fs, ydgeometry, nproma)
    type(GEOMETRY), intent(in) :: ydgeometry
    type(atlas_functionspace_BlockStructuredColumns), intent(inout) :: fs
    integer(JPIM), intent(in) :: nproma
    type(atlas_StructuredGrid) :: grid
    type(atlas_GridDistribution) :: distribution
    grid = ATLAS_REDUCEDGAUSSIANGRID( YDGEOMETRY%YRGEM%NLOENG(1:YDGEOMETRY%YRDIM%NDGLG) )
    distribution = ATLAS_GRIDDISTRIBUTION( YDGEOMETRY%YRMP%NGLOBALPROC, PART0=1 )
    fs = atlas_functionspace_BlockStructuredColumns(grid=grid, distribution=distribution, nproma=nproma)
    call distribution%final
    call grid%final
end subroutine

subroutine setup_rad_blocked_fs(fs, radgrid, nproma)
    type(RADIATION_GRID_STRUCT), intent(in) :: radgrid
    type(atlas_functionspace_BlockStructuredColumns), intent(inout) :: fs
    integer(JPIM), intent(in) :: nproma
    type(atlas_StructuredGrid) :: grid
    type(atlas_GridDistribution) :: distribution
    grid = ATLAS_REDUCEDGAUSSIANGRID( radgrid%NLOENG(1:radgrid%NDGLG) )
    distribution = ATLAS_GRIDDISTRIBUTION( radgrid%NGLOBALPROC, PART0=1 )
    fs = atlas_functionspace_BlockStructuredColumns(grid=grid, distribution=distribution, nproma=nproma)
    call distribution%final
    call grid%final
end subroutine

end program