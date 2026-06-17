module atlas_ecrad_mod
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only : jprb => real64, real32, int32

use atlas_module, only : atlas_functionspace_blockstructuredcolumns

  use radiation_io,             only : nulout
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, INO2, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type
  use radiation_save,           only : save_net_fluxes
  use radiation_setup,          only : tradiation, setup_radiation_scheme
  use radiation_constants,      only : Pi

!   use ecrad_driver_config,      only : driver_config_type
!   use ecrad_driver_read_input,  only : read_input
!   use easy_netcdf
!   use print_matrix_mod,         only : print_matrix

  ! character(len=512), parameter :: ecrad_namelist = "/Users/willem/work/radintg/atlas-bundle/configCY50R1.nam"
  ! character(len=512), parameter :: ecrad_namelist = "/Users/willem/work/radintg/atlas-bundle/source/atlas/atlas_radintg/data/configCY50R1.nam"
  !character(len=512), parameter :: ecrad_namelist = "/Users/willem/work/radintg/atlas-bundle/source/atlas/atlas_radintg/data/configCY49R1.nam"
  character(len=512), parameter :: ecrad_namelist = "/Users/willem/work/radintg/atlas-bundle/source/ecrad/test/ifs/configCY50R1.nam"

type :: ecrad_type
  type(atlas_functionspace_blockstructuredcolumns) :: fs

  ! Configuration for the radiation scheme, IFS style
  type(tradiation)          :: yradiation

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Additional arrays passed to radiation_scheme
  real(jprb), allocatable, dimension(:) :: ccn_land, ccn_sea, sin_latitude, longitude_rad, land_frac
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band

  integer :: ncol

end type

type :: ecrad_view_type
  real(jprb), pointer :: iamu0(:,:), its(:,:), iald(:,:,:), ialp(:,:,:), iemiss(:,:,:), &
                       & iccnl(:,:), iccno(:,:), igelam(:,:), igemu(:,:), islm(:,:), &
                       & ipr(:,:,:), iti(:,:,:), iaprs(:,:,:), ihti(:,:,:), iwv(:,:,:), &
                       & iico2(:,:,:), iich4(:,:,:), iin2o(:,:,:), ino2(:,:,:), ic11(:,:,:), &
                       & ic12(:,:,:), ic22(:,:,:), icl4(:,:,:), ioz(:,:,:), iclc(:,:,:), ilwa(:,:,:), &
                       & iiwa(:,:,:), irwa(:,:,:), iswa(:,:,:), &
                       & ifrso(:,:,:), ifrth(:,:,:), iswfc(:,:,:), ilwfc(:,:,:), ifrsod(:,:), &
                       & ifrted(:,:,:), ifrsodc(:,:), ifrtedc(:,:), ifdir(:,:), icdir(:,:), &
                       & isudu(:,:), iuvdf(:,:), iparf(:,:), iparcf(:,:), itincf(:,:), &
                       & iemit(:,:), ilwderivative(:,:,:), iswdiffuseband(:,:,:), iswdirectband(:,:,:)
  logical :: found_iico2, found_iich4, found_iin2o, found_ino2, found_ic11, found_ic12, found_ic22, &
            & found_icl4, found_ioz, found_iaer
  real(jprb), pointer :: iaer_packed(:,:,:)
  real(jprb), pointer :: iaero_packed(:,:,:)
contains
  procedure :: load => load_ecrad_view
end type

contains

subroutine atlas_ecrad_setup(self_c, fs_c, nlev) bind(C, name="atlas_ecrad_setup")
  use atlas_module, only : atlas_functionspace_blockstructuredcolumns, atlas_Grid
  type(c_ptr), intent(out) :: self_c
  type(c_ptr), value, intent(in) :: fs_c
  integer, value, intent(in) :: nlev
  type(ecrad_type), pointer :: self_f
  type(atlas_Grid) :: grid
  logical :: namelist_exists
  allocate(self_f)
  self_c = c_loc(self_f)

    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
    if (jprb == real32) then
        write(nulout,'(a)') 'Floating-point precision: single'
    else
        write(nulout,'(a)') 'Floating-point precision: double'
    end if

    self_f%fs = atlas_functionspace_blockstructuredcolumns(fs_c)
    self_f%ncol = self_f%fs%size()
    grid = self_f%fs%grid()
    write(nulout,'(a,a)')  "grid:   ", grid%name()
    write(nulout,'(a,i0)') "ncol:   ", self_f%ncol
    write(nulout,'(a,i0)') "nlev:   ", nlev
    write(nulout,'(a,i0)') "nproma: ", self_f%fs%nproma()
    
  ASSOCIATE(yradiation => self_f%yradiation)
  ! Read in the configuration for the radiation scheme, from an IFS-style namelist
  inquire(file=ecrad_namelist, exist=namelist_exists)
  if (.not. namelist_exists) then
    write(nulout,'(a)') 'Error: ECRAD namelist file not found: '//trim(ecrad_namelist)
    error stop
  end if
  call yradiation%rad_config%read(file_name=ecrad_namelist)
  ! ! write(0,*) "--------------------------- Finished reading ECRAD namelist --------------------------"
  ! call yradiation%rad_config%print(iverbose=2)
  ! ! write(0,*) "--------------------------- Finished printing ECRAD config --------------------------"
  
  ! write(0,*) "use_areosols:", yradiation%rad_config%use_aerosols
  ! write(0,*) "n_aerosol_types:", yradiation%rad_config%n_aerosol_types

  ! yradiation%rad_config%n_aerosol_types = 13

    ! Setup aerosols
  if (yradiation%rad_config%use_aerosols) then
    yradiation%yrerad%naermacc = 1 ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  else
    write(0,*) "Not using aerosols in radiation scheme"
    yradiation%yrerad%naermacc = 0
  endif

  ! yradiation%yrephy%NEMISSSCHEME=1 ! 6 interval scheme

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  ! yradiation%yrerad%nlwout = 6
  ! yradiation%yrerad%nlwemiss = 6
  call setup_radiation_scheme(yradiation, .true., file_name=ecrad_namelist)
  ! Or call without specifying the namelist filename, in which case
  ! the default settings are from yoerad.F90
  !call setup_radiation_scheme(yradiation, .true.)

  yradiation%yrerad%nlwout = 6
  yradiation%yrerad%nlwemiss = 6

!  call yradiation%yrephy%print(2,0)
!  call yradiation%yrerad%print(2,0)
  END ASSOCIATE

end subroutine

subroutine atlas_ecrad_teardown(self_c) bind(C, name="atlas_ecrad_teardown")
  type(c_ptr), intent(inout) :: self_c
  type(ecrad_type), pointer :: self_f

  if (c_associated(self_c)) then
    write(0,*) "Tearing down ECRAD in Fortran"
    call c_f_pointer(self_c, self_f)
    deallocate(self_f)
    self_c = c_null_ptr
  end if
end subroutine


subroutine load_ecrad_view(self, fields)
  use atlas_module, only : atlas_fieldset
  implicit none
  class(ecrad_view_type), intent(inout) :: self
  type(atlas_fieldset) :: fields
  call load_required_2d("IMU0",   self%iamu0)
  call load_required_2d("ITS",    self%its)
  call load_required_3d("IALD",   self%iald)
  call load_required_3d("IALP",   self%ialp)
  call load_required_3d("IEMISS", self%iemiss)
  call load_required_2d("ICCNL",  self%iccnl)
  call load_required_2d("ICCNO",  self%iccno)
  call load_required_2d("IGELAM", self%igelam)
  call load_required_2d("IGEMU",  self%igemu)
  call load_required_2d("ISLM",   self%islm)
  call load_required_3d("IPR",    self%ipr)
  call load_required_3d("ITI",    self%iti)
  call load_required_3d("IAPRS",  self%iaprs)
  call load_required_3d("IHTI",   self%ihti)
  call load_required_3d("IWV",    self%iwv)
  call load_optional_3d("IICO2",  self%found_iico2, self%iico2)
  call load_optional_3d("IICH4",  self%found_iich4, self%iich4)
  call load_optional_3d("IIN2O",  self%found_iin2o, self%iin2o)
  call load_optional_3d("INO2",   self%found_ino2, self%ino2)
  call load_optional_3d("IC11",   self%found_ic11, self%ic11)
  call load_optional_3d("IC12",   self%found_ic12, self%ic12)
  call load_optional_3d("IC22",   self%found_ic22, self%ic22)
  call load_optional_3d("ICL4",   self%found_icl4, self%icl4)
  call load_optional_3d("IOZ",    self%found_ioz, self%ioz)
  call load_required_3d("ICLC",   self%iclc)
  call load_required_3d("ILWA",   self%ilwa)
  call load_required_3d("IIWA",   self%iiwa)
  call load_required_3d("IRWA",   self%irwa)
  call load_required_3d("ISWA",   self%iswa)
  call load_optional_3d("IAER",   self%found_iaer, self%iaer_packed)
  call load_required_3d("IAERO",  self%iaero_packed)
  call load_required_3d("IFRSO",  self%ifrso)
  call load_required_3d("IFRTH",  self%ifrth)
  call load_required_3d("ISWFC",  self%iswfc)
  call load_required_3d("ILWFC",  self%ilwfc)
  call load_required_2d("IFRSOD", self%ifrsod)
  call load_required_3d("IFRTED", self%ifrted)
  call load_required_2d("IFRSODC",self%ifrsodc)
  call load_required_2d("IFRTEDC",self%ifrtedc)
  call load_required_2d("IFDIR",  self%ifdir)
  call load_required_2d("ICDIR",  self%icdir)
  call load_required_2d("ISUDU",  self%isudu)
  call load_required_2d("IUVDF",  self%iuvdf)
  call load_required_2d("IPARF",  self%iparf)
  call load_required_2d("IPARCF", self%iparcf)
  call load_required_2d("ITINCF", self%itincf)
  call load_required_2d("IEMIT",  self%iemit)
  call load_required_3d("ILWDERIVATIVE",  self%ilwderivative)
  call load_required_3d("ISWDIFFUSEBAND", self%iswdiffuseband)
  call load_required_3d("ISWDIRECTBAND",  self%iswdirectband)

contains

  subroutine load_required_2d(name, field)
    character(len=*), intent(in) :: name
    real(jprb), pointer :: field(:,:)
    call fields%data(name, field)
  end subroutine

  subroutine load_required_3d(name, field)
    character(len=*), intent(in) :: name
    real(jprb), pointer :: field(:,:,:)
    call fields%data(name, field)
  end subroutine

  subroutine load_optional_3d(name, found, field)
    character(len=*), intent(in) :: name
    logical, intent(out) :: found
    real(jprb), pointer :: field(:,:,:)
    found = fields%has(name)
    if (found) call fields%data(name, field)
  end subroutine
end subroutine

subroutine atlas_ecrad_execute(self_c, fieldset_c) bind(C, name="atlas_ecrad_execute")
  use atlas_module, only : atlas_fieldset
  implicit none
  type(c_ptr), value, intent(in) :: self_c
  type(c_ptr), value, intent(in) :: fieldset_c
  type(ecrad_type), pointer :: self_f
  type(atlas_fieldset) :: fields
  if (.not. c_associated(self_c)) then
    write(0,*) "Error: ECRAD not setup before execution"
    return
  end if
  call c_f_pointer(self_c, self_f)
  if (.not. c_associated(fieldset_c)) then
    write(0,*) "Error: Fieldset not setup before execution"
    return
  end if
  fields = atlas_fieldset(fieldset_c)
  call execute(self_f, fields)
end subroutine

subroutine execute(self, fields)
  use atlas_module, only : atlas_fieldset
  implicit none
  type(ecrad_type) :: self
  type(atlas_fieldset) :: fields
  type(ecrad_view_type) :: view
  real(jprb), pointer :: iaero_packed(:,:,:)
  real(jprb), pointer, contiguous :: iaero_packed_blk(:,:), iaer_packed_blk(:,:)
  real(jprb), pointer, contiguous :: iaero_blk(:,:,:), iaer_blk(:,:,:)
  integer(int32) :: naerosol, nlev, nproma, nblk, ncol, nsw, nlwemiss, nlwout
  integer(int32) :: ib, il, jrl, ibeg, iend, jrof
  real(jprb) :: sol_irrad ! solar irradiance

  integer(int32), allocatable :: iseed(:)

#include "radiation_scheme.intfb.h"
  
  ! ----------------------------------------------------------------
  ! Unpack the data from the fields
  ! ----------------------------------------------------------------
  call validate_dimensions(self, fields)
  call view%load(fields)
  ! if (.not. found_iico2) write(0,*) "Warning: IICO2 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_iich4) write(0,*) "Warning: IICH4 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_iin2o) write(0,*) "Warning: IIN2O not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_ino2)  write(0,*) "Warning: INO2 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_ic11)  write(0,*) "Warning: IC11 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_ic12)  write(0,*) "Warning: IC12 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_ic22)  write(0,*) "Warning: IC22 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_icl4)  write(0,*) "Warning: ICL4 not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_ioz)   write(0,*) "Warning: IOZ not found in fieldset (ldiagforcing condition)"
  ! if (.not. found_iaer) write(0,*) "Warning: IAER not found in fieldset (ldiagforcing condition)"

    sol_irrad = self%single_level%solar_irradiance
    ! write(0,*) "Number of aerosol types:", self%yradiation%rad_config%n_aerosol_types
    ! write(0,*) "Solar irradiance:", sol_irrad
    ! write(0,*) "nsw:", self%yradiation%yrerad%nsw
    ! write(0,*) "nlwemiss:", self%yradiation%yrerad%nlwemiss
    nproma   = size(view%ipr, 1)
    nlev     = size(view%ipr, 2)
    nblk     = size(view%ipr, 3)
    naerosol = self%yradiation%rad_config%n_aerosol_types
    nsw = self%yradiation%yrerad%nsw
    nlwemiss = self%yradiation%yrerad%nlwemiss
    nlwout = self%yradiation%yrerad%nlwout
    ! write(0,*) "naerosol:", naerosol

    ! Map iaero(nproma, naerosol * nlev, nblk) to iaero_4d(nproma, nlev, naerosol, nblk)
    ! if (size(iaero_packed, 2) /= naerosol * nlev) then
      ! write(0,*) "Error: IAERO middle dimension does not match expected nlev * naerosol"
      ! error stop
    ! end if

  if (.not. view%found_iico2) then
    allocate(view%iico2(nproma,nlev,1))
    view%iico2 = 0.0_jprb
    write(0,*) 'Warning: IICO2 missing from input; using zero fallback field'
  end if
  if (.not. view%found_iich4) then
    allocate(view%iich4(nproma,nlev,1))
    view%iich4 = 0.0_jprb
    write(0,*) 'Warning: IICH4 missing from input; using zero fallback field'
  end if
  if (.not. view%found_iin2o) then
    allocate(view%iin2o(nproma,nlev,1))
    view%iin2o = 0.0_jprb
    write(0,*) 'Warning: IIN2O missing from input; using zero fallback field'
  end if
  if (.not. view%found_ino2) then
    allocate(view%ino2(nproma,nlev,1))
    view%ino2 = 0.0_jprb
    write(0,*) 'Warning: INO2 missing from input; using zero fallback field'
  end if
  if (.not. view%found_ic11) then
    allocate(view%ic11(nproma,nlev,1))
    view%ic11 = 0.0_jprb
    write(0,*) 'Warning: IC11 missing from input; using zero fallback field'
  end if
  if (.not. view%found_ic12) then
    allocate(view%ic12(nproma,nlev,1))
    view%ic12 = 0.0_jprb
    write(0,*) 'Warning: IC12 missing from input; using zero fallback field'
  end if
  if (.not. view%found_ic22) then
    allocate(view%ic22(nproma,nlev,1))
    view%ic22 = 0.0_jprb
    write(0,*) 'Warning: IC22 missing from input; using zero fallback field'
  end if
  if (.not. view%found_icl4) then
    allocate(view%icl4(nproma,nlev,1))
    view%icl4 = 0.0_jprb
    write(0,*) 'Warning: ICL4 missing from input; using zero fallback field'
  end if
  if (.not. view%found_ioz) then
    allocate(view%ioz(nproma,nlev,1))
    view%ioz = 0.0_jprb
    write(0,*) 'Warning: IOZ missing from input; using zero fallback field'
  end if
  if (.not. view%found_iaer) then
    allocate(view%iaer_packed(nproma,6*nlev,1))
    view%iaer_packed = 0.0_jprb
  end if

    if (.not. allocated(iseed)) allocate(iseed(nproma))

    ! write(0,*) iaero_blk
    do jrl=1,self%ncol,nproma
      ibeg=jrl
      iend=min(ibeg+nproma-1,self%ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      ! Set a globally reproducible seed independent of chosen nproma
      do jrof=1,nproma
        iseed(jrof) = jrl + jrof - 1
      end do

        ! if (driver_config%iverbose >= 3) then
! #ifndef NO_OPENMP
          ! write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
              !  &  ' processing columns ', ibeg, '-', iend
! #else
        !write(nulout,'(a,i0,a,i0)')  '    Processing columns ', ibeg, '-', iend
! #endif
        ! end if

      iaero_packed_blk => view%iaero_packed(:,:,ib)
      iaero_blk(1:nproma,1:nlev,1:size(view%iaero_packed, 2)/nlev) => iaero_packed_blk

      iaer_packed_blk => view%iaer_packed(:,:,min(ib,size(view%iaer_packed,3)))
      iaer_blk(1:nproma,1:nlev,1:size(view%iaer_packed, 2)/nlev) => iaero_packed_blk

      call radiation_scheme &
             & (self%yradiation, &
        ! *** array dimensions and ranges
             &  KIDIA=1, &                       ! Start column to process
             &  KFDIA=il, &                      ! End column to process
             &  KLON=nproma, &                   ! Number of columns
             &  KLEV=nlev, &                     ! Number of levels
             &  KAEROSOL=naerosol, &             ! Number of aerosol types
         ! *** Single-level fields
             & PSOLAR_IRRADIANCE=sol_irrad, &    ! (W m-2)
             & PMU0=view%iamu0(:,ib), &               ! Cosine of solar zenith angle at surface     DIMS(KLON)
             & PTEMPERATURE_SKIN=view%its(:,ib), &    ! Skin temperature (K)                        DIMS(KLON)
             ! Diffuse and direct components of surface shortwave albedo    
             & PALBEDO_DIF=view%iald(:,:,ib), &       ! albedo_dif                                  DIMS(KLON,NSW)
             & PALBEDO_DIR=view%ialp(:,:,ib), &       ! albedo_dir                                  DIMS(KLON,NSW)
             ! Longwave spectral emissivity    
             & PSPECTRALEMISS=view%iemiss(:,:,ib), &  ! Longwave spectral emissivity                DIMS(KLON,NLWEMISS)
             ! CCN concentration, land and sea    
             & PCCN_LAND=view%iccnl(:,ib), & !                                                      DIMS(KLON)
             & PCCN_SEA=view%iccno(:,ib) ,&           ! CCN concentration, land and sea             DIMS(KLON)
             ! Longitude (radians) and sine of latitude    
             & PGELAM=view%igelam(:,ib), &            ! longitude (radians)                         DIMS(KLON)
             & PGEMU=view%igemu(:,ib), &              ! sine of latitude                            DIMS(KLON)
             ! Land-sea mask    
             & PLAND_SEA_MASK=view%islm(:,ib), &      ! land sea mask                               DIMS(KLON)
        ! *** Variables on full levels    
             ! full level pressure and temperature    
             & PPRESSURE=view%ipr(:,:,ib), &          ! full level pressure (Pa)                    DIMS(KLON,KLEV)
             & PTEMPERATURE=view%iti(:,:,ib),  &      ! full level temperature (K)                  DIMS(KLON,KLEV)
        ! *** Variables on half levels    
             & PPRESSURE_H=view%iaprs(:,:,ib), &      ! half-level pressure (Pa)                    DIMS(KLON,KLEV+1)
             & PTEMPERATURE_H=view%ihti(:,:,ib), &    ! half-level temperature (K)                  DIMS(KLON,KLEV+1)
        ! *** Gas mixing ratios on full levels    
             & PQ=view%iwv(:,:,ib), &                 ! water vapor mixing ratio (kg/kg)            DIMS(KLON,KLEV)
             & PCO2=view%iico2(:,:,min(ib,size(view%iico2,3))), &             ! carbon dioxide mixing ratio (kg/kg)         DIMS(KLON,KLEV)
             & PCH4=view%iich4(:,:,min(ib,size(view%iich4,3))), &             ! methane mixing ratio (kg/kg)                DIMS(KLON,KLEV)
             & PN2O=view%iin2o(:,:,min(ib,size(view%iin2o,3))), &             ! nitrous oxide mixing ratio (kg/kg)          DIMS(KLON,KLEV)
             & PNO2=view%ino2(:,:,min(ib,size(view%ino2,3))), &              ! nitrogen dioxide mixing ratio (kg/kg)       DIMS(KLON,KLEV)
             & PCFC11=view%ic11(:,:,min(ib,size(view%ic11,3))), &            ! CFC-11 mixing ratio (kg/kg)                 DIMS(KLON,KLEV)
             & PCFC12=view%ic12(:,:,min(ib,size(view%ic12,3))), &            ! CFC-12 mixing ratio (kg/kg)                 DIMS(KLON,KLEV)
             & PHCFC22=view%ic22(:,:,min(ib,size(view%ic22,3))), &           ! HCFC-22 mixing ratio (kg/kg)                DIMS(KLON,KLEV)
             & PCCL4=view%icl4(:,:,min(ib,size(view%icl4,3))), &             ! CCl4 mixing ratio (kg/kg)                   DIMS(KLON,KLEV)
             & PO3=view%ioz(:,:,min(ib,size(view%ioz,3))), &                ! Ozone mixing ratio (kg/kg)                  DIMS(KLON,KLEV)
        ! *** Cloud fraction and hydrometeor mass mixing ratios on full levels    
             & PCLOUD_FRAC=view%iclc(:,:,ib), &       ! cloud fraction                              DIMS(KLON,KLEV)
             & PQ_LIQUID=view%ilwa(:,:,ib), &         ! liquid water mixing ratio (kg/kg)           DIMS(KLON,KLEV)
             & PQ_ICE=view%iiwa(:,:,ib), &            ! ice water mixing ratio (kg/kg)              DIMS(KLON,KLEV)
             & PQ_RAIN=view%irwa(:,:,ib), &           ! rain water mixing ratio (kg/kg)             DIMS(KLON,KLEV)
             & PQ_SNOW=view%iswa(:,:,ib), &            ! snow mixing ratio (kg/kg)                   DIMS(KLON,KLEV)
        ! *** Aerosol mass mixing ratios    
             & PAEROSOL_OLD=iaer_blk, &          ! old aerosol mixing ratio (kg/kg)            DIMS(KLON,6,KLEV)
             & PAEROSOL=iaero_blk(:,:,1:naerosol), &   ! aerosol mixing ratio (kg/kg)                DIMS(KLON,KLEV,KAEROSOL)
    ! OUTPUT ARGUMENTS
        ! *** Net fluxes on half-levels (W m-2)
             & PFLUX_SW=view%ifrso(:,:,ib), &       
             & PFLUX_LW=view%ifrth(:,:,ib), &
             & PFLUX_SW_CLEAR=view%iswfc(:,:,ib), &
             & PFLUX_LW_CLEAR=view%ilwfc(:,:,ib),&
        ! *** Surface flux components (W m-2)
             & PFLUX_SW_DN=view%ifrsod(:,ib), &
             & PFLUX_LW_DN=view%ifrted(:,:,ib), &
             & PFLUX_SW_DN_CLEAR=view%ifrsodc(:,ib), &
             & PFLUX_LW_DN_CLEAR=view%ifrtedc(:,ib),&
             ! Direct component of surface flox into horizontal plane
             & PFLUX_DIR=view%ifdir(:,ib), &
             & PFLUX_DIR_CLEAR=view%icdir(:,ib), &
             ! As PFLUX_DIR but into a plane perpendicular to the sun
             & PFLUX_DIR_INTO_SUN=view%isudu(:,ib), &
        !*** Ultraviolet and photosynthetically active radiation fluxes (W m-2)
             & PFLUX_UV=view%iuvdf(:,ib), &
             & PFLUX_PAR=view%iparf(:,ib), &
             & PFLUX_PAR_CLEAR=view%iparcf(:,ib), &
        ! *** Other single-level diagnostics
             ! Top-of-atmosphere incident solar flux (W m-2)
             & PFLUX_SW_DN_TOA=view%itincf(:,ib), &
             ! Diagnosed longwave surface emmisivity across the whole spectrum
             & PEMIS_OUT=view%iemit(:,ib), &
             ! Partial derivative of total-sky longwave upward flux at each level
             ! with respect to upward flux at surface, used to correct heating
             ! rates at gridpoints/timesteps between calls to the full radiation
             ! scheme.  Note that this version uses the convention of level index
             ! increasing downwards, unlike the local variable ZLwDerivative that
             ! is returned from the LW radiation scheme.
             & PLWDERIVATIVE=view%ilwderivative(:,:,ib), &
             ! Surface diffuse and direct downwelling shortwave flux in each
             ! shortwave albedo band, used in RADINTG to update the surface fluxes
             ! accounting for high-resolution albedo information
             & PSWDIFFUSEBAND=view%iswdiffuseband(:,:,ib), &
             & PSWDIRECTBAND=view%iswdirectband(:,:,ib), &
             & ISEED=iseed &! seed for random number generator in stochastic cloud overlap
#ifdef BITIDENTITY_TESTING
       ! OPTIONAL ARGUMENTS for bit-identical results in tests
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             & , &
             & PRE_LIQ=zrgp(1,ifs_config%ire_liq,ib), &
             & PRE_ICE=zrgp(1,ifs_config%ire_ice,ib), &
             & PCLOUD_OVERLAP=zrgp(1,ifs_config%ioverlap,ib), &
             & ISEED=iseed(:,ib) &
#endif
      & )
    end do
    if (.not. view%found_iico2) deallocate(view%iico2)
    if (.not. view%found_iich4) deallocate(view%iich4)
    if (.not. view%found_iin2o) deallocate(view%iin2o)
    if (.not. view%found_ino2)  deallocate(view%ino2)
    if (.not. view%found_ic11)  deallocate(view%ic11)
    if (.not. view%found_ic12)  deallocate(view%ic12)
    if (.not. view%found_ic22)  deallocate(view%ic22)
    if (.not. view%found_icl4)  deallocate(view%icl4)
    if (.not. view%found_ioz)   deallocate(view%ioz)
    if (.not. view%found_iaer)  deallocate(view%iaer_packed)

end subroutine

subroutine validate_dimensions(self, fields)
  use atlas_module, only : atlas_fieldset
  implicit none
  type(ecrad_type), intent(in) :: self
  type(atlas_fieldset), intent(in) :: fields
  integer(int32) :: nproma, nblk, nlev, nsw, nlwemiss, naerosol, nlwout
  integer(int32), allocatable :: field_shape(:)

  nproma = fields%shape("IPR",1)
  nlev   = fields%shape("IPR",2)
  nblk   = fields%shape("IPR",3)
  nsw = self%yradiation%yrerad%nsw
  nlwemiss = self%yradiation%yrerad%nlwemiss
  naerosol = self%yradiation%rad_config%n_aerosol_types
  nlwout = self%yradiation%yrerad%nlwout

  write(0,*) "naerosol types (from config):", naerosol

  call assert_dims("IMU0", [nproma, nblk])
  call assert_dims("ITS",  [nproma, nblk])
  call assert_dims("IALD", [nproma, nsw, nblk])
  call assert_dims("IALP", [nproma, nsw, nblk])
  call assert_dims("IEMISS", [nproma, nlwemiss, nblk])
  call assert_dims("ICCNL", [nproma, nblk])
  call assert_dims("ICCNO", [nproma, nblk])
  call assert_dims("IGELAM", [nproma, nblk])
  call assert_dims("IGEMU", [nproma, nblk])
  call assert_dims("ISLM", [nproma, nblk])
  call assert_dims("IPR", [nproma, nlev, nblk])
  call assert_dims("ITI", [nproma, nlev, nblk])
  call assert_dims("IAPRS", [nproma, nlev + 1, nblk])
  call assert_dims("IHTI", [nproma, nlev + 1, nblk])
  call assert_dims("IWV", [nproma, nlev, nblk])
  call assert_dims("ICLC", [nproma, nlev, nblk])
  call assert_dims("ILWA", [nproma, nlev, nblk])
  call assert_dims("IIWA", [nproma, nlev, nblk])
  call assert_dims("IRWA", [nproma, nlev, nblk])
  call assert_dims("ISWA", [nproma, nlev, nblk])
  call assert_dims("IAERO", [nproma, nlev * naerosol, nblk])
  call assert_dims("IFRSO", [nproma, nlev + 1, nblk])
  call assert_dims("IFRTH", [nproma, nlev + 1, nblk])
  call assert_dims("ISWFC", [nproma, nlev + 1, nblk])
  call assert_dims("ILWFC", [nproma, nlev + 1, nblk])
  call assert_dims("IFRSOD", [nproma, nblk])
  call assert_dims("IFRTED", [nproma, nlwout, nblk])
  call assert_dims("IFRSODC", [nproma, nblk])
  call assert_dims("IFRTEDC", [nproma, nblk])
  call assert_dims("IFDIR", [nproma, nblk])
  call assert_dims("ICDIR", [nproma, nblk])
  call assert_dims("ISUDU", [nproma, nblk])
  call assert_dims("IUVDF", [nproma, nblk])
  call assert_dims("IPARF", [nproma, nblk])
  call assert_dims("IPARCF", [nproma, nblk])
  call assert_dims("ITINCF", [nproma, nblk])
  call assert_dims("IEMIT", [nproma, nblk])
  call assert_dims("ILWDERIVATIVE", [nproma, nlev+1, nblk])
  call assert_dims("ISWDIFFUSEBAND", [nproma, nsw, nblk])
  call assert_dims("ISWDIRECTBAND", [nproma, nsw, nblk])

contains
  
subroutine assert_dims(name, shape)
  implicit none
  character(len=*), intent(in) :: name
  integer(int32), intent(in) :: shape(:)
  integer :: i
  integer(int32), allocatable :: field_shape(:)

  if (.not. fields%has(name)) then
    write(0,*) "Warning: Required field not found in fieldset: ", trim(name)
    error stop
  end if

  field_shape = fields%shape(name)
  if (size(shape) /= size(field_shape)) then
    write(0,*) "Error: ", trim(name), " rank mismatch. expected=", size(shape), " got=", size(field_shape)
    error stop
  end if
  do i = 1, size(shape)
    if (field_shape(i) /= shape(i)) then
      if (shape(i) >= 0) then
        write(0,*) "Error: ", trim(name), " dim ", i, " mismatch. expected=", shape(i), " got=", field_shape(i)
        error stop
      end if
    end if
  end do
end subroutine
end subroutine

end module
