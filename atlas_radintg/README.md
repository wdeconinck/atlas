# atlas-radintg

`radintg` is a small Atlas-based driver for experimenting with interpolation between IFS-style fields and an ecRad radiation state. The project combines C++ for the Atlas workflow and Fortran for the ecRad handle setup.

## What it contains

- `atlas-radintg`: the main C++ executable that reads input fields, builds Atlas fieldsets, performs interpolation, and calls into the ecRad setup layer.
- `atlas_ecrad_mod.F90`: the Fortran module that owns the `ecrad_type` handle and exposes C bindings for setup and teardown.
- `MultiFieldCreatorRad.*`: a custom Atlas multifield creator used to pack radiation variables into a contiguous field layout.

## Dependencies

The project expects these packages to be available in your CMake environment:

- Atlas with Fortran enabled
- ecRad
- HDF5

You also need a Fortran compiler with C interoperability support and a C++ compiler compatible with the Atlas build.

## Build

From the project root:

```bash
cmake -S . -B build
cmake --build build -j
```

The top-level `CMakeLists.txt` builds the `radintg` static library and the `atlas-radintg` executable.

If CMake cannot find Atlas or ecRad, set the appropriate package prefixes first, for example:

```bash
cmake -S . -B build \
	-DCMAKE_PREFIX_PATH="/path/to/atlas;/path/to/ecrad"
```

## Run

`atlas-radintg` is configured through AtlasTool options. The most important ones are:

- `--config-file`: Atlas configuration describing the radiation variables
- `--data-file`: HDF5 input file containing the data to read
- `--rad.grid`: target grid name for the radiation grid
- `--rad.nproma`: target `nproma` value for the radiation grid
- `--interpolation`: interpolation method, defaulting to `bicubic`
- `--niter`: number of iterations to run
- `--multifield`: use the custom multifield creator for input packing
- `--on-device`: move data to device memory before interpolation when supported
- `--output.gmsh`: write Gmsh output for the input and output fieldsets

Example:

```bash
atlas-radintg \
	--config-file data/radiation_fields_config.yaml \
	--data-file ../../zrgp_dumps_from_ifs/t21/zrpg_input.hdf5 \
	--rad.grid O32 \
	--rad.nproma 64 \
	--niter 1
```

## ecRad handoff

The C++ side now passes an opaque pointer into `atlas_ecrad_setup`. The Fortran module allocates the `ecrad_type` instance and returns the handle to C++, and `atlas_ecrad_teardown` releases it at the end of the run.

This keeps the ecRad lifetime explicit and avoids creating a new handle inside the iteration loop.

## Notes

- Input variables are filtered to remove empty fields before the Atlas fieldsets are built.
- The project writes optional Gmsh output for debugging the IFS and radiation field layouts.
- `MultiFieldCreatorRad` currently implements the packed multifield creation path used by the radiation workflow.

## TODO

- Besides zrgp_in, I need the state of yradiation from the IFS experiment
- Where is the domain decomposition of the radiation grid done in IFS?
- ecrad_driver_read_input.F90 defines `read_input` to read data from a netcdf file. What is here not needed, what is essential?
- Ultimately we call in `ecrad_driver_ifs_blocked.F90`:
```fortran
        call radiation_scheme &
             & (yradiation, &
             &  1, il, nproma, &                       ! startcol, endcol, ncol
             &  nlev, yradiation%rad_config%n_aerosol_types, &    ! nlev, naerosols
             &  single_level%solar_irradiance, &                               ! solar_irrad
             ! array inputs
             &  zrgp(1,ifs_config%iamu0,ib), zrgp(1,ifs_config%its,ib), &    ! mu0, skintemp
             &  zrgp(1,ifs_config%iald,ib) , zrgp(1,ifs_config%ialp,ib), &    ! albedo_dif, albedo_dir
             &  zrgp(1,ifs_config%iemiss,ib), &                   ! spectral emissivity
             &  zrgp(1,ifs_config%iccnl,ib), zrgp(1,ifs_config%iccno,ib) ,&  ! CCN concentration, land and sea
             &  zrgp(1,ifs_config%igelam,ib),zrgp(1,ifs_config%igemu,ib), &  ! longitude, sine of latitude
             &  zrgp(1,ifs_config%islm,ib), &                     ! land sea mask
             &  zrgp(1,ifs_config%ipr,ib),   zrgp(1,ifs_config%iti,ib),  &   ! full level pressure and temperature
             &  zrgp(1,ifs_config%iaprs,ib), zrgp(1,ifs_config%ihti,ib), &   ! half-level pressure and temperature
             &  zrgp(1,ifs_config%iwv,ib),   zrgp(1,ifs_config%iico2,ib), &
             &  zrgp(1,ifs_config%iich4,ib), zrgp(1,ifs_config%iin2o,ib), &
             &  zrgp(1,ifs_config%ino2,ib),  zrgp(1,ifs_config%ic11,ib), &
             &  zrgp(1,ifs_config%ic12,ib),  zrgp(1,ifs_config%ic22,ib), &
             &  zrgp(1,ifs_config%icl4,ib),  zrgp(1,ifs_config%ioz,ib), &
             &  zrgp(1,ifs_config%iclc,ib),  zrgp(1,ifs_config%ilwa,ib), &
             &  zrgp(1,ifs_config%iiwa,ib),  zrgp(1,ifs_config%irwa,ib), &
             &  zrgp(1,ifs_config%iswa,ib), &
             &  zrgp(1,ifs_config%iaer,ib),  zrgp(1,ifs_config%iaero,ib), &
             ! flux outputs
             &  zrgp(1,ifs_config%ifrso,ib), zrgp(1,ifs_config%ifrth,ib), &
             &  zrgp(1,ifs_config%iswfc,ib), zrgp(1,ifs_config%ilwfc,ib),&
             &  zrgp(1,ifs_config%ifrsod,ib),zrgp(1,ifs_config%ifrted,ib), &
             &  zrgp(1,ifs_config%ifrsodc,ib),zrgp(1,ifs_config%ifrtedc,ib),&
             &  zrgp(1,ifs_config%ifdir,ib), zrgp(1,ifs_config%icdir,ib), &
             &  zrgp(1,ifs_config%isudu,ib), &
             &  zrgp(1,ifs_config%iuvdf,ib), zrgp(1,ifs_config%iparf,ib), &
             &  zrgp(1,ifs_config%iparcf,ib),zrgp(1,ifs_config%itincf,ib), &
             &  zrgp(1,ifs_config%iemit,ib) ,zrgp(1,ifs_config%ilwderivative,ib), &
             &  zrgp(1,ifs_config%iswdiffuseband,ib), zrgp(1,ifs_config%iswdirectband,ib)&
#ifdef BITIDENTITY_TESTING
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  ,pre_liq=zrgp(1,ifs_config%ire_liq,ib), &
             &  pre_ice=zrgp(1,ifs_config%ire_ice,ib), &
             &  pcloud_overlap=zrgp(1,ifs_config%ioverlap,ib), &
             &  iseed=iseed(:,ib) &
#endif
             & )
```
- Perhaps we only need to look at single_level%solar_irradiance
- BITIDENTITY_TESTING
  - Why the need to overwrite effective radii?
  - What does it mean?
  - Should I do it?
  - Why is it not bit-identical otherwise?
- ZRGP_IN contains coordinates. These are IFS grid points then. Do they need to be interpolated or can we use the radiation coordinates. Where are the radiation coordinates computed.
- Looking at the signature of radiation_schere above we list fields like ic11 and many others as input, but they are part of ZRGP_OUT instead of ZRGP_IN in the radiation_fields_config.yaml...