module radconstants

! This module contains constants that are specific to the radiative transfer
! code used in the RRTMGP model.

! This comment from E3SM implementation, and is entirely relevant here:
! TODO: Should this data be handled in a more robust way? Much of this contains
! explicit mappings to indices, which would probably be better handled with get_
! functions. I.e., get_nswbands() could query the kdist objects in case of
! RRTMGP, and the diag indices could look up the actual bands used in the kdist
! objects as well. On that note, this module should probably go away if
! possible in the future, and we should provide more robust access to the
! radiation interface.


use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_abortutils, only: endrun

implicit none
private
save

! SHORTWAVE DATA

! number of shorwave spectral intervals
integer, parameter, public :: nswbands = 14

! Wavenumbers of band boundaries
!
! Note: Currently rad_solar_var extends the lowest band down to
! 100 cm^-1 if it is too high to cover the far-IR. Any changes meant
! to affect IR solar variability should take note of this.

! NOTE: these follow the non-monotonic ordering used for RRTMG 
! - This is necessary because the optical properties files made for RRTMG use this order too.

! NOTE: aside from order, as noted, these values match the ones in 
! RRTMGP coefficients files. But I think we should be *setting* these
! values based on what is in that file, rather than hard-coding it here. 

real(r8),parameter :: wavenum_low(nswbands) = & ! in cm^-1
   (/2600._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
   8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,  820._r8/)
real(r8),parameter :: wavenum_high(nswbands) = & ! in cm^-1
   (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
   12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)

! Mapping from RRTMG shortwave bands to RRTMGP
integer, parameter, dimension(14), public :: rrtmg_to_rrtmgp_swbands = &
   (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)

! Solar irradiance at 1 A.U. in W/m^2 assumed by radiation code
! Rescaled so that sum is precisely 1368.22 and fractional amounts sum to 1.0
real(r8), parameter :: solar_ref_band_irradiance(nswbands) = & 
   (/                                                        &
      12.11_r8,  20.3600000000001_r8, 23.73_r8,              &
      22.43_r8,  55.63_r8, 102.93_r8, 24.29_r8,              &
      345.74_r8, 218.19_r8, 347.20_r8,                       &
      129.49_r8,  50.15_r8,   3.08_r8, 12.89_r8              &
   /)

! These are indices to the band for diagnostic output
! CHANGE: rather than make these parameters, provide subroutines that set them 
!         using the function get_band_index_by_value (which should be called on initializing radiation)
! integer, parameter, public :: idx_sw_diag = 10 ! index to sw visible band (441 - 625 nm) 
! integer, parameter, public :: idx_nir_diag = 8 ! index to sw near infrared (778-1240 nm) band
! integer, parameter, public :: idx_uv_diag = 11 ! index to sw uv (345-441 nm) band

! integer, parameter, public :: rrtmg_sw_cloudsim_band = 9  ! rrtmgp band for .67 micron
! integer, parameter, public :: rrtmgp_sw_cloudsim_band = 10 ! b/c one band moves to beginning

integer, public :: idx_sw_diag ! index to sw visible band (441 - 625 nm) 
integer, public :: idx_nir_diag! index to sw near infrared (778-1240 nm) band
integer, public :: idx_uv_diag ! index to sw uv (345-441 nm) band

! CHANGE: instead of setting rrtmg[p]_sw_cloudsim_band in radconstants, just make it in radiation
! rrtmgp_sw_cloudsim_band = get_band_index_by_value('sw', 0.67_r8, 'micron')  ! rrtmgp band for .67 micron
! same for lw:
! rrtmgp_lw_cloudsim_band = get_band_index_by_value('lw', 10.5_r8, 'micron')

! Number of evenly spaced intervals in rh
! The globality of this mesh may not be necessary
! Perhaps it could be specific to the aerosol
! But it is difficult to see how refined it must be
! for lookup.  This value was found to be sufficient
! for Sulfate and probably necessary to resolve the
! high variation near rh = 1.  Alternative methods
! were found to be too slow.
! Optimal approach would be for cam to specify size of aerosol
! based on each aerosol's characteristics.  Radiation 
! should know nothing about hygroscopic growth!
integer, parameter, public :: nrh = 1000  

! LONGWAVE DATA

! These are indices to the band for diagnostic output (see comment above about change)
! integer, parameter, public :: idx_lw_diag = 7 ! index to (H20 window) LW band
integer, public :: idx_lw_diag

! number of lw bands
integer, parameter, public :: nlwbands = 16

real(r8), parameter :: wavenumber1_longwave(nlwbands) = &! Longwave spectral band limits (cm-1)
    (/   10._r8,  350._r8, 500._r8,   630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, &
       1180._r8, 1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2380._r8, 2600._r8 /)

real(r8), parameter :: wavenumber2_longwave(nlwbands) = &! Longwave spectral band limits (cm-1)
    (/  350._r8,  500._r8,  630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, 1180._r8, &
       1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2380._r8, 2600._r8, 3250._r8 /)

! GASES TREATED BY RADIATION (line spectrae)
integer, public, parameter :: gasnamelength = 5
integer, public, parameter :: nradgas = 8
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)

! what is the minimum mass mixing ratio that can be supported by radiation implementation?
real(r8), public, parameter :: minmmr(nradgas) &
   = epsilon(1._r8)

! Length of "optics type" string specified in optics files.
integer, parameter, public :: ot_length = 32

public :: rad_gas_index

public :: get_number_sw_bands, &
          get_sw_spectral_boundaries, &
          get_lw_spectral_boundaries, &
          get_ref_solar_band_irrad, &
          get_ref_total_solar_irrad, &
          get_solar_band_fraction_irrad, &
          get_idx_sw_diag, &
          get_idx_nir_diag, &
          get_idx_uv_diag, &
          get_idx_lw_diag, &
          get_band_index_by_value

contains
!------------------------------------------------------------------------------
subroutine get_solar_band_fraction_irrad(fractional_irradiance)
   ! provide Solar Irradiance for each band in RRTMG

   ! fraction of solar irradiance in each band
   real(r8), intent(out) :: fractional_irradiance(1:nswbands)
   real(r8) :: tsi ! total solar irradiance

   tsi = sum(solar_ref_band_irradiance)
   fractional_irradiance = solar_ref_band_irradiance / tsi

end subroutine get_solar_band_fraction_irrad
!------------------------------------------------------------------------------
subroutine get_ref_total_solar_irrad(tsi)
   ! provide Total Solar Irradiance assumed by RRTMG

   real(r8), intent(out) :: tsi

   tsi = sum(solar_ref_band_irradiance)

end subroutine get_ref_total_solar_irrad
!------------------------------------------------------------------------------
subroutine get_ref_solar_band_irrad( band_irrad )

   ! solar irradiance in each band (W/m^2)
   real(r8), intent(out) :: band_irrad(nswbands)
 
   band_irrad = solar_ref_band_irradiance

end subroutine get_ref_solar_band_irrad
!------------------------------------------------------------------------------
subroutine get_number_sw_bands(number_of_bands)

   ! number of solar (shortwave) bands in the rrtmg code
   integer, intent(out) :: number_of_bands

   number_of_bands = nswbands

end subroutine get_number_sw_bands

!------------------------------------------------------------------------------
subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)
   ! provide spectral boundaries of each longwave band

   real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries  = wavenumber1_longwave
      high_boundaries = wavenumber2_longwave
   case('m','meter','meters')
      low_boundaries  = 1.e-2_r8/wavenumber2_longwave
      high_boundaries = 1.e-2_r8/wavenumber1_longwave
   case('nm','nanometer','nanometers')
      low_boundaries  = 1.e7_r8/wavenumber2_longwave
      high_boundaries = 1.e7_r8/wavenumber1_longwave
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries  = 1.e4_r8/wavenumber2_longwave
      high_boundaries = 1.e4_r8/wavenumber1_longwave
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenumber2_longwave
      high_boundaries = 1._r8/wavenumber1_longwave
   case default
      call endrun('get_lw_spectral_boundaries: spectral units not acceptable'//units)
   end select

end subroutine get_lw_spectral_boundaries

!------------------------------------------------------------------------------
subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)
   ! provide spectral boundaries of each shortwave band

   real(r8), intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
   character(*), intent(in) :: units ! requested units

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries = wavenum_low
      high_boundaries = wavenum_high
   case('m','meter','meters')
      low_boundaries = 1.e-2_r8/wavenum_high
      high_boundaries = 1.e-2_r8/wavenum_low
   case('nm','nanometer','nanometers')
      low_boundaries = 1.e7_r8/wavenum_high
      high_boundaries = 1.e7_r8/wavenum_low
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries = 1.e4_r8/wavenum_high
      high_boundaries = 1.e4_r8/wavenum_low
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenum_high
      high_boundaries = 1._r8/wavenum_low
   case default
      call endrun('rad_constants.F90: spectral units not acceptable'//units)
   end select

end subroutine get_sw_spectral_boundaries

!------------------------------------------------------------------------------
integer function rad_gas_index(gasname)

   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname
   integer :: igas

   rad_gas_index = -1
   do igas = 1, nradgas
      if (trim(gaslist(igas)).eq.trim(gasname)) then
         rad_gas_index = igas
         return
      endif
   enddo
   call endrun ("rad_gas_index: can not find gas with name "//gasname)
end function rad_gas_index
!------------------------------------------------------------------------------
subroutine get_idx_sw_diag()
   idx_sw_diag = get_band_index_by_value('sw', 500.0_r8, 'nm')
end subroutine

subroutine get_idx_nir_diag()
   idx_nir_diag = get_band_index_by_value('sw', 1000.0_r8, 'nm')
end subroutine

subroutine get_idx_uv_diag()
   idx_uv_diag = get_band_index_by_value('sw', 400._r8, 'nm')
end subroutine

subroutine get_idx_lw_diag()
   idx_lw_diag = get_band_index_by_value('lw', 1000.0_r8, 'cm^-1')
   ! value chosen to match the band used in CESM1/CESM2
end subroutine

function get_band_index_by_value(swlw, targetvalue, units) result(ans)
   character(len=*),intent(in) :: swlw  ! sw or lw bands
   real(r8),intent(in) :: targetvalue  
   character(len=*),intent(in) :: units ! units of targetvale
   integer :: ans
   ! local
   real(r8), allocatable, dimension(:) :: lowboundaries, highboundaries
   real(r8) :: tgt
   integer  :: nbnds, i

   select case (swlw)
   case ('sw','SW','shortwave')
      nbnds = nswbands
      allocate(lowboundaries(nbnds), highboundaries(nbnds))
      lowboundaries = wavenum_low
      highboundaries = wavenum_high
   case ('lw', 'LW', 'longwave')
      nbnds = nlwbands
      allocate(lowboundaries(nbnds), highboundaries(nbnds))
      lowboundaries = wavenumber1_longwave
      highboundaries = wavenumber2_longwave
   case default
      call endrun('rad_constants.F90: get_band_index_by_value: type of bands not accepted '//swlw)
   end select
   ! band info is in cm^-1 but target value may be other units,
   ! so convert targetvalue to cm^-1
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      tgt = targetvalue
   case('m','meter','meters')
      tgt = 1.0_r8 / (targetvalue * 1.e2_r8)
   case('nm','nanometer','nanometers')
      tgt = 1.0_r8 / (targetvalue * 1.e-7_r8)
   case('um','micrometer','micrometers','micron','microns')
      tgt = 1.0_r8 / (targetvalue * 1.e-4_r8)
   case('cm','centimeter','centimeters')
      tgt = 1._r8/targetvalue
   case default
      call endrun('rad_constants.F90: get_band_index_by_value: units not acceptable'//units)
   end select
   ! now just loop through the array
   do i = 1,nbnds
      if ((tgt > lowboundaries(i)) .and. (tgt <= highboundaries(i))) then
         ans = i
         exit
      end if
   end do
   ! Do something if the answer is not found? 
end function get_band_index_by_value

end module radconstants
