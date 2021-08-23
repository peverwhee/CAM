module rrtmgp_inputs

!--------------------------------------------------------------------------------
! Transform data for state inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.
!
! This code is currently set up to send RRTMGP vertical layers ordered bottom
! to top of model.  Although the RRTMGP is supposed to be agnostic about the   
! vertical ordering problems have arisen trying to use the top to bottom order
! as used by CAM's infrastructure.
!
!--------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use ppgrid,           only: pcols, pver, pverp

use physconst,        only: stebol

use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc
use camsrfexch,       only: cam_in_t

! use radconstants,     only: get_ref_solar_band_irrad, rad_gas_index  ! Not needed anymore (?)
use rad_solar_var,    only: get_variability
use rad_constituents, only: rad_cnst_get_gas

use mcica_subcol_gen, only: mcica_subcol_sw, mcica_subcol_lw

use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp 
use mo_optical_props, only: ty_optical_props, ty_optical_props_2str, ty_optical_props_1scl
! unneeded use mo_rrtmgp_util_string, only: lower_case 
use cam_logfile,         only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   rrtmgp_inputs_init,  &
   rrtmgp_set_state,    &
   rrtmgp_set_gases_lw, &
   rrtmgp_set_gases_sw, &
   rrtmgp_set_cloud_lw, &
   rrtmgp_set_cloud_sw, &
   rrtmgp_set_aer_lw,   &
   rrtmgp_set_aer_sw

real(r8), parameter :: cldmin = 1.0e-80_r8   ! min cloud fraction

real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

! Indices for copying data between cam and rrtmgp arrays
! Assume the rrtmgp vertical index goes bottom to top of atm
integer :: ktopcamm ! cam index of top layer
integer :: ktopradm ! rrtmgp index of layer corresponding to ktopcamm
integer :: ktopcami ! cam index of top interface
integer :: ktopradi ! rrtmgp index of interface corresponding to ktopcami

!==================================================================================================
contains
!==================================================================================================

subroutine rrtmgp_inputs_init(ktcamm, ktradm, ktcami, ktradi)
      
   integer, intent(in) :: ktcamm
   integer, intent(in) :: ktradm
   integer, intent(in) :: ktcami
   integer, intent(in) :: ktradi
   !--------------------------------------------------------------------------------

   ktopcamm = ktcamm
   ktopradm = ktradm
   ktopcami = ktcami
   ktopradi = ktradi

end subroutine rrtmgp_inputs_init

!==================================================================================================

subroutine rrtmgp_set_state( &
   pstate, cam_in, ncol, nlay, nlwbands, &
   nswbands, ngpt_sw, nday, idxday, coszrs, &
   kdist_sw, eccf, t_sfc, emis_sfc, t_rad, &
   pmid_rad, pint_rad, t_day, pmid_day, pint_day, &
   coszrs_day, alb_dir, alb_dif) 
   ! solin was the last (output) parameter originally (bpm) , when we remove it, don't need radconstants anymore
   ! NOTE: We should allow solin to be specified either using old way (like commented here) or with RRTMGP
   !       - add a namelist parameter.

   ! arguments
   type(physics_state), target, intent(in) :: pstate
   type(cam_in_t),              intent(in) :: cam_in
   integer,                     intent(in) :: ncol
   integer,                     intent(in) :: nlay
   integer,                     intent(in) :: nlwbands
   integer,                     intent(in) :: nswbands
   integer,                     intent(in) :: ngpt_sw
   integer,                     intent(in) :: nday
   integer,                     intent(in) :: idxday(:)
   real(r8),                    intent(in) :: coszrs(:)
   real(r8),                    intent(in) :: eccf       ! Earth orbit eccentricity factor

   class(ty_gas_optics_rrtmgp), intent(in) :: kdist_sw  ! spectral information, but is not used here.

   real(r8), intent(out) :: t_sfc(ncol)              ! surface temperature [K] 
   real(r8), intent(out) :: emis_sfc(nlwbands,ncol)  ! emissivity at surface []
   real(r8), intent(out) :: t_rad(ncol,nlay)         ! layer midpoint temperatures [K]
   real(r8), intent(out) :: pmid_rad(ncol,nlay)      ! layer midpoint pressures [Pa]
   real(r8), intent(out) :: pint_rad(ncol,nlay+1)    ! layer interface pressures [Pa]
   real(r8), intent(out) :: t_day(nday,nlay)         ! layer midpoint temperatures [K]
   real(r8), intent(out) :: pmid_day(nday,nlay)      ! layer midpoint pressure [Pa]
   real(r8), intent(out) :: pint_day(nday,nlay+1)    ! layer interface pressures [Pa]
   real(r8), intent(out) :: coszrs_day(nday)         ! cosize of solar zenith angle
   real(r8), intent(out) :: alb_dir(nswbands,nday)   ! surface albedo, direct radiation
   real(r8), intent(out) :: alb_dif(nswbands,nday)   ! surface albedo, diffuse radiation
   ! real(r8), intent(out) :: solin(pcols)             ! incident flux at domain top [W/m2] ! REMOVED -- LET RADIATION FIGURE IT OUT

   ! local variables
   integer :: k, kk, i, iband

   real(r8) :: solar_band_irrad(nswbands) ! rrtmg-assumed solar irradiance in each sw band

   ! real(r8) :: sfac(nswbands)             ! time varying scaling factors due to Solar Spectral  / REMOVED
   !                                        ! Irrad at 1 A.U. per band
   real(r8) :: bnd_irrad
   real(r8) :: solin_day(nday)
   real(r8) :: wavenumber_limits(2,nswbands)

   character(len=*), parameter :: sub='rrtmgp_set_state'

   !--------------------------------------------------------------------------------

   t_sfc = sqrt(sqrt(cam_in%lwup(:ncol)/stebol))

   ! Set surface emissivity to 1.0 here, this is treated in land surface model??
   emis_sfc = 1._r8

   t_rad(:,:ktopradm)    = pstate%t(:ncol,pver:ktopcamm:-1)     ! t_rad, pmid_rad, pint_rad
   pmid_rad(:,:ktopradm) = pstate%pmid(:ncol,pver:ktopcamm:-1)  ! are ordered bottom-to-top (assumed radiation convention)
   pint_rad(:,:ktopradi) = pstate%pint(:ncol,pverp:ktopcami:-1) ! reversed from CAM's top-to-bottom ordering.
   
   if (nlay == pverp) then

      ! add midpoint and top interface values for extra layer
      t_rad(:,nlay)      = pstate%t(:ncol,1)
      pmid_rad(:,nlay)   = 0.5_r8 * pstate%pint(:ncol,1)

      ! pint_rad(:,nlay+1) = 1.e-2_r8 ! rrtmg value
      pint_rad(:,nlay+1) = 1.01_r8

   end if

   ! properties needed at day columns
   do i = 1, nday
      t_day(i,:)    = t_rad(idxday(i),:)
      pmid_day(i,:) = pmid_rad(idxday(i),:)
      pint_day(i,:) = pint_rad(idxday(i),:)
      coszrs_day(i) = coszrs(idxday(i))
   end do
 
   ! <-- bpm -- let radiation deal with SOLIN & variability -->
   ! <-- But keep this code. We should provide an option to use this method. -->
   ! Define solar incident radiation
   ! (bpm) This is done within radiation, so we should remove this block, and not output solin.
   ! call get_ref_solar_band_irrad(solar_band_irrad)
   ! call get_variability(sfac)

   ! solin_day = 0._r8
   ! do i = 1, nday
   !    do iband = 1, nswbands
   !       bnd_irrad = sfac(iband) * solar_band_irrad(iband) * eccf * coszrs_day(i)
   !       solin_day(i) = solin_day(i) + bnd_irrad
   !    end do
   ! end do

   ! solin = 0._r8
   ! do i = 1, nday
   !    solin(idxday(i)) = solin_day(i)
   ! end do
   !<-- end remove -->

   ! <-- begin: old way of setting albedo hard-wired to 14 SW bands -->
   ! ! Surface albedo (band mapping is hardcoded for RRTMG(P) code)
   ! ! This mapping assumes nswbands=14.
   ! if (nswbands /= 14) &
   !    call endrun(sub//': ERROR: albedo band mapping assumes nswbands=14')

   ! do i = 1, nday
   !    ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
   !    alb_dir(1:8,i) = cam_in%aldir(idxday(i))
   !    alb_dif(1:8,i) = cam_in%aldif(idxday(i))
   !    alb_dir(14,i)  = cam_in%aldir(idxday(i))
   !    alb_dif(14,i)  = cam_in%aldif(idxday(i))

   !    ! Set band 24 (or, band 9 counting from 1) to use linear average of UV/visible
   !    ! and near-IR values, since this band straddles 0.7 microns:
   !    alb_dir(9,i) = 0.5_r8*(cam_in%aldir(idxday(i)) + cam_in%asdir(idxday(i)))
   !    alb_dif(9,i) = 0.5_r8*(cam_in%aldif(idxday(i)) + cam_in%asdif(idxday(i)))

   !    ! UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
   !    alb_dir(10:13,i) = cam_in%asdir(idxday(i))
   !    alb_dif(10:13,i) = cam_in%asdif(idxday(i))
   ! enddo
   ! <-- end: old way of setting albedo hard-wired to 14 SW bands -->

   ! More flexible way to assign albedo (from E3SM implementation)
   ! adapted here to loop over bands and cols b/c cam_in has all cols but albedos are daylit cols
   ! We could remove cols loop if we just set albedos for all columns separate from rrtmgp_set_state.
   ! Albedos are input as broadband (visible, and near-IR), and we need to map
   ! these to appropriate bands. Bands are categorized broadly as "visible" or
   ! "infrared" based on wavenumber, so we get the wavenumber limits here
   wavenumber_limits = kdist_sw%get_band_lims_wavenumber()
   ! Loop over bands, and determine for each band whether it is broadly in the
   ! visible or infrared part of the spectrum (visible or "not visible")
   do iband = 1,nswbands
      if (is_visible(wavenumber_limits(1,iband)) .and. &
         is_visible(wavenumber_limits(2,iband))) then

         ! Entire band is in the visible
         do i = 1, nday
            alb_dir(iband,i) = cam_in%asdir(idxday(i))
            alb_dif(iband,i) = cam_in%asdif(idxday(i))
         end do

      else if (.not.is_visible(wavenumber_limits(1,iband)) .and. &
               .not.is_visible(wavenumber_limits(2,iband))) then
         ! Entire band is in the longwave (near-infrared)
         do i = 1, nday
            alb_dir(iband,i) = cam_in%aldir(idxday(i))
            alb_dif(iband,i) = cam_in%aldif(idxday(i))
         end do
      else
         ! Band straddles the visible to near-infrared transition, so we take
         ! the albedo to be the average of the visible and near-infrared
         ! broadband albedos
         do i = 1, nday
            alb_dir(iband,i) = 0.5 * (cam_in%aldir(idxday(i)) + cam_in%asdir(idxday(i)))
            alb_dif(iband,i) = 0.5 * (cam_in%aldif(idxday(i)) + cam_in%asdif(idxday(i)))
         end do
      end if
   end do


   ! Strictly enforce albedo bounds
   where (alb_dir < 0)
       alb_dir = 0.0_r8
   end where
   where (alb_dir > 1)
       alb_dir = 1.0_r8
   end where
   where (alb_dif < 0)
       alb_dif = 0.0_r8
   end where
   where (alb_dif > 1)
       alb_dif = 1.0_r8
   end where

end subroutine rrtmgp_set_state
!

! Function to check if a wavenumber is in the visible or IR
logical function is_visible(wavenumber)

   ! wavenumber in inverse cm (cm^-1)
   real(r8), intent(in) :: wavenumber

   ! Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1
   real(r8), parameter :: visible_wavenumber_threshold = 14286._r8  ! cm^-1

   ! Wavenumber is in the visible if it is above the visible threshold
   ! wavenumber, and in the infrared if it is below the threshold
   if (wavenumber > visible_wavenumber_threshold) then
      is_visible = .true.
   else
      is_visible = .false.
   end if

end function is_visible


!==================================================================================================

subroutine rrtmgp_set_gases_lw(icall, pstate, pbuf, nlay, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2

   ! The memory management for the gas_concs object is internal.  The arrays passed to it
   ! are copied to the internally allocated memory.  Each call to the set_vmr method checks
   ! whether the gas already has memory allocated, and if it does that memory is deallocated
   ! and new memory is allocated.

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: pstate
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: ncol

   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_lw'
   integer :: i !! debugging index
   !--------------------------------------------------------------------------------

   ncol = pstate%ncol

   ! allocate array to pass to RRTMGP gas description object
   allocate(gas_vmr(ncol,nlay))

   ! Access gas mmr from CAM data structures.  Subset and convert mmr -> vmr.
   ! If an extra layer is used copy the mixing ratios from CAM's top model layer.

   ! H20
   call rad_cnst_get_gas(icall, 'H2O', pstate, pbuf, gas_mmr)
   ! water vapor represented as specific humidity in CAM 
   ! gas_vmr goes bottom to top (assumed radiation convention), 
   ! but gas_mmr goes top to bottom (CAM convention)
   ! and gas_vmr can have an extra level 
   gas_vmr(:,:ktopradm) = (gas_mmr(:ncol,pver:ktopcamm:-1) / &
                  (1._r8 - gas_mmr(:ncol,pver:ktopcamm:-1))) * amdw
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm)  ! sets gas @top layer to same as at top of radiation's VMR (e.g., layer 33 is set to same value as was given at layer 32)
   end if
   errmsg = gas_concs%set_vmr('h2o', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting H2O: '//trim(errmsg))
   end if

   ! CO2
   call rad_cnst_get_gas(icall, 'CO2', pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdc
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if

   errmsg = gas_concs%set_vmr('co2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO2: '//trim(errmsg))

   ! O3
   call rad_cnst_get_gas(icall, 'O3',  pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('o3', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O3: '//trim(errmsg))

   ! N2O
   call rad_cnst_get_gas(icall, 'N2O', pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdn
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('n2o', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2O: '//trim(errmsg))

   ! CO not available
   gas_vmr  = 1.e-7_r8

   errmsg = gas_concs%set_vmr('co', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO: '//trim(errmsg))

   ! CH4
   call rad_cnst_get_gas(icall, 'CH4', pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdm
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('ch4', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CH4: '//trim(errmsg))

   ! O2
   call rad_cnst_get_gas(icall, 'O2',  pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo2
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('o2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O2: '//trim(errmsg))

   ! N2 not available
   gas_vmr  = 0.7906_r8
   errmsg = gas_concs%set_vmr('n2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2'//trim(errmsg))


   call rad_cnst_get_gas(icall, 'CFC11', pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdc1
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   errmsg = gas_concs%set_vmr('cfc11', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CFC11: '//trim(errmsg))

   call rad_cnst_get_gas(icall, 'CFC12', pstate, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdc2
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   errmsg = gas_concs%set_vmr('cfc12', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CFC12: '//trim(errmsg))

   deallocate(gas_vmr)

end subroutine rrtmgp_set_gases_lw

!==================================================================================================

subroutine rrtmgp_set_gases_sw( &
   icall, pstate, pbuf, nlay, nday, &
   idxday, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: pstate
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   integer,                     intent(in)    :: nday
   integer,                     intent(in)    :: idxday(:)
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: ncol, i

   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_sw'
   !--------------------------------------------------------------------------------

   ncol = pstate%ncol

   ! allocate array to pass to RRTMGP gas description object
   allocate(gas_vmr(nday,nlay))

   ! Access gas mmr from CAM data structures.  Subset and convert mmr -> vmr.
   ! If an extra layer is used copy the mixing ratios from CAM's top model layer.
  
   ! This could be condensed by making a function, and more robust by using a list of gases.

   ! H20
   call rad_cnst_get_gas(icall, 'H2O', pstate, pbuf, gas_mmr)
   do i = 1, nday
      ! water vapor represented as specific humidity in CAM
      gas_vmr(i,:ktopradm) = (gas_mmr(idxday(i),pver:ktopcamm:-1) / &
                             (1._r8 - gas_mmr(idxday(i),pver:ktopcamm:-1))) * amdw
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('h2o', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting H2O: '//trim(errmsg))
   end if

   ! CO2
   call rad_cnst_get_gas(icall, 'CO2', pstate, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdc
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('co2', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting CO2: '//trim(errmsg))
   end if

   ! O3
   call rad_cnst_get_gas(icall, 'O3',  pstate, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdo
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('o3', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting O3: '//trim(errmsg))
   end if

   ! N2O
   call rad_cnst_get_gas(icall, 'N2O', pstate, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdn
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm)
   end if 
   errmsg = gas_concs%set_vmr('n2o', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting N2O: '//trim(errmsg))
   end if

   ! CO not available
   ! gas_vmr(1:nday,1:ktopradm) = 1.e-7_r8
   errmsg = gas_concs%set_vmr('co', 1.e-7_r8)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting CO: '//trim(errmsg))
   end if

   ! CH4
   call rad_cnst_get_gas(icall, 'CH4', pstate, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdm
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('ch4', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting CH4: '//trim(errmsg))
   end if

   ! O2
   call rad_cnst_get_gas(icall, 'O2',  pstate, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdo2
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('o2', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting O2: '//trim(errmsg))
   end if

   ! N2 not available
   ! gas_vmr(1:nday,1:ktopradm) = 0.7906_r8
   errmsg = gas_concs%set_vmr('n2', 0.7906_r8)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting N2'//trim(errmsg))
   end if

   call rad_cnst_get_gas(icall, 'CFC11', pstate, pbuf, gas_mmr)
   do i = 1,nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i), pver:ktopcamm:-1) * amdc1
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm)
   end if 
   errmsg = gas_concs%set_vmr('cfc11', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting CFC11: '//trim(errmsg))
   end if

   call rad_cnst_get_gas(icall, 'CFC12', pstate, pbuf, gas_mmr)
   do i = 1,nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i), pver:ktopcamm:-1) * amdc2
   end do
   if (nlay == pverp) then
      gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 
   end if
   errmsg = gas_concs%set_vmr('cfc12', gas_vmr)
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': error setting CFC12: '//trim(errmsg))
   end if

   deallocate(gas_vmr)

end subroutine rrtmgp_set_gases_sw

!==================================================================================================

subroutine rrtmgp_set_cloud_lw(state, nlwbands, cldfrac, c_cld_lw_abs, lwkDist, cloud_lw)

   ! Create MCICA stochastic arrays for cloud LW optical properties.

   ! arguments
   type(physics_state),         intent(in)    :: state
   integer,                     intent(in)    :: nlwbands
   real(r8),                    intent(in)    :: cldfrac(pcols,pver)               ! combined cloud fraction (snow plus regular)
   real(r8),                    intent(in)    :: c_cld_lw_abs(nlwbands,pcols,pver) ! combined cloud absorption optics depth (LW)
   class(ty_gas_optics_rrtmgp), intent(in)    :: lwkDist
   type(ty_optical_props_1scl), intent(inout) :: cloud_lw
   ! local vars
   integer :: i
   integer :: ncol
   integer :: ngptlw
   real(r8), allocatable :: taucmcl(:,:,:) ! cloud optical depth [mcica]
   !--------------------------------------------------------------------------------

   ncol   = state%ncol
   ngptlw = lwkDist%get_ngpt()

   allocate(taucmcl(ngptlw,ncol,pver))
   
   !***NB*** this code is currently set up to create the subcols for all model layers
   !         not just the ones where the radiation calc is being done.  Need
   !         to subset cldfrac and c_cld_lw_abs to avoid computing unneeded random numbers.
   
   call mcica_subcol_lw( &
      lwkdist, nlwbands, ngptlw, ncol, ngptlw, &
      state%pmid, cldfrac, c_cld_lw_abs, taucmcl)

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   cloud_lw%tau = 0.0_r8

   ! flip vertical ordering to match RRTMGP (from top-to-bottom in taucmcl to bottom-to-top in cloud_lw%tau)
   do i = 1, ngptlw
      cloud_lw%tau(:,:ktopradm,i) = taucmcl(i,:,pver:ktopcamm:-1)
   end do

   deallocate(taucmcl)

end subroutine rrtmgp_set_cloud_lw

!==================================================================================================

subroutine rrtmgp_set_aer_lw(ncol, nlwbands, aer_lw_abs, aer_lw)

   ! Load aerosol optical properties into the RRTMGP object.

   ! arguments
   integer,                intent(in)    :: ncol
   integer,                intent(in)    :: nlwbands
   real(r8),               intent(in)    :: aer_lw_abs(pcols,pver,nlwbands) ! aerosol absorption optics depth (LW)
   type(ty_optical_props_1scl), intent(inout) :: aer_lw

   !--------------------------------------------------------------------------------

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there.
   aer_lw%tau = 0.0_r8

   ! Subset vertical layers and flip index ordering to match RRTMGP
   aer_lw%tau(:,:ktopradm,:) = aer_lw_abs(:ncol,pver:ktopcamm:-1,:)


end subroutine rrtmgp_set_aer_lw

!==================================================================================================

subroutine rrtmgp_set_cloud_sw( &
   nswbands, nday, nlay, idxday, pmid, cldfrac, &
   c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, c_cld_tau_w_f, kdist_sw, &
   cloud_sw)

   ! Create MCICA stochastic arrays for cloud SW optical properties.

   ! arguments
   integer,                intent(in) :: nswbands
   integer,                intent(in) :: nday
   integer,                intent(in) :: nlay           ! number of layers in rad calc (may include "extra layer")
   integer,                intent(in) :: idxday(:)

   real(r8),               intent(in) :: pmid(nday,nlay)                    ! pressure at layer midpoints (Pa)
   real(r8),               intent(in) :: cldfrac(pcols,pver)                ! combined cloud fraction (snow plus regular)
   real(r8),               intent(in) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8),               intent(in) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8),               intent(in) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud assymetry parameter * w * tau
   real(r8),               intent(in) :: c_cld_tau_w_f(nswbands,pcols,pver) ! combined cloud forward scattered fraction * w * tau

   class(ty_gas_optics_rrtmgp), intent(in)    :: kdist_sw  ! shortwave gas optics object
   type(ty_optical_props_2str), intent(inout) :: cloud_sw  ! cloud optical properties object

   ! local vars
   integer, parameter :: changeseed = 1

   integer :: i, k, kk, ns
   integer :: ngptsw
   integer :: nver       ! nver is the number of cam layers in the SW calc.  It
                         ! does not include the "extra layer".

   real(r8), allocatable :: cldf(:,:)
   real(r8), allocatable :: tauc(:,:,:)
   real(r8), allocatable :: ssac(:,:,:)
   real(r8), allocatable :: asmc(:,:,:)
   real(r8), allocatable :: taucmcl(:,:,:)
   real(r8), allocatable :: ssacmcl(:,:,:)
   real(r8), allocatable :: asmcmcl(:,:,:)

   character(len=32)  :: sub = 'rrtmgp_set_cloud_sw'
   character(len=128) :: errmsg
   !--------------------------------------------------------------------------------

   ngptsw = kdist_sw%get_ngpt()
   nver   = pver - ktopcamm + 1
   
   ! Compute the input quantities needed for the 2-stream optical props
   ! object.  Also subset the vertical levels and the daylight columns
   ! here.  But don't reorder the vertical index because the mcica sub-column
   ! generator assumes the CAM vertical indexing.
   allocate( &
      cldf(nday,nver),           &
      tauc(nswbands,nday,nver),  &
      ssac(nswbands,nday,nver),  &
      asmc(nswbands,nday,nver),  &
      taucmcl(ngptsw,nday,nver), &
      ssacmcl(ngptsw,nday,nver), &
      asmcmcl(ngptsw,nday,nver)  )
   do i = 1, nday
      k = 0
      do kk = ktopcamm, pver
         k = k + 1
         cldf(i,k) = cldfrac(idxday(i),kk)
         do ns = 1, nswbands
            if (c_cld_tau_w(ns,IdxDay(i),kk) > 0._r8) then
               asmc(ns,i,k) =     c_cld_tau_w_g(ns,idxday(i),kk) / &
                              max(c_cld_tau_w  (ns,idxday(i),kk), 1.e-80_r8)
            else
               asmc(ns,i,k) = 0._r8
            endif
   
            tauc(ns,i,k) = c_cld_tau(ns,idxday(i),kk)
            if (tauc(ns,i,k) > 0._r8) then
               ssac(ns,i,k) = max(c_cld_tau_w(ns,idxday(i),kk), 1.e-80_r8) / &
                              max(tauc(ns,i,k), 1.e-80_r8)
            else
               tauc(ns,i,k) = 0._r8
               asmc(ns,i,k) = 0._r8
               ssac(ns,i,k) = 1._r8
            endif
         enddo
      enddo
   enddo

   call mcica_subcol_sw( &
      kdist_sw, nswbands, ngptsw, nday, nlay, nver, changeseed, &
      pmid, cldf, tauc, ssac, asmc,     &
      taucmcl, ssacmcl, asmcmcl)

   ! If there is an extra layer in the radiation then this initialization
   ! will provide the optical properties there.
   ! These should be shape (ncol, nlay, ngpt)
   cloud_sw%tau = 0.0_r8
   cloud_sw%ssa = 1.0_r8
   cloud_sw%g   = 0.0_r8
   ! flip vertical ordering to match RRTMGP
   do i = 1, ngptsw
      cloud_sw%g  (:,:ktopradm,i) = asmcmcl(i,:,pver:ktopcamm:-1)
      cloud_sw%ssa(:,:ktopradm,i) = ssacmcl(i,:,pver:ktopcamm:-1)
      cloud_sw%tau(:,:ktopradm,i) = taucmcl(i,:,pver:ktopcamm:-1)
   end do

   ! delta scaling adjusts for forward scattering
   errmsg = cloud_sw%delta_scale()
   if (len_trim(errmsg) > 0) then
      call endrun(sub//': ERROR: cloud_sw%delta_scale: '//trim(errmsg))
   end if
   
   deallocate( &
      cldf, tauc, ssac, asmc, &
      taucmcl, ssacmcl, asmcmcl )

end subroutine rrtmgp_set_cloud_sw

!==================================================================================================

subroutine rrtmgp_set_aer_sw( &
   nswbands, nday, idxday, aer_tau, aer_tau_w, &
   aer_tau_w_g, aer_tau_w_f, aer_sw)

   ! Load aerosol SW optical properties into the RRTMGP object.
   !
   ! *** N.B. *** The input optical arrays from CAM are dimensioned in the vertical
   !              as 0:pver.  The index 0 is for the extra layer used in the radiation
   !              calculation.  The bottom layer is index pver which corresponds to
   !              index 1 in the RRTMGP arrays.

   ! arguments
   integer,   intent(in) :: nswbands
   integer,   intent(in) :: nday
   integer,   intent(in) :: idxday(:)
   real(r8),  intent(in) :: aer_tau    (pcols,0:pver,nswbands) ! extinction optical depth
   real(r8),  intent(in) :: aer_tau_w  (pcols,0:pver,nswbands) ! single scattering albedo * tau
   real(r8),  intent(in) :: aer_tau_w_g(pcols,0:pver,nswbands) ! asymmetry parameter * w * tau
   real(r8),  intent(in) :: aer_tau_w_f(pcols,0:pver,nswbands) ! forward scattered fraction * w * tau
   type(ty_optical_props_2str), intent(inout) :: aer_sw

   ! local variables
   integer  :: ns
   integer  :: k, kk
   integer  :: i
   !--------------------------------------------------------------------------------
   ! If there is an extra layer in the radiation then this initialization
   ! will provide default values there.
   aer_sw%tau = 0.0_r8
   aer_sw%ssa = 1.0_r8
   aer_sw%g   = 0.0_r8


   ! Rearrange the aerosol optics data directly into the aer_sw object.
   ! Both the subsetting of layers used in the SW calculation, and reversing
   ! the indexing are done here.

   do ns = 1, nswbands
      do k = 1, ktopradm
         ! CAM uses opposite vertical ordering of RRTMGP.  The following
         ! code subsets the vertical layers used in the SW calculation by
         ! starting in the bottom layer and moving up.  It works for
         ! either the case where CAM uses an extra layer above the model
         ! top, or when the radiation calculation doesn't go all the way
         ! to the top of CAM (e.g. in WACCM configurations).
         kk = pver - k + 1

         do i = 1, nday

            if (aer_tau_w(idxday(i),kk,ns) > 1.e-80_r8) then

               aer_sw%g(i,k,ns) = aer_tau_w_g(idxday(i),kk,ns) / &
                                  aer_tau_w  (idxday(i),kk,ns)
            end if

            if (aer_tau(idxday(i),kk,ns) > 0._r8) then

               aer_sw%ssa(i,k,ns) = aer_tau_w(idxday(i),kk,ns) / &
                                    aer_tau  (idxday(i),kk,ns)

               aer_sw%tau(i,k,ns) = aer_tau(idxday(i),kk,ns)
            end if

         end do
      end do
   end do

end subroutine rrtmgp_set_aer_sw

!==================================================================================================

end module rrtmgp_inputs
