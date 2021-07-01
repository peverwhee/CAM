module radiation

!---------------------------------------------------------------------------------
!
! CAM interface to RRTMGP radiation parameterization.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
use spmd_utils,          only: masterproc
use ppgrid,              only: pcols, pver, pverp, begchunk, endchunk
use ref_pres,            only: pref_edge
use physics_types,       only: physics_state, physics_ptend
use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
use camsrfexch,          only: cam_out_t, cam_in_t
use physconst,           only: cappa, cpair, gravit

use time_manager,        only: get_nstep, is_first_restart_step, &
                               get_curr_calday, get_step_size

use rad_constituents,    only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info, &
                               rad_cnst_get_gas, rad_cnst_out

use radconstants,        only: nswbands, nlwbands, idx_sw_diag
                         ! not implemented in rrtmgp yet:
                         !  rrtmg_sw_cloudsim_band, rrtmg_lw_cloudsim_band, &

! use scamMod,             only: scm_crm_mode, single_column, have_cld, cldobs
                               
use cam_history,         only: addfld, add_default, horiz_only, outfld, hist_fld_active
use cam_history_support, only: fillvalue

! RRTMGP modules:
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp

use ioFileMod,           only: getfil
use cam_pio_utils,       only: cam_pio_openfile
use pio,                 only: file_desc_t, var_desc_t,               &
                               PIO_NOERR, PIO_INTERNAL_ERROR,         &
                               pio_seterrorhandling, PIO_BCAST_ERROR, &
                               pio_inq_dimlen, pio_inq_dimid,         &
                               pio_inq_varid, pio_def_var,            &
                               pio_put_var, pio_get_var,              &
                               PIO_NOWRITE, pio_closefile

use cam_abortutils,      only: endrun
use error_messages,      only: handle_err
use perf_mod,            only: t_startf, t_stopf
use cam_logfile,         only: iulog


implicit none
private
save

public :: &
   radiation_readnl,         &! read namelist variables
   radiation_register,       &! registers radiation physics buffer fields
   radiation_nextsw_cday,    &! calendar day of next radiation calculation
   radiation_do,             &! query which radiation calcs are done this timestep
   radiation_init,           &! initialization
   radiation_define_restart, &! define variables for restart
   radiation_write_restart,  &! write variables to restart
   radiation_read_restart,   &! read variables from restart
   radiation_tend,           &! compute heating rates and fluxes
   rad_out_t                  ! type for diagnostic outputs

type rad_out_t
   real(r8) :: solin(pcols)         ! Solar incident flux

   real(r8) :: qrsc(pcols,pver)

   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux

   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsutoa(pcols)        ! upwelling solar flux at TOA

   real(r8) :: fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) :: fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) :: fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns

   real(r8) :: fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) :: fsn200c(pcols)       ! fcns interpolated to 200 mb
   real(r8) :: fsnr(pcols)          ! fns interpolated to tropopause

   real(r8) :: qrlc(pcols,pver)

   real(r8) :: flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) :: flut(pcols)          ! Upward flux at top of model
   real(r8) :: flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) :: lwcf(pcols)          ! longwave cloud forcing

   real(r8) :: fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) :: fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) :: flnr(pcols)          ! net longwave flux interpolated to tropopause

   real(r8) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) :: fldsc(pcols)         ! Clear sky lw flux at srf (down)

   real(r8) :: tot_cld_vistau(pcols,pver)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: tot_icld_vistau(pcols,pver)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8) :: ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth for output on history files
end type rad_out_t

! Control variables set via namelist
character(len=cl) :: coefs_lw_file ! filepath for lw coefficients
character(len=cl) :: coefs_sw_file ! filepath for sw coefficients

integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).

integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: use_rad_dt_cosz  = .false. ! if true, use radiation dt for all cosz calculations
logical :: spectralflux     = .false. ! calculate fluxes (up and down) per band.

logical :: use_rad_uniform_angle = .false. ! if true, use the namelist rad_uniform_angle for the coszrs calculation

! Physics buffer indices
integer :: qrs_idx      = 0 
integer :: qrl_idx      = 0 
integer :: su_idx       = 0 
integer :: sd_idx       = 0 
integer :: lu_idx       = 0 
integer :: ld_idx       = 0 
integer :: fsds_idx     = 0
integer :: fsns_idx     = 0
integer :: fsnt_idx     = 0
integer :: flns_idx     = 0
integer :: flnt_idx     = 0
integer :: cldfsnow_idx = 0 
integer :: cld_idx      = 0 

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

! averaging time interval for zenith angle
real(r8) :: dt_avg = 0._r8
real(r8) :: rad_uniform_angle = -99._r8

! Number of layers in radiation calculations.
integer :: nlay

! Indices for copying data between cam and rrtmgp arrays
! Assume the rrtmgp vertical index goes bottom to top of atm
integer :: ktopcamm ! cam index of top layer
integer :: ktopradm ! rrtmgp index of layer corresponding to ktopcamm
integer :: ktopcami ! cam index of top interface
integer :: ktopradi ! rrtmgp index of interface corresponding to ktopcami

! LW coefficients
! type(ty_gas_optics_specification) :: kdist_lw
type(ty_gas_optics_rrtmgp) :: kdist_lw ! bpm changed here
integer :: ngpt_lw

! SW coefficients
! type(ty_gas_optics_specification) :: kdist_sw
type(ty_gas_optics_rrtmgp) :: kdist_sw ! bpm changed here
integer :: ngpt_sw


! Gases to use in the radiative calculations. 
! RRTMGP kdist initialization needs to know the names of the
! gases before these are available via the rad_cnst interface. 
! TODO: Move this to namelist or somewhere appropriate.
character(len=3), dimension(8) :: active_gases = (/ &
'H2O', 'CO2', 'O3 ', 'N2O', &
'CO ', 'CH4', 'O2 ', 'N2 ' &
/)


!===============================================================================
contains
!===============================================================================


subroutine radiation_readnl(nlfile)

   ! Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=*), parameter :: sub = 'radiation_readnl'

   character(len=cl) :: rrtmgp_coefs_lw_file, rrtmgp_coefs_sw_file
   integer :: rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always
   logical :: rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux
   logical :: rrtmgp_use_rad_uniform_angle
   real    :: rrtmgp_rad_uniform_angle

   namelist /rrtmgp_nl/ rrtmgp_coefs_lw_file, rrtmgp_coefs_sw_file,       &
                        rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always, &
                        rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux,      &
                        rrtmgp_use_rad_uniform_angle,                     &
                        rrtmgp_rad_uniform_angle
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rrtmgp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rrtmgp_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(rrtmgp_coefs_lw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_coefs_lw_file")
   call mpi_bcast(rrtmgp_coefs_sw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_coefs_sw_file")
   call mpi_bcast(rrtmgp_iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_iradsw")
   call mpi_bcast(rrtmgp_iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_iradlw")
   call mpi_bcast(rrtmgp_irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_irad_always")
   call mpi_bcast(rrtmgp_use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_use_rad_dt_cosz")
   call mpi_bcast(rrtmgp_spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_spectralflux")
   call mpi_bcast(rrtmgp_use_rad_uniform_angle, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_use_rad_uniform_angle")
   call mpi_bcast(rrtmgp_rad_uniform_angle, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_rad_uniform_angle")   
   ! Set module data
   coefs_lw_file   = rrtmgp_coefs_lw_file
   coefs_sw_file   = rrtmgp_coefs_sw_file
   iradsw          = rrtmgp_iradsw
   iradlw          = rrtmgp_iradlw
   irad_always     = rrtmgp_irad_always
   use_rad_dt_cosz = rrtmgp_use_rad_dt_cosz
   spectralflux    = rrtmgp_spectralflux
   use_rad_uniform_angle = rrtmgp_use_rad_uniform_angle
   rad_uniform_angle = rrtmgp_rad_uniform_angle

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   dtime  = get_step_size()
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   !----------------------------------------------------------------------- 
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'RRTMGP radiation scheme parameters:'
      write(iulog,10) trim(coefs_lw_file), trim(coefs_sw_file), iradsw, iradlw, &
                      irad_always, use_rad_dt_cosz, spectralflux
   end if

10 format('  LW coefficents file: ',                                a/, &
          '  SW coefficents file: ',                                a/, &
          '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/)

end subroutine radiation_readnl

!================================================================================================

subroutine radiation_register

   ! Register radiation fields in the physics buffer

   use physics_buffer, only: pbuf_add_field, dtype_r8
   use radiation_data, only: rad_data_register

   call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate 
   call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate 

   call pbuf_add_field('FSDS' , 'global',dtype_r8,(/pcols/), fsds_idx) ! Surface solar downward flux
   call pbuf_add_field('FSNS' , 'global',dtype_r8,(/pcols/), fsns_idx) ! Surface net shortwave flux
   call pbuf_add_field('FSNT' , 'global',dtype_r8,(/pcols/), fsnt_idx) ! Top-of-model net shortwave flux

   call pbuf_add_field('FLNS' , 'global',dtype_r8,(/pcols/), flns_idx) ! Surface net longwave flux
   call pbuf_add_field('FLNT' , 'global',dtype_r8,(/pcols/), flnt_idx) ! Top-of-model net longwave flux

   ! If the namelist has been configured for preserving the spectral fluxes, then create
   ! physics buffer variables to store the results.
   if (spectralflux) then
      call pbuf_add_field('SU'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), su_idx) ! shortwave upward flux (per band)
      call pbuf_add_field('SD'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), sd_idx) ! shortwave downward flux (per band)
      call pbuf_add_field('LU'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), lu_idx) ! longwave upward flux (per band)
      call pbuf_add_field('LD'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), ld_idx) ! longwave downward flux (per band)
   end if

   call rad_data_register()

end subroutine radiation_register

!================================================================================================

function radiation_do(op, timestep)

   ! Return true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   !-----------------------------------------------------------------------

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradsw == 1                     &
                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('lw') ! do a longwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradlw == 1                     &
                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
   ! Return calendar day of next sw radiation calculation

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

subroutine radiation_init(pbuf2d)

   ! Initialize the radiation, cloud, and aerosol optics, and solar variability
   ! parameterizations.  
   ! Add fields to the history buffer.

   use physics_buffer,  only: pbuf_get_index
   use phys_control,    only: phys_getopts
   use rad_solar_var,   only: rad_solar_var_init
   use radiation_data,  only: rad_data_init
   use cloud_rad_props, only: cloud_rad_props_init
   use modal_aer_opt,   only: modal_aer_opt_init
   use rrtmgp_inputs,   only: rrtmgp_inputs_init

   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   character(len=128) :: errmsg

   ! names of gases that are available in the model 
   ! -- needed for the kdist initialization routines 
   type(ty_gas_concs) :: available_gases

   integer :: icall, nmodes
   logical :: active_calls(0:N_DIAG)
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer :: ierr

   integer :: dtime

   character(len=*), parameter :: sub = 'radiation_init'
   !-----------------------------------------------------------------------
   
   !
   ! replacement of RRTMG's rrtmg_state_init
   !

   ! Number of layers in radiation calculation is capped by the number of
   ! pressure interfaces below 1 Pa.  When the entire model atmosphere is
   ! below 1 Pa then an extra layer is added to the top of the model for
   ! the purpose of the radiation calculation.
   nlay = count( pref_edge(:) > 1._r8 ) ! pascals (0.01 mbar)

   if (nlay == pverp) then
      ! All pver cam layers are set by radiation
      ktopcamm = 1
      ktopradm = pver
      ! All pverp cam interfaces are set by radiation
      ktopcami = 1
      ktopradi = pverp
   else ! nlay < pverp
      ! nlay layers are set by radiation
      ktopcamm = pverp - nlay
      ktopradm = nlay
      ! nlay+1 interfaces are set by radiation
      ktopcami = pverp - nlay
      ktopradi = nlay + 1
   end if


   call rad_solar_var_init()

   call rrtmgp_inputs_init(ktopcamm, ktopradm, ktopcami, ktopradi)

   call set_available_gases(active_gases, available_gases)   

   ! Initialize kdist_lw object
   call coefs_init(coefs_lw_file, kdist_lw, available_gases)
   ngpt_lw = kdist_lw%get_ngpt()

   ! Initialize kdist_sw object
   call coefs_init(coefs_sw_file, kdist_sw, available_gases)
   ngpt_sw = kdist_sw%get_ngpt()

   call rad_data_init(pbuf2d)  ! initialize output fields for offline driver

   call cloud_rad_props_init()
   
   cld_idx      = pbuf_get_index('CLD')
   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=ierr)

   ! Set the radiation timestep for cosz calculations if requested using
   ! the adjusted iradsw value from radiation
   if (use_rad_dt_cosz)  then
      dtime  = get_step_size()
      dt_avg = iradsw*dtime
   end if

   call phys_getopts(history_amwg_out   = history_amwg,    &
                     history_vdiag_out  = history_vdiag,   &
                     history_budget_out = history_budget,  &
                     history_budget_histfile_num_out = history_budget_histfile_num)

   ! Determine whether modal aerosols are affecting the climate, and if so
   ! then initialize the modal aerosol optics module
   call rad_cnst_get_info(0, nmodes=nmodes)
   if (nmodes > 0) call modal_aer_opt_init()

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if (irad_always > 0) then
      nstep       = get_nstep()
      irad_always = irad_always + nstep
   end if

   call addfld('TOT_CLD_VISTAU',  (/ 'lev' /), 'A',   '1', 'Total gbx cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Total in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Liquid in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Ice in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)

   if (cldfsnow_idx > 0) then
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', 'Snow in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

   ! get list of active radiation calls
   call rad_cnst_get_call_list(active_calls)

   ! Add shortwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('SOLIN'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar insolation', sampling_seq='rad_lwsw')

         call addfld('QRS'//diag(icall),      (/ 'lev' /),  'A', 'K/s',  'Solar heating rate', sampling_seq='rad_lwsw')
         call addfld('QRSC'//diag(icall),     (/ 'lev' /),  'A', 'K/s',  'Clearsky solar heating rate',                     &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNT'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at top of model',                  &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of model',         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net solar flux at top of atmosphere',             &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOAC'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of atmosphere',    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SWCF'//diag(icall),     horiz_only,   'A', 'W/m2', 'Shortwave cloud forcing',                         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSUTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Upwelling solar flux at top of atmosphere',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNIRTOA'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                               'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAC'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                      'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAS'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                              'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', sampling_seq='rad_lwsw')

         call addfld('FSN200'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net shortwave flux at 200 mb',                    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSN200C'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net shortwave flux at 200 mb',           &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSNR'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at tropopause',                    &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('SOLL'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward near infrared direct  to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward visible direct  to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLLD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward near infrared diffuse to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLSD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward visible diffuse to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at surface',                       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at surface',              &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSDS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Downwelling solar flux at surface',               &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSDSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky downwelling solar flux at surface',      &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FUS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave upward flux')
         call addfld('FDS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave downward flux')
         call addfld('FUSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky upward flux')
         call addfld('FDSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky downward flux')

         if (history_amwg) then
            call add_default('SOLIN'//diag(icall),   1, ' ')
            call add_default('QRS'//diag(icall),     1, ' ')
            call add_default('FSNT'//diag(icall),    1, ' ')
            call add_default('FSNTC'//diag(icall),   1, ' ')
            call add_default('FSNTOA'//diag(icall),  1, ' ')
            call add_default('FSNTOAC'//diag(icall), 1, ' ')
            call add_default('SWCF'//diag(icall),    1, ' ')
            call add_default('FSNS'//diag(icall),    1, ' ')
            call add_default('FSNSC'//diag(icall),   1, ' ')
            call add_default('FSUTOA'//diag(icall),  1, ' ')
            call add_default('FSDSC'//diag(icall),   1, ' ')
            call add_default('FSDS'//diag(icall),    1, ' ')
         endif

      end if
   end do

   ! Add longwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('QRL'//diag(icall),     (/ 'lev' /), 'A', 'K/s',  'Longwave heating rate', sampling_seq='rad_lwsw')
         call addfld('QRLC'//diag(icall),    (/ 'lev' /), 'A', 'K/s',  'Clearsky longwave heating rate',                   &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at top of model',                &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at top of model',       &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Upwelling longwave flux at top of model',          &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky upwelling longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
         call addfld('LWCF'//diag(icall),    horiz_only,  'A', 'W/m2', 'Longwave cloud forcing', sampling_seq='rad_lwsw')

         call addfld('FLN200'//diag(icall),  horiz_only,  'A', 'W/m2', 'Net longwave flux at 200 mb',                      &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLN200C'//diag(icall), horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at 200 mb',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNR'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at tropopause',                  &
                                                                           sampling_seq='rad_lwsw')

         call addfld('FLNS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at surface',                     &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at surface',            &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Downwelling longwave flux at surface',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky Downwelling longwave flux at surface',    &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FUL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave upward flux')
         call addfld('FDL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave downward flux')
         call addfld('FULC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky downward flux')

         if (history_amwg) then
            call add_default('QRL'//diag(icall),   1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
         endif

      end if
   end do

   ! Heating rate needed for d(theta)/dt computation
   call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

   if ( history_budget .and. history_budget_histfile_num > 1 ) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if

end subroutine radiation_init

!===============================================================================

subroutine radiation_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   !----------------------------------------------------------------------------


end subroutine radiation_define_restart
  
!===============================================================================

subroutine radiation_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   !----------------------------------------------------------------------------

end subroutine radiation_write_restart
  
!===============================================================================

subroutine radiation_read_restart(file)

   ! read variables from restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   !----------------------------------------------------------------------------

end subroutine radiation_read_restart
  
!===============================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   !----------------------------------------------------------------------- 
   ! 
   ! Driver for radiation computation.
   ! 
   !-----------------------------------------------------------------------
    
   use phys_grid,          only: get_rlat_all_p, get_rlon_all_p
   use cam_control_mod,    only: eccen, mvelpp, lambm0, obliqr

   use shr_orb_mod,        only: shr_orb_decl, shr_orb_cosz

   use mo_gas_concentrations, only: ty_gas_concs
   use rrtmgp_inputs,      only: rrtmgp_set_state, rrtmgp_set_gases_lw, rrtmgp_set_cloud_lw, &
                                 rrtmgp_set_aer_lw, rrtmgp_set_gases_sw, rrtmgp_set_cloud_sw, &
                                 rrtmgp_set_aer_sw

   use aer_rad_props,      only: aer_rad_props_sw, aer_rad_props_lw

   use cloud_rad_props,    only: get_ice_optics_sw,    ice_cloud_get_rad_props_lw,    &
                                 get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
                                 get_snow_optics_sw,   snow_cloud_get_rad_props_lw
   use mo_optical_props,   only: ty_optical_props, ty_optical_props_2str, ty_optical_props_1scl

   use mo_fluxes_byband,   only: ty_fluxes_byband

   use mo_rrtmgp_clr_all_sky, only: rte_lw, rte_sw

   use radheat,            only: radheat_tend

   use radiation_data,     only: rad_data_write

   use interpolate_data,   only: vertinterp
   use tropopause,         only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE

   ! Arguments
   type(physics_state), intent(in), target :: state
   type(physics_ptend), intent(out)        :: ptend
   type(physics_buffer_desc), pointer      :: pbuf(:)
   type(cam_out_t),     intent(inout)      :: cam_out
   type(cam_in_t),      intent(in)         :: cam_in
   real(r8),            intent(out)        :: net_flx(pcols)

   type(rad_out_t), target, optional, intent(out) :: rd_out


   ! Local variables
   type(rad_out_t), pointer :: rd  ! allow rd_out to be optional by allocating a local object
                                   ! if the argument is not present
   logical  :: write_output
  
   integer  :: i, k
   integer  :: lchnk, ncol
   logical  :: dosw, dolw

   real(r8) :: calday          ! current calendar day
   real(r8) :: delta           ! Solar declination angle  in radians
   real(r8) :: eccf            ! Earth orbit eccentricity factor
   real(r8) :: clat(pcols)     ! current latitudes(radians)
   real(r8) :: clon(pcols)     ! current longitudes(radians)
   real(r8) :: coszrs(pcols)   ! Cosine solar zenith angle

   ! Gathered indices of day and night columns 
   !  chunk_column_index = IdxDay(daylight_column_index)
   integer :: Nday           ! Number of daylight columns
   integer :: Nnite          ! Number of night columns
   integer :: IdxDay(pcols)  ! Indices of daylight columns
   integer :: IdxNite(pcols) ! Indices of night columns

   integer :: itim_old

   real(r8), pointer :: cld(:,:)      ! cloud fraction
   real(r8), pointer :: cldfsnow(:,:) ! cloud fraction of just "snow clouds- whatever they are"
   real(r8), pointer :: qrs(:,:)      ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:)      ! longwave  radiative heating rate 
   real(r8), pointer :: fsds(:)  ! Surface solar down flux
   real(r8), pointer :: fsns(:)  ! Surface solar absorbed flux
   real(r8), pointer :: fsnt(:)  ! Net column abs solar flux at model top
   real(r8), pointer :: flns(:)  ! Srf longwave cooling (up-down) flux
   real(r8), pointer :: flnt(:)  ! Net outgoing lw flux at model top

   real(r8), pointer, dimension(:,:,:) :: su => NULL()  ! shortwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: sd => NULL()  ! shortwave spectral flux down
   real(r8), pointer, dimension(:,:,:) :: lu => NULL()  ! longwave  spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld => NULL()  ! longwave  spectral flux down

   ! tropopause diagnostic
   integer  :: troplev(pcols)
   real(r8) :: p_trop(pcols)

   ! state data passed to radiation calc
   real(r8), allocatable :: t_sfc(:)
   real(r8), allocatable :: emis_sfc(:,:)
   real(r8), allocatable :: t_rad(:,:)
   real(r8), allocatable :: pmid_rad(:,:)
   real(r8), allocatable :: pint_rad(:,:)
   real(r8), allocatable :: t_day(:,:)
   real(r8), allocatable :: pmid_day(:,:)
   real(r8), allocatable :: pint_day(:,:)
   real(r8), allocatable :: coszrs_day(:)
   real(r8), allocatable :: alb_dir(:,:)
   real(r8), allocatable :: alb_dif(:,:)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: ice_tau    (nswbands,pcols,pver) ! ice extinction optical depth
   real(r8) :: ice_tau_w  (nswbands,pcols,pver) ! ice single scattering albedo * tau
   real(r8) :: ice_tau_w_g(nswbands,pcols,pver) ! ice assymetry parameter * tau * w
   real(r8) :: ice_tau_w_f(nswbands,pcols,pver) ! ice forward scattered fraction * tau * w
   real(r8) :: ice_lw_abs (nlwbands,pcols,pver)   ! ice absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: liq_tau    (nswbands,pcols,pver) ! liquid extinction optical depth
   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! liquid single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! liquid assymetry parameter * tau * w
   real(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! liquid forward scattered fraction * tau * w
   real(r8) :: liq_lw_abs (nlwbands,pcols,pver) ! liquid absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: cld_tau    (nswbands,pcols,pver) ! cloud extinction optical depth
   real(r8) :: cld_tau_w  (nswbands,pcols,pver) ! cloud single scattering albedo * tau
   real(r8) :: cld_tau_w_g(nswbands,pcols,pver) ! cloud assymetry parameter * w * tau
   real(r8) :: cld_tau_w_f(nswbands,pcols,pver) ! cloud forward scattered fraction * w * tau
   real(r8) :: cld_lw_abs (nlwbands,pcols,pver) ! cloud absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: snow_tau    (nswbands,pcols,pver) ! snow extinction optical depth
   real(r8) :: snow_tau_w  (nswbands,pcols,pver) ! snow single scattering albedo * tau
   real(r8) :: snow_tau_w_g(nswbands,pcols,pver) ! snow assymetry parameter * tau * w
   real(r8) :: snow_tau_w_f(nswbands,pcols,pver) ! snow forward scattered fraction * tau * w
   real(r8) :: snow_lw_abs (nlwbands,pcols,pver)! snow absorption optics depth (LW)

   ! combined cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: cldfprime(pcols,pver)              ! combined cloud fraction (snow plus regular)
   real(r8) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud assymetry parameter * w * tau
   real(r8) :: c_cld_tau_w_f(nswbands,pcols,pver) ! combined cloud forward scattered fraction * w * tau
   real(r8) :: c_cld_lw_abs (nlwbands,pcols,pver) ! combined cloud absorption optics depth (LW)

   ! RRTMGP cloud objects (McICA sampling of cloud optical properties)
   type(ty_optical_props_1scl) :: cloud_lw
   type(ty_optical_props_2str) :: cloud_sw
   real(r8), allocatable :: band_lims_wavenum(:,:)  ! should this be here, or module-level data?

   integer :: icall                 ! index through climate/diagnostic radiation calls
   logical :: active_calls(0:N_DIAG)

   ! gas vmr
   type(ty_gas_concs) :: gas_concs_lw
   type(ty_gas_concs) :: gas_concs_sw

   ! Aerosol radiative properties
   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau
   real(r8) :: aer_lw_abs (pcols,pver,nlwbands)   ! aerosol absorption optics depth (LW)

   ! RRTMGP aerosol objects
   type(ty_optical_props_1scl) :: aer_lw
   type(ty_optical_props_2str) :: aer_sw

   ! Fluxes
   type(ty_fluxes_byband) :: fsw, fswc
   real(r8), allocatable, target :: fsw_up(:,:)
   real(r8), allocatable, target :: fsw_dn(:,:)
   real(r8), allocatable, target :: fsw_up_bnd(:,:,:)
   real(r8), allocatable, target :: fsw_dn_bnd(:,:,:)
   real(r8), allocatable, target :: fswc_up(:,:)
   real(r8), allocatable, target :: fswc_dn(:,:)
   real(r8), allocatable, target :: fswc_up_bnd(:,:,:)
   real(r8), allocatable, target :: fswc_dn_bnd(:,:,:)

   type(ty_fluxes_byband) :: flw, flwc
   real(r8), allocatable, target :: flw_up(:,:)
   real(r8), allocatable, target :: flw_dn(:,:)
   real(r8), allocatable, target :: flw_up_bnd(:,:,:)
   real(r8), allocatable, target :: flw_dn_bnd(:,:,:)
   real(r8), allocatable, target :: flwc_up(:,:)
   real(r8), allocatable, target :: flwc_dn(:,:)
   real(r8), allocatable, target :: flwc_up_bnd(:,:,:)
   real(r8), allocatable, target :: flwc_dn_bnd(:,:,:)

   real(r8) :: fns(pcols,pverp)     ! net shortwave flux
   real(r8) :: fcns(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) :: fnl(pcols,pverp)     ! net longwave flux
   real(r8) :: fcnl(pcols,pverp)    ! net clear-sky longwave flux
  
   real(r8) :: ftem(pcols,pver)        ! Temporary workspace for outfld variables

   character(len=128) :: errmsg

   character(len=*), parameter :: sub = 'radiation_tend'


   real(r8), pointer     :: gas_mmr(:,:)  !!!!! debugging 
   !--------------------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol = state%ncol

   if (present(rd_out)) then
      rd => rd_out
      write_output = .false.
   else
      allocate(rd)
      write_output=.true.
   end if

   dosw = radiation_do('sw')      ! do shortwave heating calc this timestep?
   dolw = radiation_do('lw')      ! do longwave heating calc this timestep?

   ! Cosine solar zenith angle for current time step
   calday = get_curr_calday()
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

   call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, &
                     delta, eccf)

   if (use_rad_uniform_angle) then
      do i = 1, ncol
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg, uniform_angle=rad_uniform_angle)
      end do
   else
      do i = 1, ncol
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg)
      end do
   end if

   ! Gather night/day column indices.
   Nday = 0
   Nnite = 0
   do i = 1, ncol
      if ( coszrs(i) > 0.0_r8 ) then
         Nday = Nday + 1
         IdxDay(Nday) = i
      else
         Nnite = Nnite + 1
         IdxNite(Nnite) = i
      end if
   end do

   ! bpm: if Nday = 0, then we should not do shortwave:
   if (nday == 0) then
      dosw = 0
   end if


   ! Associate pointers to physics buffer fields
   itim_old = pbuf_old_tim_idx()
   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf,                            &
                          cldfsnow_idx,                    &
                          cldfsnow,                        &
                          start=(/1,1,itim_old/),          &
                          kount=(/pcols,pver,1/) )
   end if
   call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, qrs_idx, qrs)
   call pbuf_get_field(pbuf, qrl_idx, qrl)

   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsds_idx, fsds)

   call pbuf_get_field(pbuf, flns_idx, flns)
   call pbuf_get_field(pbuf, flnt_idx, flnt)

   if (spectralflux) then
      call pbuf_get_field(pbuf, su_idx, su)
      call pbuf_get_field(pbuf, sd_idx, sd)
      call pbuf_get_field(pbuf, lu_idx, lu)
      call pbuf_get_field(pbuf, ld_idx, ld)
   end if

   ! Find tropopause height if needed for diagnostic output
   if (hist_fld_active('FSNR') .or. hist_fld_active('FLNR')) then
      call tropopause_find(state, troplev, tropP=p_trop, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
   end if

   ! bpm -- debug -- check for negative water
   call rad_cnst_get_gas(0, 'H2O', state, pbuf, gas_mmr)
   if (maxval(gas_mmr) > 1.0_r8) then
      write(iulog,*) '[radiation_tend] H2O gas mmr min/max:',minval(gas_mmr),', ',maxval(gas_mmr)
      call endrun(sub//': Water vapor mass mixing ratio out of bounds.')
   end if


   if (dosw .or. dolw) then

      allocate( &
         t_sfc(ncol),             &
         emis_sfc(nlwbands,ncol), &
         t_rad(ncol,nlay),        &
         pmid_rad(ncol,nlay),     &
         pint_rad(ncol,nlay+1),   &
         t_day(nday,nlay),        &
         pmid_day(nday,nlay),     &
         pint_day(nday,nlay+1),   &
         coszrs_day(nday),        &
         alb_dir(nswbands,nday),  &
         alb_dif(nswbands,nday)   )
      
      call rrtmgp_set_state( &
         state, cam_in, ncol, nlay, nlwbands,             &
         nswbands, ngpt_sw, nday, idxday, coszrs,         &
         kdist_sw, eccf, t_sfc, emis_sfc, t_rad,          &
         pmid_rad, pint_rad, t_day, pmid_day, pint_day,   &
         coszrs_day, alb_dir, alb_dif, rd%solin)

      call clipper(t_rad, kdist_lw%get_temp_min(), kdist_lw%get_temp_max())      ! check bounds for t_rad


      call t_startf('cldoptics')

      if (cldfsnow_idx > 0) then
         do k = 1, pver
            do i = 1, ncol
               cldfprime(i,k) = max(cld(i,k), cldfsnow(i,k))
            end do
         end do
      else
         cldfprime(:ncol,:) = cld(:ncol,:)
      end if
      
      
      if (dosw) then
         call get_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
         call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
         cld_tau(:,:ncol,:)     =  liq_tau(:,:ncol,:)     + ice_tau(:,:ncol,:)
         cld_tau_w(:,:ncol,:)   =  liq_tau_w(:,:ncol,:)   + ice_tau_w(:,:ncol,:)
         cld_tau_w_g(:,:ncol,:) =  liq_tau_w_g(:,:ncol,:) + ice_tau_w_g(:,:ncol,:)
         cld_tau_w_f(:,:ncol,:) =  liq_tau_w_f(:,:ncol,:) + ice_tau_w_f(:,:ncol,:)
         if (cldfsnow_idx > 0) then
            ! add in snow
            call get_snow_optics_sw(state, pbuf, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)
            do i = 1, ncol
               do k = 1, pver
                  if (cldfprime(i,k) > 0.) then
                     c_cld_tau(:,i,k)     = ( cldfsnow(i,k)*snow_tau(:,i,k) &
                                             + cld(i,k)*cld_tau(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w(:,i,k)   = ( cldfsnow(i,k)*snow_tau_w(:,i,k)  &
                                             + cld(i,k)*cld_tau_w(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w_g(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_g(:,i,k) &
                                             + cld(i,k)*cld_tau_w_g(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w_f(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_f(:,i,k) &
                                             + cld(i,k)*cld_tau_w_f(:,i,k) )/cldfprime(i,k)
                  else
                     c_cld_tau(:,i,k)     = 0._r8
                     c_cld_tau_w(:,i,k)   = 0._r8
                     c_cld_tau_w_g(:,i,k) = 0._r8
                     c_cld_tau_w_f(:,i,k) = 0._r8
                  end if
               end do
            end do
         else
            c_cld_tau(:,:ncol,:)     = cld_tau(:,:ncol,:)
            c_cld_tau_w(:,:ncol,:)   = cld_tau_w(:,:ncol,:)
            c_cld_tau_w_g(:,:ncol,:) = cld_tau_w_g(:,:ncol,:)
            c_cld_tau_w_f(:,:ncol,:) = cld_tau_w_f(:,:ncol,:)
         end if

         call t_startf('cldoptics_sw_mcica')

         ! Create MCICA stochastic arrays for cloud optical properties.
         ! => cloud_sw 
         ! Note the init method is allocating memory, but there is no method to
         ! deallocate the memory.  Not sure whether Fortran will deallocate this
         ! memory when the local object goes out of scope.
         ! bpm / updated to `init_and_alloc_2str` method from `init_2stream`
         !     / arguments: ncol, nlay, band_lims_wvn, band_lims_gpt, name
         ! but that is private, so have to use alloc_2str .. which can do alloc_only_2str(ncol, nlay)
         ! band_lims_wvn = band_lims_wavenum from kdist, how to get it?
         ! band_lims_wavenum = kdist_sw%get_band_lims_wavenumber()
         ! ngpt_sw is not needed for alloc_only_2str
         ! errmsg = cloud_sw%alloc_2str(nday, nlay)
         write(iulog,*) 'Find out what nday and nlay are: ',nday,', ',nlay
         errmsg = cloud_sw%alloc_2str(nday, nlay, kdist_sw, name='shortwave cloud optics')
         if (len_trim(errmsg) > 0) &
            call endrun(sub//': ERROR: cloud_sw%alloc_2str: '//trim(errmsg))
         cloud_sw%tau = 0
         cloud_sw%ssa = 1
         cloud_sw%g   = 0
   
         call rrtmgp_set_cloud_sw( &
            nswbands, nday, nlay, idxday, pmid_day, cldfprime, &
            c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, c_cld_tau_w_f, kdist_sw, &
            cloud_sw)

         call t_stopf('cldoptics_sw_mcica')
         ! cloud optical depth fields for the visible band
         rd%tot_icld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)
         rd%liq_icld_vistau(:ncol,:) = liq_tau(idx_sw_diag,:ncol,:)
         rd%ice_icld_vistau(:ncol,:) = ice_tau(idx_sw_diag,:ncol,:)
         if (cldfsnow_idx > 0) then
            rd%snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
         endif

         ! multiply by total cloud fraction to get gridbox value
         rd%tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)

         ! add fillvalue for night columns
         do i = 1, Nnite
            rd%tot_cld_vistau(IdxNite(i),:)   = fillvalue
            rd%tot_icld_vistau(IdxNite(i),:)  = fillvalue
            rd%liq_icld_vistau(IdxNite(i),:)  = fillvalue
            rd%ice_icld_vistau(IdxNite(i),:)  = fillvalue
            if (cldfsnow_idx > 0) then
               rd%snow_icld_vistau(IdxNite(i),:) = fillvalue
            end if
         end do

         if (write_output) call radiation_output_cld(lchnk, ncol, rd)

      end if   ! if (dosw)

      if (dolw) then
         call ice_cloud_get_rad_props_lw(state, pbuf, ice_lw_abs)
         call liquid_cloud_get_rad_props_lw(state, pbuf, liq_lw_abs)
         cld_lw_abs(:,:ncol,:) = liq_lw_abs(:,:ncol,:) + ice_lw_abs(:,:ncol,:)
         if (cldfsnow_idx > 0) then
            ! add in snow
            call snow_cloud_get_rad_props_lw(state, pbuf, snow_lw_abs)
            do i = 1, ncol
               do k = 1, pver
                  if (cldfprime(i,k) > 0._r8) then
                     c_cld_lw_abs(:,i,k) = ( cldfsnow(i,k)*snow_lw_abs(:,i,k) &
                                            + cld(i,k)*cld_lw_abs(:,i,k) )/cldfprime(i,k)
                  else
                     c_cld_lw_abs(:,i,k) = 0._r8
                  end if
               end do
            end do
         else
            c_cld_lw_abs(:,:ncol,:) = cld_lw_abs(:,:ncol,:)
         end if

         call t_startf('cldoptics_lw_mcica')

         ! Create MCICA stochastic arrays for cloud optical properties.
         ! Note the init method is allocating memory, but there is not method to
         ! deallocate the memory.  Not sure whether Fortran will deallocate this
         ! memory when the local object goes out of scope.
         ! nit_1scalar(ncol, nlay, ngpt_lw) is replaced by alloc_1scl(ncol, nlay)
         errmsg = cloud_lw%alloc_1scl(ncol, nlay, kdist_lw, name='longwave cloud optics') 
         if (len_trim(errmsg) > 0) &
            call endrun(sub//': ERROR: cloud_lw%init_1scalar: '//trim(errmsg))
         cloud_lw%tau = 0.0
         !! reorder cloud_lw%tau -> (ncol, nlay, ngpt)
         !! but     c_cld_lw_abs -> (ngpt, ncol, nlay)
         do i = 1, ncol
            do k = 2, nlay
               cloud_lw%tau(i,k,:) = c_cld_lw_abs(:, i, k-1)
            end do
         end do
         
         call rrtmgp_set_cloud_lw( &
            state, nlwbands, cldfprime, c_cld_lw_abs, kdist_lw, &
            cloud_lw)

         call t_stopf('cldoptics_lw_mcica')

      end if   ! if (dolw)

      call t_stopf('cldoptics')

      !============================
      ! Solar radiation computation
      !============================
      if (dosw) then

         ! allocate object for aerosol optics (bpm change from init_2stream to alloc_2str)
         errmsg = aer_sw%alloc_2str(nday, nlay, kdist_sw%get_band_lims_wavenumber(), &
                                    name='shortwave aerosol optics') 
         if (len_trim(errmsg) > 0) &
            call endrun(sub//': ERROR: aer_sw%alloc_2str: '//trim(errmsg))

         ! allocate output storage and associate it with object components
         allocate(fsw_up(nday,nlay+1), fsw_dn(nday,nlay+1))
         allocate(fsw_up_bnd(nday,nlay+1,nswbands), fsw_dn_bnd(nday,nlay+1,nswbands))

         fsw%flux_up     => fsw_up
         fsw%flux_dn     => fsw_dn
         fsw%bnd_flux_up => fsw_up_bnd
         fsw%bnd_flux_dn => fsw_dn_bnd

         allocate(fswc_up(nday,nlay+1), fswc_dn(nday,nlay+1))
         allocate(fswc_up_bnd(nday,nlay+1,nswbands), fswc_dn_bnd(nday,nlay+1,nswbands))

         fswc%flux_up     => fswc_up
         fswc%flux_dn     => fswc_dn
         fswc%bnd_flux_up => fswc_up_bnd
         fswc%bnd_flux_dn => fswc_dn_bnd

         ! Get the active climate/diagnostic shortwave calculations
         call rad_cnst_get_call_list(active_calls)

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1

            if (active_calls(icall)) then

               ! set gas concentrations
               ! bpm initialize gas_concs_sw: (available_gases or active_gases?)
               errmsg = gas_concs_sw%init(active_gases)
               if (len_trim(errmsg) > 0) then
                  call endrun(sub//': ERROR code returned by gas_concs_sw%init: '//trim(errmsg))
               end if
               call rrtmgp_set_gases_sw(icall, state, pbuf, nlay, nday, idxday, gas_concs_sw)

               ! Aerosol shortwave optical properties
               call t_startf('aeroptics')
               call aer_rad_props_sw(icall, state, pbuf, nnite, idxnite, &
                                     aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)
               ! NOTE: CAM fields are products tau, tau*ssa, tau*ssa*asy, tau*ssa*asy*fsf
               !       but RRTMGP is expecting just the values per band.
               !       rrtmgp_set_aer_sw does the division and puts values into aer_sw:
               !       aer_sw%g   = aer_tau_w_g / aer_taw_w
               !       aer_sw%ssa = aer_tau_w / aer_tau
               !       aer_sw%tau = aer_tau
               call rrtmgp_set_aer_sw( &
                  nswbands, nday, idxday, aer_tau, aer_tau_w, &
                  aer_tau_w_g, aer_tau_w_f, aer_sw)
               call t_stopf('aeroptics')               
                  
               ! Compute SW fluxes

               cam_out%sols  = 0.0_r8
               cam_out%soll  = 0.0_r8
               cam_out%solsd = 0.0_r8
               cam_out%solld = 0.0_r8

               fsw%flux_up     = 0.0_r8
               fsw%flux_dn     = 0.0_r8
               fsw%bnd_flux_up = 0.0_r8
               fsw%bnd_flux_dn = 0.0_r8

               fswc%flux_up     = 0.0_r8
               fswc%flux_dn     = 0.0_r8
               fswc%bnd_flux_up = 0.0_r8
               fswc%bnd_flux_dn = 0.0_r8

               ! check that optical properties are in bounds:
               call clipper(cloud_sw%tau, 0._r8, huge(cloud_sw%tau))
               call clipper(cloud_sw%ssa, 0._r8, 1._r8)
               call clipper(cloud_sw%g,  -1._r8, 1._r8)
               call clipper(aer_sw%tau,   0._r8, huge(aer_sw%tau))
               call clipper(aer_sw%ssa,   0._r8, 1._r8)
               call clipper(aer_sw%g,    -1._r8, 1._r8)

               call t_startf('rad_sw')
               errmsg = rte_sw( kdist_sw,     &
                                gas_concs_sw, &
                                pmid_day,     & 
                                t_day,        &
                                pint_day,     &
                                coszrs_day,   &
                                alb_dir,      &
                                alb_dif,      &
                                cloud_sw,     &
                                fsw,          & 
                                fswc,         &
                                aer_props=aer_sw)
               call t_stopf('rad_sw')

               if (len_trim(errmsg) > 0) then
                  print*, 'ERROR return from rte_sw at lchnk=', lchnk
                  call endrun(sub//': ERROR code returned by rte_sw: '//trim(errmsg))
               end if

               ! Transform RRTMGP outputs to CAM outputs
               call set_sw_diags()

               if (write_output) call radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)

            end if
         end do    ! loop over diagnostic calcs

         deallocate(fsw_up, fsw_dn)
         deallocate(fsw_up_bnd, fsw_dn_bnd)
         deallocate(fswc_up, fswc_dn)
         deallocate(fswc_up_bnd, fswc_dn_bnd)

      end if  ! if (dosw)

      ! Output aerosol mmr
      call rad_cnst_out(0, state, pbuf)
                 
      !===============================
      ! Longwave radiation computation
      !===============================
      if (dolw) then
         ! initialize/allocate object for aerosol optics (note, don't just give it nlwbands b/c wrong type)
         errmsg = aer_lw%alloc_1scl(ncol,                                & 
                                    nlay,                                &
                                    kdist_lw%get_band_lims_wavenumber(), &
                                    name='longwave aerosol optics') 

         if (len_trim(errmsg) > 0) &
            call endrun(sub//': ERROR: aer_lw%init_1scalar: '//trim(errmsg))

         ! allocate output storage and associate it with object components
         allocate(flw_up(ncol,nlay+1), flw_dn(ncol,nlay+1))
         allocate(flw_up_bnd(ncol,nlay+1,nlwbands), flw_dn_bnd(ncol,nlay+1,nlwbands))

         flw%flux_up     => flw_up
         flw%flux_dn     => flw_dn
         flw%bnd_flux_up => flw_up_bnd
         flw%bnd_flux_dn => flw_dn_bnd

         allocate(flwc_up(ncol,nlay+1), flwc_dn(ncol,nlay+1))
         allocate(flwc_up_bnd(ncol,nlay+1,nlwbands), flwc_dn_bnd(ncol,nlay+1,nlwbands))

         flwc%flux_up     => flwc_up
         flwc%flux_dn     => flwc_dn
         flwc%bnd_flux_up => flwc_up_bnd
         flwc%bnd_flux_dn => flwc_dn_bnd

         ! get list of diagnostic calls
         call rad_cnst_get_call_list(active_calls)

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1

            if (active_calls(icall)) then
               ! initialize the gas concentrations (do we need to do this inside the loop?)
               errmsg = gas_concs_lw%init(active_gases)
               if (len_trim(errmsg) > 0) then
                  call endrun(sub//': ERROR code returned by gas_concs_lw%init: '//trim(errmsg))
               end if
               call rrtmgp_set_gases_lw(icall, state, pbuf, nlay, gas_concs_lw)

               call t_startf('aeroptics')
               call aer_rad_props_lw(icall, state, pbuf,  aer_lw_abs)
               call rrtmgp_set_aer_lw(ncol, nlwbands, aer_lw_abs, aer_lw)
               call t_stopf('aeroptics')
               
               ! check that optical properties are in bounds:
               call clipper(cloud_lw%tau, 0._r8, huge(cloud_lw%tau))
               call clipper(aer_lw%tau,   0._r8, huge(aer_lw%tau))
               
               ! Compute LW fluxes
               call t_startf('rad_lw')

               errmsg = rte_lw(kdist_lw,        &
                               gas_concs_lw,    &
                               pmid_rad, t_rad, &
                               pint_rad,        &
                               t_sfc,           &
                               emis_sfc,        &
                               cloud_lw,        &
                               flw,             &
                               flwc,            &
                               aer_props=aer_lw)  
               if (len_trim(errmsg) > 0) then
                  call endrun(sub//': ERROR code returned by rte_lw: '//trim(errmsg))
               end if

               call t_stopf('rad_lw')

               call set_lw_diags() ! Reverse direction of LW fluxes back to TOP-to-BOTTOM
                                   ! And derive LW dry static energy tendency (QRL, rd%QRLC (J/kg/s))

               if (write_output) then
                  ! QRL retrieved from pbuf and divided by cpair [(J/(kg s)) / (J/(K kg)) = K/s]
                  call radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out)
               end if

            end if
         end do

         ! deallocate output storage
         deallocate(flw_up, flw_dn)
         deallocate(flw_up_bnd, flw_dn_bnd)
         deallocate(flwc_up, flwc_dn)
         deallocate(flwc_up_bnd, flwc_dn_bnd)

      end if  ! if (dolw)

      deallocate( &
         t_sfc, emis_sfc, t_rad, pmid_rad, pint_rad,     &
         t_day, pmid_day, pint_day, coszrs_day, alb_dir, &
         alb_dif)

   else   !  if (dosw .or. dolw) then
      ! convert radiative heating rates from Q*dp to Q for energy conservation
      ! ** if you change qrs and qrl from J/kg/s here, then it won't be a DSE tendency,
      !    yet it is expected to be in radheat_tend to get ptend%s
      !    Does not matter if qrs and qrl are zero on these time steps

      do k =1 , pver
         do i = 1, ncol
            write(iulog,*) 'No radiation. k=',k,' qrs= ',qrs(i,k), ' qrl= ',qrl(i,k) 
            ! qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            ! qrl(i,k) = qrl(i,k)/state%pdel(i,k)
         end do
      end do

   end if   ! if (dosw .or. dolw) then

   ! output rad inputs and resulting heating rates
   call rad_data_write(pbuf, state, cam_in, coszrs)

   ! NET RADIATIVE HEATING TENDENCY
   ! INPUT: state, qrl, qrs, fsns, fsnt, flns, flnt, asdir
   ! OUTPUT: 
   !     ptend%s = (qrs + qrl)
   !     net_flx = fsnt - fsns - flnt + flns
   ! pbuf is an argument, but *is not used*
   call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                     fsnt, flns, flnt, cam_in%asdir, net_flx)

   if (write_output) then
      ! Compute heating rate for dtheta/dt
      do k = 1, pver
         do i = 1, ncol
            ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
         end do
      end do
      call outfld('HR', ftem, pcols, lchnk)
   end if

   ! [PBUF] convert radiative heating rates to Q*dp for energy conservation
   ! QR* in J/(kg s)
   do k = 1, pver
      do i = 1, ncol
         qrs(i,k) = qrs(i,k)*state%pdel(i,k)
         qrl(i,k) = qrl(i,k)*state%pdel(i,k)
      end do
   end do

   cam_out%netsw(:ncol) = fsns(:ncol)

   if (.not. present(rd_out)) then
      deallocate(rd)
   end if

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

   subroutine set_sw_diags()

      ! Transform RRTMGP output for CAM
      ! RRTMGP levels are bottom to top.  CAM levels are top to bottom.

      integer :: i
      !-------------------------------------------------------------------------

      fns  = 0._r8
      fcns = 0._r8

      do i = 1, nday
         fns(idxday(i),ktopcami:pverp)  = fsw%flux_dn(i,ktopradi:1:-1) - fsw%flux_up(i,ktopradi:1:-1)
         fcns(idxday(i),ktopcami:pverp) = fswc%flux_dn(i,ktopradi:1:-1) - fswc%flux_up(i,ktopradi:1:-1)
      end do

      call heating_rate('SW', ncol, fns, qrs)
      call heating_rate('SW', ncol, fcns, rd%qrsc)

      fsns(:ncol)       = fns(:ncol,pverp)
      fsnt(:ncol)       = fns(:ncol,1)
      rd%fsnsc(:ncol)   = fcns(:ncol,pverp)
      rd%fsntoa(:ncol)  = fsnt(:ncol)
      rd%fsntc(:ncol)   = fcns(:ncol,1)
      rd%fsntoac(:ncol) = fcns(:ncol,1)
      fsds              = 0._r8
      rd%fsdsc          = 0._r8
      rd%fsutoa         = 0._r8
      do i = 1, nday
         fsds(idxday(i))      = fsw%flux_dn(i,1)
         rd%fsdsc(idxday(i))  = fswc%flux_dn(i,1)
         rd%fsutoa(idxday(i)) = fsw%flux_up(i,nlay+1)
      end do

      ! Output fluxes at 200 mb
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns,  rd%fsn200)
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, rd%fsn200c)
      if (hist_fld_active('FSNR')) then
         do i = 1,ncol
            call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fns(i,:), rd%fsnr(i))
         end do
      end if

      if (spectralflux) then
         su  = 0._r8
         sd  = 0._r8
         do i = 1, nday
            su(idxday(i),ktopcami:pverp,:) = fsw%bnd_flux_up(i,ktopradi:1:-1,:)
            sd(idxday(i),ktopcami:pverp,:) = fsw%bnd_flux_dn(i,ktopradi:1:-1,:)
         end do
      end if

      ! Export surface fluxes
      ! The RRTMGP output isn't broken down into direct and diffuse.  For first cut
      ! put entire flux into direct and leave the diffuse set to zero.  To break the
      ! fluxes into the UV/vis and near-IR bands use the same scheme as for the albedos
      ! which is hardcoded for 14 spectral bands.
      
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)

      ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns

      ! Put half of band 9 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:

      ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron

      do i = 1, nday

         cam_out%soll(idxday(i)) = sum(fsw%bnd_flux_dn(i,1,1:8)) + &
                                   0.5_r8*fsw%bnd_flux_dn(i,1,9) + &
                                   fsw%bnd_flux_dn(i,1,14)

         cam_out%sols(idxday(i)) = 0.5_r8*fsw%bnd_flux_dn(i,1,9) + &
                                   sum(fsw%bnd_flux_dn(i,1,10:13))
                                   
      end do


   end subroutine set_sw_diags

   !-------------------------------------------------------------------------------

   subroutine set_lw_diags()

      ! Transform RRTMGP output for CAM
      ! RRTMGP levels are bottom to top.  CAM levels are top to bottom.
      !----------------------------------------------------------------------------

      fnl = 0._r8
      fcnl = 0._r8

      fnl(:ncol,ktopcami:pverp)  = flw%flux_up(:,ktopradi:1:-1)  - flw%flux_dn(:,ktopradi:1:-1)
      fcnl(:ncol,ktopcami:pverp) = flwc%flux_up(:,ktopradi:1:-1) - flwc%flux_dn(:,ktopradi:1:-1)

      call heating_rate('LW', ncol, fnl, qrl)
      call heating_rate('LW', ncol, fcnl, rd%qrlc)

      flns(:ncol) = fnl(:ncol,pverp)
      flnt(:ncol) = fnl(:ncol,1)

      rd%flnsc(:ncol) = fcnl(:ncol,pverp)
      rd%flntc(:ncol) = fcnl(:ncol,1)

      cam_out%flwds(:ncol) = flw%flux_dn(:,1)
      rd%fldsc(:ncol)      = flwc%flux_dn(:,1)

      rd%flut(:ncol)  = flw%flux_up(:,nlay+1)
      rd%flutc(:ncol) = flwc%flux_up(:,nlay+1)

      ! Output fluxes at 200 mb
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  rd%fln200)
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, rd%fln200c)
      if (hist_fld_active('FLNR')) then
         do i = 1,ncol
            call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fnl(i,:), rd%flnr(i))
         end do
      end if

      if (spectralflux) then
         lu  = 0._r8
         ld  = 0._r8
         lu(:ncol,ktopcami:pverp,:)  = flw%bnd_flux_up(:,ktopradi:1:-1,:)
         ld(:ncol,ktopcami:pverp,:)  = flw%bnd_flux_dn(:,ktopradi:1:-1,:)
      end if

   end subroutine set_lw_diags

   !-------------------------------------------------------------------------------

   subroutine heating_rate(type, ncol, flux_net, hrate)

      ! Compute heating rate as a dry static energy tendency
      
      ! arguments
      character(2), intent(in)  :: type ! either LW or SW
      integer,      intent(in)  :: ncol
      real(r8),     intent(in)  :: flux_net(pcols,pverp)  ! W/m^2
      real(r8),     intent(out) :: hrate(pcols,pver)      ! J/kg/s

      ! local vars
      integer :: k

      select case (type)
      case ('LW')

         do k = 1, pver
            ! bottom - top
            hrate(:ncol,k) = (flux_net(:ncol,k+1) - flux_net(:ncol,k)) * &
                          gravit / state%pdel(:ncol,k)
         end do

      case ('SW')

         do k = 1, pver
            ! top - bottom
            hrate(:ncol,k) = (flux_net(:ncol,k) - flux_net(:ncol,k+1)) * &
                          gravit / state%pdel(:ncol,k)
         end do

      end select

   end subroutine heating_rate

   !-------------------------------------------------------------------------------

end subroutine radiation_tend

!===============================================================================


subroutine radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)

   ! Dump shortwave radiation information to history buffer.

   integer ,               intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out

   ! local variables
   real(r8), pointer :: qrs(:,:)
   real(r8), pointer :: fsnt(:)
   real(r8), pointer :: fsns(:)
   real(r8), pointer :: fsds(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrs_idx,  qrs)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsds_idx, fsds)

   call outfld('SOLIN'//diag(icall),    rd%solin,      pcols, lchnk)

   call outfld('QRS'//diag(icall),      qrs(:ncol,:)/cpair,     ncol, lchnk)
   call outfld('QRSC'//diag(icall),     rd%qrsc(:ncol,:)/cpair, ncol, lchnk)

   call outfld('FSNT'//diag(icall),     fsnt,          pcols, lchnk)
   call outfld('FSNTC'//diag(icall),    rd%fsntc,      pcols, lchnk)
   call outfld('FSNTOA'//diag(icall),   rd%fsntoa,     pcols, lchnk)
   call outfld('FSNTOAC'//diag(icall),  rd%fsntoac,    pcols, lchnk)

   ftem(:ncol) = rd%fsntoa(:ncol) - rd%fsntoac(:ncol)
   call outfld('SWCF'//diag(icall),     ftem,          pcols, lchnk)

   call outfld('FSUTOA'//diag(icall),   rd%fsutoa,     pcols, lchnk)

   call outfld('FSNIRTOA'//diag(icall), rd%fsnirt,     pcols, lchnk)
   call outfld('FSNRTOAC'//diag(icall), rd%fsnrtc,     pcols, lchnk)
   call outfld('FSNRTOAS'//diag(icall), rd%fsnirtsq,   pcols, lchnk)

   call outfld('FSN200'//diag(icall),   rd%fsn200,     pcols, lchnk)
   call outfld('FSN200C'//diag(icall),  rd%fsn200c,    pcols, lchnk)

   call outfld('FSNR'//diag(icall),     rd%fsnr,       pcols, lchnk)

   call outfld('SOLS'//diag(icall),     cam_out%sols,  pcols, lchnk)
   call outfld('SOLL'//diag(icall),     cam_out%soll,  pcols, lchnk)
   call outfld('SOLSD'//diag(icall),    cam_out%solsd, pcols, lchnk)
   call outfld('SOLLD'//diag(icall),    cam_out%solld, pcols, lchnk)

   call outfld('FSNS'//diag(icall),     fsns,          pcols, lchnk)
   call outfld('FSNSC'//diag(icall),    rd%fsnsc,      pcols, lchnk)

   call outfld('FSDS'//diag(icall),     fsds,          pcols, lchnk)
   call outfld('FSDSC'//diag(icall),    rd%fsdsc,      pcols, lchnk)

end subroutine radiation_output_sw


!===============================================================================

subroutine radiation_output_cld(lchnk, ncol, rd)

   ! Dump shortwave cloud optics information to history buffer.

   integer ,        intent(in) :: lchnk
   integer,         intent(in) :: ncol
   type(rad_out_t), intent(in) :: rd
   !----------------------------------------------------------------------------

   call outfld('TOT_CLD_VISTAU',  rd%tot_cld_vistau,  pcols, lchnk)
   call outfld('TOT_ICLD_VISTAU', rd%tot_icld_vistau, pcols, lchnk)
   call outfld('LIQ_ICLD_VISTAU', rd%liq_icld_vistau, pcols, lchnk)
   call outfld('ICE_ICLD_VISTAU', rd%ice_icld_vistau, pcols, lchnk)
   if (cldfsnow_idx > 0) then
      call outfld('SNOW_ICLD_VISTAU', rd%snow_icld_vistau, pcols, lchnk)
   endif

end subroutine radiation_output_cld

!===============================================================================

subroutine radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out)

   ! Dump longwave radiation information to history buffer

   integer,                intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out

   ! local variables
   real(r8), pointer :: qrl(:,:)
   real(r8), pointer :: flnt(:)
   real(r8), pointer :: flns(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrl_idx,  qrl)
   call pbuf_get_field(pbuf, flnt_idx, flnt)
   call pbuf_get_field(pbuf, flns_idx, flns)

   call outfld('QRL'//diag(icall),     qrl(:ncol,:)/cpair,     ncol, lchnk)
   call outfld('QRLC'//diag(icall),    rd%qrlc(:ncol,:)/cpair, ncol, lchnk)

   call outfld('FLNT'//diag(icall),    flnt,          pcols, lchnk)
   call outfld('FLNTC'//diag(icall),   rd%flntc,      pcols, lchnk)

   call outfld('FLUT'//diag(icall),    rd%flut,       pcols, lchnk)
   call outfld('FLUTC'//diag(icall),   rd%flutc,      pcols, lchnk)
   
   ftem(:ncol) = rd%flutc(:ncol) - rd%flut(:ncol)
   call outfld('LWCF'//diag(icall),    ftem,          pcols, lchnk)

   call outfld('FLN200'//diag(icall),  rd%fln200,     pcols, lchnk)
   call outfld('FLN200C'//diag(icall), rd%fln200c,    pcols, lchnk)

   call outfld('FLNR'//diag(icall),    rd%flnr,       pcols, lchnk)

   call outfld('FLNS'//diag(icall),    flns,          pcols, lchnk)
   call outfld('FLNSC'//diag(icall),   rd%flnsc,      pcols, lchnk)

   call outfld('FLDS'//diag(icall),    cam_out%flwds, pcols, lchnk)
   call outfld('FLDSC'//diag(icall),   rd%fldsc,      pcols, lchnk)

end subroutine radiation_output_lw

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)

   ! Compute the column mean mass mixing ratio.  

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

!===============================================================================

subroutine coefs_init(coefs_file, kdist, available_gases)

   ! Read data from coefficients file.  Initialize the kdist object.

   ! arguments
   character(len=*),                  intent(in)  :: coefs_file
   class(ty_gas_optics_rrtmgp),       intent(out) :: kdist
   class(ty_gas_concs),               intent(in)  :: available_gases ! Which gases does the host model have available?

   ! local variables
   type(file_desc_t)  :: fh    ! pio file handle
   character(len=256) :: locfn ! path to actual file used

   ! File dimensions
   integer ::            &
      absorber,          &
      atmos_layer,       &
      bnd,               &
      pressure,          &
      temperature,       &
      absorber_ext,      & ! replaces `major_absorber`
      pressure_interp,   &
      mixing_fraction,   &
      gpt,               &
      temperature_Planck
   
   integer :: i, j, k
   integer :: did, vid
   integer :: ierr

   character(32), dimension(:),  allocatable :: gas_names
   integer,  dimension(:,:,:),   allocatable :: key_species
   integer,  dimension(:,:),     allocatable :: band2gpt  ! -> file : 'bnd_limits_gpt'
   real(r8), dimension(:,:),     allocatable :: band_lims_wavenum ! -> file : 'bnd_limits_wavenumber'
   real(r8), dimension(:),       allocatable :: press_ref, temp_ref
   real(r8)                                  :: press_ref_trop, temp_ref_t, temp_ref_p
   real(r8), dimension(:,:,:),   allocatable :: vmr_ref
   real(r8), dimension(:,:,:,:), allocatable :: kmajor
 ! ?  real(r8), dimension(:,:,:),   allocatable :: selfrefin, forrefin
   real(r8), dimension(:,:,:),   allocatable :: kminor_lower, kminor_upper
   real(r8), dimension(:,:),     allocatable :: totplnk
   real(r8), dimension(:,:,:,:), allocatable :: planck_frac
   real(r8), dimension(:),       allocatable :: solar_src_quiet, solar_src_facular, solar_src_sunspot  ! updated from solar_src
   real(r8), dimension(:,:,:),   allocatable :: rayl_lower, rayl_upper
   character(len=32), dimension(:),  allocatable :: gas_minor,         &
                                                    identifier_minor,  &
                                                    minor_gases_lower, &
                                                    minor_gases_upper, &
                                                    scaling_gas_lower, &
                                                    scaling_gas_upper
   integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower, &
                                                    minor_limits_gpt_upper
   ! Send these to RRTMGP as logicals,
   ! but they have to be read from the netCDF as integers
   logical, dimension(:),            allocatable :: minor_scales_with_density_lower, &
                                                    minor_scales_with_density_upper
   logical, dimension(:),            allocatable :: scale_by_complement_lower, &
                                                    scale_by_complement_upper
   integer, dimension(:), allocatable :: int2log   ! use this to convert integer-to-logical.
   integer, dimension(:),            allocatable :: kminor_start_lower, kminor_start_upper
   real(r8), dimension(:,:),         allocatable :: optimal_angle_fit
   real(r8)                                      :: tsi_default, mg_default, sb_default

   integer :: pairs, &
              minorabsorbers, &
              minor_absorber_intervals_lower, &
              minor_absorber_intervals_upper, &
              contributors_lower, &
              contributors_upper, &
              fit_coeffs

   character(len=128) :: error_msg
   character(len=*), parameter :: sub = 'coefs_init'
   !----------------------------------------------------------------------------

   ! Open file
   call getfil(coefs_file, locfn, 0)
   call cam_pio_openfile(fh, locfn, PIO_NOWRITE)

   call pio_seterrorhandling(fh, PIO_BCAST_ERROR)


   ! Get variables and validate them, then put into kdist
   
   ! Get dimensions and check for consistency with parameter values

   ierr = pio_inq_dimid(fh, 'absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorber not found')
   ierr = pio_inq_dimlen(fh, did, absorber)

   ierr = pio_inq_dimid(fh, 'atmos_layer', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': atmos_layer not found')
   ierr = pio_inq_dimlen(fh, did, atmos_layer)

   ierr = pio_inq_dimid(fh, 'bnd', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd not found')
   ierr = pio_inq_dimlen(fh, did, bnd)

   ierr = pio_inq_dimid(fh, 'pressure', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure not found')
   ierr = pio_inq_dimlen(fh, did, pressure)

   ierr = pio_inq_dimid(fh, 'temperature', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': temperature not found')
   ierr = pio_inq_dimlen(fh, did, temperature)

   ierr = pio_inq_dimid(fh, 'absorber_ext', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorber_ext not found')
   ierr = pio_inq_dimlen(fh, did, absorber_ext)

   ierr = pio_inq_dimid(fh, 'pressure_interp', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure_interp not found')
   ierr = pio_inq_dimlen(fh, did, pressure_interp)

   ierr = pio_inq_dimid(fh, 'mixing_fraction', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': mixing_fraction not found')
   ierr = pio_inq_dimlen(fh, did, mixing_fraction)

   ierr = pio_inq_dimid(fh, 'gpt', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': gpt not found')
   ierr = pio_inq_dimlen(fh, did, gpt)

   temperature_Planck = 0
   ierr = pio_inq_dimid(fh, 'temperature_Planck', did)
   if (ierr == PIO_NOERR) then
      ierr = pio_inq_dimlen(fh, did, temperature_Planck)
   end if
   ierr = pio_inq_dimid(fh, 'pair', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pair not found')
   ierr = pio_inq_dimlen(fh, did, pairs)
   ierr = pio_inq_dimid(fh, 'minor_absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber not found')
   ierr = pio_inq_dimlen(fh, did, minorabsorbers)
   ierr = pio_inq_dimid(fh, 'minor_absorber_intervals_lower', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber_intervals_lower not found')
   ierr = pio_inq_dimlen(fh, did, minor_absorber_intervals_lower)
   ierr = pio_inq_dimid(fh, 'minor_absorber_intervals_upper', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber_intervals_upper not found')
   ierr = pio_inq_dimlen(fh, did, minor_absorber_intervals_upper)
   ierr = pio_inq_dimid(fh, 'contributors_lower', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': contributors_lower not found')
   ierr = pio_inq_dimlen(fh, did, contributors_lower)
   ierr = pio_inq_dimid(fh, 'contributors_upper', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': contributors_upper not found')
   ierr = pio_inq_dimlen(fh, did, contributors_upper)

   ierr = pio_inq_dimid(fh, 'fit_coeffs', did)
   if (ierr == PIO_NOERR) then
      ierr = pio_inq_dimlen(fh, did, fit_coeffs)
   end if


   ! Get variables

   ! names of absorbing gases
   allocate(gas_names(absorber))
   ierr = pio_inq_varid(fh, 'gas_names', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': gas_names not found')
   ierr = pio_get_var(fh, vid, gas_names)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading gas_names')

   ! key species pair for each band
   allocate(key_species(2,atmos_layer,bnd))
   ierr = pio_inq_varid(fh, 'key_species', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': key_species not found')
   ierr = pio_get_var(fh, vid, key_species)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading key_species')

   ! beginning and ending gpoint for each band
   allocate(band2gpt(2,bnd))
   ierr = pio_inq_varid(fh, 'bnd_limits_gpt', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_gpt not found')
   ierr = pio_get_var(fh, vid, band2gpt)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_gpt')

   ! beginning and ending wavenumber for each band
   allocate(band_lims_wavenum(2,bnd))
   ierr = pio_inq_varid(fh, 'bnd_limits_wavenumber', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_wavenumber not found')
   ierr = pio_get_var(fh, vid, band_lims_wavenum)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_wavenumber')

   ! pressures [hPa] for reference atmosphere; press_ref(# reference layers)
   allocate(press_ref(pressure))
   ierr = pio_inq_varid(fh, 'press_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref not found')
   ierr = pio_get_var(fh, vid, press_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref')

   ! reference pressure separating the lower and upper atmosphere
   ierr = pio_inq_varid(fh, 'press_ref_trop', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref_trop not found')
   ierr = pio_get_var(fh, vid, press_ref_trop)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref_trop')

   ! temperatures [K] for reference atmosphere; temp_ref(# reference layers)
   allocate(temp_ref(temperature))
   ierr = pio_inq_varid(fh, 'temp_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': temp_ref not found')
   ierr = pio_get_var(fh, vid, temp_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading temp_ref')

   ! standard spectroscopic reference temperature [K]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_T', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_T not found')
   ierr = pio_get_var(fh, vid, temp_ref_t)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_T')

   ! standard spectroscopic reference pressure [hPa]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_P', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_P not found')
   ierr = pio_get_var(fh, vid, temp_ref_p)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_P')

   ! volume mixing ratios for reference atmosphere
   !  vmr_ref(temperature, absorber_ext, atmos_layer)
   allocate(vmr_ref(atmos_layer, absorber_ext, temperature))
   ierr = pio_inq_varid(fh, 'vmr_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': vmr_ref not found')
   ierr = pio_get_var(fh, vid, vmr_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading vmr_ref')

   ! absorption coefficients due to major absorbing gases
   allocate(kmajor(gpt,mixing_fraction,pressure_interp,temperature))
   ierr = pio_inq_varid(fh, 'kmajor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kmajor not found')
   ierr = pio_get_var(fh, vid, kmajor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kmajor')

   ! -bpm - variable wv_self & wv_for not in the newer files.
   ! ! absorption coefficients due to water vapor self continuum
   ! allocate(selfrefin(gpt,mixing_fraction,temperature))
   ! ierr = pio_inq_varid(fh, 'wv_self', vid)
   ! if (ierr /= PIO_NOERR) call endrun(sub//': wv_self not found')
   ! ierr = pio_get_var(fh, vid, selfrefin)
   ! if (ierr /= PIO_NOERR) call endrun(sub//': error reading wv_self')

   ! ! absorption coefficients due to water vapor foreign continuum
   ! allocate(forrefin(gpt,mixing_fraction,temperature))
   ! ierr = pio_inq_varid(fh, 'wv_for', vid)
   ! if (ierr /= PIO_NOERR) call endrun(sub//': wv_for not found')
   ! ierr = pio_get_var(fh, vid, forrefin)
   ! if (ierr /= PIO_NOERR) call endrun(sub//': error reading wv_for')

   ! absorption coefficients due to minor absorbing gases in lower part of atmosphere
   allocate(kminor_lower(contributors_lower, mixing_fraction, temperature))
   ierr = pio_inq_varid(fh, 'kminor_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_lower not found')
   ierr = pio_get_var(fh, vid, kminor_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_lower')

   ! absorption coefficients due to minor absorbing gases in upper part of atmosphere
   allocate(kminor_upper(contributors_upper, mixing_fraction, temperature))
   ierr = pio_inq_varid(fh, 'kminor_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_upper not found')
   ierr = pio_get_var(fh, vid, kminor_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_upper')

   ! integrated Planck function by band
   ierr = pio_inq_varid(fh, 'totplnk', vid)
   if (ierr == PIO_NOERR) then
      allocate(totplnk(temperature_Planck,bnd))
      ierr = pio_get_var(fh, vid, totplnk)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading totplnk')
   end if

   ! Planck fractions
   ierr = pio_inq_varid(fh, 'plank_fraction', vid)
   if (ierr == PIO_NOERR) then
      allocate(planck_frac(gpt,mixing_fraction,pressure_interp,temperature))
      ierr = pio_get_var(fh, vid, planck_frac)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading plank_fraction')
   end if

   ierr = pio_inq_varid(fh, 'optimal_angle_fit', vid)
   if (ierr == PIO_NOERR) then
      allocate(optimal_angle_fit(fit_coeffs, bnd))
      ierr = pio_get_var(fh, vid, optimal_angle_fit)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading optimal_angle_fit')
   end if

   ! solar_src 
   ! !bpm -- solar_source is not in file, there are solar_source_[facular, sunspot, quiet]
   ! There's a method that adds them together to get solar_source.
   ! ierr = pio_inq_varid(fh, 'solar_source', vid)
   ! if (ierr == PIO_NOERR) then
   !    allocate(solar_src(gpt))
   !    ierr = pio_get_var(fh, vid, solar_src)
   !    if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source')
   ! end if
   ierr = pio_inq_varid(fh, 'solar_source_quiet', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_quiet(gpt))
      ierr = pio_get_var(fh, vid, solar_src_quiet)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_quiet')
   end if
   ierr = pio_inq_varid(fh, 'solar_source_facular', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_facular(gpt))
      ierr = pio_get_var(fh, vid, solar_src_facular)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_facular')
   end if
   ierr = pio_inq_varid(fh, 'solar_source_sunspot', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_sunspot(gpt))
      ierr = pio_get_var(fh, vid, solar_src_sunspot)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_sunspot')
   end if

   ! +bpm also need to have tsi_default, mg_default, and sb_default
   ierr = pio_inq_varid(fh, 'tsi_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, tsi_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading tsi_default')
   end if

   ierr = pio_inq_varid(fh, 'mg_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, mg_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading mg_default')
   end if

   ierr = pio_inq_varid(fh, 'sb_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, sb_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading sb_default')
   end if

   ! rayleigh scattering contribution in lower part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_lower', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_lower(gpt,mixing_fraction,temperature))
      ierr = pio_get_var(fh, vid, rayl_lower)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_lower')
   end if

   ! rayleigh scattering contribution in upper part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_upper', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_upper(gpt,mixing_fraction,temperature))
      ierr = pio_get_var(fh, vid, rayl_upper)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_upper')
   end if

   ! +bpm the others
   allocate(gas_minor(minorabsorbers))
   ierr = pio_inq_varid(fh, 'gas_minor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': gas_minor not found')
   ierr = pio_get_var(fh, vid, gas_minor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading gas_minor')

   allocate(identifier_minor(minorabsorbers))
   ierr = pio_inq_varid(fh, 'identifier_minor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': identifier_minor not found')
   ierr = pio_get_var(fh, vid, identifier_minor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading identifier_minor')
   
   allocate(minor_gases_lower(minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'minor_gases_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_gases_lower not found')
   ierr = pio_get_var(fh, vid, minor_gases_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_gases_lower')

   allocate(minor_gases_upper(minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'minor_gases_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_gases_upper not found')
   ierr = pio_get_var(fh, vid, minor_gases_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_gases_upper')

   allocate(minor_limits_gpt_lower(pairs,minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'minor_limits_gpt_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_limits_gpt_lower not found')
   ierr = pio_get_var(fh, vid, minor_limits_gpt_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_limits_gpt_lower')

   allocate(minor_limits_gpt_upper(pairs,minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'minor_limits_gpt_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_limits_gpt_upper not found')
   ierr = pio_get_var(fh, vid, minor_limits_gpt_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_limits_gpt_upper')

   !  Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_lower))
   allocate(minor_scales_with_density_lower(minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'minor_scales_with_density_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_scales_with_density_lower not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_scales_with_density_lower')
   do i = 1,minor_absorber_intervals_lower
      if (int2log(i) .eq. 0) then
         minor_scales_with_density_lower(i) = .false.
      else
         minor_scales_with_density_lower(i) = .true.
      end if
   end do
   deallocate(int2log)

   ! Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_upper))
   allocate(minor_scales_with_density_upper(minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'minor_scales_with_density_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_scales_with_density_upper not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_scales_with_density_upper')
   do i = 1,minor_absorber_intervals_upper
      if (int2log(i) .eq. 0) then
         minor_scales_with_density_upper(i) = .false.
      else
         minor_scales_with_density_upper(i) = .true.
      end if
   end do
   deallocate(int2log)

   ! Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_lower))
   allocate(scale_by_complement_lower(minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'scale_by_complement_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scale_by_complement_lower not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scale_by_complement_lower')
   do i = 1,minor_absorber_intervals_lower
      if (int2log(i) .eq. 0) then
         scale_by_complement_lower(i) = .false.
      else
         scale_by_complement_lower(i) = .true.
      end if
   end do
   deallocate(int2log)

   ! Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_upper))
   allocate(scale_by_complement_upper(minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'scale_by_complement_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scale_by_complement_upper not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scale_by_complement_upper')
   do i = 1,minor_absorber_intervals_upper
      if (int2log(i) .eq. 0) then
         scale_by_complement_upper(i) = .false.
      else
         scale_by_complement_upper(i) = .true.
      end if
   end do
   deallocate(int2log)

   allocate(scaling_gas_lower(minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'scaling_gas_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scaling_gas_lower not found')
   ierr = pio_get_var(fh, vid, scaling_gas_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scaling_gas_lower')

   allocate(scaling_gas_upper(minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'scaling_gas_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scaling_gas_upper not found')
   ierr = pio_get_var(fh, vid, scaling_gas_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scaling_gas_upper')

   allocate(kminor_start_lower(minor_absorber_intervals_lower))
   ierr = pio_inq_varid(fh, 'kminor_start_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_start_lower not found')
   ierr = pio_get_var(fh, vid, kminor_start_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_start_lower')

   allocate(kminor_start_upper(minor_absorber_intervals_upper))
   ierr = pio_inq_varid(fh, 'kminor_start_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_start_upper not found')
   ierr = pio_get_var(fh, vid, kminor_start_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_start_upper')

   ! Close file
   call pio_closefile(fh)

   ! Initialize the gas optics class with data. The calls look slightly different depending
   !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
   ! gas_optics%load() returns a string; a non-empty string indicates an error.
   !
   if (allocated(totplnk) .and. allocated(planck_frac)) then
      error_msg = kdist%load(available_gases, gas_names, key_species, &
      band2gpt,                        &
      band_lims_wavenum,               &
      press_ref,                       &
      press_ref_trop,                  &
      temp_ref,                        &
      temp_ref_p,                      &
      temp_ref_t,                      &
      vmr_ref,                         &
      kmajor,                          &
      kminor_lower,                    &
      kminor_upper,                    &
      gas_minor,                       &
      identifier_minor,                &
      minor_gases_lower,               &
      minor_gases_upper,               &
      minor_limits_gpt_lower,          &
      minor_limits_gpt_upper,          &
      minor_scales_with_density_lower, &
      minor_scales_with_density_upper, &
      scaling_gas_lower,               &
      scaling_gas_upper,               &
      scale_by_complement_lower,       &
      scale_by_complement_upper,       &
      kminor_start_lower,              &
      kminor_start_upper,              &
      totplnk, planck_frac,            &
      rayl_lower, rayl_upper,          &
      optimal_angle_fit)
   else if (allocated(solar_src_quiet)) then
      error_msg = kdist%load(available_gases, &
      gas_names,                              &
      key_species,                            &
      band2gpt,                               &
      band_lims_wavenum,                      &
      press_ref,                              &
      press_ref_trop,                         &
      temp_ref,                               &
      temp_ref_p,                             &
      temp_ref_t,                             &
      vmr_ref,                                & 
      kmajor,                                 &
      kminor_lower,                           &
      kminor_upper,                           &
      gas_minor,                              &
      identifier_minor,                       &
      minor_gases_lower,                      &
      minor_gases_upper,                      &
      minor_limits_gpt_lower,                 &
      minor_limits_gpt_upper,                 &
      minor_scales_with_density_lower,        &
      minor_scales_with_density_upper,        &
      scaling_gas_lower,                      &
      scaling_gas_upper,                      &
      scale_by_complement_lower,              &
      scale_by_complement_upper,              &
      kminor_start_lower,                     &
      kminor_start_upper,                     &
      solar_src_quiet,                        &
      solar_src_facular,                      &
      solar_src_sunspot,                      &
      tsi_default,                            &
      mg_default,                             &
      sb_default,                             &
      rayl_lower,                             &
      rayl_upper)
   else
      error_msg = 'must supply either totplnk and planck_frac, or solar_src_[*]'
   end if

   if (error_msg /= ' ') call endrun(sub//': ERROR: '//trim(error_msg))

   deallocate( &
      gas_names, key_species,               &
      band2gpt, band_lims_wavenum,          &
      press_ref, temp_ref, vmr_ref,         &
      kmajor, kminor_lower, kminor_upper,   &
      gas_minor, identifier_minor,          &
      minor_gases_lower, minor_gases_upper, &
      scaling_gas_lower, scaling_gas_upper, & 
      minor_limits_gpt_lower,               & 
      minor_limits_gpt_upper,               &
      minor_scales_with_density_lower,      &
      minor_scales_with_density_upper,      &
      scale_by_complement_lower,            & 
      scale_by_complement_upper,            &
      kminor_start_lower, kminor_start_upper)
   if (allocated(optimal_angle_fit)) deallocate(optimal_angle_fit)
   if (allocated(totplnk))     deallocate(totplnk)
   if (allocated(planck_frac)) deallocate(planck_frac)
   if (allocated(solar_src_quiet))   deallocate(solar_src_quiet)
   if (allocated(solar_src_facular)) deallocate(solar_src_facular)
   if (allocated(solar_src_sunspot)) deallocate(solar_src_sunspot)
   if (allocated(rayl_lower))  deallocate(rayl_lower)
   if (allocated(rayl_upper))  deallocate(rayl_upper)
end subroutine coefs_init



subroutine set_available_gases(gases, gas_concentrations)
   ! This subroutine is based on the E3SM implementation. -bpm
   ! For each gas name in gases, initialize that gas in gas_concentrations.
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_rrtmgp_util_string, only: lower_case
   ! Arguments
   type(ty_gas_concs), intent(inout) :: gas_concentrations
   character(len=*),   intent(in)    :: gases(:)
   ! Local
   character(len=32), dimension(size(gases)) :: gases_lowercase
   integer :: igas
   character(len=128) :: error_msg
   ! Initialize with lowercase gas names; we should work in lowercase
   ! whenever possible because we cannot trust string comparisons in RRTMGP
   ! to be case insensitive
   do igas = 1,size(gases)
      gases_lowercase(igas) = trim(lower_case(gases(igas)))
   end do
   error_msg = gas_concentrations%init(gases_lowercase)
   if (error_msg /= '') then
      call endrun('Setting available gases. ERROR: '//trim(error_msg))
   end if
end subroutine set_available_gases


subroutine reset_fluxes(fluxes)

   use mo_rte_kind, only: wp
   use mo_fluxes_byband, only: ty_fluxes_byband
   type(ty_fluxes_byband), intent(inout) :: fluxes

   ! Reset broadband fluxes
   fluxes%flux_up(:,:) = 0._wp
   fluxes%flux_dn(:,:) = 0._wp
   fluxes%flux_net(:,:) = 0._wp
   if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) = 0._wp

   ! Reset band-by-band fluxes
   fluxes%bnd_flux_up(:,:,:) = 0._wp
   fluxes%bnd_flux_dn(:,:,:) = 0._wp
   fluxes%bnd_flux_net(:,:,:) = 0._wp
   if (associated(fluxes%bnd_flux_dn_dir)) fluxes%bnd_flux_dn_dir(:,:,:) = 0._wp

end subroutine reset_fluxes


subroutine initialize_rrtmgp_fluxes(ncol, nlevels, nbands, fluxes, do_direct)
   ! This closely follows the E3SM implementation.
   use mo_fluxes_byband, only: ty_fluxes_byband
   integer, intent(in) :: ncol, nlevels, nbands
   type(ty_fluxes_byband), intent(inout) :: fluxes
   logical, intent(in), optional :: do_direct

   logical :: do_direct_local

   if (present(do_direct)) then
      do_direct_local = .true.
   else
      do_direct_local = .false.
   end if

   ! Allocate flux arrays
   ! NOTE: fluxes defined at interfaces, so need to either pass nlevels as
   ! number of model levels plus one, or allocate as nlevels+1 if nlevels
   ! represents number of model levels rather than number of interface levels.

   ! Broadband fluxes
   allocate(fluxes%flux_up(ncol, nlevels))
   allocate(fluxes%flux_dn(ncol, nlevels))
   allocate(fluxes%flux_net(ncol, nlevels))
   if (do_direct_local) allocate(fluxes%flux_dn_dir(ncol, nlevels))

   ! Fluxes by band
   allocate(fluxes%bnd_flux_up(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_dn(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_net(ncol, nlevels, nbands))
   if (do_direct_local) allocate(fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands))

   ! Initialize
   call reset_fluxes(fluxes)

end subroutine initialize_rrtmgp_fluxes

!
! a simple clipping subroutine
!
elemental subroutine clipper(scalar, minval, maxval)
   real(r8), intent(inout) :: scalar
   real(r8), intent(in) :: minval, maxval
   if (minval < maxval) then
      if (scalar < minval) then
         scalar = minval
      end if
      if (scalar > maxval) then
         scalar = maxval
      end if
   end if
end subroutine clipper

end module radiation

