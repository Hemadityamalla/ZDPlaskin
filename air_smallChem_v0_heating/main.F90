!
! Simple air reactions test case
! ZDPLASKIN
!
!
!

program main
!
! declare variables and modules
!
  use ZDPlasKin
  implicit none
  double precision, parameter :: boltzmann_const = 1.38d-23
  double precision, parameter :: bar_to_pascal = 1.0d5  

  double precision :: gas_pressure    = 1.0d0 ! pressure, bar
  double precision :: gas_temperature = 300.0d0  ! temperature, K
  double precision :: gas_density     = 1.0d0 !density in kg/m3
  double precision :: field_amplitude    = 4d+6! Applied field in V/m
  double precision :: SI_to_Townsend = 1.0d21
  double precision :: density_ini_elec = 1.0d13     ! initial electron density, cm-3
  double precision            :: time  = 0.0d0, time_end = 2.5d-8, dtime = 1.0d-12! times, s
  integer                     :: i, t_iter, Nsteps, ifile_unit
  double precision :: field, reduced_field
  double precision, allocatable, dimension(:,:) :: outputData
!
! print
!
  !Unit conversion
  gas_density = (bar_to_pascal*gas_pressure)/(boltzmann_const*gas_temperature)*1.0d-6 !(gas number density in cm-3)
!
! initialization of ZDPlasKin
!
  call ZDPlasKin_init()

  Nsteps = int((time_end - time)/dtime)
! set the physical conditions of the calculation:
!     the gas temperature and the reduced electric field
 call get_field(field, time)
  reduced_field = field*SI_to_Townsend/(gas_density*1e6)
!
  print *, reduced_field
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=reduced_field, GAS_HEATING=.true. )
!
! set initial densities
!
  call ZDPlasKin_set_density(  'N2',0.8d0*gas_density, LDENS_CONST=.true.)
  call ZDPlasKin_set_density(  'O2',0.2d0*gas_density, LDENS_CONST=.true.)
  call ZDPlasKin_set_density(   'e',density_ini_elec)
  call ZDPlasKin_set_density('O2^+',density_ini_elec)
! print column headers and initial values
!
  write(*,'(4(A12))') 'Time_s', ( trim(species_name(i)), i = 1, species_max )
  write(*,'(4(1pe12.4))') time, density(:)
!
! time integration
  t_iter = 0
  do while(time .lt. time_end)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
    !call get_field(field, time)
    !reduced_field = field*SI_to_Townsend/(gas_density*1e6)
    reduced_field = field_amplitude*SI_to_Townsend/(gas_density*1e6)
    call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_field)
    call ZDPlasKin_set_config(BOLSIG_IGNORE_GAS_TEMPERATURE=.true.)
    call ZDPlaskin_write_qtplaskin(time)
    !write(*,'(4(1pe12.4))') time, density(:)
    !print *,time
    t_iter = t_iter+1
    enddo
contains
        subroutine get_field(field_val, t)
                implicit none
                double precision, intent(inout) :: field_val
                double precision, intent(in) :: t
                double precision :: trise = 3.0d-9
                double precision :: tfall = 6.0d-9
                double precision :: tinter = 50d-9
                double precision :: tconst = 8d-9
                if (t <= (trise)) then
                   field_val = field_amplitude*(t/trise)
                else if (t <= (trise + tconst)) then
                   field_val = field_amplitude
                else if (t <= (trise+tfall+tconst)) then
                   field_val = field_amplitude*(1 - (t - (trise+tconst))/tfall)
                else if (t <= (trise+tfall+tconst+tinter)) then
                   field_val = 0.0d0
                end if

        end subroutine get_field
end program main
