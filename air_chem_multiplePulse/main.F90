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
  double precision :: field_amplitude    = 0.8d+6! Applied field in V/m
  double precision :: SI_to_Townsend = 1.0d21
  double precision :: density_ini_elec = 1.0d13     ! initial electron density, cm-3
  double precision            :: time  = 0.0d0, time_end = 53.0d-9, dtime = 1.0d-9! times, s
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
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=reduced_field )
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
    !call ZDPlasKin_timestep(time,dtime)
    call ZDPlasKin_timestep_explicit(time,dtime, 1.0d3, 1.0d10)
    time = time + dtime
    call get_field(field, time)
    reduced_field = field*SI_to_Townsend/(gas_density*1e6)
    call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_field, GAS_HEATING=.true.)
    !call ZDPlaskin_write_qtplaskin(time)
    call ZDPlaskin_set_config(QTPLASKIN_SAVE=.true.)
    !write(*,'(4(1pe12.4))') time, density(:)
    !print *,time,dtime
    t_iter = t_iter+1
    enddo
contains
        subroutine get_field(field_val, time)
                implicit none
                double precision, intent(inout) :: field_val
                double precision, intent(in) :: time
                double precision :: trise = 5.0d-9
                double precision :: tfall = 5.0d-9
                double precision :: tinter = 150d-9
                double precision :: tconst = 5.0d-9
                double precision :: tperiod, tmp, t
                integer :: npulses = 5
                
                t = time
                tperiod = trise+tfall+tconst+tinter
                if (t < tperiod*npulses) then
                        t = modulo(t, tperiod)
                        if (t <= (trise)) then
                           field_val = field_amplitude*(t/trise)
                        else if (t <= (trise + tconst)) then
                           field_val = field_amplitude
                        else
                           tmp = max(0.01, (1 - (t - (trise+tconst))/tfall))
                           field_val = field_amplitude*tmp
                                !tmp = trise+tconst+tfall - t
                                !field_val = field_val*max(100.0, tmp/tfall)
                                !print *, "Fall region", field_val
                        !else if (t <= (trise+tfall+tconst)) then
                        !   field_val = field_amplitude*(1 - (t - (trise+tconst))/tfall)
                        !else if (t <= (trise+tfall+tconst+tinter)) then
                        !   field_val = field_amplitude/1.0d2
                        end if
                else
                        field_val = 0.0
                end if

        end subroutine get_field
end program main
