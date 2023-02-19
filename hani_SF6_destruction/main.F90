program pulse_air

  use ZDPlasKin

  implicit none

  double precision, parameter :: boltzmann_const = 1.38d-23
  double precision, parameter :: bar_to_pascal = 1.0d5  

  double precision:: gas_pressure    = 3.32e-3 ! pressure, bar
  double precision :: gas_temperature = 300.0d0  ! temperature, K
  double precision:: gas_density     = 1.0d0 !density in kg/m3
  double precision :: pulse_amplitude = 2445.0d0 !Pulse amplitude in V/m
  double precision :: e_field
  double precision :: time, dt, t_final = 7.5e-2

  double precision :: e_density = 3.46e8 !starting electron density in cm-3
  double precision :: sf6_density = 6.41e5 !starting electron density in cm-3
  !double precision :: o2plus_density = 1.0d14 !starting positive ion density in cm-3


  !Unit conversion
  gas_density = (bar_to_pascal*gas_pressure)/(boltzmann_const*gas_temperature)*1.0d-6 !(gas number density in cm-3)
  !Converting field from V/m to Td
  pulse_amplitude = 1d17*(pulse_amplitude)*(1e-2/gas_density)

  print *, gas_density, pulse_amplitude
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! initialization of ZDPlasKin; initial densities and gas temperature
  call ZDPlasKin_init()
  call ZDPlasKin_set_density('N2',0.8d0*gas_density)
  call ZDPlasKin_set_density("O2", 0.2d0*gas_density)
  call ZDPlasKin_set_density('E',e_density)
  call ZDPlasKin_set_density('SF6',sf6_density)
  !call ZDPlasKin_set_density("O2^+", o2plus_density)
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature, GAS_HEATING = .true.)


!------------------------------------------------------------------------------------------
!Time integration
  dt = 0.1d-9 !TODO: Make this based on the dielectric relaxation time
  do while (time .le. t_final)
    call ZDPlasKin_set_conditions(REDUCED_FIELD=pulse_amplitude)
    print *, "Time: ", time
    call ZDPlasKin_timestep(time, dt) 
    call ZDPlaskin_write_qtplaskin(time)
  time = time + dt
  end do

  contains
    function pulse(time, field_max) result(field_value)
    implicit none
    double precision, intent(in) :: time, field_max
    double precision ::  field_value
    double precision :: t_rise, t_fall, t_pulse

    t_rise = 3.0d-9
    t_fall = 3.0d-9
    t_pulse = 100.0d-9
    
    if (time .le. t_rise) then
      field_value = (field_max/t_rise)*time
    else if (time .le. (t_rise + t_pulse) .and. time .ge. t_rise) then
      field_value = field_max
    else if (time .le. ((t_rise + t_pulse + t_fall)) .and. time .ge. (t_rise + t_pulse)) then
      field_value = field_max - (field_max/t_fall)*(time - (t_rise + t_pulse))
    else
      field_value = 1.3d0
    end if

    
    end function pulse

end program pulse_air
