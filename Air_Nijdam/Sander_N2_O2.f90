! 
! Aram_N2_O2.f90
! N2-O2 chemistry
! 
! Created by Aram Markosyan on 7/23/13.
! Copyright (c) 2013 Aram Markosyan. All rights reserved.
!

program Aram_N2_O2
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_temperature  = 300.0d0, &  
                                 pressure         = 133.0d-3
  double precision            :: time  = 0.0d0, time_end = 1e-1, dtime = 1.0d-8 
  integer                     :: i, s

  double precision           :: oxygen_conc 
  double precision           :: density_0, external_field, reduced_field, N2_density, O2_density, e_density   

 
  oxygen_conc = 0.01;

  call ZDPlasKin_init()

  density_0 = pressure * 2.5d19  
  external_field = 0
  reduced_field = 0

  call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=external_field )

  N2_density = (1 - oxygen_conc) * density_0
  O2_density = oxygen_conc * density_0
  e_density = 1e9

  call ZDPlasKin_set_density(  'N2', N2_density)
  call ZDPlasKin_set_density(  'O2', O2_density)
  call ZDPlasKin_set_density(   'e', e_density)
  call ZDPlasKin_set_density(   'N2^+', e_density*(1 - oxygen_conc))
  call ZDPlasKin_set_density(   'O2^+', e_density*(oxygen_conc))

  do while(time < time_end)
      dtime = 1.0e-7
      reduced_field = 0
	  if (time > 1.0e-5) then
		dtime = time / 1.0e2
	  endif
	  
!	  if (time > 1.0e-4) then
!		dtime = 1.0e-6
!	  endif

!	if (time > pulse_duration) then 
!      if (density(species_max) <= 8.0e11) then
!        EXIT
!      endif 
!    endif

    call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_field)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
!    write(*,'(4(1pe12.4))') time, density(:)
  enddo

end program Aram_N2_O2
