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
                                 pressure         = 100.0d-3, & 
                                 domain_size      = 1.3d1, &    
                                 voltage          = 8.7d0, &    
                                 pulse_duration   = 250.0d-9, & 
                                 min_field_stp    = 5, &        
                                 max_field_stp    = 150        
  double precision            :: time  = 0.0d0, time_end = 1, dtime = 1.0d-8 
  integer                     :: i, s

  double precision           :: oxygen_conc, oxygen_procent, max_field, min_field, permetivity, charge, mobility
  double precision           :: density_0, external_field, reduced_field, N2_density, O2_density, e_density   
  double precision           :: N2p_density, O2p_density
  double precision, dimension(46)  :: eField, mobil 

  eField = (/ 0.1, 0.13, 0.17, 0.21, 0.27, 0.35, 0.46, 0.59, &
            & 0.77, 1.0, 1.3, 1.7, 2.1, 2.7, 3.5, 4.6, &
            & 5.9, 7.7, 10.0, 13.0, 17.0, 21.0, 27.0, 35.0, &
            & 46.0, 59.0, 77.0, 100.0, 130.0, 170.0, 210.0, 270.0, &
            & 350.0, 460.0, 590.0, 770.0, 1000.0, 1300.0, 1700.0, 2100.0,&
            & 2700.0, 3500.0, 4600.0, 5900.0, 7700.0, 10000.0 /);

  mobil = (/ 0.81671, 0.6742, 0.54323, 0.45804, 0.37889, &
            & 0.3175, 0.26759, 0.23054, 0.19731, 0.17029, 0.1484, 0.13081, 0.11972, &
            & 0.10883, 0.09944, 0.0909, 0.08397, 0.07729, 0.07151, 0.06644, 0.06208, &
            & 0.05903, 0.05583, 0.05274, 0.04968, 0.04694, 0.04426, 0.0421, 0.04031, &
            & 0.03867, 0.03741, 0.03577, 0.03401, 0.03198, 0.02994, 0.02787, 0.02587, &
            & 0.02391, 0.02204, 0.02068, 0.01924, 0.01797, 0.01689, 0.01613, 0.0155, 0.01497 /);

  permetivity = 8.85418782d-12;
  charge = 1.60217646d-19;

  oxygen_conc = 0.01
  
  oxygen_procent = 1.0d2 * oxygen_conc

  call ZDPlasKin_init()

  density_0 = pressure * 2.5d19  
  external_field = (1.0d3 * voltage) / domain_size / density_0 / 1.0d-17  
  max_field = 1.0e3 * max_field_stp * pressure / density_0 / 1.0d-17
  min_field = 1.0d3 * min_field_stp * pressure / density_0 / 1.0d-17
  reduced_field = max_field

  !call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=external_field )

  N2_density = (1.0d2 - oxygen_procent) * density_0 /1.0d2
  O2_density = oxygen_procent * density_0 / 1.0d2
  e_density = 0.4e14 * (pressure * pressure)   
  N2p_density = e_density * (1.0d2 - oxygen_procent) /1.0d2
  O2p_density = e_density * oxygen_procent / 1.0d2

  call ZDPlasKin_set_density(  'N2', N2_density)
  call ZDPlasKin_set_density(  'O2', O2_density)
  call ZDPlasKin_set_density(   'e', e_density)
  call ZDPlasKin_set_density(   'N2^+', N2p_density)
  call ZDPlasKin_set_density(   'O2^+', O2p_density)
  
  
  do while(time .lt. time_end)
    if (reduced_field > min_field) then
      dtime = 1.0e-11
      e_density = density(species_max)
      mobility = pressure/1.01325 * INTERP(reduced_field, eField, mobil)
      reduced_field = reduced_field - dtime * charge/permetivity * mobility * reduced_field * e_density * 1.0d6
    else
      dtime = 1.0e-9
      reduced_field = min_field
    endif


    if (time .ge. pulse_duration) then
      if (time .le. (pulse_duration + 50.0e-9) ) then
        dtime = 1.0e-10
        reduced_field = min_field - min_field/50.0e-9 * (time - pulse_duration)
       
        if (reduced_field .le. 0.0) then
        reduced_field = 0.0
        endif
       
      else  
		if (time > 1.0e-5) then
			dtime = time / 1.0e2
		endif
        
        reduced_field = 0.0
      endif
    endif 
 

    call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_field)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
 !   write(*,'(4(1pe12.4))') time, density(:)
  enddo

  CONTAINS

  FUNCTION INTERP(xe,xes,fes)
  double precision, dimension(46):: xes, fes
  double precision :: xe, fe
  double precision :: INTERP
  integer i

  if (xe < xes(1)) then
    INTERP = 0;  RETURN
  end if

  if (xe == xes(1)) then
    INTERP = fes(1);  RETURN
  end if

  if (xe .GE. xes(46)) then
    INTERP = fes(46);  RETURN
  end if

  i = 1;
  do while(i .LE. 46)
    if (xes(i) .GE. xe) EXIT
    i = i+1;
  end do
  INTERP = fes(i-1)+(fes(i)-fes(i-1))*(xe-xes(i-1))/(xes(i)-xes(i-1));
  END FUNCTION INTERP

end program Aram_N2_O2
