!
! Ar plasma with diffusion and external circuit
! ZDPLASKIN
!
! June 2008, by B Eismann
!


program test_extcirc

  use options
  use ZDPlasKin
  implicit none

  integer            		      :: i, icount = 0
  double precision   		      :: time, dtime, EN, &                           ! time (s), dtime (s), EN (Td)
                      				   density_old(species_max), Vdr, &             ! density (cm-3), Vdr (cm/s)
							                   t1, t2, tc, V, J, &                          ! V (V), J(A) 
							                   species_source(species_max), all_neutral, &  ! source (cm-3), all_neutral(cm-3)
							                   source_terms(reactions_max), varstep         ! source terms (cm-3/s)
  double precision, parameter :: density_min = 1.0d5, &                       ! Min density solved
                        			   absvar = 0.1, relvar = 0.1, &                ! Timestep adjusters, knowing the variation in densities
							                   max_field_err = 1.0d-3,&                     ! max allowed variation of the field
							                   dtime_min = 1.0d-19                          ! min value for the timestep
	character(len=99) 		      :: filename
  integer, parameter          :: icount_save = 100

! Print
  write(*,'(/,A,/)') 'AR PLASMA WITH DIFFUSION AND EXTERNAL CIRCUIT'

! Initialisation of the program
  call CPU_TIME(t1)         ! to get the processing time at the end
  call set_parameters()     ! read and set the parameters defined in the condition.ini file
  call ZDPlasKin_init()     ! initialize the ZDPlasKin module (species, reaction, etc...)

! open and write the name of the output file
  call file_name(filename)  ! the name is composed of every important parameter of the run
  open(1,file=filename)

! set temperature and initial densities for given pressure
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)
  call ZDPlasKin_set_density('Ar',gas_density)
  call ZDPlasKin_set_density('e',elec_density)
  call ZDPlasKin_set_density('Ar^+',elec_density)

! initialization of variables
  time  = 0.0d0
	dtime = dtime_min
  EN    = voltage / gap_length / gas_density * 1.0d17
! Tells Bolsig the value of the reduced field (to do each time EN changes)
  call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)

! write the headers of the output file
  write(*,'(99(A14))') 'Time_s', 'dtime_s', 'EN_Td', 'Elec_cm-3'
  write(1,'(99(A14))', advance='no') 'Time', 'dTime', 'E/N', ( trim(species_name(i)), i = 1, species_max )
  write(1,'(99(i14))', advance='no') (i, i=1, reactions_max)
  write(1,'(99(A14))') 'V', 'J', 'J/S', 'theta', 'P*Lambda'

! time cycling
	do while(time .lt. time_end)

!		store previous value of the densities for variation calculation
    density_old(:)  = density(:)

!		set the NEW field at each step
    call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)

!		Central routine of the program.
!		Here the master code call ZDPlasKin which will call DVODE knowing
!		the densities, the source terms, the electron temperature, etc...
!		The results returned are given at time = time + dtime
    call ZDPlasKin_timestep(time,dtime)
		time = time + dtime

!		gives the electron drift velocity
!		(from Bolsig, knowing the field, the composition and 
!		temperature of the gas and the density of electrons)	
    call ZDPlasKin_get_conditions(ELEC_DRIFT_VELOCITY=Vdr)

!------------------------------------
!		circuit part
!------------------------------------
!		Current calculation (considering that the electrons are the only current carriers)
    J  = 1.6d-19 * gap_area * density(species_electrons) * Vdr
!		Knowing the current, deduces the voltage at the plasma
!		considering a simple voltage generator/resistance circuit
    V  = voltage - resistance * J

		call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
!500		EN(t + dt) = V0 / (d + e*R*S*mu(t)*ne(t)) 
!		with d = gap_length, e = electron charge, R = resistance, S = gap_area
!		mu(t) = electron mobility at previous time,
!		ne(t) = electron density at previous time
    EN = Voltage / ( gap_length + resistance * J / ( EN*all_neutral/1.0d17 + 1.0d-99 ) ) / all_neutral * 1.0d17
    EN = 0.5d0 * ( EN + abs(EN) ) ! non-negative
!------------------------------------

!		get the exact source terms
    call ZDPlasKin_get_rates(SOURCE_TERMS=species_source(:),REACTION_RATES=source_terms(:))

!		output
    if(mod(icount,icount_save) .eq. 0) then
!		  screen print
		  write(*,'(99(1pe14.6))') time, dtime, EN, density(species_electrons)
!		  write in the output file
      write(1,'(99(1pe14.6))',advance='no') time, dtime, EN, density(:)
			write(1,'(99(1pe14.6))',advance='no') source_terms(:)
			write(1,'(99(1pe14.6))') V, J, J/gap_area, V*J/(gap_length*gap_area*density(species_electrons)), &
			                         gas_pressure/sqrt(((2.405/radius)**2)+((3.141/gap_length)**2))
		endif

!------------------------------
!           new timestep
!------------------------------
!		Here, the purpose is to adapt the timestep so that
!		it can increase if the variations are
!		small, and is also limited by physical parameters

!		First, the relative variation of the densities
		varstep = maxval(abs(density(:)-density_old(:))/(0.5d0*(density(:)+density_old(:))+ density_min))

!		if this variation is bigger than relvar (given at the beginning),
!		this will decrease the timestep
!		while if the variation of densities is much lower than
!		the chosen relative variation (relvar), then
!		the absvar (absolute variation, lower than 1)
!		will increase the timestep
		dtime = dtime/(varstep/relvar + absvar)

!		seconde variation cheking : the field
!		the field is given line 500, and this expression combined with the limit :
!		during dt,	|dE/E|<<1,
!		it gives us a limitation for the timestep :
!		dt << V0/EN/(e*R*S*mu*|Q(e)|)
!		with |Q(e)| the absolute value of the source term for the electrons
		varstep = abs( voltage / ( resistance * J + 1.0d-99 ) - 1.0d0 ) &
                      * (density(species_electrons) + 1.0d-99) / ( abs(species_source(species_electrons)) + 1.0d-99 )

!		finally, we must ensure that dt lies in the ranges we put :
!		dtime_min < dtime << V0/EN/(e*R*S*mu*|Q(e)|)
		dtime = min(dtime,varstep*max_field_err)
		dtime = max(dtime,dtime_min)		
!------------------------------

    icount = icount + 1

  enddo 	
!	end of the time loop

  close(1)
  call CPU_TIME(t2)
	tc = (t2 - t1) / 60.
	if(tc .lt. 1.5d0) then
		tc = t2 - t1
		write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc,' (s)'
	else
		tc = floor( (t2-t1) / 60. )
		write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc,' (min)'
	endif

  write(*,'(A,$)') 'PRESS ENTER TO EXIT ...'
  read(*,*)

end program test_extcirc
