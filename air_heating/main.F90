!
! ZDPLASKIN
! Air reactions with pulse and afterglow with different types of gas heating
! Author: Hemaditya Malla
! Date: 28 April 2023

program main
!
! declare variables and modules
!
    use ZDPlasKin
    implicit none
    double precision, parameter :: boltzmann_const = 1.38d-23
    double precision, parameter :: elementary_charge = 1.6d-19
    double precision, parameter :: bar_to_pascal = 1.0d5  
    double precision, parameter :: ev_to_K = 1.16e4  
    double precision, parameter :: SI_to_Townsend = 1.0d21

    
    double precision :: gas_pressure    = 1.0d0 ! pressure, bar
    double precision :: gas_temperature = 300.0d0  ! temperature, K
    double precision :: gas_density     = 1.0d0 !density in kg/m3
    double precision :: field_amplitude    = 0.8d+6! Applied field in V/m
    double precision :: density_ini_elec = 1.0d13! initial electron density, cm-3
    double precision :: time  = 0.0d0, time_end = 53.0d-9, dtime = 1.0d-10! times, s
    integer          :: i, t_iter, ifile_unit
    double precision :: field, reduced_field
    double precision :: JdotE, mob, e_power

    gas_density = (bar_to_pascal*gas_pressure)/ &
        (boltzmann_const*gas_temperature)*1.0d-6 !(gas number density in cm-3)

    ! initialization of ZDPlasKin
    call ZDPlasKin_init()
    
    call get_field(field, time)
    reduced_field = field*SI_to_Townsend/(gas_density*1e6)
    !
    print *, reduced_field
    call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=reduced_field )

    ! set initial densities
    call ZDPlasKin_set_density('N2',0.8d0*gas_density, LDENS_CONST=.true.)
    call ZDPlasKin_set_density('O2',0.2d0*gas_density, LDENS_CONST=.true.)
    call ZDPlasKin_set_density('e',density_ini_elec)
    call ZDPlasKin_set_density('O2^+',density_ini_elec)

    write(*,'(4(A12))') 'Time_s', ( trim(species_name(i)), i = 1, species_max )
    write(*,'(4(1pe12.4))') time, density(:)

    ! time integration
    t_iter = 0
    do while(time .lt. time_end)
        call ZDPlasKin_timestep(time,dtime)
        time = time + dtime
        call get_field(field, time)
        reduced_field = field*SI_to_Townsend/(gas_density*1e6)
        print *, reduced_field
        call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_field, GAS_HEATING=.true.)
        call ZDPlaskin_write_qtplaskin(time)
        call ZDPlasKin_get_conditions(ELEC_MOBILITY_N=mob)
        call ZDPlasKin_get_conditions(ELEC_POWER_ELASTIC_N=e_power)
        call ZDPlasKin_get_conditions(REDUCED_FIELD=reduced_field)
        JdotE = elementary_charge*mob*density(species_electrons)*&
            reduced_field*reduced_field
        !print *, time, elementary_charge*mob*reduced_field**2, e_power*elementary_charge*gas_density
        print *, time, mob, reduced_field
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
            end if
        else
            field_val = 0.0
        end if
    end subroutine get_field
end program main
