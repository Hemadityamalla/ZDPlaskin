program SimpleAir

  use ZDPlasKin

  implicit none

!------------------------------------------------------------------------------------------
!
! configuration
!
!------------------------------------------------------------------------------------------
!
! gas pressure, temperature and density
!
  double precision, parameter :: gas_pressure    = 1.0d0,  &  ! pressure, bar
                                 gas_temperature = 300.0d0,&  ! temperature, K
                                 gas_density = 101325.0d0 * gas_pressure &
                                                 / gas_temperature &
                                                 / 1.38d-17,  & ! gas density, cm-3
                                 O2_percent      = 0.02d0

! filename
!
  character(*), parameter     :: file_in = 'Afivo_data.dat', file_out = 'out.dat'

!------------------------------------------------------------------------------------------
!
! local variables
!
  integer :: i, j, idata_max
  double precision :: dtime, time, e_density, electric_fld, photo, &
                      N4_plus, N3_min, N2_plus, &
                      O2_min, O4_plus, O_min, O2_plus
  integer, parameter :: n_columns = 11
  double precision, allocatable :: data_in(:,:)

!------------------------------------------------------------------------------------------
!
! write
!

  write(*,'(/,A)') ' SimpleAir TEST CASE'

!------------------------------------------------------------------------------------------
!
! load profiles of field and electron density from external file
!

  open(1,file=file_in)
  read(1,*)

  idata_max = 0
  do while( .true. )
    read(1,*,end=11) time, e_density, electric_fld, photo, N4_plus, N3_min, &
                     O2_min, O4_plus, O_min, O2_plus, N2_plus
    idata_max = idata_max + 1
  enddo

11 write(*,'(x,A,i0,A)') 'found ', idata_max, ' points in ' // file_in // ' file'

  if(idata_max .le. 1) then
    write(*,*) 'wrong or missing data in the file'
    close(1)
    goto 999
  else
    allocate(data_in(idata_max,n_columns))
    rewind(1)
    read(1,*)
    do i = 1, idata_max
      read(1,*) data_in(i,1:n_columns)
    enddo
    close(1)
  endif

  time = data_in(1,1)


!------------------------------------------------------------------------------------------
!
! initialization of ZDPlasKin; initial densities and gas temperature
!

  !Converting the electric_fld to Td units
  data_in(:,3) = (data_in(:,3)*1.d-02/gas_density)/1.0d-17
  !Converting the densities to cm-3
  data_in(:,2) = data_in(:,2)*1.d-06
  data_in(:,4:n_columns) = data_in(:,4:n_columns)*1.d-06 
 
  call ZDPlasKin_init()
  call ZDPlasKin_set_density(  'N2', (1.0d0-O2_percent)*gas_density)
  call ZDPlasKin_set_density(  'O2', O2_percent*gas_density)
  call ZDPlasKin_set_density(   'e', data_in(1, 2), LDENS_CONST=.true.)
  call ZDPlasKin_set_density(   'N2^+', data_in(1, 10))
  call ZDPlasKin_set_density(   'O2^+', 0.0d0)
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)

!------------------------------------------------------------------------------------------
!
! open file for output and write headers
!

  open(1,file=file_out)
  write(1,'(99A14)') 'TIME', 'E/N', ( trim(species_name(i)), i = 1, species_max )
  write(*,'(99A14)') 'TIME'

!------------------------------------------------------------------------------------------
!
! time integration
!
!------------------------------------------------------------------------------------------

    do j = 1, idata_max - 1
      

      dtime = data_in(j+1,1) - data_in(j,1)
      electric_fld    = 0.5d0 * ( data_in(j+1,3) + data_in(j,3) )
      e_density    = 0.5d0 * ( data_in(j+1,2) + data_in(j,2) )

      call ZDPlasKin_set_conditions(REDUCED_FIELD=electric_fld)
      call ZDPlasKin_set_density('E',e_density)

      call ZDPlasKin_timestep(time,dtime)
      time = time + dtime

      write(1,'(99(1pe14.5))') time, electric_fld, density(:)

      write(*,'(99(1pe14.5))') time
    enddo

    

!------------------------------------------------------------------------------------------
!
! close file and end
!
  close(1)

999 continue
  if( allocated(data_in) ) deallocate(data_in)
end program SimpleAir
