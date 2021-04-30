!
! DBD0D TEST CASE
! ZDPLASKIN
!
! June 2008, by S Pancheshnyi
!

program DBD0D

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
                                 gas_density     = 101325.0d0 * gas_pressure &
                                                 / gas_temperature &
                                                 / 1.38d-17  ! gas density, cm-3
!
! number of pulses
!
  integer, parameter          :: ipulses = 10
!
! filename
!
  character(*), parameter     :: file_in = 'data_in.dat', file_out = 'out.dat'

!------------------------------------------------------------------------------------------
!
! local variables
!
  integer :: i, j, idata_max
  double precision :: time, dtime, EN, NE
  double precision, allocatable :: data_in(:,:)

!------------------------------------------------------------------------------------------
!
! write
!

  write(*,'(/,A)') ' DBD0D TEST CASE'

!------------------------------------------------------------------------------------------
!
! load profiles of field and electron density from external file
!

  open(1,file=file_in)
  read(1,*)

  idata_max = 0
  do while( .true. )
    read(1,*,end=11) dtime, EN, NE
    idata_max = idata_max + 1
  enddo

11 write(*,'(x,A,i0,A)') 'found ', idata_max, ' points in ' // file_in // ' file'

  if(idata_max .le. 1) then
    write(*,*) 'wrong or missing data in the file'
    close(1)
    goto 999
  else
    allocate(data_in(idata_max,3))
    rewind(1)
    read(1,*)
    do i = 1, idata_max
      read(1,*) data_in(i,1:3)
    enddo
    close(1)
  endif

  time = data_in(1,1)

!------------------------------------------------------------------------------------------
!
! initialization of ZDPlasKin; initial densities and gas temperature
!

  call ZDPlasKin_init()
  call ZDPlasKin_set_density('N2',gas_density)
  call ZDPlasKin_set_density('E',LDENS_CONST=.true.) ! fixed electron density
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

  do i = 1, ipulses

    do j = 1, idata_max - 1

      dtime = data_in(j+1,1) - data_in(j,1)
      EN    = 0.5d0 * ( data_in(j+1,2) + data_in(j,2) )
      NE    = 0.5d0 * ( data_in(j+1,3) + data_in(j,3) )

      call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
      call ZDPlasKin_set_density('E',NE)
      call ZDPlasKin_write_qtplaskin(time)
      call ZDPlasKin_timestep(time,dtime)
      time = time + dtime

      write(1,'(99(1pe14.5))') time, EN, density(:)

    enddo

    write(*,'(99(1pe14.5))') time

  enddo

!------------------------------------------------------------------------------------------
!
! close file and end
!
  close(1)

999 continue
  if( allocated(data_in) ) deallocate(data_in)
end program DBD0D
