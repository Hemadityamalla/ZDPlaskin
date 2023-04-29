!******************************************************************

    PROGRAM RUNEXAMPLE1

      USE DVODE_F90_M

      IMPLICIT NONE
      DOUBLE PRECISION ATOL, RTOL, T, TOUT, Y, RSTATS
      INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
      DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31)

      TYPE (VODE_OPTS) :: OPTIONS

      print *, "Simulation starts"
      OPEN (UNIT=6,FILE='example1.dat')
      IERROR = 0
      NEQ = 3
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      T = 0.0D0
      TOUT = 0.4D0
      RTOL = 1.D-4
      ATOL(1) = 1.D-8
      ATOL(2) = 1.D-14
      ATOL(3) = 1.D-6
      ITASK = 1
      ISTATE = 1
      OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR=RTOL, &
        USER_SUPPLIED_JACOBIAN=.TRUE.)
      DO IOUT = 1, 12
        CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
        CALL GET_STATS(RSTATS,ISTATS)
        WRITE (6,90003) T, Y(1), Y(2), Y(3)
        DO I = 1, NEQ
          IF (Y(I)<0.0D0) IERROR = 1
        END DO
        IF (ISTATE<0) THEN
          WRITE (6,90004) ISTATE
          STOP
        END IF
        TOUT = TOUT*10.0D0
      END DO
      WRITE (*,*) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(19), &
        ISTATS(20), ISTATS(21), ISTATS(22)
      IF (IERROR==1) THEN
        WRITE (*,90001)
      ELSE
        WRITE (*,90002)
      END IF
90000 FORMAT (/'  No. steps =',I4,'   No. f-s =',I4,'  No. J-s =',I4, &
        '   No. LU-s =',I4/'  No. nonlinear iterations =', &
        I4/'  No. nonlinear convergence failures =', &
        I4/'  No. error test failures =',I4/)
90001 FORMAT (/' An error occurred.')
90002 FORMAT (/' No errors occurred.')
90003 FORMAT (' At t =',D12.4,'   y =',3D14.6)
90004 FORMAT (///' Error halt: ISTATE =',I3)
      STOP
      contains

      SUBROUTINE FEX(neq,T,Y,YDOT)
        IMPLICIT NONE
        integer, intent(in) :: neq
        double precision , intent(in) :: T, Y(neq)
        double precision, intent(out) :: YDOT(neq)

        YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
        YDOT(3) = 3.E7*Y(2)*Y(2)
        YDOT(2) = -YDOT(1) - YDOT(3)
        RETURN
      END SUBROUTINE FEX

      SUBROUTINE JEX(neq,t,y,ml,mu,pd,nrpd)
        implicit none
        integer,          intent(in)  :: neq, ml, mu, nrpd
        double precision, intent(in)  :: t, y(neq)
        double precision, intent(out) :: pd(nrpd,neq)

        PD(1,1) = -.04D0
        PD(1,2) = 1.D4*Y(3)
        PD(1,3) = 1.D4*Y(2)
        PD(2,1) = .04D0
        PD(2,3) = -PD(1,3)
        PD(3,2) = 6.E7*Y(2)
        PD(2,2) = -PD(1,2) - PD(3,2)
        RETURN
      END SUBROUTINE JEX
END PROGRAM RUNEXAMPLE1
