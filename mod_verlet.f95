
Module mod_verlet
  implicit none

    REAL*8  ::  DT,Time_sim, Time_elap, M
    PARAMETER ( DT       = 5d-15 )
    PARAMETER ( Time_sim = 500d-15 )
    PARAMETER ( M        = 1.67262158d-27 )

    INTEGER     I , N


    REAL*8 ::   c4, c2
    PARAMETER ( c4 = 3.353203901e21 )
    PARAMETER ( c2 = -14.85730198 )

    REAL*8 ,dimension(3) :: R,V,A
    REAL                  E

!!!..........................................!!!
!!!..........................................!!!

contains
subroutine vloop()
  implicit none

    REAL*8        DT2, DTSQ2
    !INTEGER     I , N
    !REAL*8 ,dimension(3) :: R,V,A
    !REAL                    E

    !REAL :: c4, c2
    Time_elap = 0.d0
    R = (/3.32820446d-11,0.d0,0.d0/)
    V = (/(1280.576419e0),0.e0,0.e0/)

    print *, R
    DT2   = DT / 2.0
    DTSQ2 = DT * DT2

    N = Time_sim/DT

    ! output data into a file
    open(1, file = 'sim7.dat', status = 'new')
        DO I = 1, N
           call Acc(I)

           R = R + DT * V + DTSQ2 * A
           V = V + DT2 * A
           E = (1.0/2.0)*M*(DOT_PRODUCT(V, V)) + c4*(R(1))**4 +c2*(R(1))**2

           Time_elap = Time_elap + DT
            print *, I ,R(1),V(1), E , Time_elap
           write (1,*) I ,R(1),V(1), E , Time_elap
        ENDDO
    close(1)
end subroutine vloop


subroutine Acc(I)
    implicit none

    INTEGER     I
    !REAL*8 ,dimension(3) :: R ,A

    !REAL :: c4, c2

    !REAL*8    M
    !! F = -grad(V)
    !! F = -dV_dx(rX)   !! Where a fuction named dV_dx returns the grad value at point rX.

    REAL*8 :: FX, FY, FZ
    FX = -(4*c4*(R(1))**3 + 2*c2*R(1))    ! R(1)
    FY = 0                                ! R(2)
    FZ = 0                                ! R(3)

    A = (/FX,FY,FZ/)*(1/M)

end subroutine Acc

end Module mod_verlet

