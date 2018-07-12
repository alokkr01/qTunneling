
Module mod_verlet
  implicit none

    REAL*8               :: e ,K_b
    PARAMETER ( e   = 2.71828182846d0 )
    PARAMETER ( K_b = 1.38064852d-23 )

    REAL*8  ::  DT,Time_sim, Time_elap, M
    PARAMETER ( DT       = 5d-15 )
    PARAMETER ( Time_sim = 500d-15 )
    PARAMETER ( M        = 1.67262158d-27 )

    INTEGER     I , N


    REAL*8 ::   c4, c2
    PARAMETER ( c4 = 3.353203901e21 )
    PARAMETER ( c2 = -14.85730198 )

    REAL*8 ,dimension(3) :: R,V,A
    REAL*8               :: TE,Temp , Omg
    PARAMETER ( Temp = 298.d0)
    PARAMETER ( Omg = 18.8495559215d13 )            ! (1000cm`1)

!!!..........................................!!!
!!!..........................................!!!

contains

real*8 function PE(x)
  implicit none
    real, intent(in) :: x
    PE = c4*(x)**4 +c2*(x)**2
end function PE

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)
end subroutine gaussian_random_number

subroutine trapzInt(f,a,b,n,integral)
  implicit none
     real*8 f
     real*8, intent(in)    :: a, b
     integer, intent(in)   :: n
     real*8, intent(out)   :: integral
     integer               :: k
     real*8                :: s
     s = 0
     do k=1, n
       s = s + f(a + (b-a)*k/n)
     end do
     integral = (b-a) * (s) / n
end subroutine trapzInt

real*8 function ActionIntegrand(x)
  implicit none
    real, intent(in) :: x
    print*, TE ,PE(x)
    ActionIntegrand = sqrt(2*M*((PE(x))-TE))        !! Make sure PE(x) > TE,i.e. particle is in tunneling region
end function ActionIntegrand

subroutine init()
  implicit none
    real*8 ::r1,v1
    call gaussian_random_number(r1)
    call gaussian_random_number(v1)
    R = (/r1*sqrt(K_b*Temp/M)*(1/Omg) + (3.32820446d-11), 0.d0,0.d0/)
    V = (/v1*sqrt(K_b*Temp/M)                         , 0.d0,0.d0/)
    !generate R0 and V0                            ! ***
    !R = (/3.32820446d-11,0.d0,0.d0/)
    !V = (/(1280.576419e0),0.e0,0.e0/)

end subroutine init

subroutine vloop()
  implicit none

    REAL*8        DT2, DTSQ2

    INTEGER       vSign
    real*8               :: integral
    REAL*8               :: P_tunn
    REAL*8 ,dimension(N) :: P_left
    !INTEGER     I , N
    !REAL*8 ,dimension(3) :: R,V,A
    !REAL                    TE

    !REAL :: c4, c2
    Time_elap = 0.d0

    print *, R,V
    DT2   = DT / 2.0
    DTSQ2 = DT * DT2

    vSign = V(1)/ABS(V(1))

    N = Time_sim/DT


    ! output data into a file
    !open(1, file = 'sim7.dat', status = 'new')
        DO I = 1, N
           call Acc(I)

           if (R(1) <0) then                          ! ***
             P_left(I) = P_left(I) + 1
           end if

           R = R + DT * V + DTSQ2 * A
           V = V + DT2 * A
           TE = (1.0/2.0)*M*(DOT_PRODUCT(V, V)) + c4*(R(1))**4 +c2*(R(1))**2

           Time_elap = Time_elap + DT

           print *, I ,R(1),V(1), vsign, TE , Time_elap
           !write (1,*) I ,R(1),V(1), TE , Time_elap

           !Tunneling
           if ((V(1)*vSign) <= 0 ) then
             vSign = V(1)/ABS(V(1))
             print *, "vchanged"

             call trapzInt(ActionIntegrand,R(1),-R(1),1000,integral)     ! #Ends of tunneling
             print *, integral
             P_tunn = e**(-ABS(integral))
             if(P_tunn>RAND(0)) then
               !R(1) = -*R(1)                                            ! *** TunnelingPath ,velocity
               !print *, R
               print *, "tunneled"
             end if

             !write (1,*) I ,R(1),V(1), TE , Time_elap
           end if

        ENDDO
    !close(1)
    P_left = P_left/N
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

