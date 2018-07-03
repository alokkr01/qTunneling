
program main
  implicit none

  integer n_dim
  real, dimension(n_dim) :: x
  real, V

  ! k_boltz ;
  ! mass,freq--->omg, Temp
  ! KT=k_boltz*Temp ; MW2=mass*(omg)**2

  subroutine derv(x,V,dV_dx,dV_dx2)     !
  implicit none
     ! Defining potential energy function
      V = c4*(x)**4 +c2*(x)**2
      dV_dx  = 4*c4*(x)**3 + 2*c2*x
      dV_dx2 = 12*c4*(x)**2 + 2*c2
  end subroutine derv

  subroutine V_endpts()
  implicit none
     ! Solving for V' and V''
      x_min(1) =  (c2/(2.0*c4))**(0.5)      ! 3.32820446*(10**(-11))
      x_min(2) = -(c2/(2.0*c4))**(0.5)
      x_max(1) =  0

     ! Solving for V''(x_min) = mass*(omg)**2
      C = (1/8)*m*(omg)**2
  end subroutine V_endpts

  subroutine pot_barrier()
  implicit none
     ! Solving V(x_max) - V(x_min) = V_barr(= k_boltz*Temp)
      factr = (KT)/(MW2)
      A = (MW2)/factr*(1+(1- 128*factr)**(0.5))
  end subroutine pot_barrier

end program main
