subroutine IC1DStep(u,q,n_x,x)
  implicit none
  integer :: n_x,mid_x,i
  real :: delx,p_2,p_1,v_2,v_1,rho_1,rho_2,gamma=1.4
  real,dimension(n_x,3) :: u,q
  real,dimension(n_x) :: x,E

  p_2 = 0.571 !0.1
  p_1 = 3.528

  rho_2 = 0.5 !0.125
  rho_1 = 0.445

  v_2 = 0.0
  v_1 = 0.698

  mid_x =(n_x+1)/2
  ! Set the primitive variables
  u(1:mid_x,1) = v_1
  u(1:mid_x,2) = p_1
  u(1:mid_x,3) = rho_1

  u(mid_x+1:n_x,1) = v_2
  u(mid_x+1:n_x,2) = p_2
  u(mid_x+1:n_x,3) = rho_2

  do i = 1,n_x
    q(i,1) = u(i,3);
    q(i,2) = u(i,3) * u(i,1);
    E(i) = u(i,2) / ( (gamma-1)*u(i,3) ) + 0.5 * u(i,1)**2
    q(i,3) = E(i) * u(i,3);
  end do
  ! Check IC setup
  ! open(unit=48,file='data')
  ! do i=1,n_x
  !   write(48,*)x(i),F(i,3)
  ! end do
  ! close(48)
  ! call system('gnuplot -p plotter.plt')
  return
end subroutine IC1DStep
