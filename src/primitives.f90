subroutine primitives(q,n_x,rho,vel,E,p,a)
  integer :: n_x
  real,dimension(n_x,3) :: q
  real :: gamma = 1.4
  real,dimension(n_x) :: rho,vel,E,p,a

  rho = q(:,1)
  vel = q(:,2)/rho
  E = q(:,3)/rho
  p = (gamma-1)*rho*(E-0.5*vel*vel)
  a = SQRT(gamma*p/rho)

end subroutine primitives
