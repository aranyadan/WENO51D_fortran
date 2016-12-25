subroutine WENO51d(lambda,q,dx,n_x,hp,hn)
  use flux
  integer :: n_x
  real,dimension(n_x,3) :: q,F
  real :: lambda
  real,dimension(n_x,3) :: u,v,umm,um,up,upp,vmm,vm,vp,vpp,p0n,p1n,p2n,B0n,B1n
  real,dimension(n_x,3) :: B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,hn
  real,dimension(n_x,3) :: p0p,p1p,p2p,B0p,B1p,B2p,alpha0p,alpha1p,alpha2p
  real,dimension(n_x,3) :: alphasump,w0p,w1p,w2p,hp,res,temp,res1,res2
  real :: d0n,d0p,d1n,d1p,d2n,d2p
  ! Build the flux
  F = build_flux(q,n_x)

  ! Lax Friedrich's Flux splitting
  ! v = 0.5*(F + lambda*q)
  ! temp = 0.5 *(F - lambda*q)
  ! u = turn(temp,n_x,-1)
  v=q
  u=turn(q,n_x,-1)

  ! Right flux
  ! compute u_{i+1/2}^{+}
  vmm = turn(v,n_x,2)
  vm = turn(v,n_x,1)
  vp = turn(v,n_x,-1)
  vpp = turn(v,n_x,-2)

  ! Polynomials
  p0p = (2.0*vmm - 7.0*vm + 11.0*v)/6
  p1p = (-vm + 5.0*v + 2.0*vp)/6
  p2p = (2.0*v + 5.0*vp - vpp)/6

  ! Smoothness indicators
  B0p = 13/12*(vmm-2*vm+v  )**2 + 1/4*(vmm-4*vm+3*v)**2
  B1p = 13/12*(vm -2*v +vp )**2 + 1/4*(vm-vp)**2
  B2p = 13/12*(v  -2*vp+vpp)**2 + 1/4*(3*v-4*vp+vpp)**2

  !constants
  d0p = 1.0/10
  d1p = 6.0/10
  d2p = 3.0/10
  epsilon = 1E-6

  !alpha weights
  alpha0p = d0p/(epsilon + B0p)**2
  alpha1p = d1p/(epsilon + B1p)**2
  alpha2p = d2p/(epsilon + B2p)**2
  alphasump = alpha0p + alpha1p + alpha2p

  ! ENO weights
  w0p = alpha0p/alphasump
  w1p = alpha1p/alphasump
  w2p = alpha2p/alphasump

  ! Flux at the boundary u_{i+1/2}^{+}
  hp = w0p*p0p + w1p*p1p + w2p*p2p


  ! Left Flux
  ! compute u_{i+1/2}^{-}
  umm = turn(u,n_x,2)
  um = turn(u,n_x,1)
  up = turn(u,n_x,-1)
  upp = turn(u,n_x,-2)

  ! Polynomials
  p0n = ( -umm + 5*um + 2*u  )/6
  p1n = ( 2*um + 5*u  - up   )/6
  p2n = (11*u  - 7*up + 2*upp)/6

  ! Smoothness indicators
  B0n = 13/12*(umm-2*um+u  )**2 + 1/4*(umm-4*um+3*u)**2
  B1n = 13/12*(um -2*u +up )**2 + 1/4*(um-up)**2
  B2n = 13/12*(u  -2*up+upp)**2 + 1/4*(3*u -4*up+upp)**2

  !constants
  d0n = 3.0/10
  d1n = 6.0/10
  d2n = 1.0/10

  ! Alpha Weights
  alpha0n = d0n/(epsilon + B0n)**2
  alpha1n = d1n/(epsilon + B1n)**2
  alpha2n = d2n/(epsilon + B2n)**2
  alphasumn = alpha0n + alpha1n + alpha2n

  ! ENO weights
  w0n = alpha0n/alphasumn
  w1n = alpha1n/alphasumn
  w2n = alpha2n/alphasumn

  ! Numerical Flux at boundary u_{i+1/2}^{-}
  hn = w0n*p0n + w1n*p1n + w2n*p2n

  return
end subroutine WENO51d
