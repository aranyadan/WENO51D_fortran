subroutine WENO51d(lambda,F,q,dx,n_x,hp,hn)
  use flux
  use transform
  implicit none
  integer :: n_x,i,j
  real,dimension(n_x,3) :: F,q,hp,hn,Ftemp,qtemp
  real,dimension(-1:n_x+3,3) :: Fnew,qnew
  real :: lambda,dx
  real,dimension(1) :: rho,u,E,p,a
  real,dimension(6,3) :: F_i,q_i,g_i,v_i,F_i2,q_i2
  real,dimension(3,3) :: R_i,Rinv_i
  real,dimension(1,3) :: q_h,hpr,hnr

  Fnew(1:n_x,:) = F(1:n_x,:)
  qnew(1:n_x,:) = q(1:n_x,:)

  qnew(-1,:) = qnew(1,:); qnew(0,:) = qnew(1,:)
  qnew(n_x+1,:) = qnew(n_x,:); qnew(n_x+2,:) = qnew(n_x,:); qnew(n_x+3,:) = qnew(n_x,:)

  Fnew(-1,:) = Fnew(1,:); Fnew(0,:) = Fnew(1,:)
  Fnew(n_x+1,:) = Fnew(n_x,:); Fnew(n_x+2,:) = Fnew(n_x,:); Fnew(n_x+3,:) = Fnew(n_x,:)

  do i=1,n_x
    F_i(1:6,:) = Fnew(i-2:i+3,:)
    q_i(1:6,:) = qnew(i-2:i+3,:)
    q_h(1,:) = (qnew(i,:) + qnew(i+1,:))/2.0

    call primitives(q_h,1,rho,u,E,p,a)
    R_i = Rcalc(u(1),a(1))
    Rinv_i = Rinv(u(1),a(1))

    call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
    hp(i,:) = hpr(1,:)
    hn(i,:) = hnr(1,:)
  end do


  return
end subroutine WENO51d
