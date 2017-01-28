subroutine WENO51d(lambda,F,q,dx,n_x,hp,hn)
  use flux
  use transform
  implicit none
  integer :: n_x,i,j
  real,dimension(n_x,3) :: F,q,hp,hn,Ftemp,qtemp
  real :: lambda,dx
  real,dimension(1) :: rho,u,E,p,a
  real,dimension(6,3) :: F_i,q_i,g_i,v_i
  real,dimension(3,3) :: R_i,Rinv_i
  real,dimension(1,3) :: q_h,hpr,hnr

  do i=1,n_x
    do j = 6,1,-1
      Ftemp = turn(F,n_x,3-j)
      qtemp = turn(q,n_x,3-j)
      F_i(j,:) = Ftemp(i,:)
      q_i(j,:) = qtemp(i,:)
    end do

    if(i == n_x) then
      q_h(1,:) = (q(i,:) + q(1,:))/2.0
    else
      q_h(1,:) = (q(i,:) + q(i+1,:))/2.0
    end if

    call primitives(q_h,1,rho,u,E,p,a)
    R_i = Rcalc(u(1),a(1))
    Rinv_i = Rinv(u(1),a(1))

    call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
    hp(i,:) = hpr(1,:)
    hn(i,:) = hnr(1,:)
  end do


  return
end subroutine WENO51d
