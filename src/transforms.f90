module transform
  implicit none

contains
  function Rinv(u,a)
    real,dimension(3,3) :: Rinv,Rn
    real,dimension(3,6) :: Rtemp
    real :: determinant,u,a,ratio,temp
    real :: gamma = 1.4
    integer :: i,j,k
    Rn = Rcalc(u,a)

  Rinv(1,1) = ((gamma - 1)/4)*u*u/(a*a) + u/(2*a)
  Rinv(2,1) = 1 - ((gamma - 1)/2) * u*u/(a*a)
  Rinv(3,1) = ((gamma - 1)/4)*u*u/(a*a) - u/(2*a)

  Rinv(1,2) = -1 * ( ((gamma - 1)/2)*u/(a*a) + 1/(2*a))
  Rinv(2,2) = (gamma - 1)*u/(a*a)
  Rinv(3,2) = -1 * ( ((gamma - 1)/2)*u/(a*a) - 1/(2*a))

  Rinv(1,3) = ((gamma - 1)/(2*a*a))
  Rinv(2,3) = -1*((gamma - 1)/(a*a))
  Rinv(3,3) = ((gamma - 1)/(2*a*a))

    ! Rinv(1,1) = 1 !
    ! Rinv(2,1) = 0 !
    ! Rinv(3,1) = 0 !
    ! Rinv(1,2) = 0 !
    ! Rinv(2,2) = 1 !
    ! Rinv(3,2) = 0 !
    ! Rinv(1,3) = 0 !
    ! Rinv(2,3) = 0 !
    ! Rinv(3,3) = 1 !


  end function Rinv


  function Rcalc(u,a)
    real,dimension(3,3) :: Rcalc
    real :: gamma = 1.4
    real :: u,a

    Rcalc(1,1) = 1
    Rcalc(2,1) = u-a
    Rcalc(3,1) = a*a/(gamma - 1) + 0.5*u*u - u*a

    Rcalc(1,2) = 1
    Rcalc(2,2) = u
    Rcalc(3,2) = 0.5*u*u

    Rcalc(1,3) = 1
    Rcalc(2,3) = u+a
    Rcalc(3,3) = a*a/(gamma - 1) + 0.5*u*u + u*a

    ! Rcalc(1,1) = 1 !
    ! Rcalc(2,1) = 0 !
    ! Rcalc(3,1) = 0 !
    ! Rcalc(1,2) = 0 !
    ! Rcalc(2,2) = 1 !
    ! Rcalc(3,2) = 0 !
    ! Rcalc(1,3) = 0 !
    ! Rcalc(2,3) = 0 !
    ! Rcalc(3,3) = 1 !


  end function Rcalc


  function u_half(q,n_x)
    integer :: n_x,i
    real,dimension(n_x,3) :: q,u_half
    ! Calculate u at half steps
    do i=1,n_x-1
      u_half(i,:) = (q(i,:) + q(i+1,:))/2.0
    end do
    u_half(n_x,:) = (q(1,:) + q(n_x,:))/2.0
  end function u_half
end module transform
