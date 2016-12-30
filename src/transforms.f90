module transform
contains
  function Rinv_mult(q,n_x)
    integer :: n_x,i
    real,dimension(n_x,3) :: q,q_h,v,Rinv_mult
    real,dimension(3,3) :: R_inv
    real,dimension(n_x) :: rho,u,E,p,a
    real :: temp

    ! q at half steps
    q_h = u_half(q,n_x)
    call primitives(q_h,n_x,rho,u,E,p,a)

    ! print*,"condition",u(n_x/2),a(n_x/2)
    do i=1,n_x
      R_inv = Rinv(u(i),a(i))
      ! if(i==n_x/2) then
      !   print*,R_inv
      ! end if
      v(i,1) = R_inv(1,1) * q(i,1) + R_inv(1,2) * q(i,2) + R_inv(1,3) * q(i,3)
      v(i,2) = R_inv(2,1) * q(i,1) + R_inv(2,2) * q(i,2) + R_inv(2,3) * q(i,3)
      v(i,3) = R_inv(3,1) * q(i,1) + R_inv(3,2) * q(i,2) + R_inv(3,3) * q(i,3)
    end do
    Rinv_mult = v
  end function Rinv_mult


  function R_mult(q,v,n_x)
    integer :: n_x,i,j,k
    real,dimension(n_x,3) :: q,q_h,v,R_mult,unew
    real,dimension(3,3) :: Rn
    real,dimension(n_x) :: rho,u,E,p,a
    real :: ratio,temp

    ! q at half steps
    q_h = u_half(q,n_x)
    call primitives(q_h,n_x,rho,u,E,p,a)

    do i=1,n_x
      Rn = Rcalc(u(i),a(i))
      ! if(i==n_x/2) then
      !   print*,Rn
      ! end if
      unew(i,1) = Rn(1,1) * v(i,1) + Rn(1,2) * v(i,2) + Rn(1,3) * v(i,3)
      unew(i,2) = Rn(2,1) * v(i,1) + Rn(2,2) * v(i,2) + Rn(2,3) * v(i,3)
      unew(i,3) = Rn(3,1) * v(i,1) + Rn(3,2) * v(i,2) + Rn(3,3) * v(i,3)
    end do
    R_mult = unew
  end function R_mult


  function Rinv(u,a)
    real,dimension(3,3) :: Rinv,Rn
    real,dimension(3,6) :: Rtemp
    real :: determinant
    real :: gamma = 1.4

    Rn = Rcalc(u,a)

    determinant = Rn(1,1)*(Rn(2,2)*Rn(3,3) - Rn(3,2)*Rn(2,3)) - &
                  Rn(1,2)*(Rn(2,1)*Rn(3,3) - Rn(3,1)*Rn(2,3)) + &
                  Rn(1,3)*(Rn(2,1)*Rn(3,2) - Rn(3,1)*Rn(2,2))

    if(abs(determinant)>1e-6) then
      do i=1,3
        Rtemp(i,1:3) = Rn(i,1:3)
        Rtemp(i,4:6) = (/0,0,0/)
        Rtemp(i,i+3) = 1
      end do

      do i=1,3
        do j = 1,3
          if(i/=j) then
            ratio = Rtemp(j,i) / Rtemp(i,i)
            do k=1,6
              Rtemp(j,k) = Rtemp(j,k) - ratio*Rtemp(i,k)
            end do
          end if
        end do
      end do
      do i=1,3
        temp = Rtemp(i,i)
        do j=1,6
          Rtemp(i,j) = Rtemp(i,j)/temp
        end do
      end do
      do i=1,3
        Rinv(i,:) = Rtemp(i,4:6)
      end do
    else
      Rinv(1,1) = ((gamma - 1)/4)*u*u/(a*a) + u/(2*a)
      Rinv(2,1) = 1 - ((gamma - 1)/2) * u*u/(a*a)
      Rinv(3,1) = ((gamma - 1)/4)*u*u/(a*a) - u/(2*a)

      Rinv(1,2) = -1 * ( ((gamma - 1)/a)*u/(a*a) + 1/(2*a))
      Rinv(2,2) = (gamma - 1)*u/(a*a)
      Rinv(3,2) = -1 * ( ((gamma - 1)/a)*u/(a*a) - 1/(2*a))

      Rinv(1,3) = ((gamma - 1)/(2*a*a))
      Rinv(2,3) = -1*((gamma - 1)/(a*a))
      Rinv(3,3) = ((gamma - 1)/(2*a*a))
    end if



  end function Rinv


  function Rcalc(u,a)
    real,dimension(3,3) :: Rcalc
    real :: gamma = 1.4

    ! Rcalc(1,1) = 1 !
    ! Rcalc(2,1) = 0 !
    ! Rcalc(3,1) = 0 !
    ! Rcalc(1,2) = 0 !
    ! Rcalc(2,2) = 1 !
    ! Rcalc(3,2) = 0 !
    ! Rcalc(1,3) = 0 !
    ! Rcalc(2,3) = 0 !
    ! Rcalc(3,3) = 1 !

    Rcalc(1,1) = 1
    Rcalc(2,1) = u-a
    Rcalc(3,1) = a*a/(gamma - 1) + 0.5*u*u - u*a

    Rcalc(1,2) = 1
    Rcalc(2,2) = u
    Rcalc(3,2) = 0.5*u*u

    Rcalc(1,3) = 1
    Rcalc(2,3) = u+a
    Rcalc(3,3) = a*a/(gamma - 1) + 0.5*u*u + u*a
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
