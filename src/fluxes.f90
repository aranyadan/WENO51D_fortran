module flux
contains
  function turn(F,n_x,times)
    integer :: n_x,times
    real,dimension(n_x,3) :: F,F_new,turn
    if(times == 0) then
      turn = F
    else
      do i = 1,n_x
        if (mod(i+times,n_x)>0) then
          F_new(mod(i+times,n_x),:) = F(i,:)
        else
          F_new(n_x,:) = F(i,:)
        end if
      end do
      turn = F_new
    end if
    return
  end function turn

  function build_flux(q,n_x)
    integer :: n_x
    real :: gamma=1.4
    real, dimension(n_x,3) :: q,F,build_flux
    real, dimension(n_x) :: u,p,rho,E,a

    call primitives(q,n_x,rho,u,E,p,a)

    F(:,1) = rho*u
    F(:,2) = rho*u*u + p
    F(:,3) = rho*u*E + p*u

    build_flux=F
    return
  end function build_flux

  function get_deriv(lambda,hp,hn,q,n_x,dx)
    integer :: n_x
    real :: lambda,dx
    real, dimension(n_x,3) :: hp,hn,res1,res2,res,get_deriv,q
    res1 = 0.5 * ( build_flux(hp,n_x) + build_flux(hn,n_x) - lambda*(hn - hp) )
    res2 = 0.5 * ( build_flux((turn(hp,n_x,1)),n_x) + build_flux((turn(hn,n_x,1)),n_x) - lambda*( turn(hn,n_x,1) - turn(hp,n_x,1)) )
    res = (res1-res2)/dx

    get_deriv= res

  end function get_deriv
end module flux
