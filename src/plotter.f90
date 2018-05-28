module plotter
contains
  function save_data(q,x,n_x,t,id)
    integer :: n_x,id
    real :: t,gamma=1.4
    real,dimension(n_x,3) :: q
    real,dimension(n_x) :: u,p,rho,E,a,M,x
    character(len=5) :: charI
    character(len=1024) :: fname

    rho = q(:,1)
    u = q(:,2)/rho
    E = q(:,3)/rho
    p = (gamma-1)*rho*(E-0.5*u*u)
    a = SQRT(gamma*p/rho)
    M = u/a

    write(charI,'(I5)') id
    write(fname,'(a,i4.4,a,f6.4,a)') "./data/frame",id,"t_",t,".dat"


    open(unit = 100,file = fname)
    do i=1,n_x
      write(100,*)x(i),u(i),p(i),rho(i),M(i)
    end do
    close(100)

  end function save_data


  function plot_data(q,x,n_x,t,id,val)
    integer :: n_x,id,check,val
    real :: t
    real,dimension(n_x,3) :: q
    real,dimension(n_x) :: u,p,rho,E,a,M,x
    character(len=5) :: charI
    character(len=1024) :: property

    check = save_data(q,x,n_x,t,id)

    select case (val)
      case(1)
        property='velocity'
      case(2)
        property='pressure'
      case(3)
        property='density'
      case(4)
        property='Mach_number'
      case default
        property='y'
    end select

    ! Create the gnuplot file
    open(unit = 100,file='./plots/1plotter.plt')

    write(100,*) 'set term png'
    write(100, '(a,a,i4.4,a)' ) " set output ""./plots/",trim(property),id,".png"""
    write(100,*)'set xlabel "x"'
    write(100,*)'set ylabel "',trim(property),'"'
    write(100,'(a,i4.4,a,f6.4,a)') " m=""./data/frame",id,"t_",t,".dat"""
    ! write(100,*)'set terminal x11 0'
    write(100,*)'set nokey'
    write(100,*)'set grid'
    write(100, '(a,i4.4,a,f6.4,a)' ) " set title 'Frame ",id,", time = ",t,"'"
    write(100, '(a,i4.4,a)' ) " plot m using 1:",val+1," with linespoints"

    close(100)

    call system('gnuplot -p ./plots/1plotter.plt')
  end function plot_data



  function get_video(val)
    integer :: val
    character(len=1024) :: property,command
    select case (val)
      case(1)
        property='velocity'
      case(2)
        property='pressure'
      case(3)
        property='density'
      case(4)
        property='Mach_number'
      case default
        property='y'
    end select
    write(command,'(a,a,a,a,a)') "avconv -i ""./plots/",trim(property),"%04d.png"" -r 30 ./plots/",trim(property),".mp4"
    call system(command)
  end function get_video
end module plotter
