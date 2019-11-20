Module fonctions

  !use mod_lecture

  implicit none 

  integer, parameter :: PR = 8

Contains


  subroutine vect_to_mat(u, u_mat, Nx, Ny)
    real(PR), dimension(Nx*Ny), intent(in) :: u
    real(PR), dimension(Nx,Ny), intent(out) :: u_mat
    integer, intent(in) :: Nx, Ny
    integer :: i, j, c
    c=1
    do i = 1, Nx
       do j = 1, Ny
          u_mat(i,j) = u(c)
          c = c + 1
       end do
    end do
  end subroutine vect_to_mat

  subroutine lect_para(Nx, Ny, Lx, Ly, D, dt,tmax)
    integer, intent(out) :: Nx, Ny
    real(PR), intent(out) :: Lx, Ly, D, dt, tmax

    open(10, file="parametres.dat")

    read(10, *) Nx
    read(10, *) Ny
    read(10, *) Lx
    read(10, *) Ly
    read(10, *) D 
    read(10, *) dt
    read(10, *) tmax

    close(10)

  end subroutine lect_para

  ! f,g et h  stationnaires


  function fsta(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Lx,Ly
    integer, intent(in) :: Nx,Ny
    Real(PR), dimension(Nx,Ny) :: fsta
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          fsta(i,j) = 2*(y - y**2 + x - x**2)
       end do
    end do
  end function fsta


  function gsta(Lx,Ly,Ny,Nx)
    Real(PR), Intent(in) :: Ly,Lx
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: gsta
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          gsta(i,j) = 0
       end do
    end do
  end function gsta

  function hsta(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Lx,Ly
    integer, intent(in) :: Nx,Ny
    Real(PR), dimension(Nx,Ny) :: hsta
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          hsta(i,j) = 0
       end do
    end do

  end function hsta


  ! f,g et h periodiques

  subroutine fper(Lx,Ly, Nx, Ny,A)
    Real(PR), Intent(in) :: Lx,Ly
    integer, intent(in) :: Nx,Ny
    Real(PR), dimension(Nx,Ny), intent (inout) :: A
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          A(i,j) = sin(x) + cos(y)
       end do
    end do

  end subroutine fper

  function gper(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Ly,Lx
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: gper
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          gper(i,j) = sin(x) + cos(y)
       end do
    end do

  end function gper

  function hper(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Ly,Lx
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: hper
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          hper(i,j) = sin(x) + cos(y)
       end do
    end do
  end function hper

  ! f,g et h instationnaires periodiques

  function finstaper(Lx,Ly,Nx,Ny,t) !appeler à chaque pas de temps
    Real(PR), Intent(in) :: Ly,Lx,t
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: finstaper
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          finstaper(i,j) =  exp(-(x-(Lx)/2)**2) * exp(-(y-(Ly)/2)**2)*cos((2*acos(-1.)/2)*t)
       end do
    end do

  end function finstaper

  function ginstaper(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Ly,Lx
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: ginstaper
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          ginstaper(i,j) = 0
       end do
    end do
  end function ginstaper

  function hinstaper(Lx,Ly,Nx,Ny)
    Real(PR), Intent(in) :: Ly,Lx
    integer, intent(in) :: Ny,Nx
    Real(PR), dimension(Nx,Ny) :: hinstaper
    integer :: i,j
    real(PR) :: dx,dy,x,y

    dx = Lx / Nx 
    dy = Ly / Ny
    do i = 1,Nx
       do j = 1,Ny
          x = real(i)*dx
          y = real(j)*dy
          hinstaper(i,j) = 1
       end do
    end do
  end function hinstaper

end Module fonctions
