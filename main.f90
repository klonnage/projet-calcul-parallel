Program Main

  use fonctions
  use second_membre_sparse
  use syslin


  implicit none

  real(PR), dimension(:), allocatable :: AA,x,y
  integer, dimension(:), allocatable :: IA, JA
  real(PR), dimension(:), allocatable :: U, F, U0
  real(PR), dimension(:,:), allocatable :: fsource, g, h, U_Mat
  integer :: Nx, Ny, N, nbNN, k, i, imax
  real(PR) :: dx, dy, D, Lx, Ly, beta, dt, tmax, t0, t1

call lect_para(Nx, Ny, Lx, Ly, D, dt, tmax)


  N = Nx * Ny
  nbNN = 5*(Nx*Ny) - 2*Nx - 2*Ny !! nombre d'éléments non nuls 

  allocate(AA(nbNN), IA(nbNN), JA(nbNN), U(Nx*Ny), F(Nx*Ny), U0(Nx*Ny), &
        fsource(Nx,Ny), g(Nx,Ny), h(Nx,Ny), U_Mat(Nx,Ny), x(Nx), y(Ny))

  dx = Lx / Nx 
  dy = Ly / Ny

  imax = int(tmax/dt)

  beta = 0.01    ! GC, tolerance et nombre max d iteration
  k = 100

  !Initialisation : 
  U = 2.
  U0 = 2.! solution initiale 

  !génération de la matrice pentadiagonale
  call Mat(Nx,Ny,dx,dy,D,AA,JA,IA,dt)

  !conditions aux bords
  ! stationnaire
!!$  fsource = fsta(Lx,Ly,Nx,Ny)
!!$  g = gsta(Lx,Ly,Nx,Ny)
!!$  h = hsta(Lx,Ly,Nx,Ny)
  ! periodique
!!$  call fper(Lx,Ly,Nx,Ny,fsource)
!!$  g = gper(Lx,Ly,Nx,Ny)
!!$  h = hper(Lx,Ly,Nx,Ny)
  ! instationnaire periodique
  !fsource = finstaper(Lx,Ly,Nx,Ny,0.)
  g = ginstaper(Lx,Ly,Nx,Ny)
  h = hinstaper(Lx,Ly,Nx,Ny)




  ! génération du termes source, conditions aux bords + fsource
  call secondMembre(F,Nx,Ny,dx,dy,D,fsource,g,h)

  call cpu_time(t0)

  ! résolution du système AU = dt*F+U0, pour chaque pas de temps 
  do i = 1,imax
     fsource = finstaper(Lx,Ly,Nx,Ny,i*dt)
     call GC_SPARSE(AA,IA,JA,dt*F+U0,U0,U,beta,k,N)
     U0 = U
  end do

  call cpu_time(t1)
  print*, t1-t0, ' s'

  call vect_to_mat(U, U_Mat, Nx, Ny)

  do i = 1, Nx
     x(i) = i * dx
  end do
  
  do i = 1, Ny
     y(i) = i * dy
  end do

  call WriteBinary(x, y, U_Mat, "resultat.dat")
  call vect_to_mat(F, U_Mat, Nx, Ny)
  call WriteBinary(x, y, U_Mat, "test.dat")
  
  deallocate(AA,IA,JA,U,F,U0, g, h, U_Mat,x, y)

Contains

 subroutine WriteBinary(x, y, v, nom_fichier)
  
  implicit none
  
  ! arguments
  real(PR), dimension(:) :: x, y
  real(PR), dimension(:, :) :: v
  character(len=*) :: nom_fichier

  ! variables locales
  integer :: i, j, Nx, Ny
  
  open(11, file = nom_fichier)
  Nx = size(x)
  Ny = size(y)
!  write(11,*) real(x(1),kind=PR), real(x(Nx),kind=PR), Nx !xmin, xmax, Nx
!  write(11,*) real(y(1),kind=PR), real(y(Ny),kind=PR), Ny !ymin, ymax, Ny
  do i = 1, Nx
     do j = 1, Ny
        write(11,*) x(i), y(j), real(v(i, j), kind=4)
     end do
  end do

  close(11)
  
end subroutine WriteBinary

end  Program Main
