Program Main
  
  Use fonctions
  Use second_membre_sparse
  Use syslin
  Use mod_parallele

  Implicit None

  !include "mpif.h"
  integer :: statinfo
  Integer, Dimension(MPI_STATUS_SIZE) :: status

  Real(PR), Dimension(:), Allocatable :: AA, x, y, AAme
  Integer, Dimension(:), Allocatable :: IA, JA, tab_i1, tab_iN, IAme, JAme
  Real(PR), Dimension(:), Allocatable :: Ume, F, U0me, S, Fme, b, U, U0
  Real(PR), Dimension(:,:), Allocatable :: fsource, g, h, U_Mat
  Integer :: Nx, Ny, N, nbNN, kmax, i, imax, me, i1, iN, Np, Nme, NbLiMe, NbColMe, NbLiAv, NbLiArr
  Real(PR) :: dx, dy, D, Lx, Ly, eps, dt, tmax
  
  call MPI_INIT(statinfo)
  if (statinfo<0) then
     write (*,*) 'Error in MPI_INIT'
     STOP
  end if
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, Me, statinfo)

  Call lect_para(Nx, Ny, Lx, Ly, D, dt, tmax)

  N = Nx * Ny
  nbNN = 5*(Nx*Ny) - 2*Nx - 2*Ny !! nombre d'éléments non nuls 

  Allocate(AA(nbNN), IA(nbNN), JA(nbNN), U(Nx*Ny), F(Nx*Ny), U0(Nx*Ny), &
        fsource(Nx,Ny), g(Nx,Ny), h(Nx,Ny), U_Mat(Nx,Ny), x(Nx), y(Ny), &
        tab_i1(Np), tab_iN(Np), S(nbNN))

  dx = Lx / Nx 
  dy = Ly / Ny

  imax = Int(tmax/dt)

  !-------------------------RESOLUTION---------------------------
  eps = 0.01    ! GC, tolerance et nombre max d iteration
  kmax = 100

  !génération de la matrice pentadiagonale
  Call Mat(Nx,Ny,dx,dy,D,AA,JA,IA,dt)
  
  !on répartit la matrice A (en sparse) dans chaque proc
  call charge(nbNN,Np,me,IA,tab_i1,tab_iN,i1,iN)
  Nme = iN-i1+1 !nb d'éléments dans les vect locaux
  NbLiMe = IA(iN)-IA(i1)+1 !nb de lignes différentes des elts stockés
  NbColMe = JA(iN)-JA(i1)+1 !nb de colonnes différentes des elts stockés

  Allocate(AAme(Nme),IAme(Nme),JAme(Nme))
  AAme = AA(i1:iN)
  IAme = IA(i1:iN)
  JAme = JA(i1:iN)
  Deallocate(AA,JA,IA) !note pour CR

  !conditions aux bords
  Call fsta(Lx,Ly,Nx,Ny,fsource)
  Call gsta(Lx,Ly,Nx,Ny,g)
  Call hsta(Lx,Ly,Nx,Ny,h)

  ! génération du terme source, conditions aux bords + fsource
  Call secondMembre(F,Nx,Ny,dx,dy,D,fsource,g,h)

  allocate(Fme(NbLiMe),b(NbLiMe),U0me(NbColMe),Ume((NbColMe)))
  Fme = F(IA(i1):IA(iN))
  Deallocate(F)

  !Initialisation : 
  U0me = 2.! solution initiale ; U initialisé dans GC
  
  b = dt*Fme+U0me(IA(i1):IA(iN)) !mise à l'échelle des vecteurs
  
  !------on récupère le nb de lignes du proc d'avant et d'après
  NbLiAv = 0
  NbLiArr = 0
  if (me==0) then
     call MPI_SEND(NbLiMe,1,MPI_INTEGER,Me+1,102,MPI_COMM_WORLD,statinfo)
     call MPI_RECV(NbLiAv,1,MPI_INTEGER,Me+1,101,MPI_COMM_WORLD,status,statinfo)
  elseif (me==Np-1) then
     call MPI_SEND(NbLiMe,1,MPI_INTEGER,Me-1,102,MPI_COMM_WORLD,statinfo)
     call MPI_RECV(NbLiArr,1,MPI_INTEGER,Me-1,101,MPI_COMM_WORLD,status,statinfo)
  else
     call MPI_SEND(NbLiMe,1,MPI_INTEGER,Me-1,101,MPI_COMM_WORLD,statinfo)
     call MPI_SEND(NbLiMe,1,MPI_INTEGER,Me+1,102,MPI_COMM_WORLD,statinfo)
     call MPI_RECV(NbLiAv,1,MPI_INTEGER,Me+1,101,MPI_COMM_WORLD,status,statinfo)
     call MPI_RECV(NbLiArr,1,MPI_INTEGER,Me-1,102,MPI_COMM_WORLD,status,statinfo)
  end if

  ! résolution du système AU = dt*F+U0, pour chaque pas de temps 
  Do i = 1,imax
     Call GC_SPARSE_PARLL(AAme,IAme,JAme,b,U0me,Ume,eps,kmax,NbLiMe,NbColMe,NbLiAv,NbLiArr,Me,statinfo,status,Np)
     U0 = U
  End Do
  
  !-------------------------AFFICHAGE-----------------------------

  Call vect_to_mat(U, U_Mat, Nx, Ny)

  Do i = 1, Nx
     x(i) = i * dx
  End Do
  
  Do i = 1, Ny
     y(i) = i * dy
  End Do

  Call WriteBinary(x, y, U_Mat, "resultat.dat")
  
!!$  !Instruction
!!$  Print *,'Donner la valeur de Np'
!!$  Read *, Np
!!$  Do me=0,Np-1
!!$     Call charge(nbNN,Np,me,IA,tab_i1,tab_iN,i1,iN)
!!$     Print *, 'Les valeurs de i1 et de iN et delta N sont pour me= ',me,' :'
!!$     Print *, i1,' ',iN,' ',iN-i1+1
!!$     Print *, IA(i1), ' ', IA(iN), ' deltaI : ', IA(iN)-IA(i1)
!!$     Print *, JA(i1), ' ', JA(iN), ' deltaJ : ', JA(iN)-JA(i1)
!!$
!!$  End Do
!!$
!!$  !!Vérification de MatVect en parallèle
!!$  S = 0.
!!$  Do me=0,Np-1
!!$     Call charge(nbNN,Np,me,IA,tab_i1,tab_iN,i1,iN)
!!$     S = S + MatVectSPARSE(AA(i1:iN),IA(i1:iN),JA(i1:iN),U0)
!!$  End Do
!!$
!!$  !on affiche différence du pdt en sequentiel et parallèle (=0?)
!!$  print *, S - MatVectSPARSE(AA,IA,JA,U0)




  Deallocate(AAme,IAme,JAme,Ume,Fme,U0me, g, h, U_Mat,x, y,tab_i1,tab_iN,S,b,U,U0)
  
  call MPI_FINALIZE

Contains

 Subroutine WriteBinary(x, y, v, nom_fichier)
  
  Implicit None
  
  ! arguments
  Real(PR), Dimension(:) :: x, y
  Real(PR), Dimension(:, :) :: v
  Character(len=*) :: nom_fichier

  ! variables locales
  Integer :: i, j, Nx, Ny
  
  Open(11, file = nom_fichier)
  Nx = Size(x)
  Ny = Size(y)
  !  write(11,*) real(x(1),kind=PR), real(x(Nx),kind=PR), Nx !xmin, xmax, Nx
  !  write(11,*) real(y(1),kind=PR), real(y(Ny),kind=PR), Ny !ymin, ymax, Ny
  Do i = 1, Nx
     Do j = 1, Ny
        Write(11,*) x(i), y(j), Real(v(i, j), kind=4)
     End Do
  End Do

  Close(11)

End Subroutine WriteBinary




End  Program Main
