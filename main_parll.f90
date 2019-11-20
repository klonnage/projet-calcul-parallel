Program Main

  Use fonction
  Use second_membre_sparse
  Use syslin_parll
  Use mod_parallele

  Implicit None

  Include "mpif.h"

  Integer :: statinfo
  Integer, Dimension(MPI_STATUS_SIZE) :: status

  Character(len=100) :: nom_fichier
  Real(PR), Dimension(:), Allocatable :: AA, x, y, AAme
  Integer, Dimension(:), Allocatable :: IA, JA, tab_i1, tab_iN, IAme, JAme
  Real(PR), Dimension(:), Allocatable :: Ume, F, U0me, S, Fme, b, U, U0
  Real(PR), Dimension(:,:), Allocatable :: fsource, g, h, U_Mat
  Integer :: Nx, Ny, N, nbNN, kmax, i, imax, me, i1, iN, Np, Nme, NbLiMe, NbColMe
  Real(PR) :: dx, dy, D, Lx, Ly, eps, dt, tmax, t0, t1

  !---------------Initialisation paralelle----------------------
  Call MPI_INIT(statinfo)
  If (statinfo<0) Then
     Write (*,*) 'Error in MPI_INIT'
     Stop
  End If
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, statinfo)
  Call MPI_COMM_RANK(MPI_COMM_WORLD, Me, statinfo)

  !--------------Lecture des parametres
  Call lect_para(Nx, Ny, Lx, Ly, D, dt, tmax)

  N = Nx * Ny
  nbNN = 5*(Nx*Ny) - 2*Nx - 2*Ny ! nombre d'éléments non nuls dans la matrice A

  Allocate(AA(1:nbNN), IA(1:nbNN), JA(1:nbNN), U(N), F(N), U0(N), &
       fsource(Nx,Ny), g(Nx,Ny), h(Nx,Ny), U_Mat(Nx,Ny), x(Nx), y(Ny), &
       tab_i1(Np), tab_iN(Np), S(nbNN))

  dx = Lx / Nx
  dy = Ly / Ny
  imax = Int(tmax/dt)

  !-------------------------RESOLUTION---------------------------
  !génération de la matrice pentadiagonale
  Call Mat(Nx,Ny,dx,dy,D,AA,JA,IA,dt)

  !on répartit la matrice A (en sparse) dans chaque proc
  Call charge(nbNN,Np,me,IA,tab_i1,tab_iN,i1,iN)
  Nme = iN-i1+1 !nb d'éléments dans les vect locaux
  NbLiMe = IA(iN)-IA(i1)+1 !nb de lignes différentes des elts stockés
  NbColMe = JA(iN)-JA(i1)+1 !nb de colonnes différentes des elts stockés
  Allocate(AAme(1:Nme),IAme(1:Nme),JAme(1:Nme))
  AAme = AA(i1:iN)
  IAme = IA(i1:iN)
  JAme = JA(i1:iN)

  !Conditions aux bords
  Call fper(Lx,Ly,Nx,Ny,fsource)
  Call fper(Lx,Ly,Nx,Ny,g)
  Call fper(Lx,Ly,Nx,Ny,h)

  !Génération du terme source, conditions aux bords + fsource
  Call secondMembre(F,Nx,Ny,dx,dy,D,fsource,g,h)
  Allocate(Fme(NbLiMe),b(NbLiMe),U0me(NbColMe),Ume((NbColMe)))
  Fme = F(IA(i1):IA(iN))
  Deallocate(F)

  !Initialisation de la résolution :
  U0me = 2._PR! solution initiale ; U initialisé dans GC

  !Création du second membre du problème matriciel
  If (me==0) Then
     b = dt*Fme+U0me(1:NbColMe-Ny)
  Elseif (me==Np-1) Then
     b = dt*Fme+U0me(Ny+1:NbColMe)
  Else
     b = dt*Fme+U0me(Ny+1:NbColMe-Ny)
  End If

  ! résolution du système AU = dt*F+U0, pour chaque pas de temps, avec différentiation pour le proc 0 et Np-1

  call cpu_time(t0)


  Do i = 1,imax
     If (me==0) Then
        Call GC_SPARSE_PARLL(AAme,IAme,JAme,dt*Fme+U0me(1:NbColMe-Ny), &
             U0me,Ume,eps,kmax,NbLiMe,NbColMe,Me,statinfo,status,Np,Ny,I1,In)
     Elseif (me==Np-1) Then
        Call GC_SPARSE_PARLL(AAme,IAme,JAme,dt*Fme+U0me(Ny+1:NbColMe), &
             U0me,Ume,eps,kmax,NbLiMe,NbColMe,Me,statinfo,status,Np,Ny,I1,In)
     Else
        Call GC_SPARSE_PARLL(AAme,IAme,JAme,dt*Fme+U0me(Ny+1:NbColMe-Ny), &
             U0me,Ume,eps,kmax,NbLiMe,NbColMe,Me,statinfo,status,Np,Ny,I1,In)
     End If
     U0me = Ume
  End Do

  call cpu_time(t1)
  print*, t1-t0, ' s'

  !-------------------------AFFICHAGE----------------------------
  !Partage de la solution à proc Me=Np-1
  U = 0
  U0 = 0
  Print*,IA(I1),IA(In)
  If (me==0) Then
     U(IA(I1):IA(In))=Ume(1:NbColMe-Ny)
  Elseif (me==Np-1) Then
     U(IA(I1):IA(In))=Ume(Ny+1:NbColMe)
  Else
     U(IA(I1):IA(In))=Ume(Ny+1:NbColMe-Ny)
  End If

  If (Me==0) Then
     Call MPI_SEND(U,Nx*Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo)
  Elseif(Me==Np-1) Then
     Call MPI_RECV(U0,Nx*Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo)
     U=U+U0
  Else
     Call MPI_RECV(U0,Nx*Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo)
     U=U+U0
     Call MPI_SEND(U,Nx*Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo)
  End  If

  !Ecriture de la solution dans un fichier
  If (Me==Np-1) Then
     Call vect_to_mat(U, U_Mat, Nx, Ny)
     Do i = 1, Nx
        x(i) = i * dx
     End Do

     Do i = 1, Ny
        y(i) = i * dy
     End Do
     nom_fichier='resultat.dat'
     Call WriteBinary(x, y, U_Mat, nom_fichier)
  End  If

  Deallocate(AA,JA,IA)
  Deallocate(AAme,IAme,JAme,Ume,Fme,U0me, g, h, U_Mat,x, y,tab_i1,tab_iN,S,b,U,U0,fsource)
  Call MPI_FINALIZE


Contains

  Subroutine WriteBinary(x, y, v, nom_fichier)
    !Ecriture dans un fichier
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
    Do i = 1, Nx
       Do j = 1, Ny
          Write(11,*) x(i), y(j), Real(v(i, j), kind=4)
       End Do
    End Do
    Close(11)

  End Subroutine WriteBinary

End  Program Main
