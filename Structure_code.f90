Program Main

  Use fonctions
  Use second_membre_sparse
  Use syslin
  Use mod_parallele

  Implicit None

  !include "mpif.h"
  integer :: statinfo
  Integer, Dimension(MPI_STATUS_SIZE) :: status

  Real(PR), Dimension(:), Allocatable :: AAme, x, y
  Integer, Dimension(:), Allocatable :: IAme, JAme
  Real(PR), Dimension(:), Allocatable :: Ume, U0me, Fme,
  Real(PR), Dimension(:,:), Allocatable :: fsourceme, gme, hme
  Integer :: Nli, Ncol, N, nbNN, kmax, i, imax, me, i1, iN, Np, Nlime
  Real(PR) :: dx, dy, D, Lx, Ly, dt, tmax, beta


  call MPI_INIT(statinfo)
  if (statinfo < 0) then
     write (*,*) 'Error in MPI_INIT'
     STOP
  end if

  call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, Me, statinfo)

  Call lect_para(Nli, Ncol, Lx, Ly, D, dt, tmax)

  ! ____________________définition des paramètres globaux_____________________
  beta = 0.01    ! tolerance
  kmax = 100     ! nombre max d'itération
  dx = Lx / Nx
  dy = Ly / Ny
  imax = Int(tmax/dt)

  ! _______________________découpage du domaine_______________________________
  call charge(Ncol, Np, me, i1, iN) ! renvoie indices min et max des lignes

  Nlime = iN - i1 + 1 ! calcul du nombre de ligne pour chaque proc

  N = Nlime*Ncol ! nombre d'éléments dans le domaine me
  nbNN = 5*N - 2*Nlime - 2*Ncol !! nombre d'éléments non nuls pour chaque domaine
  !cette valeur est necessaire pour le stockage sparse
  ! 3 vecteurs
  ! AA valeur des elements non nuls
  ! IA (indice i) des elements non nuls
  ! JA (indice j) des elements non nuls

  !____________________________________________________________________________
  ! definition de la matrice Ame et des vecteurs U0me et Ume pour chaque sous-système matriciel (par domaine)
  !____________________________________________________________________________
  Allocate(AAme(nbNN), IAme(nbNN), JAme(nbNN), U0me(N), Ume(N))

  Call Mat(Nlime, Ncol, dx, dy, D, AAme, JAme, IAme, dt)
  U0me = 0. !Etat initial puis solution à l'état précédent

  ! ____________________________________________________________________________
  ! definition des conditions aux bords et du second membre pour chaque sous-domaine
  ! ____________________________________________________________________________
  Allocate(gme(2,Ncol), hme(Nlime,2), fsourceme(Nlime,Ncol))
  call g(me, Ncol, dx, dy, U0me, gme)   ! conditions aux bords horizontaux des sous-domaines
  call h(me, Ncol, dx, dy, hme)         ! conditions aux bords verticaux des sous-domaines
  call fsource(me, Ncol, dx, dy, fsourceme)   ! terme source des sous-domaines

  call secondMembre(Fme, Nlime, Ncol, dx, dy, D, gme, hme, fsourceme)    ! calcul du second membre

  ! ______________________résolution par itérations en temps___________________
  do i = 1,imax  ! boucle sur le nombre d'itérations en temps
    call GC_SPARSE(AAme, IAme, JAme, dt*Fme + U0me, U0me, Ume, beta, kmax, N)    ! resolution par GC sur chaque proc
    U0me = Ume     ! mise à jour de la solution
    ! communications MPI pour le recouvrement // créer type MPI_LIGNE
    if (me /= Np-1) then
       call MPI_SEND ! envoie dernière ligne à proc me+1
       call MPI_RECV ! reçoit première ligne de me+1 dans la deuxième ligne de gme
    end if
    if (me /= 0) then
      call MPI_SEND ! envoie première ligne à proc me-1
      call MPI_RECV ! reçoit dernière ligne de me-1 dans la première ligne de gme
    end if
    call secondMembre(Fme, Nlime, Ncol, dx, dy, D, gme, hme, fsourceme) ! mise à jour du second membre
  end do









   ! affichage : pouvoir tracer la solution (les U) avec gnuplot
   ! à refaire


  !-------------------------AFFICHAGE-----------------------------

  Call vect_to_mat(U, U_Mat, Nli, Ncol)

  Do i = 1, Nli
     x(i) = i * dx
  End Do

  Do i = 1, Ncol
     y(i) = i * dy
  End Do

  Call WriteBinary(x, y, U_Mat, "resultat.dat")

  Deallocate(AAme, IAme, JAme, Ume, Fme, U0me, gme, hme, x, y)

  call MPI_FINALIZE

Contains

 Subroutine WriteBinary(x, y, v, nom_fichier)

  Implicit None

  ! arguments
  Real(PR), Dimension(:) :: x, y
  Real(PR), Dimension(:, :) :: v
  Character(len=*) :: nom_fichier

  ! variables locales
  Integer :: i, j, Nli, Ncol

  Open(11, file = nom_fichier)
  Nli = Size(x)
  Ncol = Size(y)
  !  write(11,*) real(x(1),kind=PR), real(x(Nli),kind=PR), Nli !xmin, xmax, Nli
  !  write(11,*) real(y(1),kind=PR), real(y(Ncol),kind=PR), Ncol !ymin, ymax, Ncol
  Do i = 1, Nli
     Do j = 1, Ncol
        Write(11,*) x(i), y(j), Real(v(i, j), kind=4)
     End Do
  End Do

  Close(11)

End Subroutine WriteBinary


End  Program Main
