Module syslin_parll

  Use fonction

  Implicit None

Contains


  Subroutine GC_SPARSE_PARLL(AA,IA,JA,b,x0,x,eps,kmax,n,nbcolme,Me,statinfo,status,Np,Ny,I1,In)

    Implicit None
    Include "mpif.h"
    Integer, Intent(In) :: n,statinfo,kmax,Me,nbcolme,Np,Ny,I1,In
    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(:), Intent(In) :: b,x0
    Integer, Dimension(MPI_STATUS_SIZE),Intent(In) :: status
    Real(PR), Dimension(:), Intent(InOut) :: x
    Real(PR), Intent(In) :: eps
    Real(PR), Dimension(:), Allocatable :: p_resultat,z,r1,r2
    Real(PR), Dimension(:), Allocatable :: p_calcul
    Real(PR) :: alpha,gamma,localbeta,pr1,pr2,pzp,globalpr1,globalpr2,globalpzp,globalbeta
    Integer :: k,i

    Allocate(p_resultat(n),p_calcul(nbcolme),z(n),r1(n),r2(n))
    x = x0
    r1 =  MatVectSPARSE(AA,IA,JA,x0,n) 
    r1 = b - r1 !initialisation du résidu
    p_resultat = r1


    !pdt scalaire ---------------------------------------------------
    localbeta = Sum(r1*r1)
    globalbeta=0._PR  ! La somme de toutes les localbeta
    globalbeta = localbeta

    ! Calcul de beta, equivalent à un ALLREDUCE

    If (me==0) Then
       Do i=1,Np-1
          Call MPI_RECV(localbeta,1,MPI_REAL8,i,i,MPI_COMM_WORLD,status,statinfo)
          globalbeta = globalbeta + localbeta
       End Do
    Else
       Call MPI_SEND(localbeta,1,MPI_REAL8,0,me,MPI_COMM_WORLD,statinfo)
    End If
    Call MPI_BCAST(globalbeta,1,MPI_REAL8,0,MPI_COMM_WORLD,statinfo)

    globalbeta = Sqrt(globalbeta)
    k=0

    ! Communication permettant à chaque proc d'avoir accès à toutes les valeurs du vecteur p 

    If (me==0) Then
       p_calcul(1:NbColMe-Ny) = p_resultat
       Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0
       Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
    Elseif (me==Np-1) Then
       p_calcul(Ny+1:nbcolme) = p_resultat
       Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
       Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0         
    Else
       p_calcul(Ny+1:NbColMe-Ny) = p_resultat
       Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend la fin de p du me-1
       Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me-1 le début de p de me
       Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me+1 la fin de p de me
       Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du me+1
    End If

    ! Debut de la boucle du GC 

    Do While (globalbeta > eps .And. k<=kmax)
       !pdt mat vect ---------------------------------------------------
       z = MatVectSPARSE(AA,IA,JA,p_calcul,n)

       ! Communication permettant à chaque proc d'avoir accès à toutes les valeurs du vecteur p 
       
       If (me==0) Then
          p_calcul(1:NbColMe-Ny) = p_resultat
          Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0
          Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
       Elseif (me==Np-1) Then
          p_calcul(Ny+1:nbcolme) = p_resultat
          Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
          Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0         
       Else
          p_calcul(Ny+1:NbColMe-Ny) = p_resultat
          Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend la fin de p du me-1
          Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me-1 le début de p de me
          Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me+1 la fin de p de me
          Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du me+1
       End If

       !pdt scalaire ----------------------------------------------
       pr1 = Sum(r1*r1)
       pzp = Sum(p_resultat*z)



       ! ALLREDUCE ! BCAST 
       globalpr1 = pr1
       If (me==0) Then
          Do i=1,Np-1
             Call MPI_RECV(pr1,1,MPI_REAL8,i,i,MPI_COMM_WORLD,status,statinfo)
             globalpr1 = globalpr1 + pr1
          End Do
       Else
          Call MPI_SEND(pr1,1,MPI_REAL8,0,me,MPI_COMM_WORLD,statinfo)
       End If
       ! BCAST 
       If (Me==0) Then 
          Do i =1, Np-1
             Call MPI_SEND(globalpr1,1,MPI_REAL8,i,i,MPI_COMM_WORLD,statinfo)
          end do
       else 
          Call MPI_RECV(globalpr1,1,MPI_REAL8,0,Me,MPI_COMM_WORLD,status,statinfo)
       end If


       ! ALLREDUCE 
       globalpzp = pzp

       If (me==0) Then
          Do i=1,Np-1
             Call MPI_RECV(pzp,1,MPI_REAL8,i,i,MPI_COMM_WORLD,status,statinfo)
             globalpzp = globalpzp + pzp
          End Do
       Else
          Call MPI_SEND(pzp,1,MPI_REAL8,0,me,MPI_COMM_WORLD,statinfo)
       End If
       !BCAST
       If (Me==0) Then 
          Do i =1, Np-1
             Call MPI_SEND(globalpzp,1,MPI_REAL8,i,i,MPI_COMM_WORLD,statinfo)
          end do
       else 
          Call MPI_RECV(globalpzp,1,MPI_REAL8,0,Me,MPI_COMM_WORLD,status,statinfo)
       end If

       alpha = globalpr1/globalpzp
       x = x + alpha*p_calcul
       r2 = r1 - alpha*z
       !pdt scalaire --------------------------------------------------
       pr2 = Sum(r2*r2)
       globalpr2 = pr2

       If (me==0) Then
          Do i=1,Np-1
             Call MPI_RECV(pr2,1,MPI_REAL8,i,i,MPI_COMM_WORLD,status,statinfo)
             globalpr2 = globalpr2 + pr2
          End Do
       Else
          Call MPI_SEND(pr2,1,MPI_REAL8,0,me,MPI_COMM_WORLD,statinfo)
       End If
       If (Me==0) Then 
          Do i =1, Np-1
             Call MPI_SEND(globalpr2,1,MPI_REAL8,i,i,MPI_COMM_WORLD,statinfo)
          end do
       else 
          Call MPI_RECV(globalpr2,1,MPI_REAL8,0,Me,MPI_COMM_WORLD,status,statinfo)
       end If

       gamma = globalpr2/globalpr1
       p_resultat = gamma*p_resultat
       p_resultat = p_resultat + r2
       globalbeta = Sqrt(globalpr1)
       k = k+1
       r1 = r2
       
       !Communication des nouvelles valeurs pour chaque p_resultats 

       If (me==0) Then
          p_calcul(1:NbColMe-Ny) = p_resultat
          Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0
          Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
       Elseif (me==Np-1) Then
          p_calcul(Ny+1:nbcolme) = p_resultat
          Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1
          Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0         
       Else
          p_calcul(Ny+1:NbColMe-Ny) = p_resultat
          Call MPI_RECV(p_calcul(1:Ny),Ny,MPI_REAL8,Me-1,Me-1,MPI_COMM_WORLD,status,statinfo) !on prend la fin de p du me-1
          Call MPI_SEND(p_resultat(1:Ny),Ny,MPI_REAL8,Me-1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me-1 le début de p de me
          Call MPI_SEND(p_resultat(n-Ny+1:n),Ny,MPI_REAL8,Me+1,Me,MPI_COMM_WORLD,statinfo) ! on donne à me+1 la fin de p de me
          Call MPI_RECV(p_calcul(nbcolme-Ny+1:nbcolme),Ny,MPI_REAL8,Me+1,Me+1,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du me+1
       End If


    End Do

    If (k>kmax) Then
       Print*,"Tolérance non atteinte:",globalbeta
    End If

    Deallocate(p_resultat,p_calcul,z,r1,r2)
  End Subroutine GC_SPARSE_PARLL


  Function MatVectSPARSE(AA,IA,JA,x,n) Result(y)
    Real(PR), Dimension(:), Intent(In) :: AA,x
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(n) :: y
    Integer, Intent(In) :: n
    Integer :: i,nnz
    nnz=Size(AA)
    y=0._PR
    Do i=1,nnz
       y(IA(i)-IA(1)+1)=y(IA(i)-IA(1)+1)+AA(i)*x(JA(i)-JA(1)+1)
    End Do

  End Function MatVectSPARSE



End Module syslin_parll
