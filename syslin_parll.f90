Module syslin

  use fonctions

  Implicit None

  !include 'mpi.h'

Contains

!!$ Subroutine GC(A,b,x0,x,beta,k,n)
!!$    Implicit None
!!$    Integer, Intent(In) :: n
!!$    Real(PR), Dimension(:,:), Intent(In) :: A
!!$    Real(PR), Dimension(:), Intent(In) :: b,x0
!!$    Real(PR), Dimension(:), Intent(Out) :: x
!!$    Real(PR), Intent(Out) :: beta
!!$    Real(PR), Dimension(:), Allocatable :: p,z,r1,r2
!!$    Real(PR) :: alpha,gamma
!!$    Integer, Intent(Out) :: k
!!$
!!$    Allocate(p(n),z(n),r1(n),r2(n))
!!$    r1=b-Matmul(A,x0)
!!$    p=r1
!!$    beta=Sqrt(Sum(r1*r1))
!!$    k=0
!!$    Do While (beta>eps .And. k<=kmax)
!!$       z=Matmul(A,p)
!!$       alpha=dot_Product(r1,r1)/dot_Product(z,p)
!!$       x=x-alpha*p
!!$       r2=r1-alpha*z
!!$       gamma=dot_Product(r2,r2)/dot_Product(r1,r1)
!!$       p=r2+gamma*p
!!$       beta=Sqrt(Sum(r1*r1))
!!$       k=k+1
!!$       r1=r2
!!$    End Do
!!$    If (k>kmax) Then
!!$       Print*,"Tolérance non atteinte:",beta
!!$    End If
!!$    Deallocate(p,z,r1,r2)
!!$  End Subroutine GC


  Subroutine GC_SPARSE_PARLL(AA,IA,JA,b,x0,x,eps,kmax,n,nbcolme,NbLiAv,NbLiArr,Me,statinfo,status,Np)
    Implicit None
    Integer, Intent(In) :: n,statinfo,NbLiAv,NbLiArr,kmax,Me,nbcolme,Np
    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(:), Intent(In) :: b,x0
    Integer, Dimension(MPI_STATUS_SIZE),Intent(In) :: status
    Real(PR), Dimension(:), Intent(Out) :: x
    Real(PR), Intent(In) :: eps
    Real(PR), Dimension(:), Allocatable :: p_resultat,z,r1,r2
    Real(PR), Dimension(:), Allocatable :: p_calcul
    Real(PR) :: alpha,gamma,localbeta,pr1,pr2,pzp,globalpr1,globalpr2,globalpzp,globalbeta
    Integer :: k
    Allocate(p_resultat(n),p_calcul(nbcolme),z(n),r1(n),r2(n))

    x = x0

    r1 = b - MatVectSPARSE(AA,IA,JA,x0) !ts les vect de la bonne taille
    p_resultat=r1
    
    !pdt scalaire ---------------------------------------------------
    localbeta=sum(r1*r1)
    call MPI_ALLREDUCE(localbeta,globalbeta,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)
    globalbeta = sqrt(globalbeta)

    k=0

    Do While (globalbeta>eps .and. k<=kmax)
       !pdt mat vect ---------------------------------------------------
       if (me==0) then
          p_calcul(1:n) = p_resultat
          call MPI_SEND(p_resultat(n-NbLiArr+1:n),NbLiArr,MPI_REAL,1,102,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0
          call MPI_RECV(p_calcul(n+1:nbcolme),NbLiAv,MPI_REAL,1,101,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1

       elseif (me==Np-1) then
          p_calcul(nbcolme-n+1:nbcolme) = p_resultat
          call MPI_SEND(p_resultat(1:NbLiArr),NbLiArr,MPI_REAL,Me-1,101,MPI_COMM_WORLD,statinfo) ! on donne à proc1 la fin de p de proc0
          call MPI_RECV(p_calcul(1:n),NbLiAv,MPI_REAL,Me-1,102,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du proc1

       else
          p_calcul(NbLiArr+1:n-NbLiAv) = p_resultat
          !proc avant (me+1)
          call MPI_SEND(p_resultat(n-NbLiAv+1:n),NbLiAv,MPI_REAL,Me+1,101,MPI_COMM_WORLD,statinfo) ! on donne à me+1 la fin de p de me
          call MPI_RECV(p_calcul(nbcolme-n+1:nbcolme),n,MPI_REAL,Me+1,102,MPI_COMM_WORLD,status,statinfo) !on prend le debut de p du me+1
          !proc arriere (me-1)
          call MPI_SEND(p_resultat(1:NbLiArr),NbLiArr,MPI_REAL,Me-1,102,MPI_COMM_WORLD,statinfo) ! on donne à me-1 le début de p de me
          call MPI_RECV(p_calcul(1:n),n,MPI_REAL,Me-1,101,MPI_COMM_WORLD,status,statinfo) !on prend la fin de p du me-1
       end if
       z=MatVectSPARSE(AA,IA,JA,p_calcul)

       !pdt scalaire --------------------------------------------------
       pr1 = sum(r1*r1)
       pzp = sum(p_resultat*z)
       call MPI_ALLREDUCE(pr1,globalpr1,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)
       call MPI_ALLREDUCE(pzp,globalpzp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)
       alpha = globalpr1/globalpzp
       x = x + alpha*p_resultat

       r2 = r1 - alpha*z

       !pdt scalaire --------------------------------------------------
       pr2 = sum(r2*r2)
       call MPI_ALLREDUCE(pr2,globalpr2,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)
       gamma = globalpr2/globalpr1
       p_resultat = r2 + gamma*p_resultat
       globalbeta = sqrt(globalpr1)
       k = k+1
       r1 = r2
    End Do
    If (k>kmax) then
       Print*,"Tolérance non atteinte:",globalbeta
    End If
    Deallocate(p_resultat,p_calcul,z,r1,r2)
  End Subroutine GC_SPARSE_PARLL


  Function MatVectSPARSE(AA,IA,JA,x) Result(y)
    !Produit MatriceCSR.vecteur plein (x) retourne vecteur plein (y)
    Real(PR), Dimension(:), Intent(In) :: AA,x
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(Size(x)) :: y
    Integer :: i,n,nnz
    n=Size(x,1)
    nnz=Size(AA,1)
    y(1:n)=0._PR
    Do i=1,nnz
          y(IA(i))=y(IA(i))+AA(i)*x(JA(i))
    End Do
  End Function MatVectSPARSE

  
End Module syslin   
