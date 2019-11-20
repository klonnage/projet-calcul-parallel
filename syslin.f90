Module syslin

 use fonctions

  Implicit None

  Integer, Parameter :: kmax=1000
  Real(PR), Parameter :: eps=0.001
  

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


  Subroutine GC_SPARSE(AA,IA,JA,b,x0,x,beta,k,n)
    Implicit None
    Integer, Intent(In) :: n
    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(:), Intent(In) :: b,x0
    Real(PR), Dimension(:), Intent(Out) :: x
    Real(PR), Intent(Out) :: beta
    Real(PR), Dimension(:), Allocatable :: p,z,r1,r2
    Real(PR) :: alpha,gamma
    Integer, Intent(Out) :: k
    Allocate(p(n),z(n),r1(n),r2(n))

    x = x0

    r1=b-MatmulSPARSE(AA,IA,JA,x0)
    p=r1
    beta=sqrt(sum(r1*r1))
    k=0

    Do While (beta>eps .and. k<=kmax)

       z=MatmulSPARSE(AA,IA,JA,p)
       alpha=dot_product(r1,r1)/dot_product(z,p)
       x = x + alpha*p

       r2=r1-alpha*z
       gamma=dot_product(r2,r2)/dot_product(r1,r1)
       p=r2+gamma*p
       beta=sqrt(sum(r1*r1))
       k=k+1
       r1=r2
    End Do
    If (k>kmax) then
       Print*,"Tolérance non atteinte:",beta
    End If
    Deallocate(p,z,r1,r2)
  End Subroutine GC_SPARSE


  Function MatmulSPARSE(AA,IA,JA,x) Result(y)
    !Produit MatriceCSR.vecteur plein (x) retourne vecteur plein(y)
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
  End Function MatmulSPARSE

  


End Module syslin   
