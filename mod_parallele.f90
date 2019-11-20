Module mod_parallele

  Implicit None

  !include "mpif.h"

Contains
  
  Function dernier_terme(IA,iN)
    !Arguments
    Integer,Intent(In),Dimension(:)::IA
    Integer,Intent(In)::iN
    Integer::dernier_terme
    !Instructions
    dernier_terme = iN
    do while (IA(dernier_terme) == IA(dernier_terme+1))
       dernier_terme=dernier_terme+1
    end do
  end Function dernier_terme


  Subroutine charge(n,Np,me,IA,tab_i1,tab_iN,i1,iN)
    !Arguments
    Integer,Intent(In)::n,Np,me
    Integer,Intent(In),Dimension(n)::IA
    Integer,Intent(Out),Dimension(Np)::tab_i1,tab_iN
    Integer,Intent(Out)::i1,iN
    !Variables
    Integer::a,k
    !Instructions
    a=Int(Real(n)/Real(Np))
    tab_i1(1) = 1
    tab_iN(1) = dernier_terme(IA,tab_i1(1) + a-1)
    do k=2,Np-1
       tab_i1(k) = tab_iN(k-1) + 1
       tab_iN(k) = dernier_terme(IA,tab_i1(k) + a-1)
    end do
    tab_i1(Np) = tab_iN(Np-1) + 1
    tab_iN(Np) = n
    i1 = tab_i1(me+1)
    iN = tab_iN(me+1)
  end Subroutine charge


end Module mod_parallele
