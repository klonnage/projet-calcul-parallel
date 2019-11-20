Module second_membre_sparse

  Use fonction

  Implicit None 

Contains

  Subroutine Mat(Nx,Ny,dx,dy,D,AA,JA,IA,dt)
    Integer, Intent(in) :: Nx,Ny
    Real(PR), Intent(in) :: dx,dy,D,dt
    Real(PR), Dimension(:), Intent(out), Allocatable :: AA
    Integer, Dimension(:), Intent(out), Allocatable :: JA,IA
    Integer :: i,j,k,nbNN

    nbNN = 5*(Nx*Ny) - 2*Nx - 2*Ny !! nombre d'éléments non nuls ?
    Allocate(AA(nbNN),IA(nbNN),JA(nbNN))

    k = 1
    Do i=1,Nx*Ny
       Do j=1,Ny*Ny
          If(i==j) Then                           !! diagonale
             AA(k) = ((2/dx**2)+(2/dy**2))*D*dt + 1.   ! on code I+dt*A
             IA(k) = i
             JA(k) = j
             k = k + 1
          Elseif((i==j-1) .And. (i/=Nx*Ny) .And. (Modulo(i,Nx)/=0)) Then  !diagonale sup
             AA(k) = (-1/(dx**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          Elseif((i==j+1) .And. (i/=1) .And. (Modulo(j,Ny)/=0)) Then  !diagonale inf
             AA(k) = (-1/(dx**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          Elseif(j==Ny+i) Then                       ! diagonale tout en haut
             AA(k) = (-1/(dy**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          Elseif(i==j+Ny) Then                           !diagonale tout en bas
             AA(k) = (-1/(dy**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          End If
       End Do
    End Do
  End Subroutine Mat

  Subroutine secondMembre(F,Nx,Ny,dx,dy,D,fsource,g,h)
    Real(PR), Dimension(Nx*Ny), Intent(out) :: F
    Real(PR), Dimension(Nx,Ny), Intent(in) :: fsource
    Real(PR), Dimension(Nx,Ny), Intent(in) :: g
    Real(PR), Dimension(Nx,Ny), Intent(in) :: h
    Integer, Intent(in) :: Nx,Ny
    Real(PR), Intent(in) :: dx,dy,D
    Integer :: i,j,C
    C=1
    Do i=1,Nx
       Do j=1,Ny
          F(C) = fsource(i,j)
          If ((i==1) .Or. (i==Nx)) Then
             F(C) = F(C) + g(i,j)*D/dy**2
             If ((j==1) .Or. (j==Ny)) Then
                F(C) = F(C) + h(i,j)*D/dx**2
             End If
          Elseif ((j==1) .Or. (j==Ny)) Then
             F(C) = F(C) + h(i,j)*D/dx**2
             If ((i==1) .Or. (i==Nx)) Then
                F(C) = F(C) + g(i,j)*D/dy**2
             End If
          End If
          C = C + 1
       End Do
    End Do

  End Subroutine secondMembre

End Module second_membre_sparse
