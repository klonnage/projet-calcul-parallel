Module second_membre_sparse

  use fonctions
  implicit none 
  !  integer, parameter :: PR = 4

Contains

  subroutine Mat(Nx,Ny,dx,dy,D,AA,JA,IA,dt)
    Integer, intent(in) :: Nx,Ny
    Real(PR), intent(in) :: dx,dy,D,dt
    Real(PR), dimension(:), intent(out), allocatable :: AA
    integer, dimension(:), intent(out), allocatable :: JA,IA
    Integer :: i,j,k,nbNN

    nbNN = 5*(Nx*Ny) - 2*Nx - 2*Ny !! nombre d'éléments non nuls ?
    allocate(AA(nbNN),IA(nbNN),JA(nbNN))

    k = 1
    do i=1,Nx*Ny
       do j=1,Ny*Ny
          if(i==j) then                           !! diagonale
             AA(k) = ((2/dx**2)+(2/dy**2))*D*dt + 1.   ! on code I+dt*A
             IA(k) = i
             JA(k) = j
             k = k + 1
          elseif((i==j-1) .AND. (i/=Nx*Ny) .AND. (modulo(i,Nx)/=0)) then!diagonale sup
             AA(k) = (-1/(dx**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          elseif((i==j+1) .AND. (i/=1) .AND. (modulo(j,Ny)/=0)) then  !diagonale inf
             AA(k) = (-1/(dy**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          elseif(j==Ny+i) then                       ! diagonale tout en haut
             AA(k) = (-1/(dy**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          elseif(i==j+Nx) then                           !diagonale tout en bas
             AA(k) = (-1/(dx**2))*D*dt
             IA(k) = i
             JA(k) = j
             k = k + 1
          end if
       end do
    end do
  End subroutine Mat

  subroutine secondMembre(F,Nx,Ny,dx,dy,D,fsource,g,h)
    real(PR), dimension(Nx*Ny), intent(out) :: F
    real(PR), dimension(Nx,Ny), intent(in) :: fsource
    real(PR), dimension(Nx,Ny), intent(in) :: g
    real(PR), dimension(Nx,Ny), intent(in) :: h
    Integer, intent(in) :: Nx,Ny
    Real(PR), intent(in) :: dx,dy,D
    Integer :: i,j,C
    C=1
    do i=1,Nx
       do j=1,Ny
          F(C) = fsource(i,j)
          if ((i==1) .OR. (i==Nx)) then
             F(C) = F(C) + g(i,j)*D/dy**2
             print*, g(i,j)
             if ((j==1) .OR. (j==Ny)) then
                F(C) = F(C) + h(i,j)*D/dx**2
             end if
          elseif ((j==1) .OR. (j==Ny)) then
             F(C) = F(C) + h(i,j)*D/dx**2
             if ((i==1) .OR. (i==Nx)) then
                F(C) = F(C) + g(i,j)*D/dy**2
             end if
          end if
          C = C + 1
       end do
    end do

  end subroutine secondMembre

End Module second_membre_sparse
