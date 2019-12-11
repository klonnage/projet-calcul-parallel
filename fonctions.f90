Module fonctions

  implicit none

  include 'mpif.h'

  integer, parameter :: PR = 8

Contains

  subroutine lect_para(Nli, Ncol, Lx, Ly, D, dt, tmax)
    integer, intent(out) :: Nli, Ncol
    real(PR), intent(out) :: Lx, Ly, D, dt, tmax

    open(10, file="parametres.dat")

    read(10, *) Nli
    read(10, *) Ncol
    read(10, *) Lx
    read(10, *) Ly
    read(10, *) D
    read(10, *) dt
    read(10, *) tmax
    read(10, *) mode
    read(10, *) beta
    read(10, *) kmax

    close(10)

  end subroutine lect_para

  subroutine g(me, Ncol, dx, Ly, gme, mode)
    integer, intent(in) :: me, Ncol, mode
    real(PR), intent(in) :: dx, Ly
    Real(PR), intent(in), dimension(2,Ncol) :: gme
    integer :: i
    real(PR) :: x

    select case (mode)
    case (1) ! conditions aux bords stationnaires
      gme = 0
    case (2) ! conditions aux bords periodiques
      gme = 0
      if (me = 0) then
        do (i=1,Ncol)
            x = real(i)*dx
            gme(1,i) = sin(x) + 1 !cos(0) = 1
        enddo
      end if
      if (me = Np) then
        do (i=1,Ncol)
            x = real(i)*dx
            gme(2,i) = sin(x) + cos(Ly)
        enddo
      end if
    case(3) ! conditions aux bords instationnaires periodiques
      gme = 0
    endselect
  end subroutine g

  subroutine h(me, Nlime, i1, dy, Lx, hme, mode)
    integer, intent(in) :: me, Nlime, i1, mode
    real(PR), intent(in) :: dy, Lx
    Real(PR), intent(in), dimension(Nlime,2) :: hme
    integer :: i

    select case (mode)
    case (1) ! conditions aux bords stationnaires
      hme = 0
    case (2) ! conditions aux bords periodiques
      do (i=1,Nlime)
          y = real(i + i1 - 1)*dy
          hme(i,1) = cos(y) !sin(0) = 0
          hme(i,2) = cos(y) + sin(Lx)
      enddo
    case(3) ! conditions aux bords instationnaires periodiques
      hme = 1
    endselect
  end subroutine h

subroutine fsource(me, Ncol, Nlime, i1, dx, dy, t, fsourceme, mode)
    integer, intent(in) :: me, Ncol, Nlime, i1, mode
    real(PR), intent(in) :: dy, Lx, t
    Real(PR), intent(in), dimension(Nlime,Ncol) :: fsourceme
    integer :: i, j

    select case (mode)
    case (1) ! terme source stationnaire
      do i = 1,Nlime
        do j = 1,Ncol
          x = real(i + i1 -1)*dx
          y = real(j)*dy
          fsourceme(i,j) = 2*(y - y**2 + x - x**2)
        end do
      end do
    case (2) ! terme source periodique
      do i = 1,Nlime
       do j = 1,Ncol
          x = real(i + i1 - 1)*dx
          y = real(j)*dy
          fsourceme(i,j) = sin(x) + cos(y)
       end do
     end do
   case(3) ! terme source instationnaires periodiques
      do i = 1,Nlime
       do j = 1,Ncol
          x = real(i + i1 - 1)*dx
          y = real(j)*dy
          fsourceme(i,j) =  exp(-(x-Lx*0.5_PR)**2) * exp(-(y-Ly*0.5_PR)**2)*cos(t*pi*0.5_PR)
       end do
    end do
   endselect
 end subroutine fsource
