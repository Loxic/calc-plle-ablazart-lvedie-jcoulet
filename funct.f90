module patate

  implicit none
  integer,parameter::wp=8

  contains
  subroutine Get_F(F,Nx,Ny,dx,dy,D,Dt,time,nb_probleme)
    implicit none
    real(wp),dimension(:),allocatable,intent(out)::F
    integer,intent(in)::Nx,Ny,nb_probleme
    real(wp),intent(in)::dx,dy,D,Dt,time

    integer::i,j

    allocate(F(Nx*Ny))
    F=0
    do i=1,Nx
       do j=1,Ny
          F(Nx*(j-1)+i)=fonction_f(i*dx,j*dy,dx*(Nx+1),dy*(Ny+1),time,nb_probleme)
       end do
    end do

!!!! Division du maillage : définir sur quel domaine boucle chaque proc (charge 1D, 2D..)

! routine partionnement maillage
! => maillage de taille Nx . Ny 
! => sous maillage de taille Mx . My 
! => réecriture : distribution du vecteur U et F sur le maillage

! Petite routine de réecriture globale à faire pour obtenir le vecteur U final ?
! Ou alors, multiplot gnuplot


!!!!! Prise en compte des CL sur un sous maillage

! les Mx première coordonnées appellent le copain de GAUCHE <-> le copain passe ses Mx dernières
! les Mx dernière coordonnées appellent le copain de DROITE <-> le copain passe ses Mx premières
! chaque '1 + (i-1)*Mx' pour i allant de 1 à My appelle le voisin du BAS <-> le voison passe chaque 'i*Mx' pour i allant de 1 à My
! chaque 'i*Mx' pour i allant de 1 à Ny prend appelle le voisin du HAUT<-> le voison passe chaque 'i*Mx' pour i allant de 1 à My

! Critère de détection des bords
! Remplacer voisin par g(en 0/en 1) et copain par h(en 0 /en 1)



!!!!!! Recouvrement

! On élargit les domaines => Nouvelle routine de partionnement
! 'CL' plus larges : au lieu de prendre juste l'indice du bord à côté, on en prend une tranche, et on l'inverse

    do j=1,Ny
       F(Nx*(j-1)+1) = F(Nx*(j-1)+1) + (D/dx**2)*fonction_h(0.0_wp,j*dy,nb_probleme)
       F(Nx*(j-1)+Nx) = F(Nx*(j-1)+Nx) + (D/dx**2)*fonction_h((Nx+1)*dx,j*dy,nb_probleme)
    end do
    do i=1,Nx
       F(Nx*(1-1)+i) = F(Nx*(1-1)+i) + (D/dy**2)*fonction_g(i*dx,0.0_wp,nb_probleme)
       F(Nx*(Ny-1)+i) = F(Nx*(Ny-1)+i) + (D/dy**2)*fonction_g(i*dx,(Ny+1)*dy,nb_probleme)
    end do
    F=Dt*F

  end subroutine Get_F

  function fonction_f(x,y,Lx,Ly,t,nb_probleme)
    implicit none
    real(wp),intent(in)::x,y
    real(wp),intent(in)::Lx,Ly,t
    integer,intent(in)::nb_probleme
    real(wp)::fonction_f

    select case(nb_probleme)
    case(1)
       fonction_f=2*(x+y-x**2-y**2)
    case(2)
       fonction_f=sin(x)+cos(y)
    case(3)
       fonction_f=exp(-(x-0.5*Lx)**2)*exp(-(y-0.5*Ly)**2)*cos(0.5*3.14116*t)
    case default
       fonction_f=0
    end select
  end function fonction_f

  function fonction_g(x,y,nb_probleme)
    implicit none
    real(wp),intent(in)::x,y
    real(wp)::fonction_g
    integer,intent(in)::nb_probleme

    select case(nb_probleme)
    case(1)
       fonction_g=0
    case(2)
       fonction_g=sin(x)+cos(y)
    case(3)
       fonction_g=0
    case default
       fonction_g=0
    end select
  end function fonction_g

  function fonction_h(x,y,nb_probleme)
    implicit none
    real(wp),intent(in)::x,y
    real(wp)::fonction_h
    integer,intent(in)::nb_probleme

    select case(nb_probleme)
    case(1)
       fonction_h=0
    case(2)
       fonction_h=sin(x)+cos(y)
    case(3)
       fonction_h=1
    case default
       fonction_h=0
    end select
  end function fonction_h

  subroutine matmul_implicit(Nx,Ny,dx,dy,D,Dt,X,AX)

    implicit none

    integer, intent(in) :: Nx,Ny  ! dimensions spatiales du problème
    real(wp), intent(in)::dx,dy,D,Dt
    real(wp), dimension(:), intent(in) :: X ! donnee
    real(wp), dimension(:), intent(out) :: AX ! Sortie

    integer :: i, j, k

    AX=0
    do i = 1, Nx
       do j = 1, Ny
          k = Nx*(j-1) + i
          ! bloc M
          AX(k) = (1+2*D*Dt*(1.0_wp/dx**2 + 1.0_wp/dy**2))*X(k)
          if (i>1) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k-1)
          end if
          if (i<Nx) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k+1)
          end if

          !Bloc E inférieur
          if(j>1) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k-Ny)
          end if
          !Bloc E supérieur
          if (j<Ny) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k+Ny)
          end if

       end do
    end do

  end subroutine matmul_implicit


  subroutine save_result(U,Nx,Ny,dx,dy,fichier)
    implicit none
    real(wp),dimension(:),intent(in)::U
    integer,intent(in)::Nx,Ny
    character(len=*),intent(in)::fichier
    real(wp),intent(in)::dx,dy


    real(wp)::x,y
    integer::i,j

    open(unit=12,file=fichier)

    do i=1,Nx
       do j=1,Ny
          write(12,*) i*dx, j*dy, U(Nx*(j-1)+i)
       end do
       write(12,*)''
    end do

    close(12)


  end subroutine save_result

  end module
