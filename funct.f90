module mathfunctions

  implicit none
  integer,parameter::wp=8

  contains


!!!!! Prise en compte des CL sur sous maillage

! => maillage de taille Nx . Ny 
! => sous maillage de taille Mx . My 

! Critère de détection des bords :
! Remplacer voisin par g(en 0/en 1) et copain par h(en 0 /en 1)


!!!!!! Recouvrement

! 'CL' plus larges : au lieu de prendre juste l'indice du bord à côté, on en prend une tranche, et on l'inverse

  !subroutine Get_F(F,Nx,Ny,dx,dy,D,Dt,time,nb_probleme)
! 4 nouvelles variables ; noeud de début-fin en x-y
  subroutine Get_F(F,Mx,My,n_x1,n_xn,n_y1,n_yn,dx,dy,D,Dt,time,voisins,nb_probleme,rank,recouv)
    implicit none
    real*8,dimension(:),allocatable,intent(out)::F
    !integer,intent(in)::Nx,Ny,nb_probleme
    integer,intent(in)::Mx,My,n_x1,n_xn,n_y1,n_yn,nb_probleme,rank,recouv
    integer,dimension(4),intent(in)::voisins
    real*8,intent(in)::dx,dy,D,Dt,time
    real*8,dimension(:),allocatable::H1,B1,G1,D1
    integer::i,j

    !allocate(F(Nx*Ny))
    !F=0
   ! do i=1,Nx
    !   do j=1,Ny
   !       F(Nx*(j-1)+i)=fonction_f(i*dx,j*dy,dx*(Nx+1),dy*(Ny+1),time,nb_probleme)
   !    end do
   ! end do

! Nouvelle taille de boucle
    allocate(F(Mx*My),H1(Mx),B1(Mx),G1(My),D1(My))

    H1 = 0 
    B1 = 0 
    G1 = 0 
    D1 = 0 

if (rank==2) then
          write(6,'(I0)')n_x1,n_xn,n_y1,n_yn
end if

    F=0
    !print*,size(H1),size(B1),size(G1),size(D1)
    do j=n_y1,n_yn
       do i=n_x1,n_xn
          F( i-n_x1 +1+ Mx*(j-n_y1) )=fonction_f(i*dx,j*dy,1.,1.,time,nb_probleme)
       end do
    end do


   ! do j=1,Ny
   !    F(Nx*(j-1)+1) = F(Nx*(j-1)+1) + (D/dx**2)*fonction_h(0.0_wp,j*dy,nb_probleme)
   !    F(Nx*(j-1)+Nx) = F(Nx*(j-1)+Nx) + (D/dx**2)*fonction_h((Nx+1)*dx,j*dy,nb_probleme)
   ! end do


!!!!!!!! COMM GAUCHE-DROITE


!!! Le bloc fait toute la rangée (pas de com')

    if ((Voisins(3) < 0).and.(Voisins(1) < 0)) then

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1) = F(Mx*(j-n_y1)+1) + (D/dx**2)*fonction_h(0.0,j*dy,nb_probleme)
          F(Mx*(j+1-n_y1) -1)  = F(Mx*(j+1-n_y1) -1)  + (D/dx**2)*fonction_h(1.,j*dy,nb_probleme)
       end do


!!! Condition de bord gauche : h(0,y) + Com' à droite

    else if (Voisins(1) < 0) then
       
       ! Com' à droite : Envoyer un vecteur D1 de taille My à voisins(3)
       do j=1,My
           D1(j) = F( Mx*j -1 -recouv)
if (rank==2) then
          write(6,'(I0)')Mx*j -1 -recouv
end if
end do
       
       ! Com' à droite : Récupérer le vecteur G1 (nouvellement D1) de taille My de voisins(3)

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1)  = F(Mx*(j-n_y1)+1)  + (D/dx**2)*fonction_h(0.0,j*dy,nb_probleme)
          F(Mx*(j+1-n_y1) -1) = F(Mx*(j+1-n_y1) -1) + (D/dx**2)*D1(j-n_y1 +1)
       end do


!!! Condition de bord droit : h(1,y) + Com' à gauche

    else if (Voisins(3) < 0) then

       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       do j=1,My
          G1(j) = F(j + recouv*Mx)
       end do
       ! Com' à gauche : Récupérer le vecteur D1 (nouvellement G1) de taille My de voisins(1)


       !call MPI_SEND(U0(1:My),1,MPI_REAL,2,tag,MPI_COMM_WORLD,statinfo)
       !call MPI_SEND(Moi(1,1),1,ligne,2,tag,MPI_COMM_WORLD,statinfo)

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1) = F(Mx*(j-n_y1)+1) + (D/dx**2)*G1(j-n_y1 +1)
          F(Mx*(j+1-n_y1) -1)  = F(Mx*(j+1-n_y1) -1)  + (D/dx**2)*fonction_h(1.0,j*dy,nb_probleme)
       end do


!!! Milieu du maillage : full com' Gauche/droite

    else

       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       ! Com' à gauche : Récupérer un vecteur D1 (nouvellement G1) de taille My de voisins(1)
       ! Com' à droite : Envoyer un vecteur D1 de taille My à voisins(3)
       ! Com' à droite : Récupérer le vecteur G1 (nouvellement D1) de taille My de voisins(3)

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1)  = F(Mx*(j-n_y1)+1)  + (D/dx**2)*G1(j-n_y1 +1)
          F(Mx*(j+1-n_y1) -1) = F(Mx*(j+1-n_y1) -1) + (D/dx**2)*D1(j-n_y1 +1)
       end do
   end if


!!!!!!!! COMM HAUT-BAS


!!! Le bloc fait toute la ligne

    if ( Voisins(4) < 0 .and. Voisins(2) < 0 ) then
       do i=n_x1,n_xn
          F(i-n_x1+1 )  = F(i-n_x1+1)  + (D/dx**2)*fonction_g(i*dx,0.0,nb_probleme)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dx**2)*fonction_g(i*dx,1.,nb_probleme)
       end do
 

!!! Condition de bord bas : g(x,0) + Com' en haut

    else if (Voisins(4) < 0) then
       ! Com' en haut : Envoyer un vecteur H1 de taille Mx à voisins(2)
       ! Com' en haut : Récupérer le vecteur B1 (nouvellement H1) de taille Mx de voisins(2)

       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dx**2)*fonction_g(i*dx,0.0,nb_probleme)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dx**2)*H1(i-n_x1+1)
       end do


!!! Condition de bord haut : g(x,1) + Com' en bas

    else if (Voisins(2) < 0) then
       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       ! Com' à gauche : Récupérer le vecteur D1 (nouvellement G1) de taille My de voisins(1)

       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dx**2)*B1(i-n_x1+1)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dx**2)*fonction_g(i*dx,1.,nb_probleme)
       end do

!!! Milieu du maillage : full com' Gauche/droite

    else

       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       ! Com' à gauche : Récupérer un vecteur D1 (nouvellement G1) de taille My de voisins(1)
       ! Com' à droite : Envoyer un vecteur D1 de taille My à voisins(3)
       ! Com' à droite : Récupérer le vecteur G1 (nouvellement D1) de taille My de voisins(3)

       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dx**2)*B1(i-n_x1+1)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dx**2)*H1(i-n_x1+1)
       end do
 end if

    !do i=1,Nx
    !   F(i) = F(i) + (D/dy**2)*fonction_g(i*dx,0.0_wp,nb_probleme)
    !   F(Nx*(Ny-1)+i) = F(Nx*(Ny-1)+i) + (D/dy**2)*fonction_g(i*dx,(Ny+1)*dy,nb_probleme)
    !end do

    F=Dt*F

  end subroutine Get_F


  function fonction_f(x,y,Lx,Ly,t,nb_probleme)
    implicit none
    real*8,intent(in)::x,y
    real*8,intent(in)::Lx,Ly,t
    integer,intent(in)::nb_probleme
    real*8::fonction_f

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
    real*8,intent(in)::x,y
    real*8::fonction_g
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
    real*8,intent(in)::x,y
    real*8::fonction_h
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
    real*8, intent(in)::dx,dy,D,Dt
    real*8, dimension(:), intent(in) :: X ! donnee
    real*8, dimension(:), intent(out) :: AX ! Sortie

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
    real*8,dimension(:),intent(in)::U
    integer,intent(in)::Nx,Ny
    character(len=*),intent(in)::fichier
    real*8,intent(in)::dx,dy


    real*8::x,y
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
