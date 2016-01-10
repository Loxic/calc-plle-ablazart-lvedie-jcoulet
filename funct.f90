module mathfunctions

use display

  implicit none
  integer,parameter::wp=8


  contains



!!!!! Prise en compte des CL sur sous maillage

! => maillage de taille Nx . Ny 
! => sous maillage de taille Mx . My 

!!!!!! Recouvrement

! 'CL' plus larges : au lieu de prendre juste l'indice du bord à côté, on en prend une tranche, et on l'inverse

subroutine Get_F(F,U,Mx,My,dx,dy,D,Dt,time,voisins,nb_probleme,rank,recouv,comm_cart,map)


    implicit none

    include "mpif.h"

    integer,intent(in)::Mx,My,nb_probleme,rank,recouv,comm_cart
    real*8,dimension(Mx*My),intent(in)::U
    real*8,dimension(Mx*My),intent(inout)::F
    integer,dimension(4),intent(in)::voisins,map
    real*8,intent(in)::dx,dy,D,Dt,time
    real*8,dimension(:),allocatable::H1,B1,G1,D1,H2,B2,G2,D2
    integer::i,j,statinfo,n_x1,n_xn,n_y1,n_yn
    integer, dimension(MPI_STATUS_SIZE) :: status


! Nouvelle taille de boucle
    allocate(H1(Mx),B1(Mx),G1(My),D1(My),H2(Mx),B2(Mx),G2(My),D2(My))

    H1 = 0 ;   B1 = 0 ;    G1 = 0 ;    D1 = 0  
    H2 = 4 ;   B2 = 4 ;    G2 = 4 ;    D2 = 4

    n_x1 = map(1) ; n_xn = map(2) ; n_y1 = map(3) ; n_yn = map(4)

    F=0

    do j=n_y1,n_yn
       do i=n_x1,n_xn
          F( i-n_x1 +1+ Mx*(j-n_y1) )=fonction_f(i*dx,j*dy,1.,1.,time,nb_probleme)
       end do
    end do

!! F OK

!!!!!!!! COMM GAUCHE-DROITE


!!! Le bloc fait toute la rangée (pas de com')

    if ( Voisins(3) < 0 .and. Voisins(1) < 0 ) then
       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1) = F(Mx*(j-n_y1)+1) + (D/dx**2)*fonction_h(0.0,j*dy,nb_probleme)
          F(Mx*(j+1-n_y1))  = F(Mx*(j+1-n_y1))  + (D/dx**2)*fonction_h(1.,j*dy,nb_probleme)
       end do


!!! Condition de bord droit : h(1,y) + Com' à gauche

    else if (Voisins(3) < 0) then

       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       ! Puis récupérer le vecteur D1 (nouvellement G2) de taille My de voisins(1)

       do j=1,My
           G1(j) = U( Mx*(j-1) + 1 + recouv)
       end do

       call MPI_SEND(G1,My,MPI_REAL8,voisins(1),101,comm_cart,statinfo)
       call MPI_RECV(G2,My,MPI_REAL8,voisins(1),102,comm_cart,status,statinfo)

       !call write_data3(rank,voisins,map,nb_probleme,G2,Mx,My,dx,dy,int2char(rank))

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1) = F(Mx*(j-n_y1)+1) + (D/dx**2)*G2(j-n_y1 +1)
          F(Mx*(j+1-n_y1))  = F(Mx*(j+1-n_y1))  + (D/dx**2)*fonction_h(1.0,j*dy,nb_probleme)
       end do


!!! Condition de bord gauche : h(0,y) + Com' à droite

    else if (Voisins(1) < 0) then
       
       ! Com' à droite : Envoyer un vecteur D1 de taille My à voisins(3)
       ! Puis récupérer le vecteur G1 (nouvellement D2) de taille My de voisins(3)
       do j=1,My
           D1(j) = U( Mx*j -recouv)
       end do

       call MPI_SEND(D1,My,MPI_REAL8,voisins(3),102,comm_cart,statinfo)
       call MPI_RECV(D2,My,MPI_REAL8,voisins(3),101,comm_cart,status,statinfo)

       !call write_data2(rank,voisins,map,nb_probleme,D2,Mx,My,dx,dy,int2char(rank))

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1)  = F(Mx*(j-n_y1)+1)  + (D/dx**2)*fonction_h(0.0,j*dy,nb_probleme)
          F(Mx*(j+1-n_y1) ) = F(Mx*(j+1-n_y1) ) + (D/dx**2)*D2(j-n_y1 +1)
       end do


!!! Milieu du maillage : full com' Gauche/droite

    else

       ! Com' à gauche : Envoyer un vecteur G1 de taille My à voisins(1)
       ! Puis récupérer le vecteur D1 (nouvellement G2) de taille My de voisins(1)
       do j=1,My
           G1(j) = U( Mx*(j-1) + 1 + recouv)
       end do

       ! Com' à droite : Envoyer un vecteur D1 de taille My à voisins(3)
       ! Puis récupérer le vecteur G1 (nouvellement D2) de taille My de voisins(3)
       do j=1,My
           D1(j) = U( Mx*j -recouv)
       end do

       call MPI_SEND(G1,My,MPI_REAL8,voisins(1),101,comm_cart,statinfo)
       call MPI_RECV(G2,My,MPI_REAL8,voisins(1),102,comm_cart,status,statinfo)

       call MPI_SEND(D1,My,MPI_REAL8,voisins(3),102,comm_cart,statinfo)
       call MPI_RECV(D2,My,MPI_REAL8,voisins(3),101,comm_cart,status,statinfo)

 !call write_data3(rank,voisins,map,nb_probleme,G2,Mx,My,dx,dy,int2char(rank))
 !call write_data2(rank,voisins,map,nb_probleme,D2,Mx,My,dx,dy,int2char(rank))

       do j=n_y1,n_yn
          F(Mx*(j-n_y1)+1)  = F(Mx*(j-n_y1)+1)  + (D/dx**2)*G2(j-n_y1 +1)
          F(Mx*(j+1-n_y1) ) = F(Mx*(j+1-n_y1) ) + (D/dx**2)*D2(j-n_y1 +1)
       end do
   end if




!!!!!!!! COMM HAUT-BAS


!!! Le bloc fait toute la ligne

    if ( Voisins(4) < 0 .and. Voisins(2) < 0 ) then
       do i=n_x1,n_xn
          F(i-n_x1+1 )  = F(i-n_x1+1)  + (D/dy**2)*fonction_g(i*dx,0.0,nb_probleme)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dy**2)*fonction_g(i*dx,1.,nb_probleme)
       end do
 

!!! Condition de bord bas : g(x,0) + Com' en haut

    else if (Voisins(4) < 0) then

       ! Com' en haut : Envoyer un vecteur H1 de taille Mx à voisins(2)
       ! Puis récupérer le vecteur B1 (nouvellement H2) de taille My de voisins(2)
       do i=1,Mx
           H1(i) = U( Mx*(My-1-recouv) + i)
       end do
       call MPI_SEND(H1,Mx,MPI_REAL8,voisins(2),103,comm_cart,statinfo)
       call MPI_RECV(H2,Mx,MPI_REAL8,voisins(2),104,comm_cart,status,statinfo)


       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dy**2)*fonction_g(i*dx,0.0,nb_probleme)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dy**2)*H2(i-n_x1+1)
       end do


!!! Condition de bord haut : g(x,1) + Com' en bas

    else if (Voisins(2) < 0) then
       ! Com' EN BAS : Envoyer un vecteur B1 de taille Mx à voisins(4)
       ! Puis récupérer le vecteur H1 (nouvellement B2) de taille My de voisins(4)
       do i=1,Mx
           B1(i) = U(i + Mx*recouv)
       end do
       call MPI_SEND(B1,Mx,MPI_REAL8,voisins(4),104,comm_cart,statinfo)
       call MPI_RECV(B2,Mx,MPI_REAL8,voisins(4),103,comm_cart,status,statinfo)

       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dy**2)*B2(i-n_x1+1)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dy**2)*fonction_g(i*dx,1.,nb_probleme)
       end do

!!! Milieu du maillage : full com' Haut/bas

    else

       ! Com' en haut : Envoyer un vecteur H1 de taille Mx à voisins(2)
       ! Puis récupérer le vecteur B1 (nouvellement H2) de taille My de voisins(2)
       do i=1,Mx
           H1(i) = U( Mx*(My-1-recouv) + i)
       end do

       ! Com' EN BAS : Envoyer un vecteur B1 de taille Mx à voisins(4)
       ! Puis récupérer le vecteur B1 (nouvellement B2) de taille My de voisins(4)
       do i=1,Mx
           B1(i) = U(i + Mx*recouv)
       end do

       call MPI_SEND(H1,Mx,MPI_REAL8,voisins(2),103,comm_cart,statinfo)
       call MPI_SEND(B1,Mx,MPI_REAL8,voisins(4),104,comm_cart,statinfo)

       call MPI_RECV(H2,Mx,MPI_REAL8,voisins(2),104,comm_cart,status,statinfo)
       call MPI_RECV(B2,Mx,MPI_REAL8,voisins(4),103,comm_cart,status,statinfo)



       do i=n_x1,n_xn
          F(i-n_x1+1)  = F(i-n_x1+1)  + (D/dy**2)*B2(i-n_x1+1)
          F(Mx*(My-1) + i-n_x1 +1) = F(Mx*(My-1) + i-n_x1 +1) + (D/dy**2)*H2(i-n_x1+1)
       end do
 end if

    F=Dt*F

    deallocate(H1,B1,G1,D1,H2,B2,G2,D2)

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
