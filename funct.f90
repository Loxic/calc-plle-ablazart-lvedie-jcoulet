module mathfunctions
  
  implicit none


contains



!!!!! Prise en compte des CL sur sous maillage

  ! => maillage de taille Nx . Ny 
  ! => sous maillage de taille Mx . My 

!!!!!! Recouvrement

  ! 'CL' plus larges : au lieu de prendre juste l'indice du bord à côté, on en prend une tranche, et on l'inverse


  subroutine Get_F(F,U,Param_real,Param_int,Param_mpi,Ubord)

    implicit none
    include "mpif.h"

    real*8,dimension(5),intent(in) :: Param_real
    integer,dimension(7),intent(in) :: Param_int
    integer,dimension(10),intent(in) :: Param_mpi
    real*8,dimension(Param_int(1)*Param_int(2)),intent(in):: U
    real*8,dimension(Param_int(1)*Param_int(2)),intent(inout):: F
    real*8,dimension(2*(Param_int(1)+Param_int(2))),intent(out)::Ubord

    integer:: Mx,My,nb_probleme,n_x1,n_xn,n_y1,n_yn
    integer,dimension(4):: voisins
    integer:: rank,recouv,comm_cart,Ligne,Colonne,Ligne_Y

    real*8,dimension(Param_int(1)):: H2,B2
    real*8,dimension(Param_int(2)):: G2,D2
    real*8:: dx,dy,Dt,D,time
    integer::i,j,statinfo
    integer, dimension(MPI_STATUS_SIZE) :: status


    ! Param_real = (/dx,dy,D,Dt,time/)
    dx = Param_real(1) ; dy = Param_real(2) ; D = Param_real (3)
    Dt = Param_real (4) ; time = Param_real(5)

    !  Param_int = (/Mx,My,nb_probleme,map(1),map(2),map(3),map(4)/)
    Mx = Param_int(1) ; My = Param_int(2) ; nb_probleme = Param_int(3)
    n_x1 = Param_int(4) ; n_xn = Param_int(5) ; n_y1 = Param_int(6) ; n_yn = Param_int(7)

    !  Param_mpi = (/rank,recouv,comm_cart,voisins(1),voisins(2),voisins(3),voisins(4),Ligne,Colonne,Ligne_Y/)
    rank = Param_mpi(1) ; recouv = Param_mpi(2) ; comm_cart = Param_mpi(3) ; voisins = Param_mpi(4:7)
    Ligne = Param_mpi(8) ; Colonne = Param_mpi(9) ; Ligne_Y = Param_mpi(10) ; 

    H2 = 0 ;   B2 = 0 ;    G2 = 0 ;    D2 = 0
    F=0
    Ubord = 0

    do j=n_y1,n_yn
       do i=n_x1,n_xn
          F( i-n_x1 +1+ Mx*(j-n_y1) )=fonction_f(i*dx,j*dy,1.,1.,time,nb_probleme)
       end do
    end do



!!!!!!!! COMM GAUCHE-DROITE


!!! Le bloc fait toute la rangée (pas de com')

    if ( Voisins(3) < 0 .and. Voisins(1) < 0 ) then
       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*fonction_h(0.0,(j-1+n_y1)*dy,nb_probleme)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*fonction_h(1.,(j-1+n_y1)*dy,nb_probleme)
       end do


!!! Condition de bord droit : h(1,y) + Com' à gauche

    else if (Voisins(3) < 0) then

       ! Com' à gauche : Envoyer un vecteur de taille My à voisins(1)
       ! Puis récupérer le vecteur G2 de taille My de voisins(1)

       call MPI_SEND(U(1+recouv),1,Colonne,voisins(1),101,comm_cart,statinfo)       
       call MPI_RECV(G2,1,Ligne_Y,voisins(1),102,comm_cart,status,statinfo)
       Ubord(1:My) = G2

       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*G2(j)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*fonction_h(1.,(j-1+n_y1)*dy,nb_probleme)
       end do


!!! Condition de bord gauche : h(0,y) + Com' à droite

    else if (Voisins(1) < 0) then

       ! Com' à droite : Envoyer un vecteur de taille My à voisins(3)
       ! Puis récupérer le vecteur D2 de taille My de voisins(3)

       call MPI_SEND(U(Mx-recouv),1,Colonne,voisins(3),102,comm_cart,statinfo)
       call MPI_RECV(D2,1,Ligne_Y,voisins(3),101,comm_cart,status,statinfo)
       Ubord(Mx+My+1:Mx+2*My) = D2

       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*fonction_h(0.0,(j-1+n_y1)*dy,nb_probleme)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*D2(j)
       end do


!!! Milieu du maillage : full com' Gauche/droite

    else

       ! Com' à gauche : Envoyer un vecteur de taille My à voisins(1)
       ! Puis récupérer le vecteur G2 de taille My de voisins(1)

       ! Com' à droite : Envoyer un vecteur de taille My à voisins(3)
       ! Puis récupérer le vecteur D2 de taille My de voisins(3)

       call MPI_SEND(U(1+recouv),1,Colonne,voisins(1),101,comm_cart,statinfo)
       call MPI_SEND(U(Mx-recouv),1,Colonne,voisins(3),102,comm_cart,statinfo)

       call MPI_RECV(G2,1,Ligne_Y,voisins(1),102,comm_cart,status,statinfo)
       call MPI_RECV(D2,1,Ligne_Y,voisins(3),101,comm_cart,status,statinfo)

       Ubord(1:My) = G2
       Ubord(Mx+My+1:Mx+2*My) = D2

       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*G2(j)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*D2(j)
       end do
    end if



!!!!!!!! COMM HAUT-BAS


!!! Le bloc fait toute la ligne

    if ( Voisins(4) < 0 .and. Voisins(2) < 0 ) then
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*fonction_g((i-1+n_x1)*dx,0.0,nb_probleme)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*fonction_g((i-1+n_x1)*dx,1.,nb_probleme)
       end do


!!! Condition de bord bas : g(x,0) + Com' en haut

    else if (Voisins(4) < 0) then

       ! Com' en haut : Envoyer un vecteur de taille Mx à voisins(2)
       ! Puis récupérer le vecteur H2 de taille My de voisins(2)

       call MPI_SEND(U(Mx*(My-1-recouv) +1),1,Ligne,voisins(2),103,comm_cart,statinfo)
       call MPI_RECV(H2,1,Ligne,voisins(2),104,comm_cart,status,statinfo)

       Ubord(My+1:Mx+My) = H2 
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*fonction_g((i-1+n_x1)*dx,0.0,nb_probleme)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*H2(i)
       end do


!!! Condition de bord haut : g(x,1) + Com' en bas

    else if (Voisins(2) < 0) then
       ! Com' EN BAS : Envoyer un vecteur de taille Mx à voisins(4)
       ! Puis récupérer le vecteur B2 de taille My de voisins(4)

       call MPI_SEND(U(Mx*recouv+1),1,Ligne,voisins(4),104,comm_cart,statinfo)
       call MPI_RECV(B2,1,Ligne,voisins(4),103,comm_cart,status,statinfo)

       Ubord(Mx+2*My+1:2*Mx+2*My) = B2
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*B2(i)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*fonction_g((i-1+n_x1)*dx,1.,nb_probleme)
       end do

!!! Milieu du maillage : full com' Haut/bas

    else

       ! Com' en haut : Envoyer un vecteur taille Mx à voisins(2)
       ! Puis récupérer le vecteur H2 de taille My de voisins(2)

       ! Com' EN BAS : Envoyer un vecteur de taille Mx à voisins(4)
       ! Puis récupérer le vecteur B2 de taille My de voisins(4)

       call MPI_SEND(U(Mx*(My-1-recouv) +1),1,Ligne,voisins(2),103,comm_cart,statinfo)
       call MPI_SEND(U(Mx*recouv+1),1,Ligne,voisins(4),104,comm_cart,statinfo)

       call MPI_RECV(H2,1,Ligne,voisins(2),104,comm_cart,status,statinfo)
       call MPI_RECV(B2,1,Ligne,voisins(4),103,comm_cart,status,statinfo)

       Ubord(My+1:Mx+My) = H2 
       Ubord(Mx+2*My+1:2*Mx+2*My) = B2

       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*B2(i)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*H2(i)
       end do
    end if

    F=Dt*F

  end subroutine Get_F




  subroutine Get_FB(F,Param_real,Param_int,Param_mpi)

    implicit none
    include "mpif.h"

    real*8,dimension(5),intent(in) :: Param_real
    integer,dimension(7),intent(in) :: Param_int
    integer,dimension(10),intent(in) :: Param_mpi
    real*8,dimension(Param_int(1)*Param_int(2)),intent(inout):: F

    integer:: Mx,My,nb_probleme,n_x1,n_xn,n_y1,n_yn
    integer,dimension(4):: voisins
    integer:: rank,recouv,comm_cart,Ligne,Colonne,Ligne_Y

    real*8,dimension(Param_int(1)):: H2,B2
    real*8,dimension(Param_int(2)):: G2,D2
    real*8:: dx,dy,Dt,D,time
    integer::i,j,statinfo
    integer, dimension(MPI_STATUS_SIZE) :: status


    dx = Param_real(1) ; dy = Param_real(2) ; D = Param_real (3)
    Dt = Param_real (4) ; time = Param_real(5)

    Mx = Param_int(1) ; My = Param_int(2) ; nb_probleme = Param_int(3)
    n_x1 = Param_int(4) ; n_xn = Param_int(5) ; n_y1 = Param_int(6) ; n_yn = Param_int(7)

    rank = Param_mpi(1) ; recouv = Param_mpi(2) ; comm_cart = Param_mpi(3) ; voisins = Param_mpi(4:7)
    Ligne = Param_mpi(8) ; Colonne = Param_mpi(9) ; Ligne_Y = Param_mpi(10) ; 

    H2 = 0 ;   B2 = 0 ;    G2 = 0 ;    D2 = 0
    F=0

    do j=n_y1,n_yn
       do i=n_x1,n_xn
          F( i-n_x1 +1+ Mx*(j-n_y1) )=fonction_f(i*dx,j*dy,1.,1.,time,nb_probleme)
       end do
    end do


!!!!!!!! COMM GAUCHE-DROITE

    if ( Voisins(3) < 0 .and. Voisins(1) < 0 ) then
       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*fonction_h(0.0,(j-1+n_y1)*dy,nb_probleme)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*fonction_h(1.,(j-1+n_y1)*dy,nb_probleme)
       end do

    else if (Voisins(3) < 0) then
       call MPI_RECV(G2,1,Ligne_Y,voisins(1),102,comm_cart,status,statinfo)
       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*G2(j)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*fonction_h(1.,(j-1+n_y1)*dy,nb_probleme)
       end do

    else if (Voisins(1) < 0) then
       call MPI_RECV(D2,1,Ligne_Y,voisins(3),101,comm_cart,status,statinfo)
       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*fonction_h(0.0,(j-1+n_y1)*dy,nb_probleme)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*D2(j)
       end do

    else
       call MPI_RECV(G2,1,Ligne_Y,voisins(1),102,comm_cart,status,statinfo)
       call MPI_RECV(D2,1,Ligne_Y,voisins(3),101,comm_cart,status,statinfo)
       do j=1,My
          F(Mx*(j-1)+1) = F(Mx*(j-1)+1) + (D/dx**2)*G2(j)
          F(Mx*j)  = F(Mx*j)  + (D/dx**2)*D2(j)
       end do
    end if


!!!!!!!! COMM HAUT-BAS

    if ( Voisins(4) < 0 .and. Voisins(2) < 0 ) then
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*fonction_g((i-1+n_x1)*dx,0.0,nb_probleme)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*fonction_g((i-1+n_x1)*dx,1.,nb_probleme)
       end do

    else if (Voisins(4) < 0) then
       call MPI_RECV(H2,1,Ligne,voisins(2),104,comm_cart,status,statinfo)
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*fonction_g((i-1+n_x1)*dx,0.0,nb_probleme)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*H2(i)
       end do

    else if (Voisins(2) < 0) then

       call MPI_RECV(B2,1,Ligne,voisins(4),103,comm_cart,status,statinfo)
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*B2(i)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*fonction_g((i-1+n_x1)*dx,1.,nb_probleme)
       end do

    else
       call MPI_RECV(H2,1,Ligne,voisins(2),104,comm_cart,status,statinfo)
       call MPI_RECV(B2,1,Ligne,voisins(4),103,comm_cart,status,statinfo)
       do i=1,Mx
          F(i)  = F(i)  + (D/dy**2)*B2(i)
          F(Mx*(My-1)+i) = F(Mx*(My-1)+i) + (D/dy**2)*H2(i)
       end do
    end if

    F=Dt*F

  end subroutine Get_FB



  subroutine Get_F_CL(U,Mx,My,Param_mpi)

    implicit none
    include "mpif.h"
    integer,intent(in):: Mx,My
    real*8,dimension(Mx*My),intent(in):: U
    integer,dimension(10),intent(in):: Param_mpi

    integer,dimension(4):: voisins
    integer:: rank,recouv,comm_cart,Ligne,Colonne,Ligne_Y
    integer:: statinfo

    rank = Param_mpi(1) ; recouv = Param_mpi(2) ; comm_cart = Param_mpi(3) ; voisins = Param_mpi(4:7)
    Ligne = Param_mpi(8) ; Colonne = Param_mpi(9) ; Ligne_Y = Param_mpi(10) ; 

!!!!!!!! COMM GAUCHE-DROITE

    if (Voisins(1) >= 0) then
       call MPI_SEND(U(1+recouv),1,Colonne,voisins(1),101,comm_cart,statinfo)   
    end if

    if (Voisins(3) >= 0) then

       call MPI_SEND(U(Mx-recouv),1,Colonne,voisins(3),102,comm_cart,statinfo)    
    end if

!!!!!!!! COMM HAUT-BAS

    if (Voisins(2) >= 0) then
       call MPI_SEND(U(Mx*(My-1-recouv) +1),1,Ligne,voisins(2),103,comm_cart,statinfo)
    end if

    if (Voisins(4) >= 0) then
       call MPI_SEND(U(Mx*recouv+1),1,Ligne,voisins(4),104,comm_cart,statinfo)
    end if


  end subroutine Get_F_CL



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
    do j = 1, Ny
       do i = 1, Nx
          k = Nx*(j-1) + i
          ! bloc M
          AX(k) = (1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2))*X(k)
          if (i>1) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k-1)
          end if
          if (i<Nx) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k+1)
          end if

          !Bloc E inférieur
          if(j>1) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k-Nx)
          end if
          !Bloc E supérieur
          if (j<Ny) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k+Nx)
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

  function conv_schwartz(V,Vnext,Mx,My,tol)
    implicit none
    include "mpif.h"
    real*8,dimension(:),intent(in)::V,Vnext
    integer,intent(in)::Mx,My
    real*8,intent(in)::tol
    logical::conv_schwartz
    logical::conv_local
    integer::statinfo
    
    conv_local = ( norme2(V(1:My)-Vnext(1:My)) < tol ) .and. &
                 ( norme2(V(My+1:My+Mx)-Vnext(My+1:My+Mx)) < tol ) .and. &
                 ( norme2(V(Mx+My+1:Mx+2*My)-Vnext(My+Mx+1:Mx+2*My)) < tol ) .and. &
                 ( norme2(V(2*My+Mx+1:2*(Mx+My))-Vnext(2*My+Mx+1:2*(My+Mx))) < tol ) 
    call MPI_ALLREDUCE(conv_local,conv_schwartz,1,mpi_logical,MPI_LAND,mpi_comm_world,statinfo)

  end function

function norme2(V)
    implicit none
    real*8,dimension(:),intent(in)::V
    real*8::norme2
    norme2=sqrt(dot_product(V,V))
  end function norme2


end module
