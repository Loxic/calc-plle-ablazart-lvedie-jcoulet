program Chaleur_2D_Seqentiel

  use solve
  use mathfunctions
  use partition
  use display

  implicit none
  include "mpif.h"

  integer::Nx,Ny,Mx,My
  real*8::Lx,Ly,D,Dt,dx,dy,Tmax,time
  real*8,dimension(:,:),allocatable::A
  real*8,dimension(:),allocatable::U0,U,F,B,Ubord,Ubord2
  integer::i,j,nb_iter,nb_probleme,recouv,Ligne,Colonne,Ligne_Y,color,dumb
  integer, dimension(MPI_STATUS_SIZE)::status

  integer,dimension(4) :: map,voisins
  ! voising :  gauche haut droite bas
  integer::cart_ndims = 2
  integer,dimension(2)::cart_dims,cart_periods,coordinates
  integer:: size, rank, group, comm_cart, statinfo
  real*8::t1,t2

  real*8,dimension(5) :: Param_real
  integer,dimension(7) :: Param_int
  integer,dimension(10) :: Param_mpi

  cart_dims(1) = 0 ; cart_dims(2) = 0 ; cart_periods = 0
  coordinates(1) = 0 ;  coordinates(2) = 0

  call MPI_INIT(statinfo)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, statinfo)
  call MPI_COMM_GROUP(MPI_COMM_WORLD, group, statinfo)

  ! Creation nouvelle numérotation proc
  call MPI_Comm_create(MPI_COMM_WORLD, group, comm_cart , statinfo)

  ! Agencement en grille cartésienne, 2 dimensions, ? procs pour chaque dimensions, pas de périodicités, reorder -> NO
  call MPI_Dims_create(size, cart_ndims, cart_dims, statinfo)
  call MPI_Cart_create(MPI_COMM_WORLD, cart_ndims, cart_dims, cart_periods, 1, comm_cart, statinfo)

  ! Assignation de la numérotation globale à rank, et de sa numérotation (x,y) à coordinates
  call MPI_COMM_RANK(comm_cart, rank, statinfo) 
  call MPI_Cart_coords(comm_cart, rank, cart_ndims, coordinates, statinfo)

  ! Assignation des numéros globaux des proc voisins au tableau voisins [Ouest - Nord - Est - Sud]
  call MPI_CART_SHIFT(comm_cart,0,1,voisins(1),voisins(3),statinfo) ! gauche/droie
  call MPI_CART_SHIFT(comm_cart,1,1,voisins(4),voisins(2),statinfo) ! bas/haut



  nb_probleme=3 !Cas à résoudre
  !Lecture des paramètres
  !Paramètres spatiaux
  open (unit=11, file='parametres')
  read(11,*) Nx, Ny, Lx, Ly, D
  close(11)
  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)
  !Paramètres temporels
  Dt=1.0
  Tmax=10.0
  D = 1.
  nb_iter=ceiling(Tmax/dt) !Utile pour problème 3, sinon nb_iter = 1 car stationnaire
  !nb_iter=1
  recouv = 2

  ! mapping du domaine sur la grille des procs
  call map_rect(Nx,Ny,cart_dims(1),cart_dims(2),coordinates(1),coordinates(2),recouv,map)

  Mx = Map(2) - map(1) +1
  My = Map(4) - map(3) +1

  ! DAMIER
  color = -1
  !call damier(rank,coordinates(1),coordinates(2),size,color)

  allocate(F(Mx*My),U(Mx*My),U0(Mx*My),Ubord(2*Mx+2*My),Ubord2(2*Mx+2*My))
  F=0 ; U=0; U0=0 ; time=0 ;

  ! Création des types colonnes et lignes
  call MPI_TYPE_CONTIGUOUS(Mx,MPI_REAL8,Ligne,statinfo)
  call MPI_TYPE_COMMIT(Ligne,statinfo)

  call MPI_TYPE_VECTOR(My,1,Mx,MPI_REAL8,Colonne,statinfo)
  call MPI_TYPE_COMMIT(Colonne,statinfo)

  call MPI_TYPE_CONTIGUOUS(My,MPI_REAL8,Ligne_Y,statinfo)
  call MPI_TYPE_COMMIT(Ligne_Y,statinfo)



  Param_real = (/dx,dy,D,Dt,time/)
  Param_int = (/Mx,My,nb_probleme,map(1),map(2),map(3),map(4)/)
  Param_mpi = (/rank,recouv,comm_cart,voisins(1),voisins(2),voisins(3),voisins(4),Ligne,Colonne,Ligne_Y/)



  !Boucle principale - additif
  if (color > -1) then
     call CPU_TIME(t1)
     do i=1,1
        if (rank==0) then
           print*,'Itération',i,'sur',nb_iter
        end if
        Param_real(5) = Param_real(5) + dt
        do j = 1,1
           if (color == 0) then
              call Get_FB(F,Param_real,Param_int,Param_mpi)
              F=F+U0
              U=0
              call Sparse_solve(3,Mx,My,dx,dy,D,Dt,U,F,0.001,1000)
              U0=U
           else
              call Get_F_CL(U0,Mx,My,Param_mpi)
           end if
           color = merge(2,color,color==0) - 1

        end do
     end do
     call CPU_TIME(t2)


  else
     call CPU_TIME(t1)
     do i=1, nb_iter
        Ubord = -5; Ubord2 = 42; j=0 ! Pour passer 1er test
        Param_real(5) = Param_real(5) + dt

        do while ( .not. conv_schwartz(Ubord,Ubord2,Mx,My,0.001) )
           Ubord = Ubord2 
           call Get_F(F,U0,Param_real,Param_int,Param_mpi,Ubord2)
           F=F+U0
           U=0
           call Sparse_solve(3,Mx,My,dx,dy,D,Dt,U,F,0.001,1000)
           U0=U
           j = j+1
        end do
        if (rank==0) then
          print*,'Itération',i,'sur',nb_iter,' terminé avec ',j,' sous itérations pour Schwartz additif'
        end if
     end do
     call CPU_TIME(t2)

  end if

  call write_data(voisins,map,nb_probleme,U,Mx,My,dx,dy,int2char(rank))


  if (rank==1) then
     call script_gnu2(size,nb_probleme)
  end if
  call script_gnuplot(rank,size,nb_probleme,color,MPI_COMM_WORLD)


  deallocate(F,U0,U)



  !print*,"Temps d'execution : ",t2-t1

  call MPI_FINALIZE(statinfo)

end program Chaleur_2D_Seqentiel
