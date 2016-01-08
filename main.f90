program Chaleur_2D_Seqentiel

  use solve
  use mathfunctions
  use partition

  implicit none
  include "mpif.h"
  
  integer::Nx,Ny,Mx,My
  real*8::Lx,Ly,D,Dt,dx,dy,Tmax
  real*8,dimension(:,:),allocatable::A
  real*8,dimension(:),allocatable::U0,U,F,B
  integer::i,j,nb_iter,nb_probleme

  integer,dimension(4) :: map,voisins
! voising :  gauche haut droite bas
  integer::cart_ndims = 2
  integer,dimension(2)::cart_dims,cart_periods,coordinates
  integer size, rank, group, comm_cart, statinfo
  real*8::t1,t2

  call MPI_INIT(statinfo)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, statinfo)
  call MPI_COMM_GROUP(MPI_COMM_WORLD, group, statinfo)

  ! Creation nouveau communicateur
  call MPI_Comm_create(MPI_COMM_WORLD, group, comm_cart , statinfo)

  !2 dimensions, ? procs pour chaque dimensions, pas de périodicités, reorder -> NOPE
  cart_dims(1) = 0
  cart_dims(2) = 0
  cart_periods = 0


  call MPI_Dims_create(size, cart_ndims, cart_dims, statinfo)
  !print*,"cart_dims(1): ",cart_dims(1)
  !print*,"cart_dims(2): ",cart_dims(2)
  call MPI_Cart_create(MPI_COMM_WORLD, cart_ndims, cart_dims, cart_periods, 1, comm_cart, statinfo)

  coordinates(1) = 0
  coordinates(2) = 0
  call MPI_COMM_RANK(comm_cart, rank, statinfo) 

  call MPI_Cart_coords(comm_cart, rank, cart_ndims, coordinates, statinfo)
  !print*,"rank: ",rank,"coordinates(1): ",coordinates(1),"coordinates(2): ",coordinates(2)

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
  nb_iter=ceiling(Tmax/dt)

  call map_rect(Nx,Ny,cart_dims(1),cart_dims(2),coordinates(1),coordinates(2),2,map)
  !print*,map
  if (coordinates(2) == 0 ) then
     !print*, map(1),map(2),map(3),map(4),coordinates(1),coordinates(2)
  end if

  Mx = Map(2) - map(1) +1
  My = Map(4) - map(3) +1

  if (rank==2) then  
    print*,'coucou proc2'
    print*,Voisins
    print*,'t ou ?'
    print*,coordinates
  end if

  allocate(U(Nx*Ny),U0(Nx*Ny))
  U0=0
  !print*,'Nx',Nx,'Ny',Ny,'dx',dx,'dy',dy,'Cas',nb_probleme

  !Boucle principale

  call CPU_TIME(t1)
  do i=1, 1
     call Get_F(F,Mx,My,map(1),map(2),map(3),map(4),dx,dy,D,Dt,i*Dt,voisins,nb_probleme,rank)
     !F=F+U0
     U=0
     !call Grad_conj_implicit(U,F,0.001_wp,1000,Nx,Ny,dx,dy,D,Dt)
     U0=U
  end do
  call CPU_TIME(t2)

  !call save_result(U,Nx,Ny,dx,dy,"resultatseq.dat")

  deallocate(F)



  !print*,"Temps d'execution : ",t2-t1

  call MPI_FINALIZE(statinfo)

end program Chaleur_2D_Seqentiel
