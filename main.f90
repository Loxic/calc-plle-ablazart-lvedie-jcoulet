program Chaleur_2D_Seqentiel

  use solve
  use patate

  implicit none
  include "mpif.h"
  
  integer::Nx,Ny
  real(wp)::Lx,Ly,D,Dt,dx,dy,Tmax
  real(wp),dimension(:,:),allocatable::A
  real(wp),dimension(:),allocatable::U0,U,F,B
  integer::i,j,nb_iter,nb_probleme

  integer::cart_ndims = 2
  integer,dimension(2)::cart_dims,cart_periods 
  integer size, rank, group, cart, statinfo

  real*8::t1,t2

  call MPI_INIT(statinfo)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, statinfo) 
  call MPI_COMM_GROUP(MPI_COMM_WORLD, group, statinfo)

  ! Creation nouveau communicateur
  call MPI_Comm_create(MPI_COMM_WORLD, group, cart , statinfo)

  !2 dimensions, ? procs pour chaque dimensions, pas de périodicités, reorder -> NOPE

  call MPI_Dims_create(size, cart_ndims, cart_dims, statinfo)
  call MPI_Cart_create(MPI_COMM_WORLD, cart_ndims, cart_dims, cart_periods, 0, statinfo)


  nb_probleme=3 !Cas à résoudre
  !Lecture des paramètres
  !Paramètres spatiaux
  open (unit=11, file='parametres')
  read(11,*) Nx, Ny, Lx, Ly, D
  close(11)
  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)
  !Paramètres temporels
  Dt=1.0_wp
  Tmax=10.0_wp
  nb_iter=ceiling(Tmax/dt)


  allocate(U(Nx*Ny),U0(Nx*Ny))
  U0=0
  print*,'Nx',Nx,'Ny',Ny,'dx',dx,'dy',dy,'Cas',nb_probleme

  !Boucle principale

  call CPU_TIME(t1)
  do i=1, nb_iter
     call Get_F(F,Nx,Ny,dx,dy,D,Dt,i*dt,nb_probleme)
     F=F+U0
     U=0
     call Grad_conj_implicit(U,F,0.001_wp,1000,Nx,Ny,dx,dy,D,Dt)
     U0=U
  end do
  call CPU_TIME(t2)

  call save_result(U,Nx,Ny,dx,dy,"resultatseq.dat")

  deallocate(F)



  print*,"Temps d'execution : ",t2-t1

  call MPI_FINALIZE(statinfo)

end program Chaleur_2D_Seqentiel
