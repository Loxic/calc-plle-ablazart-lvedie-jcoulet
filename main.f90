program Chaleur_2D_Seqentiel

  use solve
  use patate

  implicit none
  integer::Nx,Ny
  real(wp)::Lx,Ly,D,Dt,dx,dy,Tmax
  real(wp),dimension(:,:),allocatable::A
  real(wp),dimension(:),allocatable::U0,U,F,B
  integer::i,j,nb_iter,nb_probleme

  real*8::t1,t2


  call CPU_TIME(t1)

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

  do i=1, nb_iter
     call Get_F(F,Nx,Ny,dx,dy,D,Dt,i*dt,nb_probleme)
     F=F+U0
     U=0
     call Grad_conj_implicit(U,F,0.001_wp,1000,Nx,Ny,dx,dy,D,Dt)
     U0=U
  end do

  call save_result(U,Nx,Ny,dx,dy,"resultatseq.dat")

  deallocate(F)

  call CPU_TIME(t2)

  print*,"Temps d'execution : ",t2-t1

end program Chaleur_2D_Seqentiel
