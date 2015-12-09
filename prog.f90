program prog

  use fonctions
  use read_write
  use algo_num

  implicit none
  include "mpif.h"

!!!!!!!!!!!!!!!!!!!!!MODIF 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: Nx, Ny, N, i, j, k, nstat, nbiter, niter, choix, Np, Me, statinfo, i1, in, q, he, i2, j2,m1,m2
  integer, dimension(MPI_STATUS_SIZE) :: status

  real :: start, finish
  real*8 :: D, dt, dx, dy, Lx, Ly, eps, d1, d2, d3
  real*8, dimension(:), allocatable :: U, F, B, U0, Urf
  character(len=10) :: Fichier
  character*14 :: tab
  character*2 :: Mn




  Call MPI_INIT(statinfo)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, Me, statinfo) 


  if (Me==0) then
     call system("rm -r Sol*")
  end if

  call read_param(Ny,Nx,Ly,Lx,D,dt,nbiter,eps)


  tab = '(100g10.4)'
  choix=0
  N = Ny*Nx

  Call charge(me, Ny, Np, i1, in)
  i1 = (i1-1)*Nx+1
  in = in*Nx

  nstat = 1

  choix=0
  if(me==0) then
     do while ((choix/=1).AND.(choix/=2).AND.(choix/=3))
        print*,'choix du problème (1/2/3)'
        read*, choix
        if ((choix/=1).AND.(choix/=2).AND.(choix/=3)) then
           print*,'Le probleme choisit doit etre 1, 2 ou 3'
        end if
     end do

  end if

  CALL MPI_BCAST(choix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, statinfo)

  CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
 ! print*,choix
  ! Voir fonctions. Les Pb 1/2 sont stationnaires, et ne s'itèrent qu'une fois. 
  ! Le 3 s'itère 100 fois (modfiable par nbiter).

  if (choix/=3) then
     dt = 1D0
     nstat = 0D0
     ! nstat = 0 => stationnaire
     nbiter = 1
  end if

!!!!!!!!!!!!!!!!!!!!!FIN MODIF 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! I_ CALCUL DE LA MATRICE

  m1 = Me/10
  m2 = Me - 10*m1
  Mn = char(m1+48)//char(m2+48)
  call system("mkdir Sol"//Mn)

  if(me==0) then
     call CPU_time(start)
  end if


  dx = Lx/(Nx+1D0)
  dy = Ly/(Ny+1D0)

  d1= - D*dt/(dx**2)
  d2= - D*dt/(dy**2)
  d3= 1D0*nstat+2D0*(-d1 -d2)
  ! nstat traduit la modification du système en stationnaire
  ! C'est tout !! Le produit Matrice/vecteur fait le reste

  ! II_ RESOLUTION

  !call CPU_time(start) 

  !  allocate (U(N), F(N), B(N), U0(N))
  allocate (U(i1:in), F(i1:in), B(i1:in), U0(i1:in))
  U0 = 0D0
  ! En stationnaire, nbiter est maintenu à 1


  do k=1,nbiter 

     U = 0D0
     F = 0D0

     Call Init(i1 , Nx, Ny, i, j)


     !Implémentation du second membre
     do q=i1,in !i1,in

        F(q) = F(q) + dt*ff(i*dx,j*dy,k*dt,choix,Ly,Lx)

        if (i/= Nx) then
           i = i+1
        else
           i = 1
           j = j + 1
        end if

     end do

     !Condition de bord en y=0 et au bord opposé

     if (Me==0) then

        do q = 1,Nx
           F(q) = F(q) + D*dt*g(q*dx,0D0,choix,Ly,Lx)/(dy**2)
        end do

     else if (Me==Np-1) then

        do q=N-Nx+1,N
!!!!!!!!!!!!!!!!!!!!!!!!! JE SAIS PAS POURQUOI, MAIS CA MARCHE
           F(q) = F(q) + D*dt*g((q-2-(Nx-Ny))*dx,(Ny+1)*dy,choix,Ly,Lx)/(dy**2)
        end do

     end if

     ! Condition de bord en y=0 et y = dy*(Nx+1)

     Call Init(i1 , Nx, Ny, i, j)

     do q=i1,in 

        !Implémentation
        if (i==1) then
           F(q) = F(q) + dt*D* h(0D0,j*dy,choix, Ly, Lx)/(dx**2)
        elseif (i== Nx) then
           F(q) = F(q) + dt*D* h(Lx,j*dy,choix, Ly, Lx)/(dx**2)
        end if

        !déroulement de i et j
        if (i/= Nx) then
           i = i+1
        else
           i = 1
           j = j + 1
        end if

     end do


     CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

!!$      print*,'F'
!!$     ! Décommenter pour print F
!!$      write(*,tab) ( F(j), j=1,N )

     ! Résolution A.U(n+1) = B(U(n))
     B = F + U0
     call Gradient_conjugue(d1,d2,d3, B, U, eps, Ny,Nx, niter,i1,iN,Me,Np)

     U0 = U

  end do


!!$  print*,'U'
!!$ ! Décommenter pour print U
!!$  write(*,tab) ( A(i,j), j=1,N )

!!!!!!!!!!!!  REDUCTION DU VECTEUR U DANS LE VECTEUR URF !!!!!!!!!!!!!!

  if (Me/=0) then

     call MPI_SEND(i1, 1, MPI_INTEGER, 0, 43, MPI_COMM_WORLD, statinfo)
     call MPI_SEND(in, 1, MPI_INTEGER, 0, 44, MPI_COMM_WORLD, statinfo)
     call MPI_SEND(U(i1:in), in-i1+1, MPI_REAL8, 0, 42, MPI_COMM_WORLD, statinfo)

  elseif (Me==0) then

     allocate(Urf(N))

     Urf(i1:in) = U(i1:in)

     do He = 1, Np-1
        call MPI_RECV(i2, 1, MPI_INTEGER, He, 43, MPI_COMM_WORLD, status, statinfo)
        call MPI_RECV(j2, 1, MPI_INTEGER, He, 44, MPI_COMM_WORLD, status, statinfo)
        call MPI_RECV(Urf(i2:j2), j2-i2+1, MPI_REAL8, He, 42, MPI_COMM_WORLD, status, statinfo)
     end do

  end if



!!!!!!!!!!!!!!!! Methode 2 : Chaqun dans un fichier différent.


  ! Calcul de vitesse de l'algo
  if (Me==0) then
     call CPU_time(finish)
     print*,"Temps de calcul : ", finish - start
     call write_gnuplot(dx, dy, Nx, Ny, Urf, choix, "TESTOWA")
  end if


  ! Ecriture et affichage du résultat



     call Rename(Fichier,Me)
     call write_gnuplot2(Me,Np,i1,in,dx, dy, Nx, Ny, U, choix, Fichier)

     CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

  if (Me ==0) then
     print*,"Appuyez sur Enter pour continuer"
     call system("cat Sol-*.dat > SolFinale.dat")
     call system("gnuplot gnuplot.gnu")
  end if

  deallocate(B,F,U0,U)

  if (Me==0) then
     deallocate(Urf)
  end if

  call MPI_FINALIZE(statinfo)

end program prog
