program TEST

  implicit none
  include "mpif.h"

  integer::Np,Me
  integer::statinfo
  integer,dimension(2)::dims_orig
  integer::comm_cart2d
  integer::rg_d,rg_g,rg_h,rg_b
  logical,dimension(2)::period=.false.
  call MPI_INIT(statinfo)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Me,statinfo)

  if (Me==0) then 
    print*,'Nb proc',Np
  end if
  dims_orig = 0
  call MPI_DIMS_CREATE(Np,2,dims_orig,statinfo)

  if(Me==0) then
   print*,"Decoupage en deux dimensions"
   print*,"Suivant x ",dims_orig(1)
   print*,"Suivant y ",dims_orig(2)
  end if

  call MPI_CART_CREATE(mpi_comm_world,2,dims_orig,&
      period,.true.,comm_cart2d,statinfo)

  if (Me==2) then  
    print*,period
    call MPI_CART_SHIFT(comm_cart2D,0,1,rg_g,rg_d,statinfo)
    print*,'a gauche',rg_g,'a droite',rg_d
    call MPI_CART_SHIFT(comm_cart2D,1,1,rg_b,rg_h,statinfo)
    print*,'en haut',rg_h,'en bas',rg_b
  end if
  
  !QQ notes
  !dans shift : comm,direction,size,Vgauche,Vdroit,statinfo
  ! direction = 0 pour x, 1 pour y
  ! Si y, Vgauche = bas et Vdroit = haut


  call MPI_FINALIZE(statinfo)

end program TEST
