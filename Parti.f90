module partition

  implicit none

  contains

! Partionnement en n² maillages carrés
subroutine map_carre(Nx,Ny,n,Mx,My,map)
  integer,intent(in)::Nx,Ny,n
  integer,intent(out)::Mx,My
  integer,dimension(n*n,4)::map
  integer :: i,j
  
 ! sur les colonnes :

  Mx = Nx/n
  My = Ny/n

  do i=1,n
    do j=1,n
       map((i-1)*n + j,1) = (i-1)*Mx +1
       map((i-1)*n + j,2) = (j-1)*My +1
    end do
  end do
  map(:,3) = Mx
  map(:,4) = My

end subroutine map_carre



! Partionnement en n.m maillages rectangulaires
! Surcharge au bord
subroutine map_rect(Nx,Ny,n,Mx,My,Mx_bord,My_bord)
  integer,intent(in)::Nx,Ny,n
  integer,intent(out)::Mx,My
  integer,dimension(n*n,4)::map
  integer :: i,j
  
  map(:,3) = Mx
  map(:,4) = Mx
  Mx = floor(Nx/n)
  My = floor(Ny/m)
  Mx_bord = Nx - Mx*n
  My_bord = Ny - My*n
  

  do i=1,m
    do j=1,n
       map((i-1)*n + j,1) = (i-1)*Mx +1
       map((i-1)*n + j,2) = (j-1)*My +1
    end do
  end do

  map(:,3) = Mx
  map(:,4) = My

  do j=1,n
    map((m-1)*n + j :,3) = Mx_bord
  end do
  do i=1,m
    map((i-1)*n + m :,3) = My_bord
  end do


end subroutine map_rect


end program Chaleur_2D_Seqentiel
