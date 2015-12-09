module partition

  implicit none

  contains

! Partionnement en n² maillages carrés
subroutine map_carre(Nx,Ny,n,recouv,Mx,My,map)
  integer,intent(in)::Nx,Ny,n
  integer,intent(out)::Mx,My
  integer,intent(in) :: recouv
  integer,dimension(n,n,4)::map
  integer :: i,j

  Mx = Nx/n
  My = Ny/n


!!! premiers indices de chaque sous-maillage

  do i=1,n
    do j=1,n
       map(i,j,1) = (i-1)*Mx +1
       map(i,j,2) = (j-1)*My +1
    end do
  end do

!!! taille des blocs

  map(:,:,3) = map(:,:,1) + Mx -1
  map(:,:,4) = map(:,:,2) + My -1

  map(1:n-1,:,3) = map(1:n-1,:,3) + recouv
  map(:,1:n-1,4) = map(:,1:n-1,4) + recouv


end subroutine map_carre



! Partionnement en n.m maillages rectangulaires
! Surcharge au bord
subroutine charge(me, N, Np, i1, in)

    implicit none
    integer, intent(in) :: me, N, Np
    integer :: C,D
    integer, intent(out) :: i1,in

    C = n/Np
    if (c == (n*1.0)/Np) then
       ! Pas de surcharge
       i1 = c*(me-1) +1
       in = c*me

    else
       D = mod(n,Np)

       if (me-1 < Np - D ) then
          i1 = c*(me-1) +1 
          in = c*me
       else
          i1 = c*(me-1) + (me - Np + D)
          ! Charge 'standard + surcharge répartie uniformément
          in = c*me + (me - Np + D)
       end if
    end if

  end subroutine charge



subroutine map_rect(Nx,Ny,n,m,px,py,recouv,map)
  integer,intent(in)::Nx,Ny,n,m,px,py
  integer::Mx,My
  integer,intent(in) :: recouv
  integer,dimension(4),intent(out)::map
  
  call charge(px, Nx, n, map(1), map(2))
  call charge(py, Ny, m, map(3), map(4))

  
! Traitement de la longueur au bord
  if (px < n) then 
    map(2) = map(2) + recouv
  end if

  if (py < m) then 
    map(4) = map(4) + recouv
  end if


end subroutine map_rect


end module
