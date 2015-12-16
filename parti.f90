module partition

  implicit none

  contains






! Partionnement en n.m maillages rectangulaires
! Surcharge au bord
  subroutine charge(me, n, Np, i1, in)

    implicit none

    integer, intent(in) :: me, n, Np
    integer :: C,D
    integer, intent(out) :: i1,in

    C = n/Np
    if (c == (n*1.0)/Np) then
       ! Pas de surcharge
       i1 = c*me +1
       in = c*(me+1)

    else
       D = mod(n,Np)

       if (me < Np - D ) then
          i1 = c*me +1 
          in = c*(me+1)
       else
          i1 = c*me + (me - Np + D +1)
          ! Charge 'standard + surcharge répartie uniformément
          in = c*(me+1) + (me - Np + D +1)
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
  if (px < n-1) then 
    map(2) = map(2) + recouv
  end if

  if (py < m-1) then 
    map(4) = map(4) + recouv
  end if


end subroutine map_rect


end module
