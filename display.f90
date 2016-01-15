module display

  use mathfunctions
  use partition

  implicit none

contains

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

  function colors()
    implicit none
    character*25,dimension(0:19) :: colors
    colors(0) = "0'blue',1'blue'"
    colors(1) = "1'#006400',2'#006400'"
    colors(2) = "2'#00ff00',3'#00ff00'"
    colors(3) = "3'cyan',4'cyan'"
    colors(4) = "4'#ff9900',5'#ff9900'"
    colors(5) = "5'#990099',6'#990099'"
    colors(6) = "6'#cc6600',7'#cc6600'"
    colors(7) = "7'red',8'red'"
    colors(8) = "8'#000000',9'#000000'"
    colors(9) = "9'#ffff00',10'#ffff00'"
    colors(10) = "10'#66b2ff',11'#66b2ff'"
    colors(11) = "11'#ffe166',12'#ffe166'"
    colors(12) = "12'#9900ff',13'#9900ff'"
    colors(13) = "13'#008B8B',14'#008B8B'"
    colors(14) = "14'#FF1493',15'#FF1493'"
    colors(15) = "15'#999999',16'#999999'"
    colors(16) = "16'#D3D3D3',17'yellow'"
    colors(17) = "17'#D3D3D3',18'yellow'"
    colors(18) = "18'#D3D3D3',19'yellow'"
    colors(19) = "19'#D3D3D3',20'yellow'"
  end function colors



  function int2char(Me)
    implicit none
    integer,intent(in)::Me
    character*11::int2char
    character*3::tn
    integer::i1,i2,i3

    i1=Me/100
    i2=(Me-100*i1)/10
    i3=Me-100*i1-10*i2

    tn=char(i1+48)//char(i2+48)//char(i3+48)
    int2char='proc'//tn//'.dat'

  end function int2char

  subroutine write_data(voisins,ind_local,nb_probleme,U,Mx,My,dx,dy,filename)
    implicit none
    integer,intent(in)::Mx,My,nb_probleme
    real*8,dimension(Mx*My),intent(in)::U
    integer,dimension(4),intent(in)::voisins,ind_local
    real*8,intent(in)::dx,dy
    character(len=*),intent(in)::filename

    integer::i1,iN,j1,jN
    integer::i,j

    i1 = ind_local(1) ; iN = ind_local(2) ; j1 = ind_local(3) ; jN = ind_local(4)
    open(unit=12,file=filename)


    ! Bord bas
    if (voisins(4) < 0) then
       do i = i1,iN
          write(12,*) i*dx, 0.0, fonction_g(i*dx,0.0,nb_probleme)
       end do
       write(12,*)''
    end if

    ! Bord droit/gauche
    if (voisins(1) < 0 .and. voisins(3) < 0) then
       do j = j1,jN
          write(12,*) 0.0, j*dy, fonction_h(0.0,j*dy,nb_probleme)
          do i = i1,iN
             write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
          end do
          write(12,*) 1.0, j*dy, fonction_h(1.0,j*dy,nb_probleme)
          write(12,*)''
       end do

       ! Bord droit
    else if (voisins(1) < 0) then
       do j = j1,jN
          write(12,*) 0.0, j*dy, fonction_h(0.0,j*dy,nb_probleme)
          do i = i1,iN
             write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
          end do
          write(12,*)''
       end do

       ! Bord droit
    else if (voisins(3) < 0) then
       do j = j1,jN
          do i = i1,iN
             write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
          end do
          write(12,*) 1.0, j*dy, fonction_h(1.0,j*dy,nb_probleme)
          write(12,*)''
       end do

       ! Milieu du domaine
    else
       do j = j1,jN
          do i = i1,iN
             write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
          end do
          write(12,*)''
       end do
    end if

    ! Bord haut
    if (voisins(2) < 0) then
       do i = i1,iN
          write(12,*) i*dx, 1.0, fonction_g(i*dx,1.0,nb_probleme)
       end do
       write(12,*)''
    end if

    close(12)

  end subroutine write_data



  subroutine script_gnu2(Np,choix)

    implicit none

    integer,intent(in)::Np,choix
    integer :: i,p,pal
    character*25,dimension(0:19) :: palette

    if (choix /= 2) then
       palette = colors()
    end if

    open(unit = 52, file = "gnuplot_additif.gnu", &
         form = 'formatted', &
         status = 'unknown', action = 'write')

    write(52,*) "set terminal wxt persist"
    write(52,*)

    if (choix /= 2) then
       write(52,*) "set style fill transparent solid 0.75"
       write(52,*) "set palette defined (\"
       do i = 0,Np-2
          write(52, "(A23,',\')") palette(i)
       end do
       write(52, "(A23,' ) ')") palette(Np-1)
    else
       write(52,*) "set style fill transparent solid 0.55"
       write(52,*) "set palette defined (1'#00ff00',2'#00ff00')"
    end if
    write(52,*)
    write(52,*) "splot 'proc000.dat' using 1:2:3:(1*$3) w pm3d title '0', \"


    do p = 1,Np-2
       pal = mod(p,20)
       if (p < 10) then
          write(52, "(' ''proc00',I1.1,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I1.1,''',\')") p,pal,p
       else
          write(52, "(' ''proc0',I2.2,'.dat'' using 1:2:3:(1*(',I2.2,'+$3)) w pm3d title ''',I2.2,''',\')") p,pal,p
       end if
    enddo

    pal = mod(p,20)
    if (Np-1 < 10) then
       write(52, "(' ''proc00',I1.1,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I1.1,''' ')") Np-1,pal,Np-1
    else
       write(52, "(' ''proc0',I2.2,'.dat'' using 1:2:3:(1*(',I2.2,'+$3)) w pm3d title ''',I2.2,''' ')") Np-1,pal,Np-1
    end if

    close(52)

  end subroutine script_gnu2



  subroutine script_gnuplot(rank,Np,choix,c1,COMM_WORLD)

    implicit none
    include "mpif.h"

    integer,intent(in)::rank,Np,choix,c1,COMM_WORLD
    integer :: i,p,statinfo
    integer,dimension(0:Np-1):: color_case
    character*25,dimension(0:19) :: palette
    integer, dimension(MPI_STATUS_SIZE) :: status

    if (rank /=0) then
       call MPI_SEND(c1,1,MPI_INT,0,rank,COMM_WORLD,statinfo)   


    else

       color_case(0) = c1
       do i=1,Np-1
          call MPI_RECV(color_case(i),1,MPI_INT,i,i,COMM_WORLD,status,statinfo)
       end do

       palette = colors()

       open(unit = 52, file = "gnuplot_multi.gnu", &
            form = 'formatted', &
            status = 'unknown', action = 'write')

       write(52,*) "set terminal wxt persist"
       write(52,*)


       write(52,*) "set style fill transparent solid 0.75"
       write(52,*) "set palette defined (\"
       write(52, "(A23,',\')") palette(0)
       write(52, "(A23,' ) ')") palette(2)

       write(52,*)
       write(52, "(' splot ''proc000.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''0'' ,\')") color_case(0)


       do p = 1,Np-2
          if (p < 10) then
             write(52, "(' ''proc00',I1.1,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I1.1,''',\')") p,color_case(p),p
          else
             write(52, "(' ''proc0',I2.2,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I2.2,''',\')") p,color_case(p),p
          end if
       enddo

       if (Np-1 < 10) then
          write(52, "(' ''proc00',I1.1,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I1.1,''' ')") Np-1,color_case(p),Np-1
       else
          write(52, "(' ''proc0',I2.2,'.dat'' using 1:2:3:(1*(',I1.1,'+$3)) w pm3d title ''',I2.2,''' ')") Np-1,color_case(p),Np-1
       end if

    end if
    close(52)

  end subroutine script_gnuplot



end module display
