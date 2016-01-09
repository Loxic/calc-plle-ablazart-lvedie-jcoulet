module display
  implicit none

  contains

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

  subroutine write_data(Me,voisins,ind_local,nb_probleme,U,Mx,My,dx,dy,filename)
    implicit none
    integer,intent(in)::Me,Mx,My,nb_probleme
    real*8,dimension(Mx*My),intent(in)::U
    integer,dimension(4),intent(in)::voisins,ind_local
    real*8,intent(in)::dx,dy
    character(len=*),intent(in)::filename

    integer::i1,iN,j1,jN
    integer::i,j

    i1 = ind_local(1) ; iN = ind_local(2) ; j1 = ind_local(3) ; jN = ind_local(4)
    open(unit=12,file=filename)


  !  if (voisins(4)==0) then ! bord y=0
  !    do i = i1,iN
  !      write(12,*) i*dx, j*dy, fonction_g(i*dx,0.0,nb_probleme)
  !    end do
  !    write(12,*)''
  !  end if


  !  if (voisins(1)==0 .and. voisins(3)==0) then ! bord x=0 et x=1
  !    do j = j1,jN
  !      write(12,*) 0.0, j*dy, fonction_h(0.0,j*dy,,nb_probleme)
  !      do i = i1,iN
  !        write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
  !      end do
  !      write(12,*) 1.0, j*dy, fonction_h(1.0,j*dy,,nb_probleme)
  !      write(12,*)''
  !    end do

  !  else if (voisins(1)==0) then ! bord x=0
  !    do j = j1,jN
  !      write(12,*) 0.0, j*dy, fonction_h(0.0,j*dy,,nb_probleme)
  !      do i = i1,iN
  !        write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
  !      end do
  !      write(12,*)''
  !    end do

  !  else if (voisins(3)==0) then ! bord x=1
  !    do j = j1,jN
  !      write(12,*) 0.0, j*dy, fonction_h(0.0,j*dy,,nb_probleme)
  !      do i = i1,iN
  !        write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
  !      end do
  !      write(12,*) 1.0, j*dy, fonction_h(1.0,j*dy,,nb_probleme)
  !      write(12,*)''
  !    end do

  !  else if (voisins(3)==0) then ! Millieu du domaine

      do j = j1,jN
        do i = i1,iN
          write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
        end do
        write(12,*)''
      end do
  !  end if


  !  if (voisins(2)==0) then ! bord y=0
  !    do i = i1,iN
  !      write(12,*) i*dx, j*dy, fonction_g(i*dx,1.0,nb_probleme)
  !    end do
  !    write(12,*)''
  !  end if


    close(12)
    
  end subroutine write_data

end module display
