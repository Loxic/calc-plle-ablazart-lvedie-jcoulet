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

  subroutine write_data(Me,ind_local,U,Mx,My,dx,dy,filename)
    implicit none
    integer,intent(in)::Me,Mx,My
    real*8,dimension(Mx*My),intent(in)::U
    integer,dimension(4),intent(in)::ind_local
    real*8,intent(in)::dx,dy
    character(len=*),intent(in)::filename

    integer::i1,iN,j1,jN
    integer::i,j

    i1 = ind_local(1) ; iN = ind_local(2) ; j1 = ind_local(3) ; jN = ind_local(4)
    open(unit=12,file=filename)
    do i = i1,iN
      do j = j1,jN
        write(12,*) i*dx, j*dy, U(i-i1 +1 + (j-j1)*Mx)
      end do
      write(12,*)''
    end do
    close(12)
    
  end subroutine write_data

end module display
