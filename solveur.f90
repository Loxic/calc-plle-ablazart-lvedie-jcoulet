module solve

  use patate
  implicit none

contains

  subroutine Gauss_seidel()
    implicit none
    !Insert code here
    
    end subroutine Gauss_seidel


  subroutine Grad_conj_implicit(K,B,eps,Nmax,Nx,Ny,dx,dy,D0,Dt)
    implicit none
    !Entrées/sorties
    real(wp),dimension(:),intent(in)::B
    real(wp),dimension(:),intent(inout)::K
    real(wp),intent(in)::eps
    integer,intent(in)::Nmax
    !Pour la multiplication
    integer,intent(in)::Nx,Ny
    real(wp),intent(in)::dx,dy,D0,Dt
    !Locales
    real(wp),dimension(size(K))::R,R_next,D,W
    real(wp)::alpha,beta
    integer::l,i

    !R=matmul(A,K)-B
    call matmul_implicit(Nx,Ny,dx,dy,D0,Dt,K,R)
    R=R-B
    D=R

    l=0
    do while ( l<Nmax .and. (dot_product(R,R))>eps**2 )

      call matmul_implicit(Nx,Ny,dx,dy,D0,Dt,D,W)
      Alpha = dot_product(D,R) / dot_product(D,W)
      K=K-alpha*D
      R_next=R-alpha*W
      beta=dot_product(R_next,R_next) / dot_product(R,R)
      D=R_next+beta*D
      l=l+1
      R=R_next
    end do
    print*,'Convergence du GC en',l
  end subroutine Grad_conj_implicit

end module
