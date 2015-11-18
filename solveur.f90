module solve

  use patate
  implicit none

contains

  function norme_2(V)
    implicit none
    real*8,dimension(:),intent(in)::V
    real*8::norme_2
    norme_2=sqrt(dot_product(V,V))
  end function norme_2

  subroutine GaussSeidel(A,b,x,tol)
    !!Méthode itérative : X doit être donné avec une valeur initiale !!
    implicit none
    ! E/S
    real*8,dimension(:,:),intent(in)::A
    real*8,dimension(:),intent(in)::B
    real*8,dimension(:),intent(inout)::X
    ! Locales
    real*8,dimension(size(X))::X_next
    real,intent(in)::tol
    integer::i,j

    do while (norme_2(matmul(A,X)-B) > tol*norme_2(B))  !Condition convergence
      do i=1,size(X)
        X_next(i) = B(i)
        do j=1,i-1
          X_next(i) = X_next(i) - A(i,j) * X_next(j)
        end do
        do j=i+1,size(X)
          X_next(i) = X_next(i) - A(i,j)*X(j)
        end do
        X_next(i) = X_next(i) / A(i,i)
      end do
      X=X_next
    end do
  end subroutine GaussSeidel


  subroutine Grad_conj_implicit(K,B,eps,Nmax,Nx,Ny,dx,dy,D0,Dt)
    !! Gradient implicite utilisé l'an dernier. La matrice n'est pas stockées !!
    implicit none
    !Entrées/sorties
    real*8,dimension(:),intent(in)::B
    real*8,dimension(:),intent(inout)::K
    real*8,intent(in)::eps
    integer,intent(in)::Nmax
    !Pour la multiplication
    integer,intent(in)::Nx,Ny
    real*8,intent(in)::dx,dy,D0,Dt
    !Locales
    real*8,dimension(size(K))::R,R_next,D,W
    real*8::alpha,beta
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
