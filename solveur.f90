module solve

  use patate
  implicit none

contains

!> Fonction norme 
!! @param [in] V Vecteur d'entrée

!> @details Fonction utilisée dans les algoirthmes de résolution de Gauss Seidel et Jacobi


  function norme_2(V)
    implicit none
    real*8,dimension(:),intent(in)::V
    real*8::norme_2
    norme_2=sqrt(dot_product(V,V))
  end function norme_2


!> Résolution du système \f$ A.X = b \f$ par l'algorithme de Jacobi
!! @param [inout] X Vecteur second membre de taille Nx.Ny
!! @param [in] A Matrice
!! @param [in] b Vecteur second membre
!! @param [in] tol critère de convergence de l'algorithme

!> @warning X doit être initialisé
!> @warning A doit être a diagonale strictement dominante

!> @todo Adapater l'algo au système 

  subroutine Jacobi(A,b,x,tol)
    !!Méthode itérative : X doit être donné avec une valeur initiale !!
    !!Attention, convergence assurée uniquement dans les cas suivants : 
    !!    * A est a diagonale strictement dominante                  !!
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
          X_next(i) = X_next(i) - A(i,j) * X(j)
        end do
        do j=i+1,size(X)
          X_next(i) = X_next(i) - A(i,j)*X(j)
        end do
        X_next(i) = X_next(i) / A(i,i)
      end do
      X=X_next
    end do
  end subroutine Jacobi

  
  subroutine Sparse_Jacobi(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)

    implicit none

    integer, intent(in) :: Mx,My  ! dimensions spatiales du problème
    real(wp), intent(in)::dx,dy,D,Dt,tol
    integer,intent(in)::itermax
    real(wp), dimension(:), intent(inout) :: X ! donnee
    real(wp), dimension(:), intent(in) :: B ! donnee

    real*8,dimension(size(X))::X_next,AX
    integer :: i, j, k,l

    AX = 500 !Big value to enter in the loop
    l = 0

    do while (l<itermax .and. norme_2(AX - B) > tol*norme_2(B))  !Condition convergence
    call matmul_implicit(Mx,My,dx,dy,D,Dt,X,AX)
    !Computes product A*X without term aii*xi
    X_next = 0
    do i = 1, Mx
       do j = 1, My
          k = Mx*(j-1) + i
          ! bloc M
          !AX(k) = (1+2*D*Dt*(1.0_wp/dx**2 + 1.0_wp/dy**2))*X(k) !Terme diagonal
          if (i>1) then
             X_next(k) = X_next(k) - (D*Dt/dx**2)*X(k-1)   !Terme sur diagonale D
          end if
          if (i<Mx) then
             X_next(k) = X_next(k) - (D*Dt/dx**2)*X(k+1)  !Terme sous diagonal de D
          end if

          !Bloc E inférieur
          if(j>1) then
             X_next(k)=X_next(k)-(Dt*D/dy**2)*X(k-My)
          end if
          !Bloc E supérieur
          if (j<My) then
             X_next(k)=X_next(k)-(Dt*D/dy**2)*X(k+My)
          end if

       end do
    end do

  !Xi = 1/aii (bi - Aij*xj )
    X_next(:) = 1./((1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2)))*(-1*X_next(:)+B(:)) !Terme diagonal
    X = X_next
    l=l+1
  end do
  end subroutine Sparse_Jacobi


!> Résolution du système \f$ A.X = b \f$ par l'algorithme de Gauss-Seidel
!! @param [inout] X Vecteur second membre de taille Nx.Ny
!! @param [in] A Matrice
!! @param [in] b Vecteur second membre
!! @param [in] tol critère de convergence de l'algorithme

!> @warning X doit être initialisé
!> @warning A doit être SDP ou a diagonale strictement dominante

!> @todo Adapater l'algo au système 

  subroutine GaussSeidel(A,b,x,tol)
    !!Méthode itérative : X doit être donné avec une valeur initiale !!
    !!Attention, convergence assurée uniquement dans les cas suivants : 
    !!    * A est SDP                                                !!
    !!    * A est a diagonale strictement dominante                  !!
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


!> Résolution du système \f$ A.X = b \f$ par l'algorithme du Gradient Conjugué
!! @param [inout] K Vecteur second membre de taille Nx.Ny
!! @param [in] B Vecteur second membre
!! @param [in] eps critère de convergence de l'algorithme
!! @param [in] Nmax Nombre d'itération max.
!! @param [in] Nx,Ny Nombre de noeuds en x, y
!! @param [in] dx,dy pas d'espace en x, y
!! @param [in] Dt pas de temps
!! @param [in] D0 Coefficient de diffusion

!> @warning A doit être SDP

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
