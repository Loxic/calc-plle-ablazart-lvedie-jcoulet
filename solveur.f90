module solve

  use mathfunctions
  implicit none

contains

  !> Fonction Sparse_solve
  !! @param[in] solver Choix du solveur utilisé 
  !! @param[in] Mx,My,dx,dy Dimensions et espacements du sous domaine
  !! @param[in] D Paramètre de diffusion
  !! @param[in] B Second membre du système AX=B
  !! @param[out] X Solution du système AX=B
  !! @param[in] tol,itermax Critères d'arrêt des méthodes itératives

  !> @details Permet de choisir le solveur utilisé dans la résolution du système
  !creux en spécifiant les critères de convergence. 
  subroutine Sparse_solve(solver,Mx,My,dx,dy,D,Dt,X,B,tol,itermax)

    implicit none

    integer,intent(in)::solver
    integer, intent(in) :: Mx,My  ! dimensions spatiales du problème
    real*8, intent(in)::dx,dy,D,Dt,tol
    integer,intent(in)::itermax
    real*8, dimension(:), intent(inout) :: X ! donnee
    real*8, dimension(:), intent(in) :: B ! donnee

    select case (solver)
    case (1)
      call Sparse_Jacobi(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)
    case (2)
      call Sparse_GaussSeidel(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)
    case (3)
      call Sparse_GradConj(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)
    case default
      print*,'Fatal error : unknow solver'
      stop
    end select

  end subroutine Sparse_solve

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


  subroutine Sparse_Jacobi(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)

    implicit none

    integer, intent(in) :: Mx,My  ! dimensions spatiales du problème
    real*8, intent(in)::dx,dy,D,Dt,tol
    integer,intent(in)::itermax
    real*8, dimension(:), intent(inout) :: X ! donnee
    real*8, dimension(:), intent(in) :: B ! donnee

    real*8,dimension(size(X))::X_next,AX
    integer :: i, j, k,l

    AX = 500 !Big value to enter in the loop
    l = 0

    do while (l<itermax .and. norme_2(AX - B) > tol*norme_2(B))  !Condition convergence
      call matmul_implicit(Mx,My,dx,dy,D,Dt,X,AX)
      !Computes product A*X without term aii*xi
      X_next = 0
      do j = 1, My
        do i = 1, Mx
          k = Mx*(j-1) + i
          ! bloc M
          !AX(k) = (1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2))*X(k) !Terme diagonal
          if (i>1) then
            X_next(k) = X_next(k) - (D*Dt/dx**2)*X(k-1)   !Terme sur diagonale D
          end if
          if (i<Mx) then
            X_next(k) = X_next(k) - (D*Dt/dx**2)*X(k+1)  !Terme sous diagonal de D
          end if

          !Bloc E inférieur
          if(j>1) then
            X_next(k)=X_next(k)-(Dt*D/dy**2)*X(k-Mx)
          end if
          !Bloc E supérieur
          if (j<My) then
            X_next(k)=X_next(k)-(Dt*D/dy**2)*X(k+Mx)
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



  subroutine Sparse_GaussSeidel(Mx,My,dx,dy,D,Dt,X,B,tol,itermax)

    implicit none

    integer, intent(in) :: Mx,My  ! dimensions spatiales du problème
    real*8, intent(in)::dx,dy,D,Dt,tol
    integer,intent(in)::itermax
    real*8, dimension(:), intent(inout) :: X ! donnee
    real*8, dimension(:), intent(in) :: B ! donnee

    real*8,dimension(size(X))::X_next,AX
    integer :: i, j, k,l

    AX = 500 !Big value to enter in the loop
    l = 0

    do while (l<itermax .and. norme_2(AX - B) > tol*norme_2(B))  !Condition convergence
      call matmul_implicit(Mx,My,dx,dy,D,Dt,X,AX)
      !Computes product A*X without term aii*xi
      X_next = 0
      do j = 1, My
        do i = 1, Mx
          k = Mx*(j-1) + i
          ! bloc M
          !AX(k) = (1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2))*X(k) !Terme diagonal
          if (i>1) then
            X_next(k) = X_next(k) - (D*Dt/dx**2)*X_next(k-1)   !Terme sur diagonale D
          end if
          if (i<Mx) then
            X_next(k) = X_next(k) - (D*Dt/dx**2)*X(k+1)  !Terme sous diagonal de D
          end if

          !Bloc E inférieur
          if(j>1) then
            X_next(k)=X_next(k)-(Dt*D/dy**2)*X_next(k-Mx)
          end if
          !Bloc E supérieur
          if (j<My) then
            X_next(k)=X_next(k)-(Dt*D/dy**2)*X(k+Mx)
          end if

          !Xi = 1/aii (bi - Aij*xj )
          X_next(k) = 1./((1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2)))*(-1*X_next(k)+B(k)) !Terme diagonal
          X(k) = X_next(k)
        end do
      end do

      l=l+1
    end do
  end subroutine Sparse_GaussSeidel

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

  subroutine Sparse_Gradconj(Nx,Ny,dx,dy,D0,Dt,K,B,eps,Nmax)
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
    integer::l

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
  end subroutine Sparse_Gradconj

  subroutine matmul_implicit(Nx,Ny,dx,dy,D,Dt,X,AX)

    implicit none

    integer, intent(in) :: Nx,Ny  ! dimensions spatiales du problème
    real*8, intent(in)::dx,dy,D,Dt
    real*8, dimension(:), intent(in) :: X ! donnee
    real*8, dimension(:), intent(out) :: AX ! Sortie

    integer :: i, j, k

    AX=0
    do j = 1, Ny
       do i = 1, Nx
          k = Nx*(j-1) + i
          ! bloc M
          AX(k) = (1+2*D*Dt*(1.0/dx**2 + 1.0/dy**2))*X(k)
          if (i>1) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k-1)
          end if
          if (i<Nx) then
             AX(k) = AX(k) - (D*Dt/dx**2)*X(k+1)
          end if

          !Bloc E inférieur
          if(j>1) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k-Nx)
          end if
          !Bloc E supérieur
          if (j<Ny) then
             AX(k)=AX(k)-(Dt*D/dy**2)*X(k+Nx)
          end if

       end do
    end do

  end subroutine matmul_implicit


end module
