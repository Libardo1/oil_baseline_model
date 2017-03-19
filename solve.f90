module solve

  use basic
  use def_types
  use tools
  implicit none
  
contains
  
  !***********************************************************************
  ! SOLVEDP
  !***********************************************************************
  subroutine solvedp(v,x,pstar, f,P,beta,alg,v_ini)
  !***********************************************************************
  ! SOLVEDP solves a dynamic programming problem with optimal stopping
  ! Usage: 
  !	call solvedp(v,x,pstar, f,n,m,P,beta,alg,v_ini)
    
  ! INPUTS
  !     f          : reward function (n X m)
  !     n          : # of states 
  !     m          : # of actions
  !     beta       : discount factor
  !     alg        : algorithm to use (1:policy & 2:value)
  !     v_ini      : initial value for value function 
    
  ! OUTPUTS
  !     v          : value function   
  !     x          : policy function
  !     pstar      : optimal transition matrix
  !************************************************************************

    real :: start, finish 
    
    !Local 
    integer :: info
    integer, allocatable ::ipiv(:)
    integer, parameter :: maxit = 5000
    real, parameter :: tol = 1e-4
    logical, parameter :: prtiters = .true.

    integer :: i, j, it, n, m
    real, allocatable :: vold(:), fstar(:)
    real :: change
    real, allocatable :: a(:,:)
    
    !Dummy 
    real, intent(out), allocatable :: v(:), x(:)
    real, intent(out), allocatable :: pstar(:,:)

    character(*)	:: alg
    real, intent(in) :: f(:,:), beta, v_ini

    type(sparse_matrix), intent(in) :: P

	!************ Allocations
	
	n = size(f,1)
	m = size(f,2)
	allocate(ipiv(n))
	allocate(vold(n))
	allocate(fstar(n))
	allocate(a(n,n))
	allocate(v(n))
	allocate(x(n))
	allocate(pstar(n,n))
    !************ Perform policy or value function iterations
    
    v = reshape1(v_ini, n)   
        
    select case (alg)
    case('policy')
       if (prtiters)then
          write(*,*) 'Solving Bellman by policy function iteration'
       end if
       do it = 1,maxit
	  	vold  = v
	  	call valmax(v,x, vold,f,P,beta)
	  	call valpol(pstar,fstar, x,f,P,beta)
	  	!**************************************
	  	!* Solving with LAPACK
	  	!* (eye(n)-beta*pstar)*v=fstar	 
 	  	a = (eye(n)-beta*pstar) 
 	  	call sgetrf ( n, n, a, n, ipiv, info )
 	  	call sgetrs ( 'n', n, 1, a, n, ipiv, fstar, n, info )
 	  	v = fstar
	  	!**************************************
 	  	change = norm2(v-vold)
	  	if (prtiters)then
              write(*,*) it, ' ' , change 
	  	end if
	  	if (change<tol)then
             exit
	  	end if
       end do
    case('value') 
       if (prtiters)then
          write(*,*) 'Solving Bellman by value function iteration'
       end if
       do it = 1,maxit
          vold = v
          call valmax(v,x, vold,f,P,beta)
          change = norm2(v-vold)
          if (prtiters)then
             write(*,*) it, ' ' , change 
          end if
          if (change<tol)then
             exit
          end if
       end do
       if (change>tol)then
          write(*,*) 'Failure to converge in solvedp'
       end if
    end select
  end subroutine solvedp
  
  !***********************************************************************
  ! VALMAX
  !***********************************************************************
  subroutine valmax(v,x, vold,f,P,beta)
        
    !Local
    real, allocatable ::  h(:,:)
    integer ::  n, m

    !Dummy
    real, intent(out), allocatable :: v(:), x(:)
    
    real, intent(in) :: beta, f(:,:), vold(:)
    
    type(sparse_matrix) :: P
    
    n = size(f,1)
    m = size(f,2)
    allocate(v(n))
    allocate(x(n))
    allocate(h(n,m))
    h = reshape(sparse_matmulvec(P,size(P%values),vold),(/n,m/))
    v = maxval(f + beta * reshape(sparse_matmulvec(P,size(P%values),vold), (/n, m/)),2)
    x = maxloc(f + beta * reshape(sparse_matmulvec(P,size(P%values),vold), (/n, m/)),2)
  end subroutine valmax
  
  !***********************************************************************
  ! VALPOL
  !***********************************************************************
  subroutine valpol(pstar,fstar, x,f,P,beta)

     implicit none
    
    !Local
    integer :: j, n, m
    real, allocatable :: a(:,:), ff(:,:)
    
    !Dummy
    real, intent(out), allocatable :: pstar(:,:),fstar(:)

    type(sparse_matrix),intent(in) :: P

    real, intent(in) :: x(:), f(:,:), beta
        
    n = size(f,1)
    m = size(f,2)
    
    allocate(a(n*m,n))
    allocate(ff(n*m,1))
    allocate(pstar(n,n))
    allocate(fstar(n))

    
    a = sparsetomatrix(P,size(P%d1),n*m,n)
    ff = reshape(f,(/n*m,1/))
    do j = 1,n	  
      fstar(j) = ff(n*(int(x(j))-1)+j,1)
      pstar(j,:) = a(n*(int(x(j))-1)+j,:)
    end do
  end subroutine valpol
end module solve