module tools

  use basic
  use def_types
  implicit none

contains

  !***********************************************************************
  ! TAUCHENHUSSEY
  !***********************************************************************
  subroutine tauchenhussey(z, prob, N, mu, rho, sigma)
    
    !Local
    integer :: i
    real :: z_m(N,1)
    
    !Dummy 
    integer, intent(in) ::  N
    real, intent(in) :: mu, rho, sigma
    
    real, intent(out) :: z(N), prob(N, N)
	
    prob = transpose(loaddata(1, 'inputData/prob_tauchenhussey.dat', N, N))

    z_m = loaddata(2, 'inputData/z_tauchenhussey.dat', N, 1)
    z = z_m(:,1)
    
  end subroutine tauchenhussey
  
  !***********************************************************************
  ! ROUWENHORST
  !***********************************************************************
  subroutine rouwenhorst(Z, PI, NN, mu, rho, sigma)
  !**********************************************************
  ! ROUWENHORST creates a grid and a transition matrix for a
  !             AR(1) process  
  ! Usage: 
  !	call rouwenhorst(z,prob,N, mu,rho,sigma)
    
  ! INPUTS
  !     N      	: size of the grid (# number of states)
  !     mu     	: mean of the AR(1) process
  !     rho    	: persistence of the process
  !     sigma  	: standart desviation of the process
    
  ! OUTPUTS	
  !		Z		: vector grid  
  !		PI   	: matrix transition  
  !***********************************************************
    
   !Local
    integer :: n
	real, allocatable :: PIp(:,:), PIm(:,:)
	real :: p, sigmaz, fi
	
	!Dummy 
	integer, intent(in) ::  NN
	real, intent(in) :: mu, rho, sigma 
	
	real, intent(out), allocatable :: Z(:), PI(:,:)
	
	sigmaz = sigma/ sqrt(1-rho**2)
	
	p = (1+rho)/2
	allocate(PI(2,2))
	PI(1,:) =(/p,1-p/)
	PI(2,:) =(/1-p,p/)
	
	do n = 3,NN
		allocate(PIp(n,n))
		allocate(PIm(n,n))
		PIm = 0. 
		PIm(1:n-1,1:n-1) = p*PI
		PIp = PIm
		
		PIm = 0. 
		PIm(:n-1,2:) = (1-p)*PI
		PIp = PIp + PIm
		
		PIm = 0. 
		PIm(2:,:n-1) = (1-p)*PI
		PIp = PIp + PIm	
		
		PIm = 0. 
		PIm(2:,2:) = p*PI
		PIp = PIp + PIm
		
		PIp(2:n-1,:) = PIp(2:n-1,:)/2
				
		deallocate(PI)
		allocate(PI(n,n))
		PI = PIp
		
		deallocate(PIp)
		deallocate(PIm)
	enddo
	
	fi = sqrt(real(NN-1))*sigmaz
	allocate(Z(NN))
	Z = linspace(-fi,fi,NN)
	Z = Z + mu
  end subroutine rouwenhorst
  
  !***********************************************************************
  ! GRIDMAKE
  !***********************************************************************
  !******** GRIDMAKE2
  function gridmake2(a, b) result(grid)
  !**********************************************************
  ! GRIDMAKE2 creates 2-dimensional state-space
  ! Usage: 
  !	grid = gridmake2(a,a_n,b,b_n)
    
  ! INPUTS
  !     a      : vector variable
  !     b      : vector variable
    
  ! OUTPUTS	
  !		grid   : grid with the 2-dimensional state-space
  !***********************************************************
    !Local
    integer :: i, j

    !Dummy
	integer :: a_n, b_n 
    real :: a(:), b(:)
    real, allocatable :: grid(:,:)
    a_n = size(a)
    b_n = size(b)
   
    allocate(grid(a_n*b_n,2))
    grid =0.
    
	grid(:,1) = reshape(spread(a,2,b_n),(/size(grid,1)/))
	grid(:,2) = reshape(spread(b,1,a_n),(/size(grid,1)/))
	
  end function gridmake2

  !******** GRIDMAKE3
  function gridmake3(a, a_n, b, b_n, c, c_n) result(grid)
  !**********************************************************
  ! GRIDMAKE3 creates 3-dimensional state-space
  ! Usage: 
  !	grid = gridmake2(a,a_n,b,b_n,c,c_n)
    
  ! INPUTS
  !     a      : vector variable
  !     a_n    : size of vector a
  !     b      : vector variable
  !     b_n    : size of vector b
  !     c      : vector variable
  !     c_n    : size of vector c
    
  ! OUTPUTS	
  !		grid   : grid with the 3-dimensional state-space
  !***********************************************************
    !Local
    integer :: h = 1, i, j, k 

    !Dummy
    integer, intent(in) :: a_n, b_n, c_n 
    real, intent(in) :: a(a_n), b(b_n), c(c_n)

    real :: grid(a_n*b_n*c_n,3)
    do k = 1,c_n 
       do j = 1,b_n
          do i = 1,a_n
             grid(h,:) = (/a(i),b(j),c(k)/)
             h = h + 1
          end do
       end do
    end do
  end function gridmake3
  
  !******** GRIDMAKE4
  function gridmake4(a, a_n, b, b_n, c, c_n, d, d_n) result(grid)
  !**********************************************************
  ! GRIDMAKE4 creates 4-dimensional state-space
  ! Usage: 
  !	grid = gridmake2(a,a_n,b,b_n,c,c_n,d,d_n)
    
  ! INPUTS
  !     a      : vector variable
  !     a_n    : size of vector a
  !     b      : vector variable
  !     b_n    : size of vector b
  !     c      : vector variable
  !     c_n    : size of vector c
  !     d      : vector variable
  !     d_n    : size of vector d
    
  ! OUTPUTS	
  !		grid   : grid with the 4-dimensional state-space
  !***********************************************************
    
    !Local
    integer :: h = 1, i, j, k, l
    
    !Dummy
    integer :: a_n, b_n, c_n, d_n 
    real :: a(a_n), b(b_n), c(c_n), d(d_n)
    
    real :: grid(a_n*b_n*c_n*d_n,4)
    
    do i = 1,a_n
       do j = 1,b_n
          do k = 1,c_n
             do l = 1,d_n
                grid(h,:) = (/a(i),b(j),c(k),d(l)/)
                h = h + 1
             end do
          end do
       end do
    end do
  end function gridmake4
     
  !**********************************************************
  ! KRONECKER SPARSE
  !**********************************************************
  function kroneckerS(b_n,c) result(a)
  !**********************************************************
  ! KRONECKERS makes a sparse kronecker product from identity 
  ! 		   matrix and a dense matrix
  ! Usage: 
  !	a = kroneckerS(b, c) 
    
  ! INPUTS
  !     b_n     : size of square-matrix b
  !     c       : real matrix
    
  ! OUTPUTS	
  !		a       : sparse matrix with the result
  !**********************************************************
    
    ! Local 
    integer :: i, j, k, l, c_n_row, c_n_col
    integer :: element = 1
    
    ! Dummy 
    integer, intent(in) :: b_n
    real, intent(in) :: c(:,:)
    
    type(sparse_matrix) :: a

	c_n_row = size(c,1)
	c_n_col = size(c,2)
    ! Allocate
    allocate(a%values(b_n*c_n_row*c_n_col))
    allocate(a%d1(b_n*c_n_row*c_n_col))
    allocate(a%d2(b_n*c_n_row*c_n_col))
    
    do i = 1,b_n  
       do k = 1,c_n_row
          do l = 1,c_n_col
             a%values(element) = c(k,l)
             a%d1(element)     = k + (i-1)*c_n_row
             a%d2(element)     = l + (i-1)*c_n_col
             element = element + 1
          end do
       end do
    end do
  end function kroneckerS
  
  !**********************************************************
  ! SPARSE
  !**********************************************************
  function sparse(matrix) result(a) 
  !**********************************************************
  ! SPARSE converts a full matrix to sparse
  ! Usage: 
  !	result = sparse(matrix)
    
  ! INPUTS
  !     matrix 	: real full matrix
    
  ! OUTPUTS	
  !		result	: sparse matrix
  !**********************************************************
    
    ! Local 
    integer :: i, j, n, m, nonzero
    integer :: element = 1
    
    ! Dummy 
    real :: matrix(:,:)   
    logical, allocatable :: mask(:,:)
    type(sparse_matrix) :: a

!	write(*,*) 'nonzero'


    ! Allocate
    n = size(matrix,1)
    m = size(matrix,1)
    allocate(mask(n,m))
    mask = .false.
	where(matrix/=0.) mask=.true.
!	write(*,*) 'mask'

    nonzero = count(mask)
    
    allocate(a%values(nonzero))
    allocate(a%d1(nonzero))
    allocate(a%d2(nonzero))
    
!    write(*,*) 'nonzero', nonzero
    do i = 1,n  
    	do j = 1,m
    		if(matrix(i,j)/=0.)then 
            	a%values(element) = matrix(i,j)
            	a%d1(element)     = i
            	a%d2(element)     = j
            	element = element + 1
			end if
    	end do
	end do
  end function sparse
  
  !**********************************************************
  ! KRONECKER 
  !**********************************************************
  function kronecker(b,c) result(a)
  !**********************************************************
  ! KRONECKER  makes kronecker product
  ! Usage: 
  !	a = kronecker(b,c)
    
  ! INPUTS
  ! 	b	: real matrix
  ! 	c   : real matrix
    
  ! OUTPUTS	
  !		a	: real matrix with the result
  !**********************************************************       
    !Local 
    integer :: i, j, k, l
    
    !Dummy 
    integer :: b_dim(2), c_dim(2)
    real   :: b(:,:), c(:,:)
    
    real, allocatable :: a(:,:)
	b_dim = shape(b)
	c_dim = shape(c)
	
    allocate(a(b_dim(1)*c_dim(1),b_dim(2)*c_dim(2)))
    
    do i = 1,b_dim(1)
    	do j = 1,b_dim(2)
			do k = 1,c_dim(1)
				do l = 1,c_dim(2)
                   a(k + (i-1)*c_dim(1), l + (j-1)*c_dim(2)) = b(i,j)*c(k,l)
                end do
            end do
       	end do
    end do
    
  end function kronecker 

  !**********************************************************
  ! SPARSE MATMUL
  !**********************************************************
  function sparse_matmulvec(b,b_n, c) result(a)
  !**********************************************************
  ! SPARSE_MATMUL makes sparse matrix vector multiplication
  ! Usage: 
  !	a = sparse_matmul(b,b_n, c)
    
  ! INPUTS
  !     b    : real sparse matrix
  !     b_n  : # of non-zero elements of sparse matrix b
  !     c    : real vector
    
  ! OUTPUTS	
  !		a    : real sparse matrix with the result
  !********************************************************** 
    ! Local
    integer :: k
    integer :: c_n
    
    ! Dummy
    type(sparse_matrix), intent(in) :: b
    integer, intent(in) :: b_n

    real, intent(in) :: c(:)

    real, allocatable :: a(:)
	
	c_n = size(c)
	
    allocate(a(b%d1(b_n)))
    a = 0. 
    do k = 1,b_n
       a(b%d1(k)) = a(b%d1(k)) + (b%values(k)*c(b%d2(k)))
    end do
    
  end function sparse_matmulvec
  
  function dot(b,c) result(a)
 	integer :: n, i
	real :: b(:), c(:)
	real :: a 
	
	if(size(c)/=size(b))then
		write(0,*) 'Error :  different shape' 
		stop
	end if
	n = size(b)
	a = 0.
	do i = 1,n
		a= a + b(i)*c(i)
	enddo 
  end function 
      
  !**********************************************************
  ! ERGODIC DISTRIBUTION
  !**********************************************************
  function ergdist(P,d) result(b)
  !**********************************************************
  ! ERGDIST finds the ergodic distribution
  ! Usage: 
  !	dist = ergodist(P,d)
    
  ! INPUTS
  ! 	P    : 
  !     d	 : 
    
  ! OUTPUTS	
  !		b    : 
  !********************************************************** 
	real :: P(:,:), d(:,:)
  	real, allocatable :: dnext(:,:), b(:,:), din(:,:)
  	integer :: maxiter, i
  	real :: tol, change
	  	
!	if (rank(d)/=2)then
!		write(0,*) 'd must be 2 dimensional array'
!		stop
!	endif

	allocate(din(1,size(d)))
	
	if(size(d,1)/=1)then
		din = transpose(d)
	else 
		din = d
	end if
	
  	allocate( dnext(1,size(d,2)))
  	allocate(     b(1,size(d,2)))
  	
! 	write(*,*) i, maxiter, tol, change
	change =100.
	maxiter = 500
	i=1
	tol = 1.e-6
	
	do while( i<maxiter .and. tol<change)
  		dnext 	= matmul(din,P)	
  		change 	= maxval(abs(dnext(1,:)-din(1,:)),1)
  		din		= dnext
!  		WRITE(*,*) i, 'k'
  		i = i+1
  	end do 
  	
  	b = reshape(din,shape(d))
!	write(*,*) b

  end function

  !********************************************************
  ! APROXIMATION MARKOV PROCESS
  !**********************************************************
  subroutine appmarkov2(prob,g,pi, cg,lambda,sigma,width,n1 ,n)
  !**********************************************************
  ! APPMARKOV2 creates grid, ergodic distribution and transition
  !            matrix for a AR(1) process  
  ! Usage: 
  !	call appmarkov2(prob,g,pi, cg,lambdaz,prob,N, mu,rho,sigma)
    
  ! INPUTS
  !     cg        : 
  !     lambda    : 
  !     sigma     :
  !     width     : 
  !     n1        : 
  !     n         :
    
  ! OUTPUTS	
  !		prob      : ergodic distribution
  !		g         : vector grid 
  !     pi        : transition matrix 
  !***********************************************************
    !Local
    real :: g_p(1,n1), pi_p(n1,1)
    
    !Dummy 
    integer, intent(in) :: n1,n
    real, intent(in) :: cg,lambda, sigma, width

    real, intent(out) :: prob(n1,n1), g(n1), pi(n1)     
    
    prob = loaddata(1, 'prob.dat', n1, n1)
    g_p = loaddata(2, 'g.dat', 1, n1)
    pi_p = loaddata(3, 'pi.dat', n1, 1)
    
    g = g_p(1,:)
    pi = pi_p(:,1)
    prob = transpose(prob)
    
  end subroutine appmarkov2

  !**********************************************************
  ! SPARSE TO GENERAL MATRIX
  !**********************************************************  
  function sparsetomatrix(b,b_n,n,m)  result(a)
  !**********************************************************
  ! SPARSETOMATRIX converts a sparse matrix in a general full matrix  
  ! Usage: 
  !	a = sparsetomatrix(b,b_n,n,m)
    
  ! INPUTS
  !     b         : sparse matrix 
  !     b_n       : size of non-zero elements of b
  !     n         : rows of the new matrix
  !     m         : columns of the new matrix

  ! OUTPUTS	
  !		a         : real general matrix
  !***********************************************************
    !Local
    integer :: i
    
    !Dummy
    type(sparse_matrix),intent(in) :: b
    integer, intent(in) :: b_n,n,m
    
    real :: a(n,m)
    a = 0.
    do i = 1,b_n
      a(b%d1(i),b%d2(i)) = b%values(i)
    
    end do
  end function sparsetomatrix
  
  !**********************************************************
  ! PIECEWISE LINEAR INTERPOLATION WITH EXTRAPOL
  !**********************************************************
  function piecewise(xgrid,ygrid,xp) result(yp)
  !**********************************************************
  ! PIECEWISE makes a linear interpolation for the function describe  
  ! 		  by xgrid & ygrid	 
  ! Usage: 
  !	yp = piecewise(xgrid,ygrid,xp)
    
  ! INPUTS
  !     xp        : point to interpolate
  !     xgrid     : grid for independent variable 
  !     ygrid     : grid for dependent variable

  ! OUTPUTS	
  !		yp        : interpolated point
  !***********************************************************
	integer :: xr
	real :: xgrid(:), ygrid(:),xp, val, yp
	
	if(size(xgrid)==size(ygrid))then
	if    (xp>xgrid(size(xgrid))) then
		xr = size(xgrid)
	elseif(xp<xgrid(1)) then 
		xr = 2
	else
		val = minval(pack(xgrid,(xp<xgrid)) )	! value 
		xr   = minloc(abs(xgrid-val),1)			! position 
	endif
	
	yp = ((xgrid(xr)-xp)*ygrid(xr-1) + (xp-xgrid(xr-1))*ygrid(xr))/(xgrid(xr)-xgrid(xr-1))
	else
		write(0,*) ' No same size in grids'
		stop
	endif
  end function piecewise	
  
  !**********************************************************
  ! MARKOV CHAIN
  !**********************************************************
  subroutine markov(chain,state, T,n,s0,V)
  !**********************************************************
  ! MARKOV generates a markovian process
  ! 		  
  ! Usage: 
  !	call markov(chain,state, T,n,s0,V)
    
  ! INPUTS
  !     T       : matrix transition
  !     n     	: # of periods to simulate
  !     s0      : initial point 
  !     V       : 

  ! OUTPUTS	
  !		chain	: nx1 chain 
  !		state   : 

  !***********************************************************
  	! Local
  	integer :: i, c, r, k, dim(2)
  	real,dimension(:), allocatable :: X
  	integer,dimension(:), allocatable :: s
  	real, allocatable :: cum(:,:), ppi(:), ones(:,:)
  	 
  	! Dummy 
  	integer,intent(in)	:: n, s0
  	real, intent(in)	:: T(:,:), V(:)
  	real, intent(out), allocatable	:: chain(:), state(:,:)
  	
  	dim = shape(T)
  	r=dim(1)
  	c=dim(2)
	allocate(ones(r,c))
	ones =1.
	allocate(chain(n-1))
  	allocate(state(r,n-1))
  	
  	allocate(X(n-1))
  	
  	do i=1,n-1
		call random_number(X(i))
    enddo
  	allocate(s(r))
  	s=0.
  	s(s0) = 1
  	allocate(cum(r,c))

  	cum = matmul(T,triu(ones))
 
  	allocate(ppi(n+1))	
  	
  	do k=1,size(X)
  		state(:,k) = s
  		ppi = (/0.,matmul(s,cum)/)
  		s=0.
  		 do i = 1,r
  		 	if((X(k)<=ppi(i+1)) .and. (X(k)>ppi(i)))then
  		 		s(i)=1.
  		 	endif
  		 enddo
  	enddo
  	chain=matmul(V,state) 
  end subroutine
  
  !**********************************************************
  ! TRIANGULAR UPPER MATRIX
  !**********************************************************
  function triu(matrix)
  !**********************************************************
  ! TRIU gets the upper triangular part of a matrix 
  ! 		  
  ! Usage: 
  !	upper = triu(matrix)
    
  ! INPUTS
  !     matrix  : real value matrix

  ! OUTPUTS	
  !		upper   : upper part of matrix
  !***********************************************************
  	integer :: dim(2),i,j
  	real, allocatable :: triu(:,:)
  	real :: matrix(:,:)
  
  	triu = matrix
  	dim=shape(matrix)
  	
  	do i=1,dim(1)
  		do j=1,dim(2)
  			if (i>j)then
  				triu(i,j)=0.
  			endif
  		enddo
  	enddo
  end function 
  
  !**********************************************************
  ! PRINT BINARY 0
  !**********************************************************
  subroutine print_binary_0(filename,var)
  !**********************************************************
  ! PRINT_BINARY_0 prints scalar
  ! 		  
  ! Usage: 
  !	call print_binary_0(filename,var)
    
  ! INPUTS
  !     filename  : name of the file 
  !		var		  : scalar
  ! OUTPUTS	
  !		
  !***********************************************************
    real, intent(in) :: var
    character (len=*) :: filename
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine 
  
  !**********************************************************
  ! PRINT BINARY'S
  !**********************************************************
  subroutine print_binary_1(filename,var)
  !**********************************************************
  ! PRINT_BINARY_1 prints scalar
  ! 		  
  ! Usage: 
  !	call print_binary_1(filename,var)
    
  ! INPUTS
  !     filename  : name of the file 
  !		var		  : vector
  ! OUTPUTS	
  !		
  !***********************************************************
    real, intent(in) :: var(:)
    character (len=*) :: filename
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine 
  
  subroutine print_binary_2(filename,var)
  !**********************************************************
  ! PRINT_BINARY_2 prints scalar
  ! 		  
  ! Usage: 
  !	call print_binary_1(filename,var)
    
  ! INPUTS
  !     filename  : name of the file 
  !		var		  : matrix
  ! OUTPUTS	
  !		
  !***********************************************************
    real, intent(in) :: var(:,:)
    character (len=*) :: filename
	integer :: o
	open(unit=1,file=filename//'.txt',status='replace', action='write')
	write(1,*) var
	close(1)
  end subroutine
  
  subroutine print_binary_3(filename,var)
  !**********************************************************
  ! PRINT_BINARY_3 prints scalar
  ! 		  
  ! Usage: 
  !	call print_binary_1(filename,var)
    
  ! INPUTS
  !     filename  : name of the file 
  !		var		  : 3rd rank tensor 
  ! OUTPUTS	
  !		
  !***********************************************************
    integer :: dim(3), i
    character(10) :: num
    real, intent(in) :: var(:,:,:)
    character (len=*) :: filename
    
    dim = shape(var)
    if(dim(3)==1)then   
		open(unit=1,file=filename//'.txt',status='replace', action='write')
		write(1,*) var(:,:,1)
		close(1)
	else 
		do i = 1,dim(3)
			open(unit=1,file=filename//'_'//trim(str(i))//'.txt',status='replace', action='write')
			write(1,*) var(:,:,i)
			close(1)
		enddo
	endif
  end subroutine 
  
  !**********************************************************
  ! LOAD BINARY FILE
  !**********************************************************
  function load_binary_0(filename) result(var)
  !**********************************************************
  ! PRINT_BINARY_0 prints scalar
  ! 		  
  ! Usage: 
  !	var = load_binary_0(filename)
    
  ! INPUTS
  !     filename  : name of the file 
  ! OUTPUTS	
  !		var		  :
  !***********************************************************
    real :: var
    character (len=*) :: filename
	OPEN (unit = 1, file = filename, status = 'old')
	READ (1,*) var
	CLOSE (1)
  end function
  
  !**********************************************************
  ! MEAN
  !**********************************************************
  function mean(var)
  !**********************************************************
  ! MEAN calculates the mean value of matrix var through columns 
  ! 		  
  ! Usage: 
  !	meanvalue = mean(var)
    
  ! INPUTS
  !     var  	: real value m x n matrix
  ! OUTPUTS	
  !		result 	: n-vector with mean values 
  !***********************************************************
    integer :: i, j
    real :: var(:,:)
    real, allocatable :: mean(:)   	
   	allocate(mean(size(var,2)))
   	do j = 1,size(var,2)
  		mean(j) = sum(var(:,j))/size(var,1)
  	enddo 		
  end function 
  
  !**********************************************************
  ! EASY MARKOV CHAIN
  !**********************************************************
  function markov_chain(T,n,sO) result(chain)
  !**********************************************************
  ! MARKOV_CHAIN simulates a markovian process with T transition matrix
  ! 		  
  ! Usage: 
  !	chain = markov_chain(T,n,sO)
    
  ! INPUTS
  !     T  		: matrix transition 
  !     n  		: number of periods to simulate
  !     s0  	: initial point

  ! OUTPUTS	
  !		chain  	: 
  !***********************************************************
    integer :: i, j, time, n, sO, chain(n)
    real :: T(:,:), x 
 	chain = 1
 	chain(1) = sO
 	
	do time = 2,n 
		! Find current state
		i = chain(time-1)		
		call random_number(x)
!		x = rand()
		do j = 1,size(T,2)
	    	if(sum(T(i,1:j)) < x) then
	    		chain(time) = chain(time) + 1
 		   	endif
 		enddo
 		if(chain(time)>size(T,2))then
			chain(time) = size(T,2)	
		endif
    enddo
  end function 

  !**********************************************************
  ! BIND
  !**********************************************************
  function bind(matrix)
  !**********************************************************
  ! BIND 
  ! 		  
  ! Usage: 
  !	result = bind(matrix)
    
  ! INPUTS
  !     matrix  	: real value matrix

  ! OUTPUTS	
  !		result  	: matrix binded through rows
  !***********************************************************
	integer :: i,j
	real :: matrix(:,:)
	real, allocatable:: bind(:,:)
	
	allocate(bind(size(matrix,1), size(matrix,2)))
	
	bind = matrix
	do i = 1,size(matrix,1)
		do j = 2,size(matrix,2)
			if(matrix(i,j-1)>0 .and. matrix(i,j)<0)then
				bind(i,j)=0.
			end if
			if(matrix(i,j-1)<0 .and. matrix(i,j)>0)then
				bind(i,j-1)=0.
			end if 
		end do 
	end do 
  end function 
  
  !**********************************************************
  ! PACKING
  !**********************************************************
  function packing(matrix,limit,replace) 
  !**********************************************************
  ! PACKING 
  ! 		  
  ! Usage: 
  !	packing(matrix,limit,replace) 
    
  ! INPUTS
  !   

  ! OUTPUTS	
  !		
  !***********************************************************
  integer :: i,j
  real:: matrix(:,:)
  real, allocatable :: packing(:,:)
  real :: limit, replace
  
  allocate(packing(size(matrix,1),size(matrix,2)))
	do i = 1,size(matrix,1)
  		do j = 2,size(matrix,2)
			if( matrix(i,j)<=limit)then
				matrix(i,j)=replace
			end if 
		end do 
	end do 
  end function
  
  !**********************************************************
  ! SUB2IND
  !**********************************************************
  function sub2ind(dim,row,col) result(ind)
  !**********************************************************
  ! SUB2IND converts sub index to linear index 
  ! 		  
  ! Usage: 
  !	index = sub2ind(dim,row,col)
    
  ! INPUTS
  !   	dim		: desired shape
  !		row		: row index
  !		col		: column index

  ! OUTPUTS	
  !		index 	: linear index
  !***********************************************************
  	integer :: dim(2), row, col, ind
  	ind = row + dim(1)*(col-1)
  end function    
  
  !**********************************************************
  ! IND2SUB
  !**********************************************************
  function ind2sub(dim,ind) result(sub)
  !**********************************************************
  ! IND2SUB 
  ! 		  
  ! Usage: 
  !	sub = ind2sub(dim,ind) 
    
  ! INPUTS
  !   	dim		: desired shape
  ! 	ind		: linear index
  ! OUTPUTS	
  !		sub		: row, col indexs
  !***********************************************************
  	integer :: dim(2), sub(2), ind
  	sub(2) = int(ind/real(dim(1))-0.001+1)
  	sub(1) = ind - dim(1)*(sub(2)-1)
  end function 
  
  !**********************************************************
  ! MASK
  !**********************************************************
  function mask(len,pos)
  !**********************************************************
  ! MASK make mask of logical values putting .false. int pos
  ! 		  
  ! Usage: 
  !	result = mask(len,pos)  
    
  ! INPUTS
  !   	len		:
  !		pos		:

  ! OUTPUTS	
  !		result 	:
  !***********************************************************
  	integer :: pos 
  	integer :: len
  	logical :: mask(len)
  	
  	mask = .true.
  	mask(pos) = .false.	
  end function
  
  !**********************************************************
  ! CUMSUM
  !**********************************************************
  function cumsum(array) 
  !**********************************************************
  ! CUMSUM makes the cumulative sum through rows 
  ! 		  
  ! Usage: 
  !	result = cumsum(array) 
    
  ! INPUTS
  ! 	array	: real valued matrix of size n x m 

  ! OUTPUTS	
  !		result	: real valued matrix of size n x m
  !***********************************************************
  	integer :: i,j 
  	real :: array(:,:)
  	real, allocatable :: cumsum(:,:)
  	allocate(cumsum(size(array,1),size(array,2)))
  	
  	cumsum(:,1) = array(:,1)
  	do i =1,size(cumsum,1)
  		do j=2,size(cumsum,2)
			cumsum(i,j) =sum(array(i,1:j-1))+array(i,j) 
  		end do 
  	end do 
  end function 
  
  !**********************************************************
  ! BETAMULTINOMIAL
  !**********************************************************
  function betamultinomial(s,n) result(Ps) 
  !**********************************************************
  ! BETAMULTINOMIAL
  ! 		  
  ! Usage: 
  !	Ps = betamultinomial(s,n)
    
  ! INPUTS
  ! 	s	: 
  ! 	n	: 

  ! OUTPUTS	
  !		Ps	: 
  !***********************************************************
  	integer :: i,j,k, ns, TT, t
  	real :: s(:), n(:,:) 
  	
  	real, allocatable :: Ps(:,:,:), Ps1(:,:,:), m(:,:,:), NN(:,:,:,:), NT(:,:,:), Ep(:,:,:,:)
  	
  	TT = size(s)
  	ns = size(n,1)
  	allocate(Ps(ns,ns,TT)) 	
  	allocate(Ps1(ns,ns,TT))
  	allocate(m(ns,ns,TT))
  	allocate(NN(ns,ns,TT,1))
  	allocate(NT(ns,TT,1))
  	allocate(Ep(ns,ns,TT,1))
  	
  	m(:,:,1) = n
  	do t = 2,TT
  		do i = 1,ns
  			do j = 1,ns
  				if(s(t) == j .and. s(t-1)==i)then 
  					m(i,j,t) = m(i,j,t-1)+1
  				else
  					m(i,j,t) = m(i,j,t-1)
  				end if
  			end do
  		end do 
  	end do 
  	
  	do i = 1,ns
  		do j = 1,ns
  			NN(i,j,:,:) = reshape(m(i,j,:),(/size(m(i,j,:)),1/))  
  			if(j==1)then 
  				NT(i,:,:) = NN(i,j,:,:)	
  			else
  				NT(i,:,:) = NT(i,:,:) + NN(i,j,:,:)	
  			end if 		
  		end do
  		do k = 1,ns
  			Ep(i,k,:,:) = NN(i,k,:,:)/NT(i,:,:)
  		end do
  	enddo 
	do i = 1,ns
		do j = 1,ns
			Ps(i,j,:) = reshape(Ep(i,j,:,:),(/size(EP(i,j,:,:))/))
		end do 
	end do 
  end function 
  
  !**********************************************************
  ! PATHFINDER
  !**********************************************************
  
  function pathfinder(dlower,dupper,slacktry, cTtry, pai, nd,ny) result(Ztry)

	integer, intent(in) :: nd, ny
	integer :: n
	integer :: i, j, dim(2)
	real :: dlower, dupper, start, finish, cumtime(5)=0.
	
	real, dimension(ny*nd,nd) :: temp, Ztry, slacktry, cTtry
	real, dimension(ny,ny) :: pai, temp1
	real, dimension(ny*nd,ny) :: Y,Y_b,X
	real, dimension(ny,nd) :: TARP
	real, dimension(ny*nd) :: Y_sum
	real, dimension(1,ny) :: temp3
		
	n=ny*nd
	
	temp=0.

	where(cTtry<=0 .or. slacktry<0) temp = 1.
	
	TARP = 0.
	do i=1,size(temp,1)
		if(sum(temp(i,:)) == nd)then
			dim = ind2sub(shape(TARP),i)
			TARP(dim(1),dim(2)) = 1.
		end if 
	end do 
	
	dim = ind2sub(shape(TARP),4)
	where(pai>0)
		temp1=1
	end where

	X = repmat(temp1,nd,1)	
	Ztry = 8e10
!	write(*,*)'a'
	
!	call print_matrix(slacktry(34000:,92:101))
!	stop
	do i = 1,nd
!		call cpu_time(start)
		Y = spread(TARP(:,i),1,n)
!		call cpu_time(finish)
!		cumtime(1)=cumtime(1)+finish-start
!		call cpu_time(start)
		Y = Y*X
!		call cpu_time(finish)
!		cumtime(2)=cumtime(2)+finish-start
!		call cpu_time(start)
		Y_sum = sum(Y,1)
		do j = 1,size(Y,1)
			if (Y_sum(j)==0.)then  
				Ztry(j,i)=1.
			end if 
		end do
!		call cpu_time(finish)
!		cumtime(3)=cumtime(3)+finish-start
!		write(*,*) i
!		if(i==nd)then 
!			call print_matrix(Ztry(1:1000,95:101))
!			stop
!		end if
	end do

!	write(*,*) cumtime(1)
!	write(*,*) cumtime(2)	
!	write(*,*) cumtime(3)

end function 

function unique(array) result(vector2)

	integer :: i,j,e
	real :: array(:)
	real, allocatable :: vector(:), vector2(:)
	
	allocate(vector(size(array))) 	
	vector = 0.
	vector(1) = array(1) 
	e = 1
	do i = 2,size(array)
		if(array(i)/= vector(e))then 
			vector(e+1) = array(i)
			e = e + 1
		end if 
	end do 
	
	allocate(vector2(size(pack(vector,vector>0))))
	vector2 = pack(vector,vector>0.)

end function 

function lagg(y,L)

	real :: y(:,:)
	integer :: L, T, j, M, i
	real, allocatable :: lagg(:,:)

	if (L==0)then 
		write(*,*) 'L can not be zero'
	end if 
    T = size(y,1)
    M = size(y,2)     

	allocate(lagg(T-L,(L+1)*M))
	
    lagg(:,(L+1)*M-(M-1):)= y(:T-L,:)
        
   	do i = 0,M-1
    do j = 1,L
   	!	write(*,*) (L+1)*M - j*(M)-i, j+1,'h',j+1+(T-L)
   	!   lagg(:,(L+1)*M - j*(M)-i) = y(j+1:j+1+(T-L),(M-i))
    end do
    end do
end function 
end module tools
