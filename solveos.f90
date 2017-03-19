module solveos

  use basic
  use def_types
  use tools
  implicit none
  
contains
  
  !***********************************************************************
  ! SOLVEDPOS
  !***********************************************************************
  subroutine solvedpos(vbad,xbad,vgood,xgood,default, f,fd,prob,beta,lambda,azero, n,m,n1,n2, iterglobal,respuesta)

    integer :: iterglobal
    logical, intent(out) :: respuesta
    
    real :: start, finish 
    
    !Local 
    integer, parameter :: maxit = 1000		! maximum number of iterations
    real, parameter :: tol = 10e-6 			! convergence tolerance, usually tol=sqrt(eps)
    logical, parameter :: prtiters = .true.	! print iterations (1) or not (0)

    integer :: it,i
    real :: change
    
    !Dummy
    integer, intent(out) :: xbad(n), xgood(n)
    real, intent(out) :: vbad(n), vgood(n) 
    
    integer, intent(in) :: n1, n2, n, m, azero
    real, intent(in) :: f(n,m), fd(n,m), prob(n1,n1), lambda, beta(n,m)

    !**********definir
    
    real :: vbadold(n), vgoodold(n)
    logical :: default(n) 
    real :: vbg_p(m,n/m), vbg_pp(n/m), vbg(m,n/m), vbg_ppp(n1), vbg_pppp(m,n1), vbg_ppppp(n)
    
    real :: Evbad(n,m), Evgood(n,m), Evbg(n,m)

    real, dimension(n,2) :: v, vold
    
    !************ Perform policy or value function iterations
    
    
    vbad = 1.
    vgood = 1.
    vbg_p = reshape(vgood,(/m,int(n/m)/))
    vbg_pp = vbg_p(azero,:)
    vbg = repmat(reshape(vbg_pp,(/1,(n/m)/)),m,1)
    
    !************
     
    if(prtiters)then
       write(*,*) 'Solving Bellman equation by Value function iteration'
    end if
    
    do it = 1,maxit
       vbadold = vbad
       vgoodold = vgood

       Evgood = kronecker(matmul(prob,(transpose(reshape(vgood,(/m,int(n/m)/))))), ones(m,1))
       Evbad  = kronecker(matmul(prob,(transpose(reshape(vbad ,(/m,int(n/m)/))))), ones(m,1))
       Evbg   = kronecker(matmul(prob,(transpose(reshape(vbg  ,(/m,int(n/m)/))))), ones(m,1))

       vbad = maxval(fd + lambda*(beta*Evbg) + (1-lambda)*beta*Evbad,2)
       xbad = maxloc(fd + lambda*(beta*Evbg) + (1-lambda)*beta*Evbad,2)

       vgood = maxval(f + beta*Evgood,2)       
       xgood = maxloc(f + beta*Evgood,2)
       
       default = vbad>vgood .or. isnan(vgood)
       vgood(find(default)) = vbad(find(default))
       vbg = reshape(vgood,(/m,n/m/)) 
       vbg_ppp = vbg(azero,:)
       vbg_pppp = repmat(reshape(vbg_ppp,(/1,(n/m)/)),m,1)
!       vbg_pppp = repmat(vbg_ppp,m,1) 
       vbg_ppppp = reshape(vbg_pppp,(/n/))
       change = norm2((/vbad,vgood/)-(/vbadold,vgoodold/))	! compute change
       if(prtiters)then
	  write(*,*) it, ' ', change 							! print progress
       end if
       if(change<tol)then 
	  exit													! convergence check
       end if
       !***************
       !Debug
       !***************
       if( it == 10000 .and. iterglobal == 1000 )then
	  respuesta = .true.
	  exit
       end if 
       !***************
    end do    
    if (change>tol)then 
       write(*,*) 'Failure to converge in solvedp'
    end if 
  end subroutine solvedpos
end module solveos
