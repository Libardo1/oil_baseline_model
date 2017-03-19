program run_constrained_V2F_OF

	use basic
	use tools
	use solve
	use def_types
	use hamann_tools
	use omp_lib
	
	implicit none 
	
	real :: start, finish
	
	integer :: i,j,k,l,time
	integer :: right_name, hh_on, yT_uncertainty
	real ::	fixc, temp(1,1)
	
	character(100) :: kof,mod_modif, learning_type, yT_nature, file, place

	real :: Pj(2), tj(2), Pp_aux(2), pr(2), p(2)
	real :: En, R, x2s, kappa, Ep, Epr, disc, Ex, Es, q
	
	real, dimension(2,2) :: Pp, Pp_RE, Pp_BL
	
	real :: rstar, betta, a, yN, xi, sigg, aT, aN, kapa

	integer, parameter :: JJ = 19, T = 25
	real :: p_t(T), n0(2,2), M_Pxstar(2,2,T),smax
	real, allocatable ::s(:), grid(:,:), SS(:)
	
	real, allocatable :: PX(:), PX_BL(:), PX_RE(:)
	 
	integer :: tt(T), ip_t(T),R1(T), R2(T), ns, n ,np
	
	integer :: Npstar, nd, NNy, warnn, NNS, ny
	real :: dlower, dupper, rho_yT, sigma_eyT
	character(1) :: eqs 
	real, allocatable :: YYT(:), YYTprob(:,:)
	real :: yt_bar

	real, allocatable, dimension(:,:,:) :: M_spr, M_Xr, M_RREVXr, M_PROFIT_Xr, M_Pxstar_HH, M_v, M_crisis
	real, allocatable, dimension(:,:,:) :: M_p, M_cT, M_c, M_la, M_mu, M_slack, M_dpix, M_dp
	 
	real, allocatable, dimension(:,:)	:: spr, RREVXr, PROFIT_Xr, Xr, pstar
	real, allocatable, dimension(:)		:: sgrid, sp, RREVX, PROFIT_X

	real, allocatable, dimension(:,:)	:: RREVXr_BL, PROFIT_Xr_BL, Xr_BL, pstar_BL
	real, allocatable, dimension(:)		:: sp_BL, RREVX_BL, PROFIT_X_BL

	real, allocatable, dimension(:,:)	:: RREVXr_RE, PROFIT_Xr_RE, Xr_RE, pstar_RE
	real, allocatable, dimension(:)		:: sp_RE, RREVX_RE, PROFIT_X_RE

	
	real, allocatable, dimension(:,:) :: c, cT, d, dp, la, pai, yT, p_p, crisis, v, mu, slack
	real, allocatable, dimension(:) :: dgrid
	integer, allocatable :: dpix(:,:)
	real :: distTol
	integer ::dist_dpix
	
	integer, allocatable :: dpix_ini(:,:)
	
	real, allocatable :: ygrid(:),yTgridlevel(:), revX(:,:),revXgridlevel(:)
	
	real, allocatable :: Pxstar(:,:), Pxstar_BL(:,:), Pxstar_RE(:,:)
	
	call cpu_time(start)

	!! Default settings to run code
	right_name = 0 				! This value must be zero always
 
	! Setting model block's
	! Kind of Firm
	
	kof = trim('Type_2')
 	
	hh_on  = 1
	
	! Model's modification
	mod_modif = trim('')
	
	! Setting the Learning Type
	learning_type =  trim('RE')
 
	! Setting uncertainty in YT
 
	yT_uncertainty = 0
	
	if (yT_uncertainty==1)then 
   		yT_nature = trim('YT_VAR')
	else
   		yT_nature = trim('YT_FIXED')
	end if
	
	! Model's Calibration
	! Oil's Firm
    			
    select case (kof)
    case('Type_1')
    	q     = 0.991436520593215 	! Implies an annual efective risk free rate of 3.5%    % 1.1^-0.25       % disc rate (0.95 -> about 5% real interest rate)
		kappa = 15.425  			! this value is to fit 6.3*4 quarters of oil reserves % 1 is low kappa  % mg cost extraction parameter
 		disc  = 0.056   			! fixed discoveries, assume 10% of possible reserves
 		fixc  = 0       			! fixed costs
 
    case('Type_2')
   
		Pj      = (/22.4787, 56.8691/)      	! [pl ph]   values of each state
		tj      = (/58.66, 62.13/)        	! [tsl tsh] periods in each state
		Pp_aux  = (/1-1/tj(1), 1-1/tj(2)/)	! [qll qhh] prob. of staying in each state
		pr	 	= (1/Pj(2))*Pj          	! price normalization ph==1
		Pp(1,:) = (/Pp_aux(1)   , 1-Pp_aux(1)/)
		Pp(2,:)	= (/1-Pp_aux(2) , Pp_aux(2)/)
		temp	= matmul(reshape(pr,(/1,size(pr)/)),ergdist(Pp,1./size(Pp,1)*ones(size(Pp,1),1))) ! expected oil price
				
		Epr     = temp(1,1)
		p		= pr/Epr						! re-normalize prices to Ep = 1
		temp	= matmul(reshape(p ,(/1,size(p )/)),ergdist(Pp,1./size(Pp,1)*ones(size(Pp,1),1))) ! re-compute expected value to 1
   		Ep  	= temp(1,1)
		
		En      = 25.
		R       = (1.035)**(.25)

		x2s     = 1/En
		kappa   = (-Ep*(1/R-1))/(2*x2s*(1-1/R)+(1/R)*x2s*x2s)
   
		disc    = 0.056         
		Ex      = disc
   	
		Es      = En*Ex
		q      	= R**(-1)
	end select
	
    
	!! Household
	rstar   =  1/q-1;     
	betta   =  0.99;         ! discount factor (in order to have nice distrubution of debt)
	
	! rstar 	= 1/0.99-1   %those values move the debt until upper bound
	! betta 	= 0.94
	
	a       = 0.493008      ! weight on traded consumption in the CES aggregator
 
	yN      = exp(-0.510734)! nontraded output
	xi      = 0.316         ! Elasticity of substitution between traded and nontraded goods
	sigg    = 2.              ! intertemporal elasticity of consumption
	aT      = 0.073025
	aN      = 0.27507
	kapa 	= 2.
 
	! Loading Path of Expected Probabilities
 
	do i =1,T
		tt(i)=i
	end do 
	R1   = 1
	R2   = 2
	where (tt<=JJ) 
		ip_t = R2
	elsewhere
		ip_t = R1
	end where
	
	do i = 1,T
		p_t(i)  = p(ip_t(i))
	enddo
	
	n0 = 0.01400  ! These initial counters no training
 	
	M_Pxstar(:,:,tt) = betamultinomial(real(ip_t),n0);
!	do i =1,T
!		write(*,*) i
!		call print_matrix(M_Pxstar(:,:,i))
!	end do 
!	stop
		
	! Construct Firm's state space
 
!   smax   = 0.056*34          
	smax = 1.8929
	allocate(s(int((smax-0.0001)/disc*5+1)))
	s      = linspace(0.0001,smax,int((smax-0.0001)/disc*5+1))
	ns     = size(s)
!	write(*,*) smax, ns
!   	write(*,*) s
!   	stop
   	
    select case(kof)
    case('Type_1')        
 
        p      = (/0.46, 1.54/) 	! to get std(PX)=0.54 from data
		np     = size(p)
		allocate(grid(ns*np,2))
		allocate(SS(ns*np))
		allocate(PX(ns*np))

		grid = gridmake2(s,p)  
		SS = grid(:,1)
		PX = grid(:,2)
		
    case('Type_2')
       
		np     = size(p)
		allocate(grid(ns*np,2))
		allocate(SS(ns*np))
		allocate(PX(ns*np))
		grid = gridmake2(s,p) 
   		SS = grid(:,1)
		PX = grid(:,2)
		
	end select
	n	= size(SS)
	
	! Construct Household's state space
	  
	if(yT_uncertainty==1) then
		NNy	= 3 	  ! To increase speed o/w it will take years in solving
		rho_yT    = 0.62
		sigma_eyT = 0.02
		allocate(YYTprob(NNy,NNy))
		allocate(YYT(NNy))
		call rouwenhorst(YYT,YYTprob, NNy, 0.0, rho_yT,sigma_eyT)
		nd     = 101   ! # of grid points for debt, d to speed up
	else
		NNy = 1
		allocate(YYTprob(NNy,NNy))
		allocate(YYT(NNy))	
		YYT = 0	
		YYTprob = 1
		nd     = 101   ! # of grid points for debt
	end if 

	yt_bar  =    exp(-1.06692)
	YYT = exp(YYT)*yt_bar
 
	Npstar = n
	NNS = NNy*Npstar
 
	dlower = 0
	dupper = 2.75  !(to have less nodes to speed up)   %3.5;
 
	!! Choosing Household's Equilibrium
	eqs   = trim('c')
	warnn =  0
	
	! Prealloacaion Matrices
	allocate(M_spr		(np,ns,T))
	allocate(M_Xr		(np,ns,T))
	allocate(M_RREVXr	(np,ns,T))
	allocate(M_PROFIT_Xr(np,ns,T))
	allocate(M_Pxstar_HH(n,n,T))
	
	allocate(M_v	(NNS,nd,T))
	allocate(M_crisis(NNS,nd,T))
	allocate(M_p	(NNS,nd,T))
	allocate(M_cT	(NNS,nd,T))
	allocate(M_c	(NNS,nd,T))
	allocate(M_la	(NNS,nd,T))
	allocate(M_mu	(NNS,nd,T))
	allocate(M_slack(NNS,nd,T))
	allocate(M_dpix	(NNS,nd,T))
	allocate(M_dp	(NNS,nd,T))
	
	allocate(spr	(np,ns))
	allocate(RREVXr	(np,ns))
	allocate(PROFIT_Xr(np,ns))
	allocate(Xr		(np,ns))

	allocate(RREVXr_BL	(np,ns))
	allocate(PROFIT_Xr_BL(np,ns))
	allocate(Xr_BL		(np,ns))

	allocate(RREVXr_RE	(np,ns))
	allocate(PROFIT_Xr_RE(np,ns))
	allocate(Xr_RE		(np,ns))
		
	allocate(sgrid(ns))
	
	allocate(sp(n))
	allocate(RREVX(n))
	allocate(PROFIT_X(n))
	allocate(pstar(n,n))
	
	allocate(sp_BL(n))
	allocate(RREVX_BL(n))
	allocate(PROFIT_X_BL(n))
	allocate(pstar_BL(n,n))

	allocate(sp_RE(n))
	allocate(RREVX_RE(n))
	allocate(PROFIT_X_RE(n))
	allocate(pstar_RE(n,n))
	
	allocate(v		(NNS,nd))
	allocate(crisis	(NNS,nd))
	allocate(p_p	(NNS,nd))
	allocate(cT		(NNS,nd))
	allocate(c		(NNS,nd))	
	allocate(la		(NNS,nd))
	allocate(mu		(NNS,nd))
	allocate(slack	(NNS,nd))	
	allocate(dpix	(NNS,nd))
	allocate(dp		(NNS,nd))
	
	allocate(d	(NNS,nd))
	allocate(yT	(NNS,nd))
	allocate(revX(NNS,nd))
	
	allocate(pai(NNS,NNS))
	allocate(dgrid(nd))
	allocate(ygrid(NNS))
	allocate(yTgridlevel(NNS))
	allocate(revXgridlevel(NNS))
	
	allocate(Pxstar(n,n))
	allocate(Pxstar_BL(n,n))
	allocate(Pxstar_RE(n,n))
	
	allocate(dpix_ini(NNS,nd))

	!!! all data is ok
!	write(*,*) 'learning'
	select case(learning_type)
	case('RE')
		right_name = 1
		select case(kof)
		case('Type_1')
		case('Type_2')		
			call oilcompany_franz_v3(PROFIT_X,RREVX,pstar,sp,Xr,PROFIT_Xr,RREVXr,sgrid,spr,  &
								q,kappa,disc,s,Pp,SS,PX,n,ns,np)
		end select
		
!	write(*,*) Xr(1,:)
!	write(*,*) 
!	write(*,*) 
!	stop
		Pxstar = pstar
		
      	M_spr(:,:,1)       = spr ! ok
     	M_Xr(:,:,1)        = Xr
      	M_RREVXr(:,:,1)    = RREVXr
      	M_PROFIT_Xr(:,:,1) = PROFIT_Xr
    
		if(hh_on==1)then
			M_Pxstar_HH(:,:,1) = Pxstar
        	if(warnn==0)then
!        		write(*,*) 'Constrained'  , ' warnn = 0'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar,rstar,betta,a,yN,xi,sigg,			&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'RE')						! input
!       		write(*,*) 'Finish Constrained'
        	end if
        	if(warnn==1)then
        		eqs = 'a'
!        		write(*,*) 'Constrained'  , ' warnn = 1'  
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar,rstar,betta,a,yN,xi,sigg,			&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'RE')						! input  		
!        		write(*,*) 'Finish Constrained'
        		warnn = 1
        	end if 

	        dpix_ini = dpix
       		call print_binary_2('outData/init_BL',real(dpix_ini))
    	    M_p(:,:,1)       = p_p
  	      	M_cT(:,:,1)      = cT
        	M_c(:,:,1)       = c
    	    M_dpix(:,:,1)    = dpix
      		M_dp(:,:,1)      = dp
      		
      		call print_matrix_dat('outData/constrained/dgrid.txt',reshape(dgrid,(/size(dgrid),1/)))      		
      		call print_matrix_dat('outData/constrained/p.txt',p_p)
			call print_matrix_dat('outData/constrained/cT.txt',cT)
			call print_matrix_dat('outData/constrained/c.txt',c)
			call print_matrix_dat('outData/constrained/dpix.txt',real(dpix))
			call print_matrix_dat('outData/constrained/dp.txt',dp)
			
      	end if
    case('BL')
    	right_name = 1
!    !$omp parallel do 
    do time = 1,T
    	write(*,*) 'T', time
    	Pp = M_pxstar(:,:,time)
		select case(kof)
		case('Type_1')
		case('Type_2')		
			call oilcompany_franz_v3(PROFIT_X,RREVX,pstar,sp,Xr,PROFIT_Xr,RREVXr,sgrid,spr,  &
								q,kappa,disc,s,Pp,SS,PX,n,ns,np)
		end select
		
		Pxstar = pstar
      	M_spr(:,:,time)       = spr 
     	M_Xr(:,:,time)        = Xr
      	M_RREVXr(:,:,time)    = RREVXr
      	M_PROFIT_Xr(:,:,time) = PROFIT_Xr
      	
      	if(hh_on==1)then
			M_Pxstar_HH(:,:,time) = Pxstar
        	if(warnn==0)then
!        		write(*,*) 'Constrained'  , ' warnn = 0'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar,rstar,betta,a,yN,xi,sigg,			&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        	end if
        	if(warnn==1)then
        		eqs = 'a'
!        		write(*,*) 'Constrained'  , ' warnn = 1'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar,rstar,betta,a,yN,xi,sigg,			&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        		warnn = 1
        	end if 
	        dpix_ini = dpix
!       		call print_matrix_dat('outData/init_BL.txt',real(dpix_ini))
       		call print_binary_2('outData/init_BL',real(dpix_ini))
      		
    	    M_p(:,:,time)       = p_p
  	      	M_cT(:,:,time)      = cT
        	M_c(:,:,time)       = c
        	M_la(:,:,time)      = la
    	    M_dpix(:,:,time)    = dpix
      		M_dp(:,:,time)      = dp			
      	end if
    end do 
!	!$omp end parallel do 
    case('REF_BLHH')
    
  	right_name = 1
    do time = 1,T
		select case(kof)
		case('Type_1')
!		   	Pp = M_pxstar(:,:,time)
		case('Type_2')		
		   	Pp_BL = M_pxstar(:,:,time)
			call oilcompany_franz_v3(PROFIT_X_BL,RREVX_BL,pstar_BL,sp_BL,Xr_BL,PROFIT_Xr_BL,RREVXr_BL,sgrid,spr,  &
								q,kappa,disc,s,Pp_BL,SS,PX_BL,n,ns,np)
			Pp_RE = Pp
			call oilcompany_franz_v3(PROFIT_X,RREVX,pstar,sp,Xr,PROFIT_Xr,RREVXr,sgrid,spr,  &
								q,kappa,disc,s,Pp_RE,SS,PX,n,ns,np)
		end select
		
		Pxstar_BL = pstar_BL
		
      	M_spr(:,:,1)       = spr 
     	M_Xr(:,:,1)        = Xr
      	M_RREVXr(:,:,1)    = RREVXr
      	M_PROFIT_Xr(:,:,1) = PROFIT_Xr
      	
      	if(hh_on==1)then
			M_Pxstar_HH(:,:,time) = Pxstar
        	if(warnn==0)then
!        		write(*,*) 'Constrained'  , ' warnn = 0'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar_BL,rstar,betta,a,yN,xi,sigg,		&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        	end if
        	if(warnn==1)then
        		eqs = 'a'
!        		write(*,*) 'Constrained'  , ' warnn = 1'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar_BL,rstar,betta,a,yN,xi,sigg,		&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        		warnn = 1
        	end if 
    	    M_v(:,:,time)       = v
    	    M_crisis(:,:,time)  = crisis
    	    M_p(:,:,time)       = p_p
  	      	M_cT(:,:,time)      = cT
        	M_c(:,:,time)       = c
        	M_la(:,:,time)      = la
    	    M_mu(:,:,time)      = mu
       	    M_slack(:,:,time)   = slack
    	    M_dpix(:,:,time)    = dpix
      		M_dp(:,:,time)      = dp			
      	end if
    end do 
  
    case('BLF_REHH')
    
  	right_name = 1
    do time = 1,T
		select case(kof)
		case('Type_1')
		case('Type_2')		
		   	Pp_BL = M_pxstar(:,:,time)
			call oilcompany_franz_v3(PROFIT_X,RREVX,pstar,sp,Xr,PROFIT_Xr,RREVXr,sgrid,spr,  &
								q,kappa,disc,s,Pp_BL,SS,PX,n,ns,np)
			Pp_RE = Pp
			call oilcompany_franz_v3(PROFIT_X_RE,RREVX_RE,pstar_RE,sp_RE,Xr_RE,PROFIT_Xr_RE,RREVXr_RE,sgrid,spr,  &
								q,kappa,disc,s,Pp_RE,SS,PX_RE,n,ns,np)
			Pp_RE = Pp
		end select
		
		Pxstar_RE = pstar_RE
		
      	M_spr(:,:,1)       = spr 
     	M_Xr(:,:,1)        = Xr
      	M_RREVXr(:,:,1)    = RREVXr
      	M_PROFIT_Xr(:,:,1) = PROFIT_Xr
      	
      	if(hh_on==1)then
			M_Pxstar_HH(:,:,time) = Pxstar_RE
        	if(warnn==0)then
!        		write(*,*) 'Constrained'  , ' warnn = 0'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar_RE,rstar,betta,a,yN,xi,sigg,		&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        	end if
        	if(warnn==1)then
        		eqs = 'a'
!        		write(*,*) 'Constrained'  , ' warnn = 1'    		
				call constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
									la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p_p,  				&	! output 
									YYT,YYTprob,RREVX,Pxstar_RE,rstar,betta,a,yN,xi,sigg,		&	! input
									aT,aN,kapa,nd,NNy*n,dlower,dupper,eqs,'BL')						! input
        		warnn = 1
        	end if 
    	    M_v(:,:,time)       = v
    	    M_crisis(:,:,time)  = crisis
    	    M_p(:,:,time)       = p_p
  	      	M_cT(:,:,time)      = cT
        	M_c(:,:,time)       = c
        	M_la(:,:,time)      = la
    	    M_mu(:,:,time)      = mu
       	    M_slack(:,:,time)   = slack
    	    M_dpix(:,:,time)    = dpix
      		M_dp(:,:,time)      = dp			
      	end if
    end do 
	end select
    call cpu_time(finish)
    file = trim(trim('constrained_V2F_OF_')//trim(mod_modif)//trim('_')//trim(kof)//'_'//&
    			trim(learning_type)//trim('_')//trim(yT_nature)//trim('_Nnod')//trim(str(nd)))//'/'

    write(*,*) 'Model :', file
	place = '/outData/'//file
	place = 'outData/example/'
	!**********************************************
	! Saving options 

!	call print_binary_2(trim(place)//'c',c)
!	call print_binary_2(trim(place)//'cT',cT)
!	call print_binary_2(trim(place)//'d',d)
!	call print_binary_2(trim(place)//'dp',dp)
!	call print_binary_2(trim(place)//'dpix.txt',int(dpix))
!	call print_binary_2(trim(place)//'dpix_ini.txt',int(dpix_ini))
!	call print_binary_2(trim(place)//'la',la)
!	call print_binary_2(trim(place)//'yT',yT)
	
	if(1==1)then
	call print_binary_3(trim(place)//'M_v',M_v)
	call print_binary_3(trim(place)//'M_crisis',M_crisis)
	call print_binary_3(trim(place)//'M_p',M_p)
	call print_binary_3(trim(place)//'M_cT',M_cT)
	call print_binary_3(trim(place)//'M_c',M_c)
	call print_binary_3(trim(place)//'M_la',M_la)
	call print_binary_3(trim(place)//'M_mu',M_mu)
	call print_binary_3(trim(place)//'M_slack',M_slack)
	call print_binary_3(trim(place)//'M_dpix',real(M_dpix))
	call print_binary_3(trim(place)//'M_dp',M_dp)
	
	call print_matrix_dat(trim(place)//'M_dp.txt',M_dp(:,:,1))
	call print_matrix_dat(trim(place)//'M_p.txt',M_p(:,:,1))
	call print_matrix_dat(trim(place)//'M_slack.txt',M_slack(:,:,1))


	call print_binary_3(trim(place)//'M_spr',M_spr)
	call print_binary_3(trim(place)//'M_Xr',M_Xr)
	call print_binary_3(trim(place)//'M_RREVXr',M_RREVXr)
	call print_binary_3(trim(place)//'M_PROFIT_Xr',M_PROFIT_Xr)
	call print_binary_3(trim(place)//'M_Pxstar_HH',M_Pxstar_HH)
	end if
	write(*,*) 'total time of execution: ', finish-start
	
	select case(learning_type)
	case('RE')
		call ergodic(100,int(0.01*100),1,M_Pxstar,M_RREVXr,M_Xr,M_spr,M_p,M_cT,M_c,M_dpix,M_dp,Pp,yN,&
					dgrid,sgrid,YYT,YYTprob,RREVX,PX,revXgridlevel,yTgridlevel,hh_on,nd,ny,ns,np,NNy,learning_type_for)
	case('BL')
		call ergodic(100,int(0.01*100),size(M_Pxstar,3),M_Pxstar,M_RREVXr,M_Xr,M_spr,M_p,M_cT,M_c,M_dpix,M_dp,Pp,yN,&
					dgrid,sgrid,YYT,YYTprob,RREVX,PX,revXgridlevel,yTgridlevel,hh_on,nd,ny,ns,np,NNy,learning_type_for)
	end select
	

end program 