program run_model

	use basic
	use tools
	use model_tools
	use omp_lib
	
	implicit none 
	
	real :: start, finish
	
	integer :: i,j,k,l,time
	integer :: hh_on
	real ::	fixc, temp(1,1)
	
	! It's important to define a correct dimensions for discretized variables
	! before run the code in order to work propertly. This avoid the use of 
	! dynamic memory 
	integer :: np = 2, ns = 70, nd  = 50
	
	character(100) :: mod_modif, learning_type, file, place

	real, dimension(np) :: Pj, tj, Pp_aux, pr, p
	real :: En, R, x2s, kappa, Ep, Epr, disc, Ex, Es, q
	
	real, dimension(np,np) :: Pp, Pp_RE, Pp_BL
	
	real :: rstar, betta, a, yN, xi, sigg, aT, aN, kapa

	integer, parameter :: JJ = 19, T = 25
	real :: p_t(T), n0(np,np), M_Pxstar(np,np,T),smax
	real ::s(ns)
		 
	integer :: tt(T), ip_t(T),R1(T), R2(T)
	
	real :: dlower, dupper, rho_yT, sigma_eyT
	real :: YYT, yt_bar

	! Firm's structures 
	real, dimension(np,ns,T) :: M_spr, M_Xr, M_RREVXr, M_PROFIT_Xr, M_EE
	! Household's structures 
	real, dimension(np*ns,nd,T) :: M_p, M_cT, M_c, M_la, M_mu, M_slack, M_dpix, M_dp,  M_v, M_crisis
	 
	real, dimension(np,ns)	:: spr, RREVXr, PROFIT_Xr, Xr
	real, dimension(ns)		:: sgrid
	
	real, dimension(np*ns,nd) :: c, cT, d, dp, la, pai, yT, p_p, crisis, v, mu, slack
	real, dimension(nd) :: dgrid
	real :: distTol
	real :: BB(nd)
		
	real :: Pxstar(np,np)
	 
	! Setting model block's 	
	hh_on  = 1
	
	! Model's modification
	mod_modif = trim('')
	
	! Setting the Learning Type
	learning_type =  trim('RE')
 	
	! Model's Calibration
	! Oil's Firm
    			
	Pj      = (/22.4787, 56.8691/)      	! [pl ph]   values of each state
	tj      = (/58.66, 62.13/)        		! [tsl tsh] periods in each state
	Pp_aux  = (/1-1/tj(1), 1-1/tj(2)/)		! [qll qhh] prob. of staying in each state
	pr	 	= (1/Pj(2))*Pj          		! price normalization ph==1
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
    
	!! Household
	rstar   =  1/q-1;     
	betta   =  0.99;         ! discount factor (in order to have nice distrubution of debt)	
	a       =  0.493008      ! weight on traded consumption in the CES aggregator
 
	yN      = exp(-0.510734)	! nontraded output
	xi      = 0.316         	! Elasticity of substitution between traded and nontraded goods
	sigg    = 2.           		! intertemporal elasticity of consumption
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
 
	smax = 1.8929
	s      = linspace(0.0001,smax,ns)
   	
	! Construct Household's state space
	  
	YYT = 0	
	yt_bar  =    exp(-1.06692)
	YYT = exp(YYT)*yt_bar
 	blower = 0
	bupper = 2.75

	BB(1) = blower;
	do i = 2,nd
    	B(i)=B(i-1)+(bupper-B(i-1))/((nd-i+1)**1.05); 
	end do  
	
	!! Choosing Household's Equilibrium
	eqs   = trim('c')

	select case(learning_type)
	case('RE')
	
		call oilcompany(M_spr(:,:1),M_Xr(:,:,1),M_RREVXr(:,:,1),M_PROFIT_Xr(:,:,1),M_EE(:,:,1), &
		 	 			Pp,p,s,smin,smax,disc,kappa,q,ns,np)
		
!      	M_spr(:,:,1)       = spr ! ok
!     	M_Xr(:,:,1)        = Xr
!      	M_RREVXr(:,:,1)    = RREVXr
!      	M_PROFIT_Xr(:,:,1) = PROFIT_Xr
    
		if(hh_on==1)then
			call household(pnn,cT,c,bp,bpy,gdp,  bp,b,yt,cT,c,pn,B,yN,aT,r,bpmax,revx,Pp, &
						cbind,sigg,omega,gamma,kapa,blower,bupper,cn,betta,np,ns,nd)
			
		    M_p(:,:,1)       = pnn
  	      	M_cT(:,:,1)      = cT
        	M_c(:,:,1)       = c
        	M_bp(:,:,1)       = bp
        	M_bpy(:,:,1)       = bpy
        	M_gdp(:,:,1)       = gdp
      					
      	end if
    case('BL')
    	right_name = 1
!    !$omp parallel do 
    do time = 1,T
    	write(*,*) 'T', time
    	Pp = M_pxstar(:,:,time)
		call oilcompany(M_spr(:,:time),M_Xr(:,:,time),M_RREVXr(:,:,time),M_PROFIT_Xr(:,:,time),M_EE(:,:,time), &
		 	 			Pp,p,s,smin,smax,disc,kappa,q,ns,np)

!		M_spr(:,:,time)       = spr 
!		M_Xr(:,:,time)        = Xr
!		M_RREVXr(:,:,time)    = RREVXr
!		M_PROFIT_Xr(:,:,time) = PROFIT_Xr
      	
      	if(hh_on==1)then
    		call household(pnn,cT,c,bp,bpy,gdp,  bp,b,yt,cT,c,pn,B,yN,aT,r,bpmax,revx,Pp, &
						cbind,sigg,omega,gamma,kapa,blower,bupper,cn,betta,np,ns,nd)
		    M_p(:,:,time)       = pnn
  	      	M_cT(:,:,time)      = cT
        	M_c(:,:,time)       = c
        	M_bp(:,:,time)      = bp
        	M_bpy(:,:,time)     = bpy
        	M_gdp(:,:,time)     = gdp
      	end if
    end do 
!	!$omp end parallel do 
	end select
	
    write(*,*) 'Model :', learning_type
	place = '/outData/'//learning_type
	
	!**********************************************
	! Saving options 

	if(1==1)then
	call print_binary_3(trim(place)//'M_p',M_p)
	call print_binary_3(trim(place)//'M_cT',M_cT)
	call print_binary_3(trim(place)//'M_c',M_c)
	call print_binary_3(trim(place)//'M_bp',M_dbp)
	call print_binary_3(trim(place)//'M_bpy',M_dbpy)
	call print_binary_3(trim(place)//'M_gdp',M_gdp)
	

	call print_binary_3(trim(place)//'M_spr',M_spr)
	call print_binary_3(trim(place)//'M_Xr',M_Xr)
	call print_binary_3(trim(place)//'M_RREVXr',M_RREVXr)
	call print_binary_3(trim(place)//'M_PROFIT_Xr',M_PROFIT_Xr)
	end if
	
end program 