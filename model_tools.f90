module model_tools
use basic 
use tools

implicit none 
contains 

subroutine oilcompany(spr,Xr,REVXr,PROFIT_Xr,EE & 			! output
					Pxprob,P,SS,smin,smax,d,kappa,q,np,ns)  	! input

	!! OILCOMPANY  Solves optimal oil extraction problem by private firm 

	! Dummy 
	integer, intent(in) :: np, ns
	real, intent(in) 	:: Pxprob(:,:),P(:),SS(:),smin,smax,q,d,kappa,q
	
	real, intent(out), dimension(np,ns)	:: spr,Xr,REVXr,PROFIT_Xr,EE 
		
	! Local	
	integer :: i, j, k, l
		
	! Construct the firm's state space

	s  = repmat(S,1,np);
	sp = repmat(S,1,np);
	px = repmat(PX,ns,1);

	! Time Iteration Loop 
	!%%%%%%%%%%%%%%  Technical parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iter            = 0;            ! 
	d2              = 100;          !
	uptd_of         = 0.1;          ! Weight on new policy function in update to next iteration
	iter_tol_of     = 3000;         !
	tol_of          = 1e-6;         !
	outfreq_of      = 10;           ! Display frequency (shows in screen each 10th iteration)

	disp('DE Iter      Norm');

	!%%%%%%%%%%%%%%%%%%%%%%%%%% Start of iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	do while(d2>tol_of .and. iter < iter_tol_of)
    
    ! Retrieve updated policy functions
    oldsp=sp;

    ! Marginal profits
    mr = px - (2*kappa*(s+d-sp)*s - kappa.*(s+d-sp)**2 )/(s**2);
       
    ! Interpolation 
    do i=1,ns
        do j=1,np    
        	!************** Piecewise interpolation **************  
			do k = 1,np
				inter(k) =  piecewise(SS,mr(:,k),sp(i,j))
			end do	
			emr(i,j) = beta*sum(inter*Pxprob(j,:))		! Expected marginal profits
			inter=0.
			!*****************************************************

            x(i,j) = (px(i,j) - emr(i,j))*s(i,j)/(2*kappa)  
            x(i,j) = maxval(0.0001,x(i,j))
               
            sp(i,j) = maxval(s(i,j) - x(i,j) + d,smin)
            EE(i,j) = x(i,j) - (px(i,j) - emr(i,j))*s(i,j)/(2*kappa)
        end do 
    end do
    
    iter=iter+1 ! Update iteration counter
    
    ! Calculate difference between new and old policies
    d2 = max(max(abs(sp-oldsp)));
    
    ! Print results once every (outfreq) iterations
    if mod(iter, outfreq_of) == 0
    	fprintf('%d          %1.7f \n',iter,d2);
    end
    
    ! =====================Updating rule for next iteration=================
    sp = uptd_of*sp + (1-uptd_of)*oldsp;
    ! ======================================================================
    end do 
 
	write(*,*) iter, d2    ! Print last iteration result

	! Other equations of the firm 
	revx = px*x
	profit = px*x - kappa*(x^2)/s
end subroutine 

subroutine household(pnn,cT,c,bp,bpy,gdp,  bp,b,yt,cT,c,pn,B,yN,aT,r,bpmax,revx,Pp, &
						cbind,sigg,omega,gamma,kapa,blower,bupper,cn,betta,np,ns,nd)

	! Dummy	
	
	!Local
	integer :: i,j,k,l,t

	thereshold_ct = 1e-4;  ! Minimum value of consumption

	! Time Iteration Loop 

	! %%%%%%%%%%%%%%%%%%%%%%  Technical parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	uptd     = 0.9;     ! Weight on new policy function in update to next iteration
	outfreq  = 5;       ! Display frequency (shows in screen each 20th iteration)
	tol      = 3E-3;    ! Numerical tolerance (convergence criterion) 

	iter     = 0;
	d2       = 100;    

	write(*,*) 'DE Iter      Norm'
	! %%%%%%%%%%%%%%%%%%%%%%%%%% Start of iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	while( d2>tol )
    
    	! Retrieve updated policy functions
    	oldct = ct
        
    	! Marginal utility of consumption given current consumption policy 
    	mu  = real(c**-sigma*gamma**(1/omega)**(ct/c)**(-1/omega));  
        
    	! Interpolation  of expected marginal utility
    	! We interpolate the expected marginal utility for each level of
    	! reserves (k)
    	do k = 1,ns
        	mu_aux = mu((k-1)*nd+1:k*nd,:) 
        	do i = 1,nd
            	do j = 1,Ntotal 
                	!************** Piecewise interpolation **************  
					do l = 1,np
						inter(l) =  piecewise(BB,mu_aux(:,l),bp(i,j))
					end do	
					emup_aux(i,j) = sum(inter*Prob(j,:))		! Expected marginal profits
					inter=0.
					!*****************************************************
            	end do
        	end do
    	end do
    
    	do t = 1,ns
    	   emup((t-1)*nd+1:t*nd,:) = emup_aux(:,:,t)
    	end do 
   	  
    	!----------------- Find new constrained values -----------------------
    	gdp = (yt+pn*yn+revx)
    	bpmax = kappa*gdp 
    	where (bpmax>bmax) bpmax = bmax
    	where (bpmax<bmin) bpmax = bmin
    	ctbind = -(1+r)*b + revx + yt - at + bpmax
    	where (ctbind<thereshold_ct) ctbind = thereshold_ct
    	cbind  = (gamma**(1/omega)*ctbind**((omega-1)/omega)+(1-gamma)**(1/omega)*cn**((omega-1)/omega))**(omega/(omega-1))
    	!----------------------------------------------------------------------
    
    	do i=1:nd*ns
        do j=1:Ntotal
       
            ! Calculate Euler Equation Error at maximum feasible consumption cbind
            EE(i,j) = cbind(i,j)**-sigma*gamma**(1/omega)*(ctbind(i,j)/cbind(i,j))**(-1/omega)-(1+r)*beta*emup(i,j);
        
            if(EE(i,j)>tol*1e-5)then ! If positive, constraint will be binding then:
               
                ! Debt will be as big as possible                
                bp(i,j) = bpmax(i,j);      

                ! Consumption is solved from budget constraint
                c(i,j)  = cbind(i,j);                            
                ct(i,j) = ctbind(i,j);
                pn(i,j) = real(((1-gamma)/gamma*(ct(i,j)/cn))**(1/omega))
                gdp(i,j) = yt(i,j)+ pn(i,j)*yn + revx(i,j)
                bpy(i,j) = bp(i,j)/gdp(i,j);
               
            else ! Constraint not binding
            
            	! Define function that calculates the absolute value of Euler
                ! Equation Error for a given consumption
                
             	!********************* Bisection Algorithm ***********************
    			cl = thereshold_ct
				ch = 1.        
				error_EE = 1. 
	             
				cm = (cl+ch)/2
				if ( (gamma**(1/omega)*cl**((omega-1)/omega)+(1-gamma)**(1/omega)* &
					cn**((omega-1)/omega))**(omega/(omega-1)))**(-sigma+(1/omega))* &
					gamma**(1/omega)*(cl)**(-1/omega))-(1+r)*beta*emup(i,j)   * 
					 (gamma**(1/omega)*ch**((omega-1)/omega)+(1-gamma)**(1/omega)* &
					cn**((omega-1)/omega))**(omega/(omega-1)))**(-sigma+(1/omega))* &
					gamma**(1/omega)*(ch)**(-1/omega))-(1+r)*beta*emup(i,j)) > 0 ) then 
					write (*,*) 'Warning: Bisection can not find a solution in EE'
				end if
				do while (abs(error_EE)>delta) 
					cm = (cl+ch)/2
					if ( (gamma**(1/omega)*cl**((omega-1)/omega)+(1-gamma)**(1/omega)* &
						cn**((omega-1)/omega))**(omega/(omega-1)))**(-sigma+(1/omega))* &
						gamma**(1/omega)*(cl)**(-1/omega))-(1+r)*beta*emup(i,j)   * 
					 	(gamma**(1/omega)*ch**((omega-1)/omega)+(1-gamma)**(1/omega)* &
						cn**((omega-1)/omega))**(omega/(omega-1)))**(-sigma+(1/omega))* &
						gamma**(1/omega)*(ch)**(-1/omega))-(1+r)*beta*emup(i,j)) > 0 ) then
						ch = cm
					else 
						cl = cm
					end if
					
					error_EE = (gamma**(1/omega)*cm**((omega-1)/omega)+(1-gamma)**(1/omega)* &
						cn**((omega-1)/omega))**(omega/(omega-1)))**(-sigma+(1/omega))* &
						gamma**(1/omega)*(cm)**(-1/omega))-(1+r)*beta*emup(i,j)
				end do

				ct(i,j)  = cm
                EE(i,j) = error_EE
               	!********************* Bisection Algorithm ***********************
                       
                ct(i,j) = maxval(ct(i,j),thereshold_ct)
                
                c(i,j) = ((gamma**(1/omega)*ct(i,j)**((omega-1)/omega)+(1-gamma)**(1/omega)*cn**((omega-1)/omega))**(omega/(omega-1)));
    			!!!! get real part from a complex number 
                pn(i,j) = real(((1-gamma)/gamma*(ct(i,j)/cn))**(1/omega))
                gdp(i,j) = yt(i,j)+ pn(i,j)*yn + revx(i,j)
                 
                ! Solve debt from budget constraint, check if it is within grid bounds

                bp(i,j) = maxval(-revx(i,j) - yt(i,j)+(1+r)*b(i,j)+ ct(i,j)+ at, bmin)
                bp(i,j) = min(bp(i,j),bmax)
                bpy(i,j) = bp(i,j)/gdp(i,j)
            end if
        end do
    end do
  	iter=iter+1 ! Update iteration counter
   
  	! Calculate difference between new and old policies
    
  	d2 = maxval(maxval(abs(ct-oldct),1),1)

  	! Print results once every (outfreq) iterations
  	if (mod(iter, outfreq) == 0)then
  		write (*,*) iter, '     ', d2
  	end if 
    
    !=====================Updating rules for next iteration================
    ct = uptd*ct+(1-uptd)*oldct
    pn = real(((1-gamma)/gamma*(ct/cn))**(1/omega))
    c = (gamma**(1/omega)*ct**((omega-1)/omega)+(1-gamma)**(1/omega)*cn**((omega-1)/omega))**(omega/(omega-1))
    where (-revx - yt+(1+r)*b + ct+ at<bmin) bp = bmin 
    where (bp>bmax) bp = bmax
    !======================================================================
	end do
  	write (*,*) iter, '     ', d2
end subroutine 

subroutine ergodic(TT,Tburn,LE,M_Pxstar,M_REVXr,M_Xr,M_spr,M_p,M_cT,M_c,M_dpix,M_dp,Pp,yN,&
				dgrid,sgrid,YYT,YYTprob,RREVX,PX,revXgridlevel,yTgridlevel,hh_on,nd,ny,ns,np,NNy,learning_type_for)

	!Dummy
	character(*) :: learning_type_for
	integer, intent(in) :: nd, ny, ns, np, NNy
	integer, intent(in) :: TT, Tburn, LE, hh_on
	real, dimension(:,:,:), intent(in) 	:: M_Pxstar,M_REVXr,M_Xr,M_spr,M_p,M_cT,M_c,M_dpix,M_dp
	real, dimension(:,:), 	intent(in) 	:: Pp, YYTprob
	real, dimension(:), 	intent(in) 	:: sgrid, dgrid, YYT, RREVX, PX,revXgridlevel, yTgridlevel!,pXgridlevel
	real, intent(in) :: yN
	
	!Local
	integer :: i, time, h, g, t, l, m, n, u, j, dpix1, LH, LL, z, k

	character(100) :: mod_modif, kof, yT_nature, fileplace

	real, dimension(ny,nd) :: tau, crisis, mu 
	real :: v_sX0, sX0, grid(NNy*ny,2), CPp(np,np)
	real :: revX0, yT0, v_d0, d0, pX0
	real :: CYYTprob(NNy,NNy)
	real :: pXgridlevel(np), ran_num
	
	real, dimension(TT,LE) :: PX_SIM, SX_SIM, SXp_SIM, X_SIM, REVX_SIM
	real, dimension(TT,LE) :: LA_SIM, YT_SIM, D_SIM, Dp_SIM, V_SIM, PN_SIM
	real, dimension(TT,LE) :: CT_SIM, C_SIM, SLACK_SIM, MU_SIM, Y_SIM, CAY_SIM
	real, dimension(TT,LE) :: CA_SIM, DY_SIM, TB_SIM, TBY_SIM

	real, dimension(TT,LE) :: HH, II, JJ, UU
	real :: MPX, SPX, MSX, MX, MD, MCT, MPN, MCA, MY, MDY, MCAY, MPOR, MREVX  
	
	
	mod_modif           = trim('')
	kof                 = trim('Type_2')
!	learning_type_for   = trim('RE')
	yT_nature           = trim('YT_FIXED')
 
	! read policy functions

	fileplace = 'outData/ergodic/'//learning_type_for
		   
	! Initial Conditions
	CPp = cumsum(Pp)
 
	v_sX0 =  1.0823
	h = minloc(abs(sgrid-v_sX0),1)
	sX0 = sgrid(h)
 
	if(hh_on==1)then 
		revX0 = revXgridlevel(1)
		yT0   = yTgridlevel(1)
		v_d0 =  0.6440     				!	Initial debt
		g = minloc(abs(dgrid-v_d0),1)
		d0 = dgrid(g)
	end if
	
	! Initializations' Firm
	PX_SIM      = 0.
	SX_SIM      = 0.
	SXp_SIM     = 0.
	X_SIM       = 0.
	REVX_SIM    = 0.
 
	! Initializations' Household
	LA_SIM      = 0. 	!marginal utility of cT
	YT_SIM      = 0.  	!yT
	D_SIM       = 0. 	!current debt
	Dp_SIM      = 0.  	!next-period debt
	V_SIM       = 0. 	!value function
	PN_SIM      = 0. 	!relative price nontradables
	CT_SIM     = 0. 	!consumption of tradables
	C_SIM      = 0.   	!composite consumption
	SLACK_SIM  = 0.  	!slackness of collateral constraint
	MU_SIM     = 0.  	!collateral constraint multiplier
	Y_SIM      = 0.
	CAY_SIM    = 0.
	CA_SIM     = 0.
	DY_SIM     = 0.
	TB_SIM     = 0.
	TBY_SIM    = 0.
 
	pX0 = PX(size(PX)) 

 	pXgridlevel = unique(PX)

	CYYTprob = cumsum(YYTprob)		
		
	do l = 1,LE
	do t = 1,TT
		! Firm
		PX_SIM(t,l)    = pX0
		SX_SIM(t,l)    = sX0
		m        	      = minloc(abs(pXgridlevel-pX0),1)
		n           	  = minloc(abs(sgrid-sX0),1)
		REVX_SIM(t,l)  = M_REVXr(m,n,l)
		X_SIM(t,l)     = M_Xr(m,n,l)
		SXp_SIM(t,l)   = M_spr(m,n,l)
		revX0             = M_REVXr(m,n,l)
		sX0           	  = M_spr(m,n,l)
	
		! Household
		
		if(hh_on==1)then
			YT_SIM(t,l)    = yT0
			D_SIM(t,l)     = d0
			h              = minloc(abs(revX0-revXgridlevel),1)
			i              = minloc(abs(yT0-unique(yTgridlevel)),1)				! DIFERENCIA
			do z = 1,nd
!				write(*,*) dgrid(z), d0, z
				if(dgrid(z)==d0)then
					j = z
!					!z=nd
				end if 
			end do 
			
!		write(*,*) 'l', l, 't', t, 'm ', m, 'n ', n, 'j', j
!		if(l == 1 .and. t==1000)then
!			stop
!		end if 
		
			HH(t,l)         = h
			II(t,l)         = i
			JJ(t,l)         = j
			u               = NNy*(h-1) + i
			UU(t,l)         = u
		
			PN_SIM(t,l)        = M_p(u,j,l);
			CT_SIM(t,l)        = M_cT(u,j,l);
			C_SIM(t,l)         = M_c(u,j,l);
			dpix1              = M_dpix(u,j,l);
			Dp_SIM(t,l)        = dgrid(int(dpix1))
 	
!		WRITE(*,*) 'yT0', yT0, 'd0', d0, 'h',h,'i',i,'j',j, 'u',u,'l',l,'PN_SIM', PN_SIM(t,l)
!		do i =1,340
!		write(*,*) i, M_p(i,j,l)
!		end do
		
		
			d0                 = M_dp(u,j,l);
		
		 	call random_number(ran_num)
		 	LH=1
 			do j = 1,size(CYYTprob,2)	
				if(CYYTprob(i,j) < ran_num) then
   					LH = LH+1
   				end if 
			enddo
			yT0 = yTgridlevel(LH);
		end if 
	
	 	call random_number(ran_num)
	 	LL=1
		do j = 1,size(CPp,2)	
			if(CPp(i,j) < ran_num) then
   				LL = LL+1
   			endif
   		enddo
		pX0 = pXgridlevel(LL);	
	end do
		! Building Offcore Variables
		if(hh_on==1)then
   
		! Construct the current account
			do k = 1,LE
				CA_SIM = 0.
				CA_SIM(2:,k) = -(D_SIM(2:,k)-D_SIM(:size(D_SIM(:,k),1)-1,k))
			end do 		
			
			! Produce simulations of other   variables of the model
			TB_SIM(:,l)  = YT_SIM(:,1) + REVX_SIM(:,l) - CT_SIM(:,l)		! trade balance
			Y_SIM(:,l)   = YT_SIM(:,1) + yN*PN_SIM(:,l) + REVX_SIM(:,l)		! output measured in terms of tradables
 
			DY_SIM(:,l)  = Dp_SIM(:,l)/Y_SIM(:,l)							! debt-to-output ratio (annual)
			TBY_SIM(:,l) = TB_SIM(:,l)/Y_SIM(:,l)							! trade balance to output ratio
			if(LE==1)then 
				CAY_SIM(:,1) = CA_SIM(:,1) /Y_SIM(:,1)  					! trade balance to output ratio
			else 
   			 	CAY_SIM(:,l) = 0.*Y_SIM(:,l) 								! We do not need this for BL exercise
			end if 
		end if 
	end do
	
	
	! Computing Ergodic Means
	 
	MPX   = sum(mean(PX_SIM(Tburn:,:)))/Tburn
!	SPX   = std(PX_SIM(Tburn:))
  	MSX   = sum(mean(SX_SIM(Tburn:,:)))/Tburn
	MX    = sum(mean(X_SIM(Tburn:,:)))/Tburn
	MPOR  = MSX/MX
	MD    = sum(mean(D_SIM(Tburn:,:)))/Tburn
	MCT   = sum(mean(CT_SIM(Tburn:,:)))/Tburn
	MPN   = sum(mean(PN_SIM(Tburn:,:)))/Tburn
	HH    = sum(mean(CA_SIM(Tburn:,:)))/Tburn
	MY    = sum(mean(Y_SIM(Tburn:,:)))/Tburn
	MDY   = sum(mean(DY_SIM(Tburn:,:)))/Tburn
	MCAY  = sum(mean(CAY_SIM(Tburn:,:)))/Tburn
	MREVX = (MPX*MX)
  
  write(*,*) 'Rational Expectation    -  # of Nodes: ',nd    
  write(*,*) 
  write(*,*) '   Oil price      ' 
  write(*,*) '     Exp Value            ' , MPX
  write(*,*) '     Std Dev                ' , SPX
  write(*,*) '     qll                  	' , Pp(1,1)
  write(*,*) '     qhh                  	' , Pp(2,2)
! write(*,*) '     Low state persist.   '  ,Erhol(t),q(1)-(1-q(1)) 
! write(*,*) '     High state persist.  '  ,Erhoh(t),q(1)-(1-q(1))
  write(*,*) '   Firm                   '  
  write(*,*) '     Stock                ' , MSX
  write(*,*) '     Extraction             ' , MX
  write(*,*) '     Periods of Reserves (Q) ' , MPOR
  write(*,*) '     Oil  Revenues           ' , MREVX
  write(*,*) '   Household                 '  
  write(*,*) '     Debt                    ' , MD
  write(*,*) '     T Consumption           ' , MCT
  write(*,*) '     NT Price                ' , MPN
  write(*,*) '     Current Account         ' , MCA
  write(*,*) '     GDP                     ' , MY
  write(*,*) '     Debt to GDP             ' , MDY
  write(*,*) '     Current Account to GDP  ' , MCAY
  
  ! data
!  call print_binary_2('outData/PX_SIM',PX_SIM)
!  call print_binary_1('outData/M_p', M_p(:,j,l))

end subroutine 
end module 