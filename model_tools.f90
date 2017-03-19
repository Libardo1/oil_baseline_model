module hamann_tools
use basic 
use tools
use solve
use def_types

implicit none 
contains 

subroutine oilcompany_franz_v3(PROFIT_X,REVX,pstar,sp,Xr,PROFIT_Xr,REVXr,sgrid,spr,  & 	! output
								q,kappa,disc,s,Pp,SS,PX,n,ns,np)						! input
 
	!! OILCOMPANY  Solves optimal oil extraction problem by private firm 
	! Unlike oilcompany_franz.m in this code we need to rewrite S,PX and Pr in order to have s changing first and then prices. We need to do
	! in this way in order to implement this outcomes in Time Iteration codes.

	! Dummy 
	integer, intent(in) :: n
	real, intent(in) 	:: q, kappa, disc, s(:), SS(:)  
	
	real, intent(out), dimension(np,ns)	:: spr, REVXr, PROFIT_Xr, Xr
	real, intent(out)	:: sgrid(ns), sp(n), REVX(n), PROFIT_X(n)
	real, allocatable, intent(out) :: pstar(:,:) 	
	
	integer, intent(inout):: ns, np
	real, intent(inout)	:: PX(n), Pp(np,np)
	
	! Local	
	integer :: i, j, k, l
	real	:: x(n,ns), pie(n,ns)
	
	integer, allocatable :: indx(:), jndx(:), val(:)
	integer :: element, rep
	real :: EST(ns*n,ns), EST_FINAL(ns*n,ns*2), Pr_aux(ns*np*ns,ns*np)
	real :: temp(ns*n,n)
	real, allocatable :: v(:), x_p(:)
	type(sparse_matrix) :: Pr
	
	real :: XX(n)
	
	real :: f(n), Eo, Ex, Erevx, profit, value
	
	real ::	vr(np,ns), fr(ns,np), FFr(ns,np)

	x = 0.
	pie = 0.
	!! Profit function and feasible extraction x<=o
	
	do i = 1,ns*np
    	do j = 1,ns
        	x(i,j) = (SS(i)-s(j)+disc)        	
        	if(x(i,j)>=0 .and. x(i,j)<=SS(i)+disc) then        
          		pie(i,j) = PX(i)*x(i,j)-kappa*((x(i,j))**2/(SS(i)))
        	else
          		pie(i,j) = -99999999999.
        	end if 
        	!write(*,*) 'i', i , 'j', j ,x(i,j),'>=',0, x(i,j)>=0 , '----',x(i,j),'<=',SS(i)+disc,  x(i,j)<=SS(i)+disc
    	end do
!    	if (i==175)then
    		!stop
!    	end if
	end do 
	
	
!	call print_matrix(pie(:,1:10))
	!! Transition probability matrix (see Sargent and Ljundqvist) 
!	write(*,*) 'sparse part'
	! we include the functionalities of Sparse BLAS Library to handle the following lines   
	! create the VAL, INDX, JINDX entries 
	
	! EST = kron(eye(ns),ones(n,1))
	! EST_FINAL = repmat(EST,1,2)
	rep = 2
	allocate(indx(ns*n*rep)) 
	allocate(jndx(ns*n*rep))
	allocate(val(ns*n*rep))
	
    allocate(Pr%values(ns*n*rep))
    allocate(Pr%d1(ns*n*rep))
    allocate(Pr%d2(ns*n*rep))

	element = 1
!	val = 1.	
	do k = 1,rep
  		do i = 1,ns
  			do j = 1,n
  			 	Pr%values(element) = 1.
				Pr%d1(element) = (i-1)*n + j
				Pr%d2(element) = i + (k-1)*ns
				!write(*,*) Pr%values(i), Pr%d1(element), Pr%d2(element), element
				element = element + 1
  			end do
  		end do 
  	end do 
  
	!	Pr_aux = EST_FINAL*(repmat(kronecker(Pp,ones(ns,ns)),ns,1))
	temp  = repmat(kronecker(Pp,ones(ns,ns)),ns,1)
!	write(*,*)shape(temp), shape(repmat(kronecker(Pp,ones(ns,ns)),ns,1))

!	call print_matrix(temp(:,1:10))
!	stop
	
	do i = 1,ns*n*rep
		Pr%values(i) = temp(Pr%d1(i),Pr%d2(i))
!		write(*,*) Pr%values(i), temp(Pr%d1(i),Pr%d2(i))!, Pr%d1(i), Pr%d2(i)
	end do 
!	stop
	!! Solve Bellman equation for oil-firm

!	call print_matrix(pie(:,1:10))

	allocate(v(n))
	allocate(x_p(n))
	allocate(pstar(n,n))

!	write(*,*) 'Solvedp'	
	call solvedp(v,x_p,pstar, pie,Pr,q,'policy',0.)
!	call print_matrix_dat('outData/oilcompany/x.txt',reshape(x_p,(/size(x_p),1/)))
!	call print_matrix_dat('outData/oilcompany/pstar.txt',pstar)

!	write(*,*) pstar(169,:)
!	stop
	XX = SS + disc - s(int(x_p))
	where(XX<1E-6) XX =0.
!	do i =1,ns
!		write(*,*)	XX(i),  SS(i),  disc, '-',s(int(x_p(i))+1)
!	enddo
	REVX = PX*XX
	
	PROFIT_X = PX*XX - (kappa/2)*(XX**2)/SS
	f = reshape(ergdist(pstar,1./size(pstar,1)*ones(n,1)),(/n/))
	!! Compute optimal steady-state oil reserves desired by oil company
!	call print_matrix_dat('outData/oilcompany/f.txt',reshape(f,(/size(f),1/)))
 
	Eo  = dot(f,s(int(x_p)))
	Ex  = dot(f,XX)
	Erevx  = dot(f,REVX)
	profit = dot(f,(PX*XX-kappa*(XX**2)/SS)) 
	value  = dot(f,PX*SS)
		
	!! Plots
	! Reshaping
	vr = transpose(reshape(v,(/ns,np/)))
	Xr = transpose(reshape(XX,(/ns,np/)))
	fr = reshape(f,(/ns,np/))
	FFr = transpose(cumsum(transpose(fr)))
	REVXr = transpose(reshape(REVX,(/ns,np/)))
	PROFIT_Xr = transpose(reshape(PROFIT_X,(/ns,np/)))
	
	sp = s(int(x_p))
	spr = transpose(reshape(sp,(/ns,np/)))
	sgrid = s	
!	do i = 1,n
!		write(*,*) i, XX(i), ' = ', SS(i), '+', disc, '-', s(int(x_p(i))), '/////' , int(x_p(i))
!		write(*,*) i, f(i)
!	end do 
!	write(*,*) spr(1,:)
!	stop

	write(*,*)	
	write(*,*)	'Steady State Means'
	write(*,*)	
	write(*,*)	'	Stock					=', Eo 
	write(*,*)	'	Value of Stock				=', value
	write(*,*)	' 	Extraction				=', Ex
	write(*,*)	'	Reserves (in periods)			=', Eo/Ex
	write(*,*)	' 	Profits					=', profit
	write(*,*)	'	Revenues				=', Erevx
	write(*,*)
	
	!*************************************************
	! PLOTS
	!*************************************************
	! We export data to plot into a python notebook in this same folder
	!
	
	call print_matrix_dat('outData/oilcompany/s.txt',reshape(s,(/size(s),1/)))
	call print_matrix_dat('outData/oilcompany/Xr.txt',Xr)
	call print_matrix_dat('outData/oilcompany/vr.txt',vr)
	call print_matrix_dat('outData/oilcompany/PROFIT_Xr.txt',PROFIT_Xr)
	call print_matrix_dat('outData/oilcompany/REVXr.txt',REVXr)
	call print_matrix_dat('outData/oilcompany/fr.txt',fr)
	call print_matrix_dat('outData/oilcompany/Fr.txt',Fr)

end subroutine 

subroutine constrained_V2F(warnn,r,pai,dgrid,ygrid,yTgridlevel,revX,revXgridlevel, 	&	! output
						   la,dpix,dp,yT,cT,c,distTol,dist_dpix,d,p,  				&	! output 
							YYT,YTprob,RREVX,Pxstar,rstar,betta,a,yN,xi,sigg,		&	! input
							aT,aN,kapa,nd,ny,dlower,dupper,eqs,solution)				! input
 
	!Use a policy-function iteration procedure to approximate the equilibrium of  a small open economy with a flow
	!collateral constraint as presented in
	!''Multiple Equilibria in Open Economy Models with Collateral Constraints: Overborrowing Revisited'' 
	! by Stephanie Schmitt-Grohé and Martín Uribe.
	!The algorithm accounts for  the possibility that for
	!some current states there is no choice of next-period debt
	!for which consumption of tradables is positive and the collateral constraint is satisfied. This is done by calling pathfinder.m
	!©  Stephanie Schmitt-Grohé and Martín Uribe, May 2016 for matlab

	!Choose an equilibrium selection criterion
	!'b' favors the least severe of the two constrained outcomes (point B in the graphs shown in the references cited in the preamble of this program)
	!'c' favor the most severe of the two constrained outcomes (point C in the graphs)
	!'a' favor unconstrained outcomes (point A in the graphs). [We leave the possibility of this equilibrium as a conjecture that the present model delivers an equilibrium displaying overborrowing. We  encourage readers to work on this conjecture. At the present stage, the equilibrium under criterion (a) delivers a solution (i.e., artificial time series generated by simu.m, see simu_constrained_a.mat)  satisfying all equilibrium conditions except the nonnegativity of consumption (i.e., min(CT) fails to be positive).
	
	! Dummy	
	integer, intent(out):: warnn
	real, intent(out)	:: r,pai(ny,ny),dgrid(nd),ygrid(ny),yTgridlevel(ny)
	real, intent(out)	:: revX(ny,nd),revXgridlevel(ny)
	
	integer, intent(out):: dpix(ny,nd), dist_dpix
	real, intent(out), dimension(ny,nd)	:: la,dp,yT,cT,c,d,p
	real, intent(out) 	:: distTol	
							
	integer, intent(in)	:: nd,ny
	real, intent(in)	:: YYT(:),YTprob(:,:),RREVX(ny),Pxstar(ny,ny),rstar,betta,a,yN,xi,sigg,aT,aN,kapa,dlower,dupper
	character(*), intent(in) ::	 eqs, solution 
	
	!Local
	integer :: i, j, k, cont, dim(2), xR
	real, dimension(ny*nd,nd) :: cTtry, ctry, ptry, collateraltry, slacktry, Ztry, latry, LHSEulertry, AMUtry, MUtry, temp3
	real, dimension(ny*nd,nd) :: dytry, temptry
	real :: dd, revxgrid(ny), start, finish, temp
	real :: temp1(ny,nd)
	integer,dimension(ny,nd) :: dix, dpixnew, dpix_init
	real :: dY(ny,nd)
	integer :: dpix_ind(ny*nd,2), n
	logical :: selec = .true.

	real :: NDL
	
	real :: cN
	real, dimension(ny*nd,nd) :: yy, dtry 
	
	integer, dimension(ny*nd):: bind1,bind2
	real :: temp6(ny,nd)
	
	real :: cumtime(10), aa(50)

	character(1) :: eqm_selection_criterion 

	dist_dpix = 100

!	write(*,*) 'Running Constrained'
	eqm_selection_criterion = eqs
	!load tpm.mat pai ygrid rgrid; %produced by running tpm.m
	revxgrid = reshape(kronecker(reshape(RREVX,(/size(RREVX),1/)),ones(size(YYT),1)),(/ny/))
	ygrid = reshape(repmat(reshape(YYT,(/1,size(YYT)/)),size(RREVX),1),(/ny/))
	pai = kronecker(Pxstar,YTprob)
	
!	write(*,*) shape(Pxstar)
!	call print_matrix(Pxstar(:,1:10))
!	stop
	yTgridlevel = ygrid
	revXgridlevel  = revxgrid 

	NDL = minval(yTgridlevel+revXgridlevel,1)*(1+rstar)/(rstar)	!Natural debt limit 

	! Lower and upper bounds for debt grid

	! debt grid
	dgrid = linspace(dlower,dupper,nd)
	dd = dgrid(2) - dgrid(1)

	n = ny*nd
	! Let yT, r, and d  be 3 ny-by-nd matrices containing values  of the current state
	yT = transpose(spread(yTgridlevel,1,nd))
	revX  = transpose(spread(revXgridlevel,1,nd))
	r = rstar
	d  = spread(dgrid,1,ny)		

	dix = spread((/(i,i=1,nd)/),1,ny)	! ny-by-nd matrix of indices of the current debt state
!	write(*,*) shape(dix)
!	call print_matrix(real(dix(:,:)))
	! n-by-nd matrix of tradable consumption as a function of the current state (n in total) and  possible choices of next-period debt (nd in total)

	cTtry = repmat(reshape(revX,(/n,1/)),1,nd) + repmat(reshape(dgrid,(/1,nd/)),n,1) + &
			repmat(reshape(yt,(/n,1/)),1,nd) - repmat(reshape(d,(/n,1/)),1,nd)*(1+r) - aT*ones(n,nd)
	
	if (minval(maxval(cTtry,1),1)<=0)then 
    	write(0,*) 'Natural debt limit violated'
	end if 
	
	cN = yN -aN
	
	! n-by-nd matrix of relative price of nontradables as a function of the current state (n in total) and possible values for next-period debt (nd in total)

	ptry  = (1-a)/a*(cTtry/cN)**(1/xi)
	where(cTtry<0) ptry = 1E20
	
	yy = repmat(reshape(yt,(/n,1/)),1,nd) + repmat(reshape(revX,(/n,1/)),1,nd) + yN*ptry
	dytry = repmat(reshape(dgrid,(/1,nd/)),n,1)/yy
	
	collateraltry = kapa*(spread(reshape(yT,(/n/)),2,nd) + repmat(reshape(revX,(/n,1/)),1,nd) + yN*ptry)

	! n-by-nd matrix of the slackness of the collateral constraint  as a function of the current state (n in total) and possible values for next-period debt (nd in total).
	slacktry = collateraltry - spread(dgrid,1,n)

	! Gridize  values of  slacktry matrix for   which the collateral constraint is binding  by setting them equal to 0
			
	slacktry = bind(slacktry)
!	call print_matrix(slacktry(34000:,92:101))
!	stop
	
	!Identify the indices at which slack is null (binding collateral constraint).
	!If for the current state there exists exactly 1 value of next-period debt for which slack is 0, 
	!collect this index in bind1ix. If for the current state there exist  two binding points, the larger is collected in bind2ix,
	!and the other in bind1ix. 
	!The vectors (bind1i,bind1j) are the subindex representation of bind1ix, and similarly for (bind2i,bind2j).
	!Note: the reason why there can be exactly 1 binding point for a given state (as opposed to either 0 or 2)
	!is that we have already set to NaN all entries of slack associated with cT<=0 (which could have contained binding points);
	!thus  bind2i is a subset of bind1i.

	bind1 = 0
	bind2 = 0
		
	do i = 1,n
		if(minval(abs(slacktry(i,:)),1)<1e-5)then
			bind1(i) = minloc(abs(slacktry(i,:)),1)
		end if 
		if(minval(abs(pack(slacktry(i,:),mask(nd,bind1(i)))),1)<1e-5)then
			bind2(i) = minloc(abs(pack(slacktry(i,:),mask(nd,bind1(i)))),1)+1
		end if 
	enddo
	
	!n-by-nd matrix of consumption of the composite good  as a function of the current state (n in total) and possible values for next-period debt (nd in total)
	
	ctry = ((a**(1./xi)) * cTtry**(1-1/xi) + ((1-a)**(1./xi)) * cN**(1-1/xi))**(1/(1-1/xi))
	
	!n-by-nd matrix of marginal utility of tradable consumption  as a function of the current state (n in total) and possible values for next-period debt (nd in total)
	
	latry = (a**(1./xi))*ctry**(-sigg) * (cTtry/ctry)**(-1/xi)
	where (cTtry<=0) latry=1e6
	
	!n-by-nd matrix of lambda_t/(1+r_t)  as a function of the current state (n in total) and possible values for next-period debt (nd in total)
	
	LHSEulertry = latry!/spread(1+reshape(r,(/n/)),2,nd)

	!initializations
	if(solution .eq. 'BL')then 
		OPEN (unit = 1, file = 'outData/init_BL.txt', status = 'old')
		READ (1,*) temp1
		CLOSE (1)
		dpix = int(temp1)
	else if(solution .eq.'RE')then
		dpix = dix
	end if 
	

!	call cpu_time(finish)
!	write(*,*) 'time ', finish-start
!	call cpu_time(start)
	Ztry = pathfinder(dlower,dupper,slacktry, cTtry, pai, nd,ny) !This step rules out picking values for d_{t+1} that put the economy with positive probability into a state tomorrow in which either consumption is negative or the collateral constraint is violated.
!	call cpu_time(finish)
!	write(*,*) 'time ', finish-start
	cont = 0
	dpix_init = dpix	
	
!	call print_matrix(Ztry(34000:,92:101))
!	stop
	! Euler Iteration 
	
	if(selec)then 
		write(*,*) 'Equilibrium Selection : ', eqm_selection_criterion   
	else 
		write(*,*) 'Equilibrium Selection : None'   
	end if 
	write(*,*) '    dist_dpix', '        xR'
	write(*,*) 
	
!	call cpu_time(start)
	do i = 1,n
		dim = ind2sub(shape(dpix),i)
		dpix_ind(i,1) = dim(1)
		dpix_ind(i,2) = dim(2)
		!write(*,*) 'indice' , i, dim(1), dim(2)
	enddo
!	stop
!	call cpu_time(finish)
!	write(*,*) 'time ', finish-start


	do while(dist_dpix>0)
	
!		call cpu_time(start)
		do i=1,size(dpix)
			la(dpix_ind(i,1),dpix_ind(i,2)) = latry(i,dpix(dpix_ind(i,1),dpix_ind(i,2)))
		enddo
		
!		call cpu_time(finish)
!		cumtime(1)= cumtime(1) + finish-start
!		write(*,*) 'time picker', cumtime(1)

		! n-by-nd matrix of mu_t*la_t as a function of the current state (n in total) and possible values for next-period debt (nd in total)
		
!		call cpu_time(start)
		MUtry = LHSEulertry - (1+r)*betta*repmat(matmul(pai,la),nd,1)
!		MUtry = (1+r)*betta*repmat(matmul(pai,la),nd,1)		

!		call cpu_time(finish)
		cumtime(2)= cumtime(2) + finish-start
!		write(*,*) 'time repmat', finish-start

		! Find MUtry =0, by identifying sign change along rows
!		write(*,*)dist_dpix, sum(aa)
		
!		call cpu_time(start)
		MUtry = bind(MUtry)
!		call cpu_time(finish)
!		cumtime(3)= cumtime(3) + finish-start
!		write(*,*) 'time bind', finish-start
		
!		call cpu_time(start)	
		AMUtry = abs(MUtry)
		where (slacktry<0) AMUtry= 8e10
		AMUtry = AMUtry*Ztry 
!		call cpu_time(finish)
!		cumtime(4)= cumtime(4) + finish-start
	
!		call print_matrix(Ztry(1:100,95:101))
!		stop 
	
!		call cpu_time(start)	
		dpixnew = reshape(minloc(AMUtry,2),shape(dpixnew))		
!		call cpu_time(finish)
!		cumtime(5)= cumtime(5) + finish-start
		
		
		! Equilibrium selection  criterion (a)
		if(eqm_selection_criterion=='a')then
			do i = 1,n				
				temp = MUtry(i,dpixnew(dpix_ind(i,1),dpix_ind(i,2)))
				if(bind1(i)/=0 .and. temp<0 )then 
					dpixnew(dpix_ind(i,1),dpix_ind(i,2)) = bind1(i)
				end if 
			end do 
		end if 
		
!		write(*,*) ' test'
!		write(*,*) dpixnew(1,:)
!			call print_matrix(real(dpixnew(1,1:10)))
!		if(cont==10)then
!			stop
!		end if
		
		call cpu_time(start)	
		! Equilibrium selection  criterion (c)
		if(eqm_selection_criterion=='c')then 		
			do i =1,n
				if (bind1(i)/=0 .and. MUtry(i,bind1(i))>=0)then
					dpixnew(dpix_ind(i,1),dpix_ind(i,2)) = bind1(i)
				end if  
			enddo
		end	if 
		call cpu_time(finish)
		
		! Equilibrium selection  criterion (b)
		if(eqm_selection_criterion=='b')then 
			do i =1,n
				if (bind1(i)/=0 .and. MUtry(i,bind1(i))>=0)then
					dpixnew(dpix_ind(i,1),dpix_ind(i,2)) = bind1(i)
				end if
				if (bind2(i)/=0 .and. MUtry(i,bind2(i))>=0)then
					dpixnew(dpix_ind(i,1),dpix_ind(i,2)) = bind2(i)
				end if
			enddo
		end	if 
		cumtime(6)= cumtime(6) + finish-start
!		write(*,*) 'time eq', finish-start

		
!		call cpu_time(start)
		dist_dpix = maxval(maxval(abs(dpixnew-dpix),1),1)
!		call cpu_time(finish)
!		cumtime(7)= cumtime(7) + finish-start

!		write(*,*) 'time max', finish-start

!		call cpu_time(start)	
		xR = count(dpix/=dpixnew) !number of states that have not yet converged (for monitoring only, not used in approximation)
!		call cpu_time(finish)
!		write(*,*) 'time count', finish-start
		! Apply a one-step-at-a-time  updating criterion for the policy funciton
		
		where(dpixnew-dpix>0) dpix=dpix+1
		where(dpixnew-dpix<0) dpix=dpix-1
		
		write(*,*) dist_dpix, xR
		cont = cont + 1
		aa = 0.
		if(dist_dpix==1)then
			aa(cont) = 1.
		end if
!		write(*,*)sum(aa)
		if(sum(aa)==50.)then
			dist_dpix = 0
			warnn = 1
		else
			warnn = 0
		end if 
	end do
	
	p=0.
	do i=1,size(dpix)
		p(dpix_ind(i,1),dpix_ind(i,2)) = ptry(i,dpix(dpix_ind(i,1),dpix_ind(i,2)))
		dp(dpix_ind(i,1),dpix_ind(i,2)) = dgrid(dpix(dpix_ind(i,1),dpix_ind(i,2)))		
	enddo
	!write(*,*) p(:,92:101)
	!stop

	do i=1,size(dpix)
		c(dpix_ind(i,1),dpix_ind(i,2)) 	= ctry(i,dpix(dpix_ind(i,1),dpix_ind(i,2)))
   		cT(dpix_ind(i,1),dpix_ind(i,2)) = cTtry(i,dpix(dpix_ind(i,1),dpix_ind(i,2)))
		dY(dpix_ind(i,1),dpix_ind(i,2)) = dYtry(i,dpix(dpix_ind(i,1),dpix_ind(i,2)))		
	enddo
	write(*,*) 'Finish'
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