program ergodic_constrainedV2F_OF

	use basic 
	use tools
	implicit none 

	character(100) :: mod_modif, kof, learning_type_for, yT_nature, fileplace

	integer, parameter :: nd =  101, ny = 340, ns = 170, np = 2, n =ns*np , NNy = 1
	integer, parameter :: T = 10000, LE = 1
	
	integer :: hh_on

	integer :: Tburn, i, time 
 
	real, dimension(ny,nd) :: tau, crisis, mu 
	real :: Pp(np,np), CPp(np,np), v_sX0, sX0, sgrid(ns), grid(NNy*ny,2)
	real, dimension(ny*NNy) :: yTgridlevel, revXgridlevel
	real :: dgrid(nd)
	real :: revX0, yT0, v_d0, d0
	integer :: h, g 

	real :: YYT(NNy),REVX(Ny)	
	real :: CYTTprob(NNy,NNy), YYTprob(NNy,NNy)

	real :: PX(ny), pX0
	
	real, dimension(T,LE) :: PX_SIM = 0., SX_SIM = 0., SXp_SIM = 0., X_SIM = 0., REVX_SIM = 0.
	real, dimension(T,LE) :: LA_SIM = 0., YT_SIM = 0., D_SIM = 0., Dp_SIM = 0., V_SIM = 0., PN_SIM = 0.
	real, dimension(T,LE) :: CT_SIM = 0., C_SIM = 0., SLACK_SIM = 0., MU_SIM = 0., Y_SIM = 0., CAY_SIM = 0.
	real, dimension(T,LE) :: CA_SIM = 0., DY_SIM = 0., TB_SIM = 0., TBY_SIM = 0.


	mod_modif           = trim('')
	kof                 = trim('Type_2')
	learning_type_for   = trim('RE')
	yT_nature           = trim('YT_FIXED')
 
	! read policy functions
!	sgrid , Pp	

	fileplace = 'outData/example'
	
    OPEN (unit = 1, file = trim(fileplace)//'sgrid.txt', status = 'old')
	READ (1,*) sgrid
	CLOSE (1)
	
    OPEN (unit = 1, file = trim(fileplace)//'RREVX.txt', status = 'old')
	READ (1,*) RREVX
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'YYT.txt', status = 'old')
	READ (1,*) YYT
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'dgrid.txt', status = 'old')
	READ (1,*) dgrid
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'PX.txt', status = 'old')
	READ (1,*) PX
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'YYTprob.txt', status = 'old')
	READ (1,*) YYTprob
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'M_REVXr.txt', status = 'old')
	READ (1,*) M_REVXr
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'M_Xr.txt', status = 'old')
	READ (1,*) M_Xr
	CLOSE (1)
	
	OPEN (unit = 1, file = trim(fileplace)//'Pp.txt', status = 'old')
	READ (1,*) Pp
	CLOSE (1)
	
	if(learning_type_for .eq. 'RE' .and. LE/= 1)then
		write(0,*) 'Error : Bad LE initialization'
		stop
	else if(learning_type_for .eq. 'BL' .and. LE/= 25)then
		write(0,*) 'Error : Bad LE initialization'
		stop
	end if 

	Tburn    = int(0.01*T)
   
	! Initial Conditions
	CPp = cumsum(Pp)
 
	v_sX0 =  1.0823
	h = minloc(sgrid-v_sX0,1)
	sX0 = sgrid(h)
 
	if(hh_on==1)then 
		grid = gridmake2(YYT,REVX)
		yTgridlevel =grid(:,1)
		revXgridlevel = grid(:,2)
		revX0 = revXgridlevel(1)
		yT0   = yTgridlevel(1)
		v_d0 =  0.6440     				!	Initial debt
		g = minloc(dgrid-v_d0,1)
		d0 = dgrid(g)
	end if

	! Initializations' Firm
!	PX_SIM      = 0.
!	SX_SIM      = 0.
!	SXp_SIM     = 0.
!	X_SIM       = 0.
!	REVX_SIM    = 0.
 
	! Initializations' Household
!	LA_SIM      = 0. !marginal utility of cT
!	YT_SIM      = 0.  !yT
!	D_SIM       = 0. !current debt
!	Dp_SIM      = 0.  !next-period debt
!	V_SIM       = 0. !value function
!	PN_SIM      = 0. !relative price nontradables
!	CT_SIM     = 0. !consumption of tradables
!	C_SIM      = 0.   !composite consumption
!	SLACK_SIM  = 0.  !slackness of collateral constraint
!	MU_SIM     = 0.  !collateral constraint multiplier
!	Y_SIM      = 0.
!	CAY_SIM    = 0.
!	CA_SIM     = 0.
!	DY_SIM     = 0.
!	TB_SIM     = 0.
!	TBY_SIM    = 0.
 
!	pX0 = PX(2) 
 
!	CYTTprob = cumsum(YYTprob)

!	do l = 1,LE
!	do time = 1,T
!		! Firm
!		PX_SIM(time,l)    = pX0
!		SX_SIM(time,l)    = sX0
!		m        	      = minloc(pXgridlevel-pX0,1)
!		n           	  = minloc(sgrid-sX0,1)
!		REVX_SIM(time,l)  = M_REVXr(m,n,l)
!		X_SIM(time,l)     = M_Xr(m,n,l)
!		SXp_SIM(time,l)   = M_spr(m,n,l)
!		revX0             = M_REVXr(m,n,l)
!		sX0           	  = M_spr(m,n,l)
		
		! Household
!	end do 
!	end do 
end program 