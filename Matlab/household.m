function [pn,ct,c,bp,bpy,gdp] = household(bp,b,yt,ct,c,pn,B,yn,at,r,bpmax,revx,Prob, Ntotal, nd,ns,cbind,sigma,omega,gamma,kappa,bmin,bmax,cn,beta)
%OILCOMPANY  Household solves the maximazation of the houselhold using 
%
% Input: 
%
%   bp        :  Next period큦 debt
%   b         :  Debt 
%   yt        :  Tradable GDP
%   ct        :  Tradable consumption
%   c         :  Consumption
%   pn        :  Non-tradable큦 price 
%   B         :  Debt큦 grid
%   yn        :  Non tradable GDP
%   at        :  Assets
%   r         :  Interest rate
%   bpmax     :  Maximum level of tomorrow큦 debt
%   revx      :  Revenues
%   Prob      :  Transition matrix
%   Ntotal    :  Total number of grid points for the exogenous variables
%   nd        :  Size of debt큦 grid
%   ns        :  Size of reserves큦 grid
%   cbind     :  Maximum consumption level when the restriction is binding
%   sigma     :  Interteporal consumption elasticity
%   omega     :  Elasticity of substitution between types of consumption
%   gamma     :  Non tradable and tradable consumption participation on
%                total consumption bundle
%   kappa     :  Parameter od debt
%   bmin      :  Minimum of debt level
%   bmax      :  Maximum of debt level
%   cn        :  Non tradable consumption
%   beta      :  Discount factor
%
% Output:
%
%   pn        :  Non trables price
%   ct        :  Tradable consumption
%   c         :  Total consumption
%   bp        :  Tomorrow큦 debt
%   bpy       :  Debt to GDP
%   gdp       :  GDP

%%%%%%%%%%%%%%%%%%%%%%
thereshold_ct = 1e-4;  % Minimum value of consumption

%% Time Iteration Loop 

%%%%%%%%%%%%%%%%%%%%%%  Technical parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uptd     = 0.9;     % Weight on new policy function in update to next iteration
outfreq  = 5;       % Display frequency (shows in screen each 20th iteration)
tol      = 3E-3;    % Numerical tolerance (convergence criterion) 

iter     = 0;
d2       = 100;     

disp('DE Iter      Norm');


%%%%%%%%%%%%%%%%%%%%%%%%%% Start of iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while( d2>tol )
    
    % Retrieve updated policy functions
    oldct = ct;
        
    % Marginal utility of consumption given current consumption policy 
    mu  = real(c.^-sigma.*gamma.^(1/omega).*(ct./c).^(-1/omega));  
        
    % Interpolation  of expected marginal utility
    % We interpolate the expected marginal utility for each level of
    % reserves (k)
    for k=1:ns
        mu_aux = mu((k-1)*nd+1:k*nd,:); 
        for i =1:nd
            for j = 1:Ntotal 
                emup_aux(i,j,k)  = interp1(B,mu_aux,bp(i,j),'linear')*Prob(j,:)';
            end
        end
    end
    
    for t =1:ns
       emup((t-1)*nd+1:t*nd,:)=emup_aux(:,:,t);
    end
     
    %------------------ Find new constrained values -----------------------
    gdp = (yt+pn*yn+revx);
    bpmax = kappa.*gdp; 
    bpmax(bpmax>bmax) = bmax;
    bpmax(bpmax<bmin) = bmin;
    ctbind = -(1+r)*b + revx + yt - at + bpmax;
    ctbind = max(ctbind,thereshold_ct);
    cbind  = (gamma^(1/omega)*ctbind.^((omega-1)/omega)+(1-gamma)^(1/omega)*cn.^((omega-1)/omega)).^(omega/(omega-1));
    %----------------------------------------------------------------------
    
      for i=1:nd*ns
        for j=1:Ntotal
       
            % Calculate Euler Equation Error at maximum feasible consumption cbind
            EE(i,j) = cbind(i,j)^-sigma*gamma^(1/omega)*(ctbind(i,j)/cbind(i,j))^(-1/omega)-(1+r)*beta*emup(i,j);
        
            if EE(i,j)>tol*1e-5            % If positive, constraint will be binding then:
               
                % Debt will be as big as possible                
                bp(i,j) = bpmax(i,j);      

                % Consumption is solved from budget constraint
                c(i,j)  = cbind(i,j);                            
                ct(i,j) = ctbind(i,j);
                pn(i,j) = real(((1-gamma)/gamma.*(ct(i,j)/cn)).^(1/omega));
                gdp(i,j) = yt(i,j)+ pn(i,j)*yn + revx(i,j);
                bpy(i,j) = bp(i,j)/gdp(i,j);
               
            else % Constraint not binding
                
                % Define function that calculates the absolute value of Euler
                % Equation Error for a given consumption
               
                f = @(cct) abs((((gamma^(1/omega)*cct.^((omega-1)/omega)+(1-gamma)^(1/omega)*cn.^((omega-1)/omega)).^(omega/(omega-1))).^(-sigma+(1/omega)).*gamma.^(1/omega).*(cct).^(-1/omega))-(1+r)*beta*emup(i,j));
                
                [ct(i,j),EE(i,j)] = fminbnd(f,thereshold_ct,1);
                ct(i,j) = max(ct(i,j),thereshold_ct);
                
                c(i,j) = ((gamma^(1/omega)*ct(i,j).^((omega-1)/omega)+(1-gamma)^(1/omega)*cn.^((omega-1)/omega)).^(omega/(omega-1)));
    
                pn(i,j) = real(((1-gamma)/gamma.*(ct(i,j)/cn)).^(1/omega));
                gdp(i,j) = yt(i,j)+ pn(i,j)*yn + revx(i,j);
                 
                % Solve debt from budget constraint, check if it is within grid bounds

                bp(i,j)=max(-revx(i,j) - yt(i,j)+(1+r)*b(i,j)+ ct(i,j)+ at, bmin);
                bp(i,j)=min(bp(i,j),bmax); 
                bpy(i,j) = bp(i,j)/gdp(i,j);
            end
        end
    end    
  iter=iter+1; % Update iteration counter
   
  % Calculate difference between new and old policies
    
  d2 = max(max(abs(ct-oldct)));

  % Print results once every (outfreq) iterations
  if mod(iter, outfreq) == 0;
    fprintf('%d          %1.7f \n',iter,d2);
  end
    
    %=====================Updating rules for next iteration================
    ct = uptd*ct+(1-uptd)*oldct;
    pn = real(((1-gamma)/gamma.*(ct./cn)).^(1/omega));
    c = (gamma^(1/omega)*ct.^((omega-1)/omega)+(1-gamma)^(1/omega)*cn.^((omega-1)/omega)).^(omega/(omega-1));
    bp = max(-revx - yt+(1+r)*b + ct+ at, bmin);
    bp =min(bp,bmax); 
    %======================================================================
end 

fprintf('%d          %1.7f \n',iter,d2);    % Print last iteration result

