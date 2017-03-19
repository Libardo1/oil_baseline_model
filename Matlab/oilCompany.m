function [sp,x,revx,profit,EE] = oilCompany(Pxprob,PX,S,smin,smax,d,kappa,beta)
% OILCOMPANY  oilCompany solves the maximazation of the firm using 
%             Time Iteration algorithm.
%
% Usage:      [sp,x,revx,profit,EE] = oilCompany(Pxprob,PX,S,smin,smax,d,kappa,beta)
%
% Input:
%
%       Pxprob  : Probability transition matrix for oil prices
%       PX      : Oil prices grid
%       S       : Reserves grid 
%       smin    : Minimum reserves level
%       smax    : Maximum reserves level
%       d       : Discovery rate
%       kappa   : Cost parameter of the firm
%       beta    : Intertemporal discount factor
%
% Output:
%
%       sp      : Tomorrow´s reserves
%       x       : Extraction
%       revx    : Revenues
%       profit  : Profits
%       EE      : Euler errors of firm´s maximization problem
%
%% Construct the firm's state space

NPX = length(PX);
NS  = length(S);

s  = repmat(S,1,NPX);
sp = repmat(S,1,NPX);
px = repmat(PX,NS,1);

%% Time Iteration Loop 

%%%%%%%%%%%%%%  Technical parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter            = 0;            % 
d2              = 100;          %
uptd_of         = 0.1;          % Weight on new policy function in update to next iteration
iter_tol_of     = 3000;         %
tol_of          = 1e-6;         %
outfreq_of      = 10;           % Display frequency (shows in screen each 10th iteration)

disp('DE Iter      Norm');

%%%%%%%%%%%%%%%%%%%%%%%%%% Start of iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while d2>tol_of && iter < iter_tol_of
    
    % Retrieve updated policy functions
    oldsp=sp;

    % Marginal profits
    mr = px - (2*kappa*(s+d-sp).*s - kappa.*(s+d-sp).^2 )./(s.^2);
       
    % Interpolation 
    for i=1:NS
        for j=1:NPX    
            emr(i,j)   = beta*interp1(S,mr,sp(i,j),'linear','extrap')*Pxprob(j,:)';  %Expected marginal profits

            x(i,j) = (px(i,j) - emr(i,j))*s(i,j)/(2*kappa);   
            x(i,j) = max(0.0001,x(i,j));
               
            sp(i,j) = max(s(i,j) - x(i,j) + d,smin); 
            sp(i,j) = max(s(i,j) - x(i,j) + d,smin);    
            EE(i,j) = x(i,j) - (px(i,j) - emr(i,j))*s(i,j)/(2*kappa);
        end
    end
    
  
    iter=iter+1; % Update iteration counter
    
    % Calculate difference between new and old policies
    d2=max(max(abs(sp-oldsp)));
    
    % Print results once every (outfreq) iterations
    if mod(iter, outfreq_of) == 0
        fprintf('%d          %1.7f \n',iter,d2);
    end
    
    %=====================Updating rule for next iteration=================
    sp=uptd_of*sp+(1-uptd_of)*oldsp;
    %======================================================================
    
end
fprintf('%d          %1.7f \n',iter,d2);    % Print last iteration result

%% Other equations of the firm 
revx = px.*x;
profit = px.*x - kappa.*(x.^2)./s;

end