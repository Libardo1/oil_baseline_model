clear;
%% Setting model block's

%% Model's Calibration
%% Oil's Firm Parameters
   
Pj      = [22.4787 56.8691];      % [pl ph]   values of each state
tj      = [58.66   62.13];        % [tsl tsh] periods in each state
Pp_aux  = [1-1/tj(1) 1-1/tj(2)];  % [qll qhh] prob. of staying in each state
pr      = (1/Pj(2))*Pj;           % price normalization ph==1
Pp      = [ Pp_aux(1)  1-Pp_aux(1); 1-Pp_aux(2)  Pp_aux(2) ]; 
Epr     = pr*ergdist(Pp);         % expected oil price 
p       = pr./Epr;                % re-normalize prices to Ep = 1 
Ep      = p*ergdist(Pp);          % re-compute expected value to 1
En      = 25; 
r       = (1.035)^(.25)-1;
R       = r+1;
x2s     = 1/En;
kappa   = (-Ep*(1/R-1))/(2*x2s*(1-1/R)+(1/R)*x2s*x2s);   
disc    = 0.056;          
  
%% Household
 
betta   =  0.9825;         % discount factor (in order to have nice distrubution of debt)
yN      =  exp(-0.510734);  % nontraded output
sigg    =  2;               % intertemporal elasticity of consumption
gamma   =  0.493009;        % weight on traded consumption in the CES aggregator 
omega   =  0.316;
aT      = 0.073025;
aN     	= 0.27507;
cn      = yN - aN;
kapa    = 3.2;
                         
%% Construct Firm's state space 

smax   = 0.056*34; 
smin   = 0.0001; 
s      = (smin:disc/2:smax)'; 
ns     = length(s);           
np     = length(p);
[S,PXc] = gridmake(s,p');     
        
%% Construct Household's state space 

YT      = 0;
nd      = 40;               % Number of grid points for debt
yt_bar 	= exp(-1.06692);
YT      = exp(YT)*yt_bar;
blower  = 0;                % Minimum debt level
bupper  = 2.75;             % Maximum debt level

B=zeros(nd,1); 
B(1) = blower;

for i=2:nd
    B(i)=B(i-1)+(bupper-B(i-1))/((nd-i+1)^1.05); 
end

b       = repmat(B,ns,np);                % Initial debt
yt      = repmat(YT',nd*ns,np);               % Tradable GDP
cbind   = ones(nd*ns,np);                 % Initial constrained consumption values

%% Load initial guess
load('guess_RE.mat');

%% Solving the model
%% FIRM
bp = M_dpguess;   %Initial debt (obtained from Uribe)

[spr, Xr,REVXr,PROFIT_Xr,EE] = oilCompany(Pp,p,s,smin,smax,disc,kappa,betta);

%% HOUSEHOLD                   
% Preparing data of Firm solution to use in te Household problem 
% The state space of the household has on the columns the exogenous variables
% on the rows are the endogenous variables (revenues and debt)
% Each nd points represent a policy function (debt) associated with each level of reserves

revx    = kron(REVXr,ones(nd,1));
ct      = max(0.00001,-(1+r)*b + revx + yt - aT + bp);       %Budget constraint 
pn      = real(((1-gamma)/gamma.*(ct./cn)).^(1/omega));      %Price of non tradables 
bpmax   = kapa.*(yt + pn * yN + revx);                       % Debt Limit
c       = ((gamma^(1/omega)*ct.^((omega-1)/omega)+(1-gamma)^(1/omega)*cn.^((omega-1)/omega)).^(omega/(omega-1)));  %Consumption bundle

[pnn,cT,c,bp,bpy,gdp] = household(bp,b,yt,ct,c,pn,B,yN,aT,r,bpmax,revx,Pp,np,nd,ns,cbind,sigg,omega,gamma,kapa,blower,bupper,cn,betta);
                

%% Saving Outcomes
save oil_model_RE;