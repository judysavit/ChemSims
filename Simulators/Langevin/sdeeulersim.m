function sdeeulersim
% Figs. 9.2, 9.3 Euler-Maruyama Simulations (3/2007): 
%   Linear, Time-Dep. Coeff. SDE,
%   dX(t) = X(t)(mu(t)dt+sigma(t)dW(t)), X(0) = x0, t0 < t < tf,
%   Given  Initial data: x0, t0, tf, Nt; functions: f, g
clc
x0 = 1; t0 = 0; tf = 5; Nt = 2^10; 
randn('state',8);
DT = tf/Nt; sqrtdt = sqrt(DT); 
t = t0:DT:tf;
DW = randn(1,Nt)*sqrtdt;
W = cumsum(DW); % Note: omits initial zero value; count is off by 1;
%
Xexact = zeros(1,Nt+1); Xexact(1) = x0; % initialize
for k = 1:Nt % Exact formula to fine precision for exact consistency:
    Xexact(k+1) = xexact(x0,t(k+1),W(k));
end
% Lumped coarse sample from fine sample:
L = 2^3; NL = Nt/L; KL = 0:L:Nt; DTL = L*DT; tL = t0:DTL:tf;
fprintf('(N_t,NL)=(%i,%i); Size(t,KL,tL)=[(%i,%i);(%i,%i);(%i,%i)];'...
    ,Nt,NL,size(t),size(KL),size(tL));
Xeul = zeros(1,NL+1); Xeul(1) = x0;
Xdiff = zeros(1,NL+1); Xdiff(1) = Xeul(1) - Xexact(1);
for k = 1:NL % Euler-Maruyama formula to coarse precision:
    DWL = sum(DW(1,KL(k)+1:KL(k+1))); 
    Xeul(k+1) = Xeul(k) + f(Xeul(k),tL(k))*DTL+g(Xeul(k),tL(k))*DWL;
    Xdiff(k+1) = Xeul(k+1) - Xexact(KL(k+1));
end
%
scrsize = get(0,'ScreenSize');
ss = [3.0,2.8,2.6,2.4,2.2,2.0];
%
nfig = 1;
figure(nfig);
plot(tL,Xeul,'k--','linewidth',3); hold on
plot(t,Xexact,'k-','linewidth',3); hold off
axis([t0 tf 0 max(max(Xeul),max(Xexact))]);
title('Euler-Maruyama and Exact Linear SDE Simulations'...
    ,'Fontsize',44,'FontWeight','Bold');
xlabel('t, Time'...
    ,'Fontsize',44,'FontWeight','Bold');
ylabel('X(t), State'...
    ,'Fontsize',44,'FontWeight','Bold');
legend('Xeul(t): Euler','Xexact(t): Exact','Location','NorthWest');
set(gca,'Fontsize',36,'FontWeight','Bold','linewidth',3);
set(gcf,'Color','White','Position' ...
    ,[scrsize(3)/ss(nfig) 70 scrsize(3)*0.60 scrsize(4)*0.80]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Xdiffmax = max(abs(Xdiff));
fprintf('\nMaximal Euler-Exact Absolute Error:');
fprintf('\n     max(abs(Xeul(TL)-Xexact(TL)))=%8.2e=%8.2e*DTL;\n'...
    ,Xdiffmax,Xdiffmax/DTL);
% (N_t,NL) = (1024,128); Size(t,KL,tL) = [(1,1025);(1,129);(1,129)];
% Maximal Euler-Exact Abs. Error:
%      max(abs(Xeul(TL)-Xexact(TL))) = 1.31e+00 = 3.36e+01*DTL;
%
nfig = nfig+1;
figure(nfig);
plot(tL,Xdiff,'k-','linewidth',3);
axis tight;
title('Euler and Exact Linear SDE Simulations Error'...
    ,'Fontsize',44,'FontWeight','Bold');
xlabel('t, Time'...
    ,'Fontsize',44,'FontWeight','Bold');
ylabel('Xeul(t)-Xexact(t), Error'...
    ,'Fontsize',44,'FontWeight','Bold');
set(gca,'Fontsize',36,'FontWeight','Bold','linewidth',3);
set(gcf,'Color','White','Position' ...
    ,[scrsize(3)/ss(nfig) 70 scrsize(3)*0.60 scrsize(4)*0.80]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = f(x,t)
mu = 1/(1+0.5*t)^2;  
y = mu*x;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = g(x,t)  % for general case, ignore "t not used" warning.
sig = 0.5;
y = sig*x;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = xexact(x0,t,w)
% exact solution if available for general linear SDE:
mubar = 2-2/(1+0.5*t); sig = 0.5; sig2bar = sig^2*t/2; 
y = x0*exp(mubar-sig2bar + sig*w);
%
% end sdeeulersim.m