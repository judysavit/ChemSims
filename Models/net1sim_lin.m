function [T, X] = net1sim_lin(p, X0, mode, tspan)
% [T, X] = net1sim_lin(X0, p, mode)
% Input: p  = vector of parameter values, including input stim (length = 13)
%        mode = simulation type, values are 'CME' or 'langevin'
%        (opt) X0 = vector of initial X values (length = 3)
% Output: T = vector of event times
%         X = matrix of species counts (#species x length(T))
% Species:
%   AA - active A
%   BA - active B
%   CA - active C
% Reaction network (deterministic):
%   dA/dt = k_1*u(t)(AT-AA)-k_2*AA
%   dB/dt = k_3*AA*(BT-BA)-k_4*BA
%   dC/dt = k_5*AA*(CT-CA)-k_6*BA*CA./(K_7+CA)-k_8
% Reaction network (stochastic, 7 rxns):
%   c1: u + AI -> u + AA
%   c2: AA -> AI
%   c3: AA + BI -> AA + BA
%   c4: BA -> BI
%   c5: AA + CI -> AA + CA
%   c6,c7: BA + CA <-> BACA -> BA + CI
%   c8: CA -> CI
% Parameters:
%   8 stochastic rate constants (c values) above
%   1 input constant
%   3 total A, B, C constraints (AT, BT, CT)

% Check input parameters
if (length(p)~= 11); disp('Parameter vector should be 11x1'); end;
if ~ismember(mode, {'CME','Noise', 'langevin','langevinExact','ODE'}); disp('Illegal mode'); end;

stoich_matrix = [ 1 0 0 ; % prod. A
    -1 0 0 ; % deg. A
    0 1 0 ; % prod. B
    0 -1 0; % deg. B
    0 0 1 ; % prod. C
    0 0 -1 ; %inh. C by B
    0 0 -1]; %deg. C

% Run simulation
switch mode
    case 'CME'
        simMode = @directMethod;
        simFun = @prop_fcn_net1_lin;
    case 'langevin'
        simMode = @langevinCustom;
        simFun = @prop_fcn_net1_lin;
    case 'langevinExact'
        simMode = @langevinExact;
        simFun = @prop_fcn_net1_lin;
    case 'ODE'
        simMode = @deterministicODE;
        simFun = @ODE_fcn_net1_lin;
    case 'Noise'
        T = 0;
        X = noise_fun_net1sim_lin(X0,p, stoich_matrix);
        return;
    otherwise
        error('Incorrect simulation mode');
end
 
%Get to steady state
[T1,X1] = simMode(stoich_matrix, simFun, [1 100], X0, p);
%Pulse
p(1) = p(1).*5;
[T2,X2] = simMode(stoich_matrix, simFun, tspan, X1(end,:), p);

T = [T1(1:end-1); T2+T1(end)];
X = [X1(1:end-1,:); X2];
end


function a = prop_fcn_net1_lin(X, p)
u = p(1);
c1 = p(2);
c2= p(3);
c3 = p(4);
c4 = p(5);
c5 = p(6);
c6 = p(7);
c7 = p(8);
AT = p(9);
BT = p(10);
CT = p(11);

% S = [1 X]';
% m = [c1.*u.*AT -c1.*u 0 0;
%     0 c2 0 0;
%     0 c3*BT -c3.*X(2) 0;
%     0 0 c4 0;
%     0 c5.*CT -c5.*X(2) 0;
%     0 0 0 c6.*X(3);
%     0 0 0 c7];
% disp(S);
% a = m*S;
    
% Alternative implementation:
A = X(1);
B = X(2);
C = X(3);
a(1) = c1.*u.*(AT-A); 
a(2) = c2.*A; 
a(3) = c3.*A.*(BT-B);
a(4) = c4.*B;
a(5) = c5.*A.*(CT-C);
a(6) = c6.*B.*C;
a(7) = c7.*C;
a = a';

% Implementation improvement: do matrix multiplication by the stoich
% matrix.

end

function noise = noise_fun_net1sim_lin(X, p, stoich_matrix)
n = size(X,2);
a = prop_fcn_net1_lin(X,p);
noise = sum(stoich_matrix.*repmat(sqrt(a),1,n));
end

function dydt = ODE_fcn_net1_lin(t, X, p)
u = p(1);
c1 = p(2);
c2= p(3);
c3 = p(4);
c4 = p(5);
c5 = p(6);
c6 = p(7);
c7 = p(8);
AT = p(9);
BT = p(10);
CT = p(11);

A = X(1);
B = X(2);
C = X(3);

dydt = [c1*u.*(AT-A) - c2.*A;
    c3.*A.*(BT-B) - c4.*B;
    c5.*A.*(CT-C) - c6.*B.*C - c7.*C];
    
% Alternative implementation:
% A = X(1);
% B = X(2);
% C = X(3);
% a(1) = c1.*u.*(AT-A); 
% a(2) = c2.*A; 
% a(3) = c3.*A.*(BT-B);
% a(4) = c4.*B;
% a(5) = c5.*A.*(CT-C);
% a(6) = c6.*B.*C;
% a(7) = c7.*C;
% a = a';

% Implementation improvement: do matrix multiplication by the stoich
% matrix.

end