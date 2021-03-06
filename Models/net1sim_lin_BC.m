function [T, X] = net1sim_lin_BC(p, X0, mode, tspan)
% [T, X] = net1sim_lin_BC(X0, p, mode)
% Input: p  = vector of parameter values, including input stim (length = 13)
%        mode = simulation type, values are 'CME' or 'langevin'
%        (opt) X0 = vector of initial X values (length = 3)
% Output: T = vector of event times
%         X = matrix of species counts (#species x length(T))
% Species:
%   AA - active A
%   BA - active B
%   CA - active C
%   BC - BC complex
% Reaction network (deterministic):
%   dA/dt = k_1*u(t)(AT-AA)-k_2*AA
%   dB/dt = k_3*AA*(BT-BA)-k_4*BA
%   dBC/dt = k5*B*C - (k6+k7)*BC
%   dC/dt = k_8*AA*(CT-CA)-k_5*B*C + k_6*BC - k_9*C

% Reaction network (stochastic, 7 rxns):
%   c1: u + AI -> u + AA
%   c2: AA -> AI
%   c3: AA + BI -> AA + BA
%   c4: BA -> BI
%   c5: AA + CI -> AA + CA
%   c6: BA + CA -> BACA
%   c7: BACA -> BA + CA
%   c8: BACA -> BA + CI
%   c9: CA -> CI

% Parameters:
%   9 stochastic rate constants (c values) above
%   1 input constant
%   3 total A, B, C constraints (AT, BT, CT)

% Check input parameters
if (length(p)~= 13); disp('Parameter vector should be 13x1'); end;
if ~ismember(mode, {'CME', 'langevin'}); disp('Illegal mode'); end;

stoich_matrix = [ 1 0 0 0; % prod. A
    -1 0 0 0; % deg. A
    0 1 0 0; % prod. B
    0 -1 0 0; % deg. B
    0 0 0 1; % prod. C
    0 -1 1 -1; %BC goes into complex
    0 1 -1 1; %BC goes out of complex
    0 1 -1 0; %BC goes to inactive C product
    0 0 0 -1]; %C is degraded

% Run simulation
switch mode
    case 'CME'
        simfun = @directMethod;
    case 'langevin'
        simfun = @langevinApprox;
end

%Get to steady state
[T1,X1] = directMethod(stoich_matrix, @prop_fcn_net1_lin, [1 100], X0, p);
%Pulse
p(1) = p(1).*5;
[T2,X2] = directMethod(stoich_matrix, @prop_fcn_net1_lin, tspan, X1(end,:), p);

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
c8 = p(9);
c9 = p(10);
AT = p(11);
BT = p(12);
CT = p(13);

% This needs to be modified:
% S = [1 X]';
% m = [c1.*u.*AT -c1.*u 0 0;
%     0 c2 0 0;
%     0 c3*BT -c3.*X(2) 0;
%     0 0 c4 0;
%     0 c5.*CT -c5.*X(2) 0;
%     0 0 0 c6.*X(3);
%     0 0 0 c7];
% a = m*S;
    
%Alternative implementation:
A = X(1);
B = X(2);
BC = X(3);
C = X(4);
a(1) = c1.*u.*(AT-A); 
a(2) = c2.*A; 
a(3) = c3.*A.*(BT-B);
a(4) = c4.*B;
a(5) = c5.*A.*(CT-C);
a(6) = c6.*B.*C;
a(7) = c7.*BC;
a(8) = c8.*BC;
a(9) = c9.*C;
a = a';

% Implementation improvement: do matrix multiplication by the stoich
% matrix.

end

