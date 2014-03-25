function [T, X] = net1sim(p, X0, mode, tspan, scale)
% [T, X] = net1sim(X0, p, mode)
% Input: p  = vector of parameter values, including input stim (length = 13)
%        mode = simulation type, values are 'CME' or 'langevin'
%        (opt) X0 = vector of initial X values (length = 3)
% Output: T = vector of event times
%         X = matrix of species counts (#species x length(T))
% Reaction network:
%   dA/dt = k_uA*u(t)(1-A)/(1-A+K_uA)-k_A*A/(A+K_A)
%   dB/dt = k_AB*A(1-B)/(1-B+K_AB)-k_B*B/(B+K_B)
%   dC/dt = k_AC*A(1-C)/(1-C+K_AC)-k_BC*C/(C+K_BC)

% Check input parameters
if (length(p)~= 13); disp('Parameter vector should be 13x1'); end;
if ~ismember(mode, {'CME', 'langevin'}); disp('Illegal mode'); end;


stoich_matrix = [ 1 0 0 ; % prod. A
    -1 0 0 ; % deg. A
    0 1 0 ; % prod. B
    0 -1 0; % deg. B
    0 0 1 ; % prod. C
    0 0 -1]; %deg. C


% Assume system size = 1000 molecules of each type
stoich_matrix = stoich_matrix./1000;


% Run simulation
switch mode
    case 'CME'
        [T,X] = firstReactionMethod(stoich_matrix, @prop_fcn_net1, tspan, X0, p,scale);
    case 'langevin'
        tau = 1; %seconds
        try
            p_stoch = convertParams(p);
        catch
            disp('MUST IMPLEMENT PARAMETER CONVERSION FROM SSA TO LANGEVIN');
        end
        [T,X] = langevinApprox(stoich_matrix, @prop_fcn_net1, tspan, X0, p, scale);
end

end


function a = prop_fcn_net1(X, p, scale)
u = p(1);
k_uA = p(2);
K_uA = p(3);
k_A = p(4);
K_A = p(5);
k_AB = p(6);
K_AB = p(7);
k_B = p(8);
K_B = p(9);
k_AC = p(10);
K_AC = p(11);
k_BC = p(12);
K_BC = p(13);

A = X(1);
B = X(2);
C = X(3);

a(1) = k_uA.*u.*(scale-A)./(scale-A+K_uA);
a(2) = k_A.*A./(A+K_A);
a(3) = k_AB.*A.*(scale-B)/(scale-B+K_AB);
a(4) = k_B.*B./(B+K_B);
a(5) = k_AC.*A.*(scale-C)./(scale-C+K_AC);
a(6) = k_BC.*C./(C+K_BC);
a = a';
end


