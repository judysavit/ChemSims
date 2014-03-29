function [t,x] = langevinApprox(stoich_matrix, propensity_fcn, tspan, x0, rate_params, u)
%   Returns:
%       t:              time vector          (Nreaction_events x 1)
%       x:              species amounts      (Nreaction_events x Nspecies)    
%
%   Required:
%       tspan:          Initial and final times, [t_initial, t_final].
%       x0:             Initial species values, [S1_0, S2_0, ... ].
%       stoich_matrix:  Matrix of stoichiometries (Nreactions x Nspecies).
%                       Each row gives the stoichiometry of a reaction.
%       prop_fcn:       Function handle to function that calculates
%                       reaction propensities.
%                       Target function should be of the form
%                           ac = f( xc, rate_params ),
%                       where xc is the current state [S1, S2, ...], and
%                       rate_params is the user-defined rate parameters.
%                       The function should return vector ac of 
%                       propensities (Nreactions x 1) in the same order as
%                       the reactions given in stoich_matrix.
%       rate_params:    These will be passed to the propensity function

%% MAIN LOOP
mode = 'exact';
switch mode
    case 'exact';
    %Use exact langevin (no noise term)
odefun = @langevinODE_exact;
noise = [];
    case 'noisy-RV';
%Use noisy langevin with RV-based noise propagation
odefun = @langevinODE_RVnoise;
noise = [];
    case 'noisy-spline';
%Use noisy langevin with noise determined as a spline
odefun = @langevinODE_splineNoise;
tVect = linspace(tspan(1), tspan(2), 10^6);
rnoise = randn(1, length(tVect));
cumnoise = cumsum(rnoise);
noise = spline(tVect,cumnoise);
end

[T,X] = ode45(odefun, tspan, x0, [], stoich_matrix, propensity_fcn, noise, rate_params, u );

% Record output
t = T; % later can do this with a solver and deval, unsure what mesh to use now
x = X;  

end

function dydt = langevinODE_exact(T,X,stoich_matrix,prop_fcn,noise,p,u)
    a = prop_fcn(X',p,u);
    dydt = sum(stoich_matrix.*repmat(a,1,3))';   
end

function dydt = langevinODE_RVnoise(T,X,stoich_matrix,prop_fcn,noise,p)
error('THE LANGEVINODE_RVNOISE FCN DOESNT WORK YET');
D = X(1:3); % these RV's refer to the chem species
N = X(4); % the last entry propagates the noise term
a = prop_fcn(D',p);
dydt = sum(stoich_matrix.*repmat(a,1,3))'+ sum(stoich_matrix.*repmat(sqrt(a),1,3).*(N+randn(size(stoich_matrix))))';
end

function dydt = langevinODE_splineNoise(T,X,stoich_matrix,prop_fcn,noise,p)
    a = prop_fcn(X',p);
    dydt = sum(stoich_matrix.*repmat(a,1,3))'+ sum(stoich_matrix.*repmat(sqrt(a),1,3).*ppval(noise, T))';   
end
