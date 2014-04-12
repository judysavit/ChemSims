function [t,x] = deterministicODE(ode_fcn, tspan, x0, outputFcn)
% function [t,x] = deterministicODE(ode_fcn, tspan, x0, rate_params, input,
% outputFcn)
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


%% MAIN LOOP
% ASSUME UNIT VOLUME FOR SIMPLICITY - no need to change parameters

if ~exist('outputFcn','var'); outputFcn = []; end;

opts  = odeset('RelTol',1e-4, 'AbsTol', 1e-6, 'MaxStep',1, 'OutputFcn',outputFcn);
[T,X] = ode15s(ode_fcn, tspan, x0, opts);

% Record output
t = T; % later can do this with a solver and deval, unsure what mesh to use now
x = X;  

end


