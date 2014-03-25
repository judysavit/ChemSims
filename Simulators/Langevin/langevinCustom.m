function [t,x] = langevinCustom(stoich_matrix, propensity_fcn, tspan, x0, rate_params)
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


% Choice of tau is critical.
tau = 0.001; % in sec.
tVect = tspan(1):tau:tspan(2); tVect = tVect';
tLen = length(tVect);

n = length(x0); % n is number of species
X = zeros(tLen, n);
X(1,:)=x0;

% Create correlated noise
rNoise = randn(size(stoich_matrix,1),tLen);
corrNoise = cumsum(rNoise,2);

%% MAIN LOOP

for i = 1:tLen-1
    a = propensity_fcn(X(i,:),rate_params);
    exactTerm = sum(stoich_matrix.*repmat(a,1,n).*tau);
    noiseTerm = sum(stoich_matrix.*repmat(sqrt(a),1,n).*randn(size(stoich_matrix)).*sqrt(tau));%repmat(corrNoise(:,i),1,n));
    X(i+1,:) = X(i,:)+exactTerm + noiseTerm;
end

% Record output
t = tVect; % later can do this with a solver and deval, unsure what mesh to use now
x = X;  

end