function [T,X] = langevinCustom(stoich_matrix, propensity_fcn, tspan, x0, rate_params, input)%   Returns:
%function [T,X] = langevinCustom(stoich_matrix, propensity_fcn, tspan, x0, rate_params, input)%   Returns:
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
tau = 10^(-2); % in sec.
tVect = tspan(1):tau:tspan(2); tVect = tVect';
tLen = length(tVect); 

n = length(x0); % n is number of species
m = size(stoich_matrix,1); % m is the number of rxns
X = zeros(tLen, n);
X(1,:)=x0;

% Create correlated noise
rNoise = permute(repmat(randn(m,tLen)',[1,1,n]),[2 3 1]);
stochTerm = repmat(stoich_matrix, [1,1,tLen]).*rNoise;

for i = 1:tLen-1
    a = repmat(propensity_fcn(X(i,:)).*tau,1,n);
    if (any(a)<0)
        disp(['PROP < 0 at ' num2str(i)]);
    end
    exactTerm = stoich_matrix.*a;
    noiseTerm = stochTerm(:,:,i).*sqrt(a);%repmat(corrNoise(:,i),1,n));
    tempX = X(i,:)+sum(exactTerm + noiseTerm,1);
    X(i+1,:)= tempX;
    %X(i+1,:) = ternif( any(tempX<0), X(i,:), tempX );
end

% Record output
T = tVect; % later can do this with a solver and deval, unsure what mesh to use now
  

end