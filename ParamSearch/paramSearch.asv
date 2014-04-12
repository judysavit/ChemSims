addpath(genpath('/home/jsav/Documents/ChemSims'));
% Script to search for parameters that give adaptation



% Linear, implicit BC system
model = @net1sim_lin; 
N = 3; 
X0 = [100 100 100];
mode = 'ODE_afterSS';
tspan = [0 100];


paramSetSize = 10^4;
paramRange = 6; % number of orders or mag above and below 1.

% Structure of the output file: 
% 10, params, 2 input metrics (input to ss, deltaInput), 4 output metrics (Cstart, Cpeak, Cpeaktime, Cend) 
out = zeros(paramSetSize, 14); 
telapsed = zeros(paramSetSize,1);
parfor iter = 1:paramSetSize

% Choose rand params
p = zeros(10,1);
p(1:7) = 2.^(paramRange.*rand(7,1)-paramRange/2);
p(8:10) = [400 400 400]; % in the future make these variable as well

% Simulate
start = tic;
[T, X] = model(p, X0, mode, tspan);
telapsed(iter) = toc(start);
% Store output statistics
Cstart = X(1,3);
[Cmax maxInd] = max(X(:,3));
CpeakTime = T(maxInd);
Cend = X(end, 3);

out(iter, :) = [p', Cstart, Cmax, CpeakTime, Cend]; 
end

hist(telapsed)
save('ParamSearchOutSet_1_2.mat', 'out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

