addpath(genpath('C:\Users\Judy\Documents\'));
% Script to search for parameters that give adaptation



% Linear, implicit BC system
model = @net1sim_lin; 
N = 3; 
X0 = [100 100 100];
mode = 'ODE_afterSS';
tspan = [0 100];

p = zeros(10,1);
paramSetSize = 200;
paramRange = 5; % number of orders or mag above and below 1.

out = zeros(paramSetSize, length(p)+ 4); % 4 output metrics, start, peak, peaktime, end 
telapsed = zeros(paramSetSize,1);
for iter = 1:paramSetSize
    disp(iter);
% Choose rand params
%p(1:7) = 2.^(10.*rand(7,1)-5);
p(1:7) = 2.^(2.*rand(7,1)-1);
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
save('ParamSearchTest.mat', 'out');