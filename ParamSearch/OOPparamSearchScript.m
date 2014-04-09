addpath(genpath('/home/jsav/Documents/ChemSims'));
addpath(genpath('C:\Users\Judy\Documents\Noise Exploration\ChemSims'));

paramSetSize = 10;%^4;
paramRange = 5; % number of orders or mag above and below 1.
out = zeros(paramSetSize, 14); % 4 output metrics, start, peak, peaktime, end 
telapsed = zeros(paramSetSize,1);

p = zeros(10,1); 
p(8:10) = [400 400 400]; % in the future make these variable as well

input = 1;
delInput = 0.2;

nl = Net1LinModel(p, input);
X0 = [100 100 100];
odeSim = ODESimulator();
for iter = 1:paramSetSize

% Choose rand params
p(1:7) = 2.^(paramRange.*rand(7,1)-paramRange./2);
nl.params = p;

% Simulate
start = tic;
X_ss = odeSim.initialize(nl, X0);
[T,X] = odeSim.simulate(nl, [1, 100], delInput, X_ss);
telapsed(iter) = toc(start);
% Store output statistics
Cstart = X(1,3);
[Cmax maxInd] = max(X(:,3));
CpeakTime = T(maxInd);
Cend = X(end, 3);

out(iter, :) = [p', Cstart, Cmax, CpeakTime, Cend]; 
end
