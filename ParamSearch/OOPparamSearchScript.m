
function OOPparamSearchScript(homePath)
% Find parameters that do adaptivity
% Only input arg is the path for ChemSims

addpath(genpath(homePath));

netType = @Net1LinModel;
folder = 'Data Lin Impl BC Data';

clear out telapsed;
input = 1;
deltaInput = 0.2;
paramSetSize = 3;% 10^4;
paramRange = 5; % number of orders or mag above and below 1.

% Structure of the output file: 
% 10, params, 2 input metrics (input to ss, deltaInput), 4 output metrics
% (Cstart, Cpeak, Cpeaktime, Cend) 

% First keep c1=1
out = zeros(paramSetSize, 16); 
telapsed = zeros(paramSetSize,1);

parfor iter = 1:paramSetSize
    p = zeros(1,10);
    p(1) = 1;
    p(2:7) = 10.^(paramRange.*rand(1,6)-(paramRange-2));
    p(8:10) = [1000 1000 1000];
    net = Net1LinModel(p);
    net.initialize(input, deltaInput);
    start = tic;
    [T X] = net.simulate('ODE');
    telapsed(iter) = toc(start);
    
    Cstart = X(1,3);
    [Cmax maxInd] = max(X(:,3));
    CpeakTime = T(maxInd);
    Cend = X(end, 3);
    out(iter, :) = [net.params, net.input, net.deltaInput, Cstart, Cmax, CpeakTime, Cend];
end
% Save data
D.paramSet = out;
D.headings = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'input', 'deltaInput','Cstart','Cmax','CpeakTime','Cend'};
D.runTimes = telapsed;
save([homePath 'Models\' folder '\ParamSearch_c1_constant.mat'],'D')  


% Let c1 vary
clear out telapsed;

out = zeros(paramSetSize, 16); 
telapsed = zeros(paramSetSize,1);

parfor iter = 1:paramSetSize
    p = zeros(1,10);
    p(1:7) = 10.^(paramRange.*rand(1,7)-(paramRange-2));
    p(8:10) = [1000 1000 1000];
    net = Net1LinModel(p);
    net.initialize(input, deltaInput);
    start = tic;
    [T X] = net.simulate('ODE');
    telapsed(iter) = toc(start);
    
    Cstart = X(1,3);
    [Cmax maxInd] = max(X(:,3));
    CpeakTime = T(maxInd);
    Cend = X(end, 3);
    out(iter, :) = [net.params, net.input, net.deltaInput, Cstart, Cmax, CpeakTime, Cend];
end
% Save data
D.paramSet = out;
D.headings = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'input', 'deltaInput','Cstart','Cmax','CpeakTime','Cend'};
D.runTimes = telapsed;
save([homePath 'Models\' folder '\ParamSearch_c1_variable.mat'],'D')  

end


