%function quantifyNoiseSSA(D)
nP = size(D.p,1);
IO = 1; D.IO = IO;
deltaIO = 1; D.deltaIO = deltaIO;
net = Net2LinModel(D.p(1,:), IO, deltaIO);
if isfield(D,'outputStats'); D = rmfield(D,'outputStats');end;

% 1st dimension is the param set
% 2nd dimension is the measure of the distribution (mean, std dev,
% kurtosis)
% 3rd dimension is the species A, B, C
D.initNoise = zeros(nP,3,3);
D.peakNoise = zeros(nP,3,3);
D.ssNoise = zeros(nP,3,3);

nEns =  8000;
%blockSize = 10; IN UNITS OF TIME, NOT SIMULATION STEPS. MIGHT ADD THIS, RAJ SUGGESTED LOOKING AT NOISE AROUND NOT
%JUST AT PEAK, NOT SURE ITS KOSHER SINCE THESE DISTRIBUTIONS ARE TIME
%CORRELATED.
D.telapsed = zeros(1,nP);
for i = 1:nP
    tic;
    net.changeParams(D.p(i,:));
    net.initializeInput(IO, repmat(D.p(i,end),1,3));
    net.delayBeforeStim = 0.025;
    
    % Save output stats for new deltaIO = 1
    [TO,XO]= net.simulate('ODE',deltaIO, 0.8);
    C0 = net.X0(3);
    Css = XO(end,3);
    [Cpeak, maxCInd] = max(XO(:,3));
    maxCTime = TO(maxCInd);
    
    Z = deltaIO./IO;
    S = log10((abs(Cpeak-C0)./C0)./Z);
    P = log10(Z.*C0)-log10(abs(Css-C0));
    
    D.outputStats(i,:) = [C0 Cpeak maxCTime Css S P];
    
    % Simulate with ssa
    peak = -1*ones(nEns,3); 
    ss = -1*ones(nEns,3); 
    parfor j = 1:nEns
        
        [T,X] = net.simulate('langevin',deltaIO, 0.8);
        peakInd = find(T<maxCTime,1,'last');
        peak(j,:) = X(peakInd,:);
        ss(j,:) = X(end,:);
    end

    D.peakEnsP1 = peak;
    D.ssEnsP1 = ss;
    
    D.peakNoise(i,1,:) = mean(peak);
    D.peakNoise(i,2,:) = std(peak);
    D.peakNoise(i,3,:) = kurtosis(peak);
    
    D.ssNoise(i,1,:) = mean(ss);
    D.ssNoise(i,2,:) = std(ss); 
    D.ssNoise(i,3,:) = kurtosis(ss);
    D.telapsed(i)= toc;
    save('net2_langevin_full.mat','D');
    
end




