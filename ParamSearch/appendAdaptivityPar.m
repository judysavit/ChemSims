function [Dout, A] = appendAdaptivityPar(D)

%% Find sensitivity and precision and add to Data File
if ~ismember(D.headingsDataProc,'Sensitivity')
    D.headingsDataProc{end+1}='Sensitivity';
    D.headingsDataProc{end+2}='Precision';
end
nSims = size(D.paramSet,1);
D.dataProc = [D.dataProc(:,1:6) zeros(nSims,2)];

adaptiveInds = zeros(nSims,1);
sens = zeros(nSims,1);
prec = zeros(nSims,1);

CstartInd = D.dataProc(:,strcmp(D.headingsDataProc,'Cstart'));
CmaxInd = D.dataProc(:,strcmp(D.headingsDataProc,'Cmax'));
CendInd = D.dataProc(:,strcmp(D.headingsDataProc,'Cend'));
deltaIOInd = D.dataProc(:,strcmp(D.headingsDataProc,'deltaInput'));
IOInd = D.dataProc(:,strcmp(D.headingsDataProc, 'input'));

%for i = 1:5
for i = 1:nSims
   Cstart = CstartInd(i);
   Cmax = CmaxInd(i);
   Cend = CendInd(i);
   IO = deltaIOInd(i);
   deltaIO = IOInd(i);
   [adaptiveInds(i) S P ]=isAdaptive(Cstart,Cmax,Cend,IO,deltaIO);
   sens(i) = S;
   prec(i) = P;
end

D.dataProc(:,end-1)=sens;
D.dataProc(:,end)=prec;
A.adaptiveInds=adaptiveInds;
A.headings = D.headings;
A.headingsDataProc = D.headingsDataProc;

Dout = D;
end


