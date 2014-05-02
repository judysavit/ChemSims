dataPath = 'Data Lin Impl BC Data\';
homePath = 'C:\Users\Judy\Documents\Noise Exploration\ChemSims\Models\';
cd([homePath dataPath]);

D=importdata('ParamSearch_c1_constant.mat');

%% Find sensitivity and precision and add to Data File
if ~ismember(D.headings,'Sensitivity')
    D.headings{end+1}='Sensitivity';
    D.headings{end+2}='Precision';
end
nSims = size(D.paramSet,1);
D.dataProc = [D.dataProc(:,1:6) zeros(nSims,2)];

A.headings = D.headings;
A.headingsDataProc = D.headingsDataProc;
A.paramSet = [];
A.dataProc = [];
A.adaptiveInds = zeros(nSims,1);

for i = 1:nSims
   Cstart = D.dataProc(i,strcmp(D.headingsDataProc,'Cstart'));
   Cmax = D.dataProc(i, strcmp(D.headingsDataProc,'Cmax'));
   Cend = D.dataProc(i, strcmp(D.headingsDataProc,'Cend'));
   IO = D.dataProc(i, strcmp(D.headingsDataProc,'input'));
   deltaIO = D.dataProc(i, strcmp(D.headingsDataProc,'deltaInput'));
   [A.adaptiveInds(i) S P ]=isAdaptive(Cstart,Cmax,Cend,IO,deltaIO);
   D.dataProc(i,end+1) = S;
   D.dataProc(i,end+2) = P;
end

%% Save as substruct with todays date and data generation methodology
D.A_041414.inds = A.adaptiveInds;
A.fxn = 'Z = deltaIO./IO; S = log10((abs(Cpeak-C0)./C0)./Z); P = log10(Z.*C0)-log10(abs(Css-C0)); if S>-0.5 && P>1; bool = true; else bool =false; end varargout{1} = S; varargout{2} = P;';
D.A_041414.fxn = A.fxn;
D.A_041414.IO = [IO deltaIO]; 
save('ParamSearch_c1_constant.mat','D');
clear D;

%% Explore relationship between parameters, sensitivity, and precision
A.paramSet = D.paramSet(A.adaptiveInds);
A.dataProc = D.dataProc(A.adaptiveInds);
figure();
outLabels = {'Sensitivity','Precision','PeakTime'};
outInds = [7 8 5];
for out = 1:3
    for c = 7:-1:1
        subplot(3,7,(out-1)+c)
        xvals = A.paramSet(:,c);
        yvals = A.dataProc(:,outInds(out));
        scatter(xvals,yvals);
        if (out==3); title(A.headings{c});end;
    end
    ylabel(outLabels{out});
end

A.date = '041414';
save('Adaptive_c1_constant.mat','A');