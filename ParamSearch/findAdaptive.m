
D=importdata('ParamSearchOutSet_1_12.mat');
IO  = 1;
deltaIO = 0.2;
nSims = size(D,1);
outAdapt = [];
for i = 1:nSims
   C0 = D(i, 11);
   Cpeak = D(i, 12);
   Css = D(i, 14);
    if (isAdaptive(C0,Cpeak,Css,IO,deltaIO))
        outAdapt = D(i, :);
    end
end