D=importdata('ParamSearchOutSet_1_12.mat');
nSims = size(D,1);
out = [D zeros(nSims,2)];

for i = 1:nSims
   C0 = D(i, 11);
   Cpeak = D(i, 12);
   Css = D(i, 14);
   [B S P ]=isAdaptive(C0,Cpeak,Css,IO,deltaIO);
   D(i,15) = S;
   D(i,16) = P;
end