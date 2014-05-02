load('net1_langevin_full.mat');
nP = 243;
k1 = D.p(1:nP,1);
k2 = D.p(1:nP,2);
k3 = D.p(1:nP,3);
k4 = D.p(1:nP,4);
k5 = D.p(1:nP,5);
k6 = D.p(1:nP,6);
k7 = D.p(1:nP,7);
AT = D.p(1:nP,8);
BT = D.p(1:nP,9);
CT = D.p(1:nP,10);

Ass = k1.*D.IO.*AT./(k2+k1.*D.IO);
Bss = k3.*Ass.*BT./(k4 + k3.*Ass);
Css = k5.*Ass.*CT./(k6.*Bss + k7 + k5.*Ass);

fCA = k5.*(CT-Css);
fBB = -k3.*Ass-k4;
fCB = -k6.*Css;
fBA = k3.*(BT-Bss);
fAA = -k1.*D.IO - k2;
fCC = -k5.*Ass -k6.*Bss-k7;

RajParam = fCA.*fBB-fCB.*fBA;
FreqParam = 1./sqrt(fCC.*fAA);

figure()
scatter(log(RajParam), D.outputStats(:,6));
title('Precision versus Rajs Param Combo')
xlabel('LOG f_{CA}f_{BB}-f_{CB}f_{BA}')
ylabel('Precision')

figure()
scatter(RajParam, D.outputStats(:,5));
title('Sensitivity versus Rajs Param Combo')
xlabel('f_{CA}f_{BB}-f_{CB}f_{BA}')
ylabel('Senstivity')


figure()
scatter(FreqParam, D.outputStats(:,3));
title('Precision versus Rajs Param Combo')
xlabel('\frac{1}{sqrt{f_{CC}f_{AA}}}')
ylabel('Peak Time')
