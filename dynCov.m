
function dY = dynCov(y, p, u, m)
C_AA = y(1);
C_AB = y(2);
C_AC = y(3); 
C_BB = y(4);
C_BC = y(5);
C_CC = y(6);

mA = m(1);
mB = m(2);
mC = m(3);

c1 = p(1);
c2 = p(2);
c3 = p(3);
c4 = p(4);
c5 = p(5);
c6 = p(6);
c7 = p(7);
AT = p(8);
BT = p(9);
CT = p(10);

dY(1) = 2.*C_AA.*(c2+c1.*u)+mA.*(c2-c1.*u)+c1.*u.*AT;
dY(2) = C_AB.*(c1.*u+c2+c3.*mA+c4)-c3.*(BT-mB).*C_AA;
dY(3) = C_AC.*(c1.*u+c2+c5.*mA+c6.*mB+c7)-c5.*(CT-mC).*C_AA+c6.*mC.*C_AB;
dY(4) = -2.*c3.*C_AB.*(1+BT-mB)+2.*c3.*mA.*C_BB+c3.*mA.*(BT-mB)+c4.*mB;
dY(5) = C_BC.*(c3.*mA+c4+c5.*mA+c6.*mB+c7)-C_AC.*c3.*(BT-mB)-c5.*(CT-mC).*C_AB+c6.*mC.*C_BB;
dY(6) =2.*C_CC.*(c5.*mA+c6.*mB+3/4.*c7)+C_AC.*(2.*c5.*(CT-mC)-c5) +C_BC.*(2.*c6.*mC+c6)+c5.*mA.*(CT-mC)+c6.*mB.*mC+c7.*mC;



