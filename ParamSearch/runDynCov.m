i = 1;
p = D.p(i,:);
net = Net1LinModel(p, D.IO, D.deltaIO);
net.initializeInput();
tspan = [0 6];

u = net.input;

%C_AA C_AB C_AC C_BB C_BC C_CC constant
matforSS = [2.*(c2+c1.*u), 0, 0, 0, 0, 0;
-c3.*(BT-mB), (c1.*u+c2+c3.*mA+c4), 0, 0, 0, 0;
-c5.*(CT-mC), c6.*mC, (c1.*u+c2+c5.*mA+c6.*mB+c7), 0, 0, 0;
0, -2.*c3.*(1+BT-mB), 0, 2.*c3.*mA, 0,0;
0, -c5.*(CT-mC), -c3.*(BT-mB), c6.*mC, (c3.*mA+c4+c5.*mA+c6.*mB+c7), 0;
0, 0, 2.*c5.*(CT-mC)-c5, 0,(2.*c6.*mC+c6),2*(c5.*mA+c6.*mB+3/4.*c7)];
%Find steady state covariances after initialization

B = [mA.*(c2-c1.*u)+c1.*u.*AT;
    0;
    0;
    c3.*mA.*(BT-mB)+c4.*mB;
    0;
    c5.*mA.*(CT-mC)+c6.*mB.*mC+c7.*mC];
     
X = linsolve(matforSS,B);

[T,MX] = net.simulate('ODE', net.deltaInput, 6);
CY = zeros(length(T),6); 
CY(1,:) = COVINIT;
for t = 2:length(T)
    CY(t,:) = CY(t-1,:)+ dynCov(CY(t-1,:), p, D.IO+D.deltaIO, MX(t,:));
end