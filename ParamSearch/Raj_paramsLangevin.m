% SIMULATING RAJ'S PARAMETERS

D=importdata('psets-2014-04-02.csv')';
tspan=[0 10];
% Implicit BC
for i = 1:20%size(D,1)
figure();
p = pnew(i,:);
model = @net1sim_lin; 
X0 = [100 100 100];
leg= {'A','B','C'};

sims = {'ODE'};
for s = 1:length(sims)
    
    mode = sims{s};
    [T, X] = model(p, X0, mode, tspan);
    % Plot time course
    subplot(1,length(sims),s)
    title(mode);
    plot(T,X); hold on;
    xlabel('Time (s)');
    if s ==1
        ylabel('Molecules');
    end
    
end
end

c1 = 1;
c2 = p(:,2);
c3 = p(:,3);
c4 = p(:,4);
c5 = p(:,5);
c6 = p(:,6);
c7 = p(:,7);
AT = 400;
BT = 400;
CT = 400;
Ass = AT./(c2+1);
Bss = c3.*Ass.*BT./(c4+c3.*Ass);
Css = c5.*Ass.*CT./(c6.*Bss+c7+c5.*Ass);




pnew = repmat(D(2,:),20,1);
pnew(:,2)=linspace(0.01,1,20);




