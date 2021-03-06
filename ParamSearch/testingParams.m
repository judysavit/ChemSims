
netType = @Net1LinModel;
paramRange = 4; % number of orders or mag above and below 1.
net = Net1LinModel();
p = zeros(1,10);
p(1)=1;
p(8:10)=[1000 1000 1000];
CH = figure();
AH = figure();

rajParams = importdata('paramSets_higherTotals.csv')';
for j = 1:15
    %p(2:7) = 10.^(paramRange.*rand(1,6)-(paramRange-2));
    p = rajParams(j,:);

    net.changeParams(p);

    for i = -1.5:0.5:1.5

        input = 10;
        deltaInput = input;
        net.initializeInput(input);
        [T,X] = net.simulate('ODE', deltaInput);
        A = Net1Lin_analyticalSolnA(p, input, deltaInput);
        disp([i max(X(:,3))]);
                  
        cI = str2double(sprintf('%.2f',map(i, -1.5, 1.5, 0, 1)));

        colorA = [cI, 0, 0];
        colorC = [0, ternif(cI>1,1,cI), 0 ];
        
        figure()
        %subplot(3,5,j);
        plot(T,X(:,1), 'Color',colorA); hold on;
        %[AX,H1,H2]=plotyy(T,X(:,3),T,A(T)); hold on;
        %set(H1,'Color',colorA)
        %set(H2,'Color',colorC)
        xlim([0 50])
        ylim([0 1000])
        
        figure()
        %subplot(3,5,j);
        plot(T,X(:,3), 'Color',colorC); hold on;
        xlim([0 50])
        ylim([0 1000])
    end
end

for j = 1:15
    subplot(3,5,j);
            ylim([0 300])
end

j=13;
p = rajParams(j,:);
net.changeParams(p);
net.initializeInput(input);
for i = -1.5:0.5:1.5
    deltaInput = 10^i;
    [T,X] = net.simulate('ODE', deltaInput);
    disp([i max(X(:,3))]);
end
%%
p = pADAPT;
figure()
for i = -1.4:0.2:1.4
    cI = str2double(sprintf('%.2f',map(i, -1.4, 1.4, 0, 1)));
    input =10^i;
    deltaInput = input;

    net.initializeInput(input, deltaInput);
    [T,X] = net.simulate('ODE');
    
    A = Net1Lin_analyticalSolnA(p,input, deltaInput);
    subplot(3,6,1+floor(cI*15));
    colorA = [cI, 0, 0];
    % plot(T,A(T),'Color',colorA);hold on;
    colorC = [0, ternif(cI>1,1,cI), 0 ];
    plot(T,X(:,1),'Color',colorC,'Linestyle', '-'); hold on;
    xlim([0 30])
end
