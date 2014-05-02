function quantifyNoise(D)

figure();
outputs = {'sensitivity','precision','noise'};

tau = 0.0001;

net = Net1LinModel();
IO = 1;
deltaIO = 0.2;

D.initMeanNoise = zeros(921,3);
D.peakMeanNoise = zeros(921,3);
D.ssMeanNoise = zeros(921,3);
for i = 1:921
    net.changeParams(D.p(i,:));
    net.initializeInput(IO);
    [T,X] = net.simulate('ODE',deltaIO);
    C = X(:,3);
    [Cpeak, CpeakTime] = max(C);

    D.initMeanNoise(i,:) = sum(net.stoich_matrix.*repmat(sqrt(net.prop_fcn(X(1,:))*tau),[1,3]));
    D.peakMeanNoise(i,:) = sum(net.stoich_matrix.*repmat(sqrt(net.prop_fcn(X(CpeakTime,:))*tau),[1,3]));
    D.ssMeanNoise(i,:) = sum(net.stoich_matrix.*repmat(sqrt(net.prop_fcn(X(end,:))*tau),[1,3]));

end

figure();
outputs = {'Init Noise','Peak Noise','SS Noise'};
CNoise =abs([D.initMeanNoise(:,3) D.peakMeanNoise(:,3) D.ssMeanNoise(:,3)]);
for i = 1:1:6
    for j = 1:1
        subplot(1,6,(j-1)*6+i);
        
        figure();
        scatter(D.outputStats(:,1), CNoise(:,2)./D.outputStats(:,1))
        title('Peak Noise')
        
             figure();
        scatter(CNoise(:,1), CNoise(:,2))
        title('init Noise')   
        
        %scatter(log10(D.p(:,i+1)), CNoise(:,j))
        %xlim([-3,3])
        if i==1; ylabel(outputs{j}); end
        if j==1; title(D.headings{i+1}); end;
    end   
end


