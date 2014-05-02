function sensPrecVsParamsFigs(D)
figure();
outputs = {'sensitivity','precision'};
for i = 1:1:6
    for j = 1:2
        subplot(2,6,(j-1)*6+i);
        scatter(log10(D.p(:,i+1)), D.outputStats(:,j+4))
        xlim([-3,3])
        if i==1; ylabel(outputs{j}); end
        if j==1; title(D.headings{i+1}); end;
    end   
end