figure()
Cpeak = D.peakEnsP1(:,3);
repeats = 10;
z = zeros(length(1:10:10000),repeats);
sampSize = 1:10:10000;
for s=1:length(sampSize)
    for r=1:10
        z(s,r) = mean(randsample(Cpeak,sampSize(s)));
    end
end
xvals = repmat([1:10:10000]', 1, repeats);
scatter(xvals(:), z(:));

hold on;

Cpeak = L.peakEnsP1(:,3);
repeats = 1;
z = zeros(length(1:10:5000),repeats);
sampSize = 1:10:5000;
for s=1:length(sampSize)
    for r=1:1
        z(s,r) = mean(randsample(Cpeak,sampSize(s)));
    end
end
xvals = repmat([1:10:5000]', 1, repeats);
scatter(xvals(:), z(:),'ro');