
%     D.peakNoise(i,1,:) = mean(peak);
%     D.peakNoise(i,2,:) = std(peak);
%     D.peakNoise(i,3,:) = kurtosis(peak);
%     
%     D.ssNoise(i,1,:) = mean(ss);
%     D.ssNoise(i,2,:) = std(ss); 
%     D.ssNoise(i,3,:) = kurtosis(ss);

nP = size(D.outputStats,1);

% SCATTER THE MEAN AND VARIANCE OF C PEAK AS A FUNCTION OF SENSITIVITY AND
% PRECISION
figure()

COV = D.peakNoise(1:nP,2,3)./D.peakNoise(1:nP,1,3);

subplot(2,2,1);
scatter(D.peakNoise(1:nP,1,3)-D.outputStats(:,1), D.outputStats(:,5));
title('Mean C at peak versus Sensitivity')
xlabel('Mean C at Peak');
ylabel('Sensitivity');

subplot(2,2,2);
scatter(COV, D.outputStats(:,5));
title('Noise at peak versus Sensitivity')
xlabel('COV (Noise)');
ylabel('Sensitivity');

subplot(2,2,3);
scatter(D.peakNoise(1:nP,1,3)-D.outputStats(:,1), D.outputStats(:,6));
title('Mean C at peak versus Precision')
xlabel('Mean C at Peak');
ylabel('Precision');

subplot(2,2,4);
scatter(COV, D.outputStats(:,6));
title('Noise at peak versus Precision')
xlabel('COV');
ylabel('Precision');

figure()
scatter(COV, D.peakNoise(1:nP,1,3)-D.outputStats(:,1));
title('Noise at peak versus Mean C at peak')
xlabel('COV');
ylabel('Mean C at Peak');

figure()
scatter(D.outputStats(:,5),(D.peakNoise(1:nP,1,3)-D.outputStats(:,1))./D.peakNoise(1:nP,2,3));
title('Signal:Noise at peak versus Sensitivity')
xlabel('Signal-to-noise ratio = (peak-C0)/stddev');
ylabel('Sensitivity');


figure()
for i = 1:6
    subplot(1,6,i);
    scatter(D.p(1:nP,i+1),D.peakNoise(1:nP,2,3));
   xlabel(D.headings(i+1));
      title(D.headings(i+1)); 
    if i==1
        ylabel('Noise (Std. Dev. of C peak)')
    end
end

