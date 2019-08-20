

%% Additional testing of individual changepoint/oddball coefficients -- as 
%  suggested by elife reviewers:

% This figure works with all data as long as we use new method for doing
% stats (weighted regression). 

xLabs={'Changepoint', 'Oddball'}
% Model to test learning in the individual conditions:
allRelCoeffs=[relROI.ROI_params_sepConditions(:,selROIs,5), relROI.ROI_params_sepConditions(:,selROIs,6)]
allRelCoeffs=allRelCoeffs(goodSub,:);

ms=12

% Columns:
% 1) Early CP, 
% 2) Late  CP, 
% 3) Early ODD,
% 4) Late  ODD,

[h,p,b,stats]=ttest(allRelCoeffs)
        meanBeta=nanmean(allRelCoeffs)
        semBeta =nanstd(allRelCoeffs)./sqrt(length(allRelCoeffs));

for i = 1:length(p)
    disp(sprintf('Coefficient stats for %s: p=%g, t=%g, mean=%g, SEM=%g', xLabs{i},p(i), stats.tstat(i), meanBeta(i), semBeta(i)))
end

relROI.peakTime(selROIs)
Scale=max(abs(allRelCoeffs(:)));
%xJit = smartJitter(allRelCoeffs,.03,.1);
xJit = smartJitter(allRelCoeffs,.075,.035);
% is there any covariance among parameter estimates?
[r, p]=corr(allRelCoeffs)
ll=size(allRelCoeffs, 1);
plot([0, size(allRelCoeffs, 2)], [0 0], '--k')
nROI=length(selROIs);
hold on
for i = 1:size(allRelCoeffs, 2)
    if i ==1
        ff(i)=plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors(1,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
    else
        fff(i)= plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors((2),:), 'markerEdgeColor', 'k', 'lineWidth', 1);
    end
end
ylim([-Scale, Scale]);
xlim([.5, size(allRelCoeffs, 2)+.5]);

[H,P] = kstest(zscore(allRelCoeffs(:,1)))

hold on
set(gca, 'fontsize', 16)
plot([1:length(b); 1:length(b)], b, 'k', 'lineWidth', 3, 'color', [.5, .5, .5])
plot([(1:length(b)) - .4;(1:length(b)) + .4 ], [nanmean(b); nanmean(b)], 'k', 'color', [.5 .5 .5], 'lineWidth', 7)
ylabel('Coefficient')
set(gca, 'xtick', 1:length(xLabs)*2, 'xticklabel', xLabs, 'box', 'off')
set(gca, 'box','off')
saveas(gcf, 'sepConditionEEG_regs.eps', 'epsc2')
close all




