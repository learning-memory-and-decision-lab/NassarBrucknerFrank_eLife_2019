% make_p300_learningFig


num        = 1;
wid        = 17.6; % total width
hts        = [6];
cols       = {2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [3], [3], [], '');
set(axs,'Units','normalized');

% draw in each panel, one at a time

% a = schematic for analysis
% b = coefficients for subjects in phys*pe model
% c = yHat from that model across conditions and phys strength (roi 1)
% d = yHat from that model across conditions and phys strength (roi 2)

% e = schematic for residual analysis:
% f = correlation with lr-hat by condition
% g = coefficients for phys*pe residual model
% h = individual differences in ROI 1


%% second level regression:
% alternative -analysis second level regression?

        
allRelCoeffs=[relROI.ROI_params(:,selROIs,5), relROI.ROI_params(:,selROIs,6)];
allRelCoeffs=allRelCoeffs(goodSub,:);
allRelSTDs=[relROI.ROI_params_wModTerms_bootStd(:,selROIs,5), relROI.ROI_params_wModTerms_bootStd(:,selROIs,6)];
allRelSTDs=allRelSTDs(goodSub,:)

xMat=[ones(sum(goodSub),1), zscore(contextualSurpriseEffectOnLR(goodSub)), zscore(surpriseEffectOnLR(goodSub))];

% Do this for stats:
for i =1:4
   [B,BINT,R,RINT,STATS] = regress(allRelCoeffs(:,i), xMat);
   expSEM=allRelSTDs(:,i);
   [b, bint] = regressW_mrn(allRelCoeffs(:,i),expSEM,xMat);
    BINT;
    bint
end

% Redo, without z-scoring for trendlines:

% Do this for stats:
xMat=[ones(sum(goodSub),1), contextualSurpriseEffectOnLR(goodSub), surpriseEffectOnLR(goodSub)];
for i =1:4
   expSEM=allRelSTDs(:,i);
   [b, bint] = regressW_mrn(allRelCoeffs(:,i),expSEM,xMat);
   allB(i,:)=b;
end



ec=3;
lc=2;
ms=10;
trendLineColor=[.7, .7, .7];

for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1

        Scale=max(abs(surpriseEffectOnLR(goodSub)))+.25
        xVals=[-Scale, Scale];
        tLine=allB(1,1)+xVals.*allB(1,3)
        hold on
        plot([-Scale, Scale], [0, 0], '--k')
        plot(xVals, tLine, 'color', trendLineColor, 'lineWidth', 2)
        plot([surpriseEffectOnLR(goodSub), surpriseEffectOnLR(goodSub)]', [allRelCoeffs(:,1)+ allRelSTDs(:,1), allRelCoeffs(:,1)- allRelSTDs(:,1)]', '-k',  'lineWidth', 1)
        plot(surpriseEffectOnLR(goodSub), allRelCoeffs(:,1), 'o', 'markerFaceColor', cbColors(ec,:), 'markerEdgeColor', 'k', 'markerSize', ms, 'lineWidth', 1)
        %lsline
        ylabel('Direct learning EEG coefficient')
        xlabel('Behavioral surprise coefficient')
        set(gca, 'box', 'off')
        
        
    elseif xx==2
        
        Scale=max(abs(surpriseEffectOnLR(goodSub)))+.25
        xVals=[-Scale, Scale];
        tLine=allB(2,1)+xVals.*allB(2,3)
        hold on
        plot([-Scale, Scale], [0, 0], '--k')
        plot(xVals, tLine, 'color', trendLineColor, 'lineWidth', 2)
        plot([surpriseEffectOnLR(goodSub), surpriseEffectOnLR(goodSub)]', [allRelCoeffs(:,2)+ allRelSTDs(:,2), allRelCoeffs(:,2)- allRelSTDs(:,2)]', '-k',  'lineWidth', 1)
        plot(surpriseEffectOnLR(goodSub), allRelCoeffs(:,2), 'o', 'markerFaceColor', cbColors(lc,:), 'markerEdgeColor', 'k', 'markerSize', ms, 'lineWidth', 1)
        %lsline
        ylabel('Direct learning EEG coefficient')
        xlabel('Behavioral surprise coefficient')
        set(gca, 'box', 'off')
        

    end
    
    
    %setPLOT_panelLabel(gca, xx);
end

% kk=annotation('textbox')
% set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')





saveas(gcf,  'p300_indDiff.fig', 'fig')
saveas(gcf,  'p300_indDiff.eps', 'epsc2')
close(gcf)
