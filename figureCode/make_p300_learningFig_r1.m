% make_p300_learningFig


num        = 1;
wid        = 17.6; % total width
hts        = [6, 6];
cols       = {3, 3};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [], [], [], '');
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

selROIs=1;
allRelCoeffs=[relROI.ROI_params(:,selROIs,5), relROI.ROI_params(:,selROIs,6)];
allRelCoeffs=allRelCoeffs(goodSub,:);
if isfield(relROI, 'ROI_params_bootStd')
allRelSTDs=[relROI.ROI_params_bootStd(:,selROIs,5), relROI.ROI_params_bootStd(:,selROIs,6)];
allRelSTDs=allRelSTDs(goodSub,:)
end

xMat=[ones(sum(goodSub),1), zscore(contextualSurpriseEffectOnLR(goodSub))];

% Do this for stats:
clear B BINT allStats
for i =1:size(allRelCoeffs, 2)
    [B(:,i),BINT(:,:,i),R,RINT,STATS] = regress(allRelCoeffs(:,i), xMat);
    allStats(i)=regstats(allRelCoeffs(:,i), xMat(:,2));
end

disp('Main Effect of EEG signal -- intercept:');
disp(sprintf('Mean Beta: %g \n95 conf: %g / %g', B(1,1), BINT(1,1,1), BINT(1,2,1)));
disp(sprintf('t stat: %g ', allStats(1).tstat.t(1)));
disp(sprintf('pvalue: %g ', allStats(1).tstat.pval(1)));

disp('Main Effect of EEG signal -- ind diffs:');
disp(sprintf('Mean Beta: %g \n95 conf: %g / %g', B(2,1), BINT(2,1,1), BINT(2,2,1)));
disp(sprintf('t stat: %g ', allStats(1).tstat.t(2)));
disp(sprintf('pvalue: %g ', allStats(1).tstat.pval(2)));

disp('Interaction Effect of EEG signal -- intercept:');
disp(sprintf('Mean Beta: %g \n95 conf: %g / %g', B(1,2), BINT(1,1,2), BINT(1,2,2)));
disp(sprintf('t stat: %g ', allStats(2).tstat.t(1)));
disp(sprintf('pvalue: %g ', allStats(2).tstat.pval(1)));

disp('Interaction Effect of EEG signal -- ind diffs:');
disp(sprintf('Mean Beta: %g \n95 conf: %g / %g', B(2,1), BINT(2,1,2), BINT(2,2,2)));
disp(sprintf('t stat: %g ', allStats(2).tstat.t(2)));
disp(sprintf('pvalue: %g ', allStats(2).tstat.pval(2)));


% Redo, without z-scoring for trendlines:

% Do this for stats:
xMat=[ones(sum(goodSub),1), contextualSurpriseEffectOnLR(goodSub)];
clear allB
for i =1:size(allRelCoeffs, 2)
    [B, BINT] = regress(allRelCoeffs(:,i), xMat);
    allB(i,:)=B;
end

doErrBars=false
ms=8
ms2=6
xLims=[-100, 700];
erpYLims=[-3, 5];
erdYLims=[-2, 2];
trendLineColor=[.7,.7,.7];
xLims=[-100, 700];
erpYLims=[-3, 5];
erdYLims=[-2, 2];
ms=8;
ec=3;
lc=2;


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        
        hold on
        plot([-1, 1], [-1, 1], '--k')
        plot([-1, 1], [-.4, .4], 'color', cbColors(2,:))
        plot([-1, 1], [-.6, .6], 'color', cbColors(1,:))
        plot([-1, 1], [-0, 0], '--k')
        set(gca, 'box','off', 'xtick', [], 'ytick', [])
        ylabel('Update')
        xlabel('Prediction error')
        
    elseif xx==2
        
        
        hold on
        
        %   plot([-1, 1], [-1, 1], '--k')
        plot([-1, 1], [-.6, .6]+.05, 'color', cbColors(2,:), 'lineWidth', 3)
        plot([-1, 1], [-.6, .6], 'color', cbColors(1,:), 'lineWidth', 3)
        %   plot([-1, 1], [-0, 0], '--k')
        
        set(gca, 'box','off', 'xtick', [], 'ytick', [])
        ylabel('Learning rate')
        xlabel('Signal strength')
        title('Positive PE*EEG Signal effect')
        
        
    elseif xx==3
        
        hold on
        plot([-1, 1], [0, .6]+.05, 'color', cbColors(1,:), 'lineWidth', 3)
        plot([-1, 1], [0, -.6], 'color', cbColors(2,:), 'lineWidth', 3)
        set(gca, 'box','off', 'xtick', [], 'ytick', [])
        ylabel('Learning rate')
        xlabel('Signal strength')
        title('Positive PE*EEG Signal*condition effect')
        
        
        
    elseif xx==4
        
        f=1
        disp(['getting ' ROIs_forSingleTrialAnalysis{f} ])
        eval(['relROI=allROIs.' ROIs_forSingleTrialAnalysis{f} ';']) ;
        xLabs={'Direct', 'Conditional'};
        %xLabs=(repmat(relROI.peakTime(selROIs), 1, 2));
        
        
        allCoeffsForYHat=squeeze(relROI.ROI_params(:,selROIs,:));
        
        % MAKE learning rate versus EEG signal figure:
        exSigs=-4:.1:4;
        conds=[0, 1]-.5;  % -.5 = oddball, +.5 =CP
        intCol=3;                 % PE term from model.
        condCol=4;
        totLearnCol=5;
        intLearnCol=6;
        
        
        % Now we are just plotting one coefficient:
        clear predLR
        relCoeffs=allRelCoeffs;
        % loop through subjects:
        for i = 1:length(allRelCoeffs)
            
            subCoefs=allCoeffsForYHat(i,:);
            % loop through oddball/CP conditions:
            for j = 1:length(conds)
                
                % predLR  =
                predLR(i,j,:)=subCoefs(intCol) +  ...  intercept
                    conds(j)*subCoefs(condCol)   +  ...  condition
                    exSigs.*subCoefs(totLearnCol)+  ...  phys signal
                    exSigs.*subCoefs(intLearnCol)*conds(j); % phys signal* cond
            end
        end
        
        
        
        % 1) direct effect on learning, 2) interaction learning*condition
        
        
        [h, p, b, stats]=ttest(allRelCoeffs)
        meanBeta=nanmean(allRelCoeffs)
        semBeta =nanstd(allRelCoeffs)./sqrt(length(allRelCoeffs));
        for i = 1:length(p)
            disp(sprintf('Coefficient stats for %s: p=%g, t=%g, mean=%g, SEM=%s', xLabs{i},p(i), stats.tstat(i), meanBeta(i), semBeta(i)))
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
                    if i <=nROI
                        plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors(4,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                    else
                        plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                    end
                end
                %ylim([-Scale, Scale]);
                xlim([.5, size(allRelCoeffs, 2)+.5]);
          end
        %ylim([-Scale, Scale]);
        xlim([.5, size(allRelCoeffs, 2)+.5]);
        
        
        hold on
        plot([1:length(p); 1:length(p)], b, 'k', 'lineWidth', 5, 'color', [.5, .5, .5])
        plot([(1:length(p)) - .4; (1:length(p)) + .4 ], [meanBeta; meanBeta], 'k', 'color', [.5 .5 .5], 'lineWidth', 7)
        ylabel('Coefficient')
        set(gca, 'xtick', 1:length(xLabs), 'xticklabel', xLabs, 'box', 'off')
        set(gca, 'box','off')
        ylim([-.4, .4])
        
        % June 18, 2019 -- Get individual CP/Oddball stats to address
        % reviewer concerns:
        
    elseif xx==5
        
        
        bScale=[min(contextualSurpriseEffectOnLR(goodSub))-.25, max(contextualSurpriseEffectOnLR(goodSub))+.25]
        xVals=bScale;
        tLine=allB(2,1)+xVals.*allB(2,2)
        hold on
        plot(bScale, [0, 0], '--k')
        plot(xVals, tLine, 'color', trendLineColor, 'lineWidth', 2)
        
        if doErrBars
            plot([contextualSurpriseEffectOnLR(goodSub), contextualSurpriseEffectOnLR(goodSub)]', [allRelCoeffs(:,3)+ allRelSTDs(:,3), allRelCoeffs(:,3)- allRelSTDs(:,3)]', '-k',  'lineWidth', 1)
        end
        
        [r, p]=corr(contextualSurpriseEffectOnLR(goodSub), allRelCoeffs(:,2))
        
        
        plot(contextualSurpriseEffectOnLR(goodSub), allRelCoeffs(:,2), 'o', 'markerFaceColor', cbColors(ec,:), 'markerEdgeColor', 'k', 'markerSize', ms, 'lineWidth', 1)
        ylabel('Conditional learning EEG coefficient')
        xlabel('Behavioral surprise*condition coefficient')
        ylim([-.25, .5])
        xlim([-.5, 2])
        set(gca, 'box', 'off')
        
        
        
        
    elseif xx==6
        
        % MRN adding "raw" learning rates for binned data:
        %medValOfBin
        if exist('ROI_binnedLR')
            rawRegLRs_odd=squeeze(squeeze(ROI_binnedLR(:,selROIs, :, 1)));
            rawRegLRs_CP=squeeze(squeeze(ROI_binnedLR(:,selROIs, :, 2)));
            rawRegLRs_CP_ste_CP=nanstd(rawRegLRs_CP)./sqrt(size(rawRegLRs_CP, 1));
            rawRegLRs_odd_ste_odd=nanstd(rawRegLRs_odd)./sqrt(size(rawRegLRs_odd, 1));
        end
        
        
        relROI_peaks=relROI.peakTime(selROIs);
        cpCol=1;
        odCol=2;
        predToPlotCP=squeeze(predLR(:,2,:));
        predToPlotOdd=squeeze(predLR(:,1,:));
        SEM_cp=nanstd(predToPlotCP)./sqrt(length(predToPlotCP));
        SEM_odd=nanstd(predToPlotOdd)./sqrt(length(predToPlotOdd));
        hold on
        
        H1=shadedErrorBar(exSigs,nanmean(predToPlotCP), SEM_cp, {'-', 'color', cbColors(cpCol,:)}, true)
        H2=shadedErrorBar(exSigs,nanmean(predToPlotOdd), SEM_odd, {'-', 'color', cbColors(odCol,:)}, true)
        
        % plot it if we've got it!
        if exist('rawRegLRs_odd')
            % plot "rawer" CP data
            plot([nanmean(medValOfBin); nanmean(medValOfBin)], [nanmean(rawRegLRs_CP(:,:,1)) + rawRegLRs_CP_ste_CP(:,:,1); nanmean(rawRegLRs_CP(:,:,1)) - rawRegLRs_CP_ste_CP(:,:,1)], 'k')
            plot([nanmean(medValOfBin)], [nanmean(rawRegLRs_CP(:,:,1))], 'o', 'markerFaceColor', cbColors(cpCol,:), 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms2)
            % plot "rawer" oddball data
            plot([nanmean(medValOfBin); nanmean(medValOfBin)], [nanmean(rawRegLRs_odd(:,:,1)) + rawRegLRs_odd_ste_odd(:,:,1); nanmean(rawRegLRs_odd(:,:,1)) - rawRegLRs_odd_ste_odd(:,:,1)], 'k')
            plot([nanmean(medValOfBin)], [nanmean(rawRegLRs_odd(:,:,1))], 'o', 'markerFaceColor', cbColors(odCol,:), 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms2)
        end
        
        
        ylabel('Learning rate')
        xlabel('Signal strength')
        ff=legend([H1.mainLine, H2.mainLine], 'Changepoint', 'Oddball')
        set(ff, 'location', 'southwest', 'box', 'off')
       
        set(gca, 'box','off')
    end
    
    
    %setPLOT_panelLabel(gca, xx);
end

% kk=annotation('textbox')
% set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')

saveas(gcf,  'p300_learningFig_noExclude_wIndDiffs_altHP_small.fig', 'fig')
saveas(gcf,  'p300_learningFig_noExclude_wIndDiffs_altHP_small.eps', 'epsc2')
close(gcf)


