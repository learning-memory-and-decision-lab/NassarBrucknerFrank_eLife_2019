% make_p300_resLearningFig_d2


num        = 1;
wid        = 17.6; % total width
hts        = [6];
cols       = {3};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [2], [2], [], '');
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

ms=8
ms2=8
xLims=[-100, 700];
erpYLims=[-3, 5];
erdYLims=[-2, 2];
trendLineCol=[.7, .7, .7];
relRange=[300, 700];


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        set(gca, 'Visible', 'off')
        
    elseif xx==2
        
        imagesc(downSampTimes(timesToLookAt), 1:64, allROIs.basicCoeffs_cpPlusOdd.fullTMap)
        xlim([0, 1000])
        ylim([0, 64])
        plot([relRange' relRange']', [0, 64; 0, 64]', '--k', 'lineWidth', 1) ;
        set(gca, 'ytick', [], 'yticklabels', '')
        xlabel('Time')
        ylabel('Channel')
        %colorbar
        title('CP+Oddball T map')
        set(gca, 'box', 'off')
        
        
    elseif xx==3
        
        
         
        hold on
        plot([meanTime(1), meanTime(end)], [0, 0], '--k')
        plot([relRange' relRange']', [-4, 4; -4, 4]', '--k', 'lineWidth', 1) ;

        H1=shadedErrorBar(meanTime(isDat),meanBeta1(isDat),semBeta1(isDat),{'color', cbColors(4,:)},false);
        ylim([-.06, .06]);
        set(gca, 'box', 'off');
        xlim([0, 1000]);
        
        hold on
        plot([meanTime(1), meanTime(end)], [0, 0], '--k');
        H2=shadedErrorBar(meanTime(isDat),meanBeta2(isDat),semBeta2(isDat),{'color', cbColors(3,:)},false);
        ylim([-.06, .06]);
        set(gca, 'box', 'off');
        xlim([0, 1000]);
        
        
        % Get stats:
        tStat=meanBeta2(isDat)./semBeta2(isDat);
        selTime=meanTime(isDat);
        peakTime=selTime( find(tStat==max(tStat))) % Peak Time
        meanBeta2(meanTime==peakTime)              % Peak Coeff
        semBeta2(meanTime==peakTime)               % Peak SEM
        
        
        
        
        plot(maskedTime(clustCorrSig2), maskedBeta2(clustCorrSig2), '.r');
        plot(maskedTime(clustCorrSig1), maskedBeta1(clustCorrSig1), '.r');

        
        
        %     hold on
        %     plot(meanTime(isDat), stats.tstat(isDat))
        %     ylim([-4, 4])
        %     plot([min(meanTime), max(meanTime)], [0, 0], '--k')
        ylabel('Residual learning coefficient');
        xlabel('Time');
        %title('Interaction')
        xlim([0, 1000]);
        
        
        
        ff=legend([H1.mainLine, H2.mainLine],  'Direct', 'Conditional');
        set(ff, 'box', 'off');
        set(gca, 'box', 'off')
    end
    
    
    %setPLOT_panelLabel(gca, xx);
end



saveas(gcf,  'p300_learningResFig_d2.fig', 'fig')
saveas(gcf,  'p300_learningResFig_d2.eps', 'epsc2')
close(gcf)

