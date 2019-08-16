% make_p300_supplearningFig


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



xLims=[-100, 700];
erpYLims=[-3, 5];
erdYLims=[-2, 2];
ms=8


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
        
        
        
    elseif xx==4
        
        f=2
        disp(['getting ' ROIs_forSingleTrialAnalysis{f} ])
        eval(['relROI=allROIs.' ROIs_forSingleTrialAnalysis{f} ';']) ;
        selROIs=[1:2]; % if f = 2
        
        
        
        xLabs=(repmat(relROI.peakTime(selROIs), 1, 2))

        
        % MAKE learning rate versus EEG signal figure:
        exSigs=-4:.1:4;
        conds=[0, 1]-.5;  % -.5 = oddball, +.5 =CP
        
        intCol=3;                 % PE term from model.
        condCol=4;
        totLearnCol=5;
        intLearnCol=6;
        
        
        allRelCoeffs=[relROI.ROI_params(:,selROIs,:)]
        allRelCoeffs=allRelCoeffs(goodSub,:,:);
        
        
        % loop through brain signals:
        clear predLR
        for z=1:size(allRelCoeffs, 2)
            relCoeffs=squeeze(allRelCoeffs(:,z,:));
            % loop through subjects:
            for i = 1:length(allRelCoeffs)
                
                subCoefs=relCoeffs(i,:);
                % loop through oddball/CP conditions:
                for j = 1:length(conds)
                    
                    % predLR  =
                    predLR(z,i,j,:)=subCoefs(intCol) +  ...  intercept
                        conds(j)*subCoefs(condCol)   +  ...  condition
                        exSigs.*subCoefs(totLearnCol)+  ...  phys signal
                        exSigs.*subCoefs(intLearnCol)*conds(j); % phys signal* cond
                end
            end
        end

        
        % 1) direct effect on learning, 2) interaction learning*condition

        allRelCoeffs=[relROI.ROI_params(:,selROIs,5), relROI.ROI_params(:,selROIs,6)]
        allRelCoeffs=allRelCoeffs(goodSub,:);
          
        [h, p, b, stats]=ttest(allRelCoeffs)
        meanBeta=nanmean(allRelCoeffs)
        semBeta =nanstd(allRelCoeffs)./sqrt(length(allRelCoeffs));

        
        for i = 1:length(p)
            disp(sprintf('Coefficient stats for %g: p=%g, t=%g, mean=%g, SEM=%s', xLabs(i),p(i), stats.tstat(i), meanBeta(i), semBeta(i)))
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
            if i <=nROI
                plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors(i+2,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
            else
                plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', ms, 'markerFaceColor', cbColors((i+2)-nROI,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
            end
        end
        ylim([-Scale, Scale]);
        xlim([.5, size(allRelCoeffs, 2)+.5]);
        
        
        hold on
        plot([1:4; 1:4], b, 'k', 'lineWidth', 2, 'color', [.7, .7, .7])
        plot([(1:4) - .4;(1:4) + .4 ], [meanBeta; meanBeta], 'k', 'color', [.7 .7 .7], 'lineWidth', 7)

        
        ylabel('Coefficient')
        set(gca, 'xtick', 1:length(xLabs)*2, 'xticklabel', xLabs, 'box', 'off')
        set(gca, 'box','off')
        

    elseif xx==5
        
          
        
        sigToPlot=1;
        relROI_peaks=relROI.peakTime(selROIs);
        cpCol=1;
        odCol=2;
        predToPlotCP=squeeze(predLR(sigToPlot, :,2,:));
        predToPlotOdd=squeeze(predLR(sigToPlot, :,1,:));
        SEM_cp=nanstd(predToPlotCP)./sqrt(length(predToPlotCP));
        SEM_odd=nanstd(predToPlotOdd)./sqrt(length(predToPlotOdd));
        hold on
        H1=shadedErrorBar(exSigs,nanmean(predToPlotCP), SEM_cp, {'-', 'color', cbColors(cpCol,:)}, true)
        H2=shadedErrorBar(exSigs,nanmean(predToPlotOdd), SEM_odd, {'-', 'color', cbColors(odCol,:)}, true)
        ylabel('Learning rate')
        xlabel('Signal strength')
        ff=legend([H1.mainLine, H2.mainLine], 'Changepoint', 'Oddball')
        set(ff, 'location', 'east', 'box', 'off')
        title(num2str(relROI_peaks(sigToPlot)))
        set(gca, 'box','off')
        

    elseif xx==6
        sigToPlot=2;
        relROI_peaks=relROI.peakTime(selROIs);
        cpCol=1;
        odCol=2; 
        predToPlotCP=squeeze(predLR(sigToPlot, :,2,:));
        predToPlotOdd=squeeze(predLR(sigToPlot, :,1,:));
        SEM_cp=nanstd(predToPlotCP)./sqrt(length(predToPlotCP));
        SEM_odd=nanstd(predToPlotOdd)./sqrt(length(predToPlotOdd));
        hold on
        H1=shadedErrorBar(exSigs,nanmean(predToPlotCP), SEM_cp, {'-', 'color', cbColors(cpCol,:)}, true)
        H2=shadedErrorBar(exSigs,nanmean(predToPlotOdd), SEM_odd, {'-', 'color', cbColors(odCol,:)}, true)
        ylabel('Learning rate')
        xlabel('Signal strength')
        ff=legend([H1.mainLine, H2.mainLine], 'Changepoint', 'Oddball')
        set(ff, 'location', 'east', 'box', 'off')
        title(num2str(relROI_peaks(sigToPlot)))
        set(gca, 'box','off')
        
        
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

    end
    
    
    %setPLOT_panelLabel(gca, xx);
end

% kk=annotation('textbox')
% set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')


saveas(gcf,  'p300_learning_suppFig.fig', 'fig')
saveas(gcf,  'p300_learning_suppFig.eps', 'epsc2')
close(gcf)