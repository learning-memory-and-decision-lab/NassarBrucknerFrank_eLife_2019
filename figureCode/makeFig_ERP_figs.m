% makeFig1


num        = 1;
wid        = 17.6; % total width
hts        = [3, 3];
cols       = {2, 2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, [], [], [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time


xLims=[-100, 700]
erpYLims=[-3, 5];
erdYLims=[-2.7, 2.7];


cpColorInd=3;
oddColorInd=2;
%neutralColorInd=7; Just Black!





set(0, 'defaultlinelinewidth', 1);



for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        ch=1;
        erpChan=erpChans{ch};
        hold on
        plot(xLims, [0,0], '--k', 'linewidth', 1)
        
        meanERP=nanmean(storedERPs(ch).normalERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).normalERP(goodSub,:))./sqrt(sum(goodSub))';
        H1=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', [0,0,0]},false);
        xlim(xLims)

        meanERP=nanmean(storedERPs(ch).cpERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).cpERP(goodSub,:))./sqrt(sum(goodSub))';
        H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(cpColorInd,:)},false);
        meanERP=nanmean(storedERPs(ch).oddERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).oddERP(goodSub,:))./sqrt(sum(goodSub))';
        H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(oddColorInd,:)},false);
        ylabel(sprintf('Event related potential (%s)', erpChan))
        %xlabel('time relative to outcome')
        ff=legend([H1.mainLine,H2.mainLine, H3.mainLine], {'expected', 'change point', 'oddball'})
        set(ff, 'location', 'northwest', 'box', 'off')
        fig_fn=fullfile(figDir, ['ERP_fig_', erpChan, '_' date, '.eps']);
        set(gca, 'box','off')
        ylim(erpYLims)

        
        
    elseif xx==2
         
         ch=1;
         erpChan=erpChans{ch};
         cpDiff=storedERPs(ch).cpERP-storedERPs(ch).normalERP;
         oddDiff=storedERPs(ch).oddERP-storedERPs(ch).normalERP; 
         hold on
         plot(xLims, [0, 0], '--k')
         meanERP=nanmean(cpDiff(goodSub,:));
         semERP=nanstd(cpDiff)./sqrt(sum(goodSub))';
         H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(cpColorInd,:)},false);
         meanERP=nanmean(oddDiff(goodSub,:));
         semERP=nanstd(oddDiff(goodSub,:))./sqrt(sum(goodSub))';
         H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(oddColorInd,:)},false);
         xlim(xLims)
         ylabel(sprintf('Event related difference (%s)', erpChan))
         %xlabel('time relative to outcome')
         ff=legend([H2.mainLine, H3.mainLine], { 'change point', 'oddball'})
         set(ff, 'location', 'northwest', 'box', 'off')
         set(gca, 'box','off')
         fig_fn=fullfile(figDir, ['ERD_fig_', erpChan, '_', date, '.eps']);
         ylim(erdYLims)
         set(gca, 'box', 'off')

    elseif xx==3
        
        ch=2;
         erpChan=erpChans{ch};
        hold on
        plot(xLims, [0,0], '--k', 'linewidth', 1)
        meanERP=nanmean(storedERPs(ch).normalERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).normalERP(goodSub,:))./sqrt(sum(goodSub))';
        H1=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', [0,0,0]},false);
        
        meanERP=nanmean(storedERPs(ch).cpERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).cpERP(goodSub,:))./sqrt(sum(goodSub))';
        
        H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(cpColorInd,:)},false);
        meanERP=nanmean(storedERPs(ch).oddERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).oddERP(goodSub,:))./sqrt(sum(goodSub))';
        H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(oddColorInd,:)},false);
        xlim(xLims)
        ylabel(sprintf('Event related potential (%s)', erpChan))
        xlabel('time relative to outcome')
%       ff=legend([H1.mainLine,H2.mainLine, H3.mainLine], {'expected', 'change point', 'oddball'})
%       set(ff, 'location', 'northwest', 'box', 'off')
        fig_fn=fullfile(figDir, ['ERP_fig_', erpChan, '_' date, '.eps']);
        set(gca, 'box','off')
        ylim(erpYLims)

    elseif xx==4
        
         ch=2;
          erpChan=erpChans{ch};
         cpDiff=storedERPs(ch).cpERP-storedERPs(ch).normalERP;
         oddDiff=storedERPs(ch).oddERP-storedERPs(ch).normalERP; 
         hold on
         plot(xLims, [0, 0], '--k')
         meanERP=nanmean(cpDiff(goodSub,:));
         semERP=nanstd(cpDiff)./sqrt(sum(goodSub))';
         H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(cpColorInd,:)},false);
         meanERP=nanmean(oddDiff(goodSub,:));
         semERP=nanstd(oddDiff(goodSub,:))./sqrt(sum(goodSub))';
         H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(oddColorInd,:)},false);
         xlim(xLims)
         ylabel(sprintf('Event related difference (%s)', erpChan))
         xlabel('time relative to outcome')
%        ff=legend([H2.mainLine, H3.mainLine], { 'change point', 'oddball'})
%        set(ff, 'location', 'northwest', 'box', 'off')
         set(gca, 'box','off')
         fig_fn=fullfile(figDir, ['ERD_fig_', erpChan, '_', date, '.eps']);
         ylim(erdYLims)
         set(gca, 'box', 'off')
    end
    
    
    %setPLOT_panelLabel(gca, xx);
end

% kk=annotation('textbox')
% set(kk, 'string', 'Nassar et al 2018 Figure 1', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')



saveas(gcf,  'erpFig_noExclude_altH.fig', 'fig')
saveas(gcf,  'erpFig_noExclude_altH.eps', 'epsc2')
close(gcf)
defaultPlotParameters
close all
