function [h]=viewEEG_timeseries(eegDat, times, chanlocs, delay, timeInc,...
    eventTimes, manualScroll, displayIndex)

% eegDat is the data (num electrodes X timepoints)
% times is the timing of each timepoint
% chanlocs is the channel locations necessary for topoplot
% delay is the amout of pausing on each loop for the movie
% timeInc is the amount of time you want to jump through on each loop
% event times is an array of event times that you want marked in the video
% and manual scroll is a flag that allows you to move forward and backward
% using the keyboard.

h=figure
allTs=min(times):25:max(times);
close all

t=1;
while t<=length(allTs)
    
    ind=find(times>=allTs(t), 1);
    cLim=max(abs(eegDat(:)));
    subplot(2, 1, 1)
    cla
    topoplot(eegDat(:,ind), chanlocs)
    set(gca, 'clim', [-cLim, cLim])
    
    
    
    subplot(2, 1, 2)
    hold on
    imagesc(times, 1:length(chanlocs), eegDat(displayIndex,:))
    set(gca, 'clim', [-cLim, cLim])
    colorbar
    
    
    for event=1:length(eventTimes)
        plot([eventTimes(event), eventTimes(event)], [1, length(chanlocs)], '-k')
    end
        
    plot([allTs(t), allTs(t)], [1, 64], '-r')
    
    set(gca, 'box', 'off')
    ylabel('channel')
    xlabel('time relative to outcome')

    if manualScroll
        
        switch input('h to go forward, g to go back, enter to quit!', 's')
            case 'h'
                t=t+1
            case 'g'
                t=t-1
            case 'save'
                saveas(gcf, 'tmp_figure.eps', 'epsc2')
            case ''
                close all
                break
            otherwise
                disp('unknown method')
        end
    else
        pause(delay);
        t=t+1;
        
    end
end
