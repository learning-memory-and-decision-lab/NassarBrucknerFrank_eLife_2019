function timeFrequencyStruct=TF_Script(EEG, params)

% NOTE: implement linear spacing of cycles.

% Written by MRN to do time frequency analysis ala Jim Cavanagh.
% This function should basically run a nice clean form of the script that
% MRN got from SM got from JC.

% Outputs:  timeFrequencyStruct is a data structure containing the
% following fields:

%   power
%   phase angle
%   frequencies
%   times

% params is a structure containing parameters that are used in the
% analysis, EEG is a data structure in the eeglab format.

% params should contain the following fields:

% minTime   -- minimum time of frequency window (careful of edges)
% maxTime   -- maximum time of frequency window (careful of edges)
% sites     -- list of channels that you want in the analysis
% numcycles -- number of cycles (Jim's code used 4)
% num_freqs -- number of frequencies
% baseTime  -- [minBase maxBase] -- time to use as baseline
% doBaseline-- remove baseline power from each trial. Don't do this. 

% output structure will contain the following fields:

% baseSubtractTrialPower -- logarithm of baseline subtracted trial by trial
%                           power spectrum 
% A_TF_dB                -- logarithm of baseline subtracted average power
%                           spectrum
% baseline               -- logarithm of baseline power
% trialPhaseMap          -- phase map for each freq, time, trial, electrode (4d)
% trialPowerMap          -- power map for each freq, time, trial, electrode
% times                  -- times corresponding to columns of output
% freqs                  -- frequencies corresponding to rows of output
% sites                  -- just in case you forgot what sites you wanted.

% Power/phase info shoudl come out as 4d: 
% 1 = freq
% 2 = time within trial
% 3 = trial
% 4 = channel

% Baseline power info (for normalization) should come out 3d:
% 1 = freq
% 2 = trial
% 3 = chanel


                                                                 
% what timepoints are we going to compute power for?
startIndex=find(EEG.times>params.minTime,1);
endIndex=find(EEG.times<params.maxTime, 1, 'last');
storeTimes=EEG.times(startIndex:endIndex);
% baseline timing info:
b1=find(EEG.times>=params.baseTime(1),   1);
b2=find(EEG.times<params.baseTime(2),   1, 'last');


% lets store some useful dimensional info;
inSize=size(EEG.data)
tt=length(storeTimes)

% Jim's setup
frex=logspace(.01,1.7,params.num_freqs); % frequencies to probe
s=params.numcycles./(2*pi.*frex);    % Time required per cycle
t=-2:1/EEG.srate:2; % ok... this is a very large window of times... maybe he is going to 0 pad all the way out to these edges?
time2display=find(EEG.times==min(EEG.times)):find(EEG.times==max(EEG.times)); %Initial Time Frequency plot times to display


%% preallocating memory for the data structure:
% dims: freq, time, trial, electrode.
baseSubtractTrialPower= nan([params.num_freqs, tt, inSize(3), length(params.sites)]);
% A_TF_dB               = nan([params.num_freqs, tt]);
% A_TF                  = nan([params.num_freqs, tt, 2]);
trialPhaseMap         = nan([params.num_freqs, tt, inSize(3), length(params.sites)]);
trialPowerMap         = nan([params.num_freqs, tt, inSize(3), length(params.sites)]);
baselinePower         = nan([params.num_freqs, inSize(3), length(params.sites)]);



disp('space preallocated: proceeding to analysis')


%% Create Morlet wavelet for each frequency

% note, don't need to do this for each channel.. 

% Morlet wave is convolution of complex exponential (eg. e^ix --
% and remember eulers law -- so also cos x + i sin x)
% AND a gaussian envelope.

clear i
w=nan(length(frex), length(t));
for fi=1:length(frex)
    w(fi,:) = exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*s(fi)^2));
    % time res: 2sigma_t *1000 for ms
    timeres(fi)=2*s(fi) * 1000;
    % freq res: 2sigma_f
    frexres(fi)=2* (frex(fi) / params.numcycles);
end

dims=size(EEG.data);

% loop through vector of channels:
for ii = 1:length(params.sites)
    site=params.sites(ii)  %ch 40 is FCz - see chanlocs for other channels
    tic
    % convolve EEG signal with morlet wavelets:
  
    for fi=1:params.num_freqs
        
        % Jim takes an interesting approach to the edge effect issue...
        % instead of computing the convolution of wavelet and signal for
        % each trial separately with zero padding... he puts all of the
        % trials in a row as if they are a single timeseries... then
        % convolves, then removes the windows from the edges.
        
        % it seems like an obvious danger to this strategy is that
        % something that shows up in the end of the trial might also affect
        % the beginning of the trial... but assuming we start from large
        % segments fo data i guess this isn't a problem.
        
        % Anyway, here is the convolution
        stuffnjunk=fconv_JFC(reshape(EEG.data(site,:,:),1,dims(2)*dims(3)),w(fi,:));
        % And now we remove a chunk of convolved data so that our input and
        % output dimensions match.
        stuffnjunk=stuffnjunk((size(w,2)-1)/2:end-1-(size(w,2)-1)/2); %cut off 1/2 the length of the w from beg, and 1/2 from the end
        % and we put it back in the same dimensions as our original data.
        stuffnjunk=reshape(stuffnjunk,dims(2),dims(3));
        
      %  keyboard
        
%         % Power info:
%         A_TF(fi,:,1, ii) = mean(abs(stuffnjunk(startIndex:endIndex,:)).^2,2); % POWER
%         
%         % Phase information:
%         A_TF(fi,:,2, ii) = abs(mean(exp(1i*(angle(stuffnjunk(startIndex:endIndex,:)))),2)); % PHASE consistency?        
       
        baselinePower(fi,:, ii) = (mean(abs(stuffnjunk(b1:b2,:)).^2,1)); % POWER
        
        trialPowerMap(fi,:,:, ii)= abs(stuffnjunk(startIndex:endIndex,:)).^2;  % store power for each trial.
        trialPhaseMap(fi,:,:, ii)= angle(stuffnjunk(startIndex:endIndex,:));  % store phase for each trial.

       
    end % end frequency loop
    
%     if params.doBaseline
%    % keyboard
%     baseline(:,ii)=log10(mean(squeeze(A_TF(:,b1:b2,1)),2)); % OK... this baseline is subtracting the baseline from ALL trials... not baselinefrom a single trial. 
%     basecue=repmat(baseline(:,ii) ,1, size(A_TF,2));
%     % subtract baseline power from average power;
%     A_TF_dB(:,:,ii) = 10*( log10(squeeze(A_TF(:,:,1, ii)))   -  basecue);
%     
%     else
%     A_TF_dB(:,:,ii) = 10*( log10(squeeze(A_TF(:,:,1, ii)))   -  basecue);    
%     end   

    
    % subtract baseline power from trial power estimates;
    % DIM: freq, time, trials, electrode 
%    baseSubtractTrialPower(:,:,:,ii)=log10(trialPowerMap(:,:,:, ii))- repmat(basecue, [1 1 size(trialPowerMap, 3)]);
 toc
end


% create an output data structure

timeFrequencyStruct=struct;
%timeFrequencyStruct.baseSubtractTrialPower=baseSubtractTrialPower;  % baselined trial power estimates... first and last trials are probably problematic
%timeFrequencyStruct.averagePowerAndPhase  =A_TF_dB ;
timeFrequencyStruct.baselinePower         =baselinePower;
timeFrequencyStruct.trialPhaseMap         =trialPhaseMap;
timeFrequencyStruct.trialPowerMap         =trialPowerMap;
timeFrequencyStruct.times                 =storeTimes;
timeFrequencyStruct.freqs                 =frex;
timeFrequencyStruct.sites                 =params.sites;
    