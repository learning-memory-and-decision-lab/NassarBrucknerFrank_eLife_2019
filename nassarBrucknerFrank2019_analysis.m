%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for Nassar, Bruckner and Frank (2019) eLife
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% In order to run this code, you will first need download:

% 1) the associated EEG data on Dryad: HERE 
% 2) the associated behavioral data here:  https://sites.brown.edu/mattlab/resources/

% Then you will need to change the paths on lines 44-48 to match locations on your machine. 

% Then you should be good to go. 
% To reproduce the figures in the paper you should:

% 1) set forceRunRegs and doSigTesting to TRUE. (lines 58-59)
% 2) run entire script to generate spatiotemporal clusters (referred to as ROIs)       
% 3) set forceRunRegs to FALSE and runRoiAnalysis to TRUE. (line 60)
% 4) rerun script to examine the trial-to-trial activations on learning. 






% EEG TRIGGERS:
% Trigger 1: Trial Onset (subject presses buttons to indicate prediction)
% Trigger 2: Prediction/Fixation Cross 1: When subject presses Space until 500 mSec
% Trigger 3: Outcome: 500 - 1000 mSec
% Trigger 4: Fixation Cross 2: 1000 - 2000 mSec
% Trigger 5: Shield: 2000 - 2500 mSec
% Trigger 6: Fixation Cross 3: 2500 - 3500 mSec
% Trigger 7: Trial Summary


% Trial centered times:
% Trigger 1: Trial Onset (subject presses buttons to indicate prediction)
% Trigger 2: fixCross:        -500-->0
% Trigger 3: Outcome:          0   - 500 mSec
% Trigger 4: Fixation Cross 2: 500 - 1500 mSec
% Trigger 5: Shield:          1500 - 2000 mSec
% Trigger 6: Fixation Cross 3:2000 - 3000 mSec

%% SETUP ANALYSIS

clear classes
behavDataDir='~/Dropbox/helieeg/cannonBehavData'; % should point to the behavioral data
scriptDir='~/Documents/GitHub/NassarBrucknerFrank_eLife_2019'; % should point to the folder with this script
figDir=fullfile(scriptDir, 'figures');
dataDir='~/Dropbox/cannonEegData/alt_preprocess';  % should point to the EEG data


%I don't think eeglab is required... 
%eegLabDir='~/matt_work_stuff/Matt/matlabPackages/EEG/eeglab/';

% Flags:

forceRunRegs             =true; % do you want to run regression? -- if not it will load the most recently saved reg data
doSigTesting             =false;
runRoiAnalysis           =false;
showVM_approximation     =false;
byHand                   =false; % scroll through movies by hand?
baselineCategory         = 2 ;  % 0 =  no baseline term included in model,
% 1 =  mean voltage included in model,
% 2 =  -50-->0 msec baseline included in model
useWeightedReg           =false;
makeSFN_figs             =false;

makefigs                    = false; % make diagnostic figures
storeERPs                   = true;
trialEEGMeasure             = 3;     % 1 =weighted average across cluster voxels, 2 = dot product w cluster t values, 3 = dot product of residual w cluster t values (just regresses out baseline)
tSeriesRes                  = true;  % compute trial-by-trial analyses in sliding windows
recreateConnectionMat       = false;


% This stuff was all to satisfy reviewer requests:
doReviewerRequest           = true;  % in ROI analysis so this needs to be on too...
zROIdata                    = true;  % z-score ROI data before using in prediction
bootstrapSEM                = false;
useAllTrialsForAddingEEG    = false;  % Now we do stuff over time -- so this doesn't matter
unifMix                     = 0;     % uniform mixture component of circular regression fits. 
makeEEGHist                 = false;
doBinnedReg                 = true;


% Analysis Parameters
inc=4;  % downsampling rate (sample every inc trial).
nTrials=480;
goodTrialThresh=0   % exclude subjects missing more than this fraction of data to bad epochs
timeLimit=[-300 2300];
connectThresh=40; % threshold distance for "connected" electrodes
thresh=.0005; % cluster forming threshold for EEG regression
% NOTE: the threshold is set for one-tailed tests in each
% direction, so the "true" threshold is actually thresh*2;

erpChans={'FCz', 'Pz'};  % for main figure: 'FCz', 'Pz'  % At timepoint 158 the peak positive effect is on channel AFz and peak negative on channel P5
numBootstraps = 200;

sortChannel='Fz';
numPerm=1000;  % number of s for permutation testing
regMatsForROI={'allModBasedCoeffs', 'basicCoeffs', 'allLRCoeffs'};
%ROIs_forSingleTrialAnalysis={'allModBasedCoeffs_surprise', 'allModBasedCoeffs_uncertainty', 'allModBasedCoeffs_hit', 'allModBasedCoeffs_CP_cond', 'allLRCoeffs_LR_hat', 'allModBasedCoeffs_surprise_context', 'basicCoeffs_oddball', 'basicCoeffs_changePoint'};
ROIs_forSingleTrialAnalysis={'basicCoeffs_cpPlusOdd', 'basicCoeffs_cpMinusOdd'};


%% RUN ANALYSIS:



addpath(dataDir);
addpath(behavDataDir);
addpath(genpath(scriptDir));
%addpath(genpath(eegLabDir));

% datasets available:
a=dir(fullfile(dataDir, '*_Cannon*'));
fns={a.name};
cd(scriptDir);

for i = 1:length(fns)
   disp(fns{i})  
end

% set plotting defaults:
defaultPlotParameters; close all;
cbColors=[cbColors; cbColors; cbColors];
% make a figure to show the von-mises normal approximation
if showVM_approximation
    alpha=-pi : .01: pi
    p =circ_vmpdf(alpha, 0, 10)
    vmVar=1./10
    gaussStd=vmVar.^.5
    gaussP =normpdf(alpha, 0, gaussStd)
    
    hold on
    plot(alpha, p, 'b', 'linewidth', 10);
    plot(alpha, gaussP, 'r');
    
    fn=fullfile(figDir, 'vm_approx.eps');
    saveas(gcf, 'vm_approx.eps', 'epsc2');
    close all
end


% OK, this is bit of a time sink, runs through all subjects,
% Loops through all subjects and gets behavior and EEG data.
% selects behavioral data for good trials and creates 3 regression matrices
% (one for model-free and one for model-based analysis & pure LR).

% types of things we'll want to do:
% 1) behavioral analysis.
% 2) EEG regressions (requires behavioral data and EEG data)
% 3) beahvioral analyses using EEG trial data from ROIS (requires
% behavioral data, EEG data, and MUST ALREADY HAVE ROIS.

% it seems that we will always want to run the behavioral analysis.
% we will generally want to load the EEG data, but we may not want to run
% the regressions (which take quite some time).

% we will sometimes want to run additional behavioral analyses using ROI
% data.

% Check subject number 33 -- seems to be an issue. 


clear similarityCoeffs similarityCoeffs_CP similarityCoeffs_Odd allOddCoeffs allCPCoeffs allModBasedCoeffs subBehavCoeffs anyUpBehavCoeffs meanAnyUp subBehavCoeffs_anyUp
for j =1:length(fns)
    tic
    fn=fns{j};
    subID=fn(1:3);
    disp(sprintf('starting subject: %s', num2str(subID)));
    
    % load behavioral data
    fn=sprintf('ADL_*_%s.mat', subID);
    a=dir(fullfile(behavDataDir, fn));
    
    % If a is empty, keep adding zeros till you get something
    loopCount=0;
    while isempty(a)
        subID=['0', subID];
        fn=sprintf('ADL_*_%s.mat', subID);
        a=dir(fullfile(behavDataDir, fn));
        
        loopCount=loopCount+1;
        if loopCount>10
            keyboard
        end
    end
    
    
    behavFN=fullfile(behavDataDir, a.name);
    behavData=load(behavFN);
    
    %% Now we need the behavioral data...
    
    % if cBal == 1, load oddball data first (as it came first in the task)
    % if cBal == 2, load cpt data first (as it came first in the task)
    
    %  G in the name indicates that green was the high value color
    %  B in the name indicates that blue was the high value color
    
    clear allBehavData
    %if iscell( unique(eval(sprintf('behavData.DataCP_%s.cBal', subID))));
    %         Order=((unique(eval(sprintf('behavData.DataCP_%s.cBal', subID)))));
    %         Order=str2num(Order{1});
    %else
    Order=unique(eval(sprintf('behavData.DataCP_%s.cBal', subID)));
    %end
    
    subBlockOrder(j)=Order; % 1= oddball first, 2= CP first;
    
    
    if Order==1
        eval(sprintf('allBehavData=behavData.DataOddball_%s;', subID));
        allBehavData=straightStruct(allBehavData);
        eval(sprintf('allBehavData=catBehav(behavData.DataCP_%s, allBehavData);', subID));
    elseif Order==2
        eval(sprintf('allBehavData=behavData.DataCP_%s;', subID));
        allBehavData=straightStruct(allBehavData);
        eval(sprintf('allBehavData=catBehav(behavData.DataOddball_%s, allBehavData);', subID));
    else
        'wtf?';
        keyboard
    end
    allBehavData.totTrial=(1:length(allBehavData.outcome))';
    
    
    %% GET SOME other useful variables from optimal model:
    % Identify blocks
    newBlock=[true; diff(allBehavData.block)~=0];
    begBlock=find(newBlock);
    endBlock=[begBlock(2:end)-1; length(newBlock)];
    
    % create some other variables:
    blockCond=strcmp('main', allBehavData.cond);        % block condition
    trialVal=mod(allBehavData.actRew,2);                % outcome value
    allBehavData.outcome(allBehavData.outcome==360)=0;
    
    % Set parameters for generative model:
    % get average number of change-points and oddballs
    simHaz=mean([nanmean(allBehavData.cp), nanmean(allBehavData.oddBall)]);
    simDrift=unique(allBehavData.driftConc); % get drift
    simNoise=unique(allBehavData.sigma);     % get noise
    
    % preallocate space for model variable estimates
    modPred=nan(size(allBehavData.outcome)); %model prediction (regarding cannon aim)
    modSurp=nan(size(allBehavData.outcome)); %model surprise estimate
    modUnc =nan(size(allBehavData.outcome)); %model uncertainty estimate
    modRU=nan(size(allBehavData.outcome)); %model uncertainty estimate
    
    
    % compute circular PE:
    outcomeRad=deg2rad(allBehavData.outcome);
    predRad=deg2rad(allBehavData.pred);
    subsPredRad=nan(size(predRad));
    subsPredRad(2:end)=predRad(1:end-1);
    subsPredRad(newBlock)=nan;
    
    
    PE=circ_dist(outcomeRad, predRad);
    UP=[circ_dist(predRad(2:end), predRad(1:end-1)); nan];
    newB=find(newBlock(1:end));
    UP(newB(2:end)-1)=nan;
    zeroUp=abs(UP) < abs(deg2rad(.5));
    allBehavData.zeroUp=zeroUp;
    
    
    % OK, what is the best way to compute LR:
    vmVar=1./allBehavData.sigma;
    gaussStd=vmVar.^.5;
    driftVar=1./allBehavData.driftConc;
    gaussDriftStd=driftVar.^.5;
    gaussDriftStd(blockCond==1)=0;
    
    % BEHAVIORAL ANALYSIS:
    % questions...
    % 1) what factors affect whether a subject updates at all?
    % 2) what factors affect how much they update (when they do)?
    % 3) can we get a predicted LR for each trial?
    
    
    [modSurp, modRU, errBased_LR, errBased_UP]=getTrialVarsFromPEs_cannon...
        (gaussStd, PE, simHaz,  newBlock', false(size(PE)),...
        .99, 0, 1, 1, gaussDriftStd, ~blockCond, 2*pi);
 
    
    
    
    
    behavLabels={'intercept', 'PE', 'PE*condition', 'PE*RU', 'PE* CPP/ODD', 'PE*CPP/ODD*condition', 'PE*HIT'};
    
    
    rawX=[ones(size(PE)), PE, PE.* meanCentX(blockCond), PE.*meanCentX(modRU), ...
        PE.* meanCentX(modSurp), PE.* meanCentX(modSurp).* meanCentX(blockCond), ...
        PE.* meanCentX(allBehavData.hit)];

    sel=isfinite(UP) & all(isfinite(rawX), 2);
    
    
    binX=[meanCentX(blockCond), meanCentX(modRU), ...
        meanCentX(modSurp), meanCentX(modSurp).* meanCentX(blockCond), ...
        meanCentX(allBehavData.hit)];
    
    
    % Model updates versus non updates (eg. binary go/nogo)
    anyUp=  UP~=0;
    [B,DEV,STATS] = glmfit(binX(sel,:),anyUp(sel),'binomial');
    meanAnyUp(j)=nanmean(anyUp);
    anyUpBehavCoeffs(j,:)=B;

    
    
    % does regression model make big mistakes on post-oddball trials?
    % MRN added this july 1, 2019 -- it seems that model is failing for
    % post-oddball trials (higher error variance) -- and when we remove
    % these trials we get larger fit values to the behavioral regression
    % coefficients of interest. 
    
    postOdd=find(allBehavData.oddBall==1) +1;
    postOdd(postOdd>length(allBehavData.oddBall))=[];
    postOddLogical=false(size(allBehavData.oddBall));
    postOddLogical(postOdd)=true;
    postOddLogical(newBlock)=false;
    sel=~postOddLogical&isfinite(UP) & all(isfinite(rawX), 2);

    
    % Model continuous updates:
    data.Y=UP(sel);
    data.X=(rawX(sel,:));
    data.priorWidth=ones(1, size((rawX(sel,:)), 2)).*5;
    data.priorMean =[0, .5, zeros(1, 5)];
    data.startPoint=[5, 0, 1, 0, 0, 0, 0, 0, unifMix];
    data.includeUniform=true;
    data.whichParams=logical([1, 1, 1, 1, 1, 1, 1, 1, 0]);
   % keyboard
    [params, ~]=fitLinearModWCircErrs(data);
    intCoefTerms=[ones(size(PE)), meanCentX(blockCond), meanCentX(modRU), meanCentX(modSurp), ...
        meanCentX(modSurp).* meanCentX(blockCond), meanCentX(allBehavData.hit)];
     
    % lets compute how much variance our model explains:
    updateVar(j)=circ_var(data.Y);
    residualVar(j)=circ_var(circ_dist(data.Y, data.X*params(2:end)'));
    fracVarExp(j)=1-(residualVar(j)./updateVar(j));

    % 1 = conc, 2 = intercept, 3= raw PE term.
    LR_hat= intCoefTerms*params(3:sum(data.whichParams))';
    subBehavConc(j)=params(1);
    subBehavCoeffs(j,:)=params(2:end);  % MRN changed this 1/21/19 -- i'm not sure why it was end-1. 


    % June 23... Same model-- except code surprise*PE separately for each condition --
    % to address reviewer concerns:
    % Model continuous updates:
    
    rawX_sepTerms=[ones(size(PE)), PE, PE.*meanCentX(blockCond), PE.*meanCentX(modRU), ...
        PE.* meanCentX(modSurp).*(blockCond==0), PE.* meanCentX(modSurp).*(blockCond==1), ...
        PE.* meanCentX(allBehavData.hit)];
    
    data2.Y=UP(sel);
    data2.X=(rawX_sepTerms(sel,:));
    data2.priorWidth=ones(1, size((rawX_sepTerms(sel,:)), 2)).*5;
    data2.priorMean =[0, .5, zeros(1, 5)];
    data2.startPoint=[5, 0, 1, 0, 0, 0, 0, 0, 0];
    data2.includeUniform=true;
    data2.whichParams=logical([1, 1, 1, 1, 1, 1, 1, 1, 0]);

    % keyboard
    [params_sepTerms, ~]=fitLinearModWCircErrs(data2);
    subBehavConc_sepTerms(j)=params_sepTerms(1);
    subBehavCoeffs_sepTerms(j,:)=params_sepTerms(2:end);  

    
    % once we've got everything right, we'll want the Y-Hat from this model:
    %     LR_hat=nan(size(PE));
    %     LR_hat(sel)=rawX(:,3:end)*params(4:end)';
    %
    
    
    
    % OK, now lets look at updating behavior ONLY for trials that included
    % some updating!
    sel=sel& anyUp;
    data.Y=UP(sel);
    data.X=(rawX(sel,:));
    data.includeUniform=false;
    data.startPoint=[5, 0, 1, 0, 0, 0, 0, unifMix];
    data=rmfield(data, 'whichParams');
    [params2, negLogLike]=fitLinearModWCircErrs(data);
    subBehavConc_anyUp(j)=params2(1);
    subBehavCoeffs_anyUp(j,:)=params2(2:(end-1));

    
    % plug new variables back into data structure:
    allBehavData.modPred=predRad+errBased_UP;
    allBehavData.modSurp=modSurp;
    allBehavData.modRU=modRU;
    allBehavData.LR_hat=LR_hat;
    allBehavData.cpCond=strcmp('main', allBehavData.cond);
    allBehavData.subNum=repmat(j, length(allBehavData.actRew), 1);
    allBehavData.UP=UP;
    allBehavData.PE=PE;
    allBehavData.outcomeRad=outcomeRad;
    allBehavData.predRad=predRad;
    allBehavData.subsPredRad=subsPredRad;
    allBehavData.newBlock=newBlock;
    allBehavData.postOdd=postOddLogical;
    allBehavData=straightStruct(allBehavData);
    
    
    if j ==1
        allSubBehavData=allBehavData;
    else
        allSubBehavData=catBehav(allBehavData, allSubBehavData);
    end
    
    % load preprocessed EEG data for a fixed time interval
    % MRN -- need to re-acquire subject id because for some subjects it
    % differs from behavioral to EEG data.
    fn=fns{j};
    subID=fn(1:3);
    
    
        fn=fullfile(dataDir, sprintf('%s_Cannon_FILT_altLow_STIM.mat', subID));
    
    
    
    disp(sprintf('loading eeg data: %s', fn'));
    eegDat=load(fn);
    
    fracGood=length(eegDat.epochNumbers)./nTrials;
    tNums=eegDat.epochNumbers';
    
    if subID==233
        tNums=tNums+62
    end
    
    
    isGood=false(480,1);
    isGood(tNums)=true;
    isBad=find(~isGood);
    
    EEG=eegDat.EEG;
    
    % select only behavior for which we have EEG data:
    gT=false(size(allBehavData.outcome));
    gT(tNums)=true;
    
    relBehav=selBehav(allBehavData, gT);
    trialVal=mod(relBehav.actRew,2);
    
    sub_goodEEGfrac(j)=nanmean(gT);
    
    
    
    
    downSampData=EEG.data(:,inc:inc:end,:);
    downSampTimes=EEG.times(inc:inc:end);
    
    
    %% MUTUAL INFO ANALYsIS
    
    
    
    
    %% STORE ERPs
    
    if storeERPs
        baselineTime=[-200, 0];
        
        for ch = 1:length(erpChans)
        erpChan=erpChans{ch};
        useChannel=strmatch(erpChan, {EEG.chanlocs.labels});
        totBase=nanmean(nanmean(downSampData(useChannel, downSampTimes>baselineTime(1)&downSampTimes<baselineTime(2),:), 3), 2);
        sel=relBehav.oddBall~=1&relBehav.cp~=1&mod(relBehav.trial-1, 60)>1;
        normalERP=nanmean(downSampData(useChannel, :, sel), 3)-totBase;
        sel=relBehav.oddBall~=1&relBehav.cp==1&mod(relBehav.trial-1, 60)>1;
        cpERP=nanmean(downSampData(useChannel, :, sel), 3)-totBase;
        sel=relBehav.oddBall==1&relBehav.cp~=1&mod(relBehav.trial-1, 60)>1;
        oddERP=nanmean(downSampData(useChannel, :, sel), 3)-totBase;

        
        %         hold on
        %         plot(downSampTimes, normalERP)
        %         plot(downSampTimes, cpERP, 'r')
        %         plot(downSampTimes, oddERP, 'g')
        storedERPs(ch).normalERP(j,:)=normalERP;
        storedERPs(ch).cpERP(j,:)=cpERP;
        storedERPs(ch).oddERP(j,:)=oddERP;
        end
        
    end
    
    
    
    clear eegDat
    timesToLookAt=find(downSampTimes>timeLimit(1) &downSampTimes < timeLimit(2));
    subtractor=min(timesToLookAt)-1;
    
    
    
    if forceRunRegs % ONLY RUN THIS REGRESSION IF YOU MEAN TO!!!
        
        % create explanatory matrix of trial by trial variabls:
        % OK, first model: change-point and oddball regressors
        xMat=meanCentX([ones(sum(gT),1), relBehav.oddBall, relBehav.cp, relBehav.cpCond, relBehav.hit, trialVal]);
        % OK, secondly: model based analysis.
        xMat_mod=meanCentX([ones(sum(gT),1), relBehav.modSurp, ...
            relBehav.modRU, relBehav.hit, relBehav.actRew, relBehav.cpCond, ...
            meanCentX(relBehav.modSurp).*meanCentX(relBehav.cpCond)]);
        xMat_LR=meanCentX([ones(sum(gT),1), relBehav.LR_hat]);
        
        basicCoeffs.labels={'int', 'oddball', 'changePoint', 'CP_block', 'hit', 'trialValue'};
        allModBasedCoeffs.labels={'intercept', 'surprise', 'uncertainty', 'hit', 'rew', 'CP_cond', 'surprise_context'};
        allLRCoeffs.labels={'intercept', 'LR_hat'};
        extraParams=0;
        
        
        if baselineCategory==1
            extraParams=extraParams+1;
            basicCoeffs.labels=[basicCoeffs.labels, 'meanVoltage'];
            allModBasedCoeffs.labels= [allModBasedCoeffs.labels, 'meanVoltage'];
            allLRCoeffs.labels=[allLRCoeffs.labels, 'meanVoltage'];
        elseif baselineCategory==2
            extraParams=extraParams+1;
            basicCoeffs.labels=[basicCoeffs.labels, 'baseline'];
            allModBasedCoeffs.labels= [allModBasedCoeffs.labels, 'baseline'];
            allLRCoeffs.labels=[allLRCoeffs.labels, 'baseline'];
        end
        
        % preallocate space for regression coefficients:
        coeffs=nan(64, length(timesToLookAt), size(xMat,2)+extraParams);
        modCoeffs=nan(64, length(timesToLookAt), size(xMat_mod,2)+extraParams);
        LRCoeffs=nan(64, length(timesToLookAt), size(xMat_LR,2)+extraParams);
        
        for channel = 1:64
            disp(sprintf('starting channel: %s', num2str(channel)));
            
            % compute variance over time as a proxy for expected
            % noise in the EEG signal:
            if useWeightedReg
                chanSignal=(squeeze(downSampData(channel, timesToLookAt,:)));
                chanStd=std(chanSignal)';
            end
            
            switch baselineCategory
                case 0
                    base=[];
                case 1
                    base=meanCentX(squeeze(nanmean(downSampData(channel, :,:), 2)));
                case 2
                    baseTimes= downSampTimes(timesToLookAt)> -50 & downSampTimes(timesToLookAt) <=0;
                    base=meanCentX(squeeze(nanmean(downSampData(channel, baseTimes,:), 2)));
            end
            
            for t=timesToLookAt
                Y=(squeeze(downSampData(channel, t,:)));
                if useWeightedReg
                    
                    coeffs(channel,t-subtractor,:)    = regressW_mike(Y,chanStd, [xMat, base]);
                    modCoeffs(channel,t-subtractor,:) = regressW_mike(Y, chanStd, [xMat_mod, base]);
                    LRCoeffs(channel,t-subtractor,:) = regressW_mike(Y, chanStd, [xMat_LR, base]);
                    
                else
                    
                    % MRN added z-scoring to standard regression:
                    coeffs(channel,t-subtractor,:)    = regress(zscore(Y),[xMat, base]);
                    modCoeffs(channel,t-subtractor,:) = regress(zscore(Y), [xMat_mod, base]);
                    LRCoeffs(channel,t-subtractor,:) = regress(zscore(Y), [xMat_LR, base]);
                    
                end 
            end
        end
        
        basicCoeffs.data(:,:,:,j)=coeffs;
        allModBasedCoeffs.data(:,:,:,j)=modCoeffs;
        allLRCoeffs.data(:,:,:,j)=LRCoeffs;
        
    end

    %keyboard
    % load ROI data and see how it relates to updating behavior
    if runRoiAnalysis
        
        % We should z-score the EEG data across trials... 
        
        relData=downSampData(:,timesToLookAt,:);
         if ~exist('allROIs', 'var')
            %fn=fullfile(scriptDir, 'ROI_data.mat')
            fn=fullfile(scriptDir, 'ROI_data_7-24-19.mat') 
            load(fn);
         end
         
         % Huh, so baseline really doesn't pull out that much variance
         % here. A bit surprising -- but i guess this is because we filter
         % out low frequencies...
         
         if trialEEGMeasure==3 % pull out variance that can be explained by baseline before taking dot product: 
             resData=nan(size(relData)); % preallocate space for residual data.
             baseTimes_roi= downSampTimes(timesToLookAt)> -50 & downSampTimes(timesToLookAt) <=0; % figure out what times are in "baseline" period
             base_roi=(squeeze(nanmean(downSampData(:, baseTimes_roi,:), 2)));
             for CH=1:size(base_roi, 1)
                 ch_base=base_roi(CH,:)';
                 ch_dat =squeeze(relData(CH,:,:)); % zscore original...
                 for TT=1:size(ch_dat, 1) % loop through time within trial
                 xMat=[ones(size(ch_base)), ch_base];    
                 [~,~,resData(CH,TT,:)] = regress(zscore(ch_dat(TT,:)'),xMat); 
                 end
             end
         end

         
         if zROIdata % For original submission we did not z-score across trials here... but we should have. 
             relData= zscore(relData, [], 3);
         end

         
        for f = 1:length(ROIs_forSingleTrialAnalysis) 
            eval(['relROI=allROIs.' ROIs_forSingleTrialAnalysis{f} ';']) ;
            if length(relROI.sign)>0;
                clear meanTrialEffect
                 
                for rr=1:size(relROI.sign);
                    
                    
                    for t=1:size(relData,3)
                        tData=relData(:,:,t); % loop through trials. 
                        if trialEEGMeasure==1
                        meanTrialEffect(t,rr)=nanmean(tData(relROI.maps(:,:,rr)));
                        elseif trialEEGMeasure==2   
                        % use Anne Collins dot product method:
                        meanTrialEffect(t,rr)=tData(relROI.maps(:,:,rr))'*relROI.fullTMap(relROI.maps(:,:,rr));
                        elseif trialEEGMeasure==3  
                        % use dot product method on baseline-regressed residuals:
                        tData=resData(:,:,t); % loop through trials. 
                        meanTrialEffect(t,rr)=tData(relROI.maps(:,:,rr))'*relROI.fullTMap(relROI.maps(:,:,rr));
                        end
                    end
                end
                
                clear ROI_params
                for rr=1:size(relROI.sign);
                    physData=zScoreX(meanTrialEffect(:,rr));
                    
                    rawX=[ones(size(relBehav.PE)), relBehav.PE, relBehav.PE.*meanCentX(relBehav.cpCond),...
                        relBehav.PE.* meanCentX(physData), relBehav.PE.* meanCentX(physData).*meanCentX(relBehav.cpCond)];
                    
                    sel=isfinite(relBehav.UP) & all(isfinite(rawX),2)&~relBehav.postOdd;
                    data.Y=relBehav.UP(sel);
                    data.X=(rawX(sel,:));
                    
                    data.includeUniform=true;   % Include uniform, but don't fit it -- just assume that we occasionally get a uniform pick
                    data.priorWidth=[2, 2, 1, .1, .1];  % strong regularization on phys terms, weak regularization on everything else.
                    data.priorMean= [0, .5, 0, 0, 0];
                    data.startPoint=[5, 0, 1, 0, 0, 0, unifMix];
                    data.whichParams=logical([1, 1, 1, 1,1, 1, 0]);
                    
                    [ROI_params(rr,:), negLogLike]=fitLinearModWCircErrs(data);
                    
                    if bootstrapSEM
                        permData=data;
                        
                        % to estimate variance, lets make priors on terms
                        % of interest as weak as possible:
                        permData.priorWidth=[2, 2, 1, 100, 100];
                        
                        for bb=1:numBootstraps
                            toKeep=randi(length(data.Y), length(data.Y), 1);
                            permData.Y=data.Y(toKeep);
                            permData.X=data.X(toKeep,:);
                            [permBoot(bb,:)]=fitLinearModWCircErrs(permData);
                        end
                        ROI_params_bootStd(rr,:)=  nanstd(permBoot);
                    end
                    
                    
                    if doBinnedReg
                        
                        % SORT data on physiological signal and condition
                        % and compute actual learning slope for each bin.
                        % Reviewers asked for this...
                        
                        nPtiles=5;
                        pTiles=100.*(1:(nPtiles-1))./nPtiles;
                        Bins=[-inf, prctile(physData,pTiles), inf];
                        data.priorWidth=[2, 2];  % strong regularization on phys terms, weak regularization on everything else.
                        data.priorMean= [0, .5];
                        data.startPoint=[5, 0, .5, unifMix];
                        data.whichParams=logical([1, 1, 1, 0]);

                        for Bin=1:(length(Bins)-1)
                            
                            inBin=physData>Bins(Bin)&physData<Bins(Bin+1);
                            
                            
                            medValOfBin(j,Bin)=nanmedian(physData(inBin));
                            
                            % in bin, oddball:
                            rawX=[ones(size(relBehav.PE)), relBehav.PE];
                            sel=relBehav.cpCond==0 & inBin & isfinite(relBehav.UP) & all(isfinite(rawX),2)&~relBehav.postOdd;
                            data.Y=relBehav.UP(sel);
                            data.X=(rawX(sel,:));
                            coeffs=fitLinearModWCircErrs(data);
                            ROI_binnedLR(j,rr,Bin,1)=coeffs(3); % just store PE slope as learning rate
                            
                            % in bin, CP:
                            rawX=[ones(size(relBehav.PE)), relBehav.PE];
                            sel=relBehav.cpCond==1 & inBin & isfinite(relBehav.UP) & all(isfinite(rawX),2)&~relBehav.postOdd;
                            data.Y=relBehav.UP(sel);
                            data.X=(rawX(sel,:));
                            coeffs=fitLinearModWCircErrs(data);
                            ROI_binnedLR(j,rr,Bin,2)=coeffs(3); % just store PE slope as learning rate
                            
                            
                        end
                    end
                    
                end
 
                clear ROI_params_sepConditions
                for rr=1:size(relROI.sign);
                    physData=zScoreX(meanTrialEffect(:,rr));
                    
                    % This model separately models the learning in the
                    % changepoint and oddball blocks: (
                    rawX=[ones(size(relBehav.PE)), relBehav.PE, relBehav.PE.*meanCentX(relBehav.cpCond),...
                        relBehav.PE.* meanCentX(physData).*(relBehav.cpCond==1), relBehav.PE.* meanCentX(physData).*(relBehav.cpCond==0)];
                    
                    sel=isfinite(relBehav.UP) & all(isfinite(rawX),2) &~relBehav.postOdd;
                    data.Y=relBehav.UP(sel);
                    data.X=(rawX(sel,:));
                    
                    data.includeUniform=true;   % Include uniform, but don't fit it -- just assume that we occasionally get a uniform pick... as in how HDDM fits RTs
                    data.priorWidth=[2, 2, 1, .1, .1];  % strong regularization on phys terms, weak regularization on everything else.
                    data.priorMean= [0, .5, 0, 0, 0];
                    data.startPoint=[5, 0, 1, 0, 0, 0, unifMix];
                    data.whichParams=logical([1, 1, 1, 1,1, 1, 0]);
                    
                    [ROI_params_sepConditions(rr,:), negLogLike]=fitLinearModWCircErrs(data);

                    if bootstrapSEM
                        permData=data;
                        permData.priorWidth=[2, 2, 1, 100, 100];
                       
                        for bb=1:numBootstraps
                            toKeep=randi(length(data.Y), length(data.Y), 1);
                            permData.Y=data.Y(toKeep);
                            permData.X=data.X(toKeep,:);
                            [permBoot(bb,:)]=fitLinearModWCircErrs(permData);
                        end
                        ROI_params_sepConditions_bootStd(rr,:)=  nanstd(permBoot);
                    end 
                end

                
                % Goal -- lets do the ROI-type analysis on learning -- but
                % do it in sliding windows across time using all channels. 
                % we can use this to get at whether the differences in ROI
                % effects are due to duration (eg. number of ch/timepoints
                % avged) or qualitative difference in what the signal is
                % representing. 
                if f ==1
                    
                    %keyboard
                     
                    % Step 1) loop through time
                    % Step 2) get thresholded t mask for relevant time
                    % Step 3) get "meanTrialSignal" from mask*trialData
                    % Step 4) fit update data with model that includes meanTrialSignal
                    % step 5) store betas
                    % step 6) finish time loop. 
                    actTimes=downSampTimes(timesToLookAt);
                    windWidth=5;
                    coeffByTime=nan(length(actTimes), 6);                 
                    inBaseWind=actTimes<0&actTimes>-50;
                    resData=nan(size(relData));
                    if trialEEGMeasure ==3
                        for baseElec=1:size(relData,1)
                            for baseTime=1:size(relData,2)
                                xmat=[ones(size(relData, 3), 1), zscore(squeeze(nanmean(relData(baseElec, inBaseWind,:),2)))];
                                [B,BINT,R]=regress(squeeze(relData(baseElec,baseTime, :)), xmat);
                                %coeffMap(baseElec, baseTime)=B(2);                  
                                resData(baseElec, baseTime,:)=R;
                            end
                        end
                    end
%                     imagesc(actTimes, 1:64,  coeffMap)
%                     colorbar
    
                        
                    % Loop through time:
                    for ttt = (windWidth+1):(length(actTimes)-windWidth)
                      
                        minTime(ttt)=actTimes(ttt-windWidth);
                        maxTime(ttt)=actTimes(ttt+windWidth);
                        mask=relROI.fullTMap;
                        % set stuff before relevant time to zero:
                        mask (:, 1:ttt-windWidth-1)=0;
                        % set stuff after relevant time to zero:
                        mask (:, ttt+windWidth+1:end)=0;
                        
                        
                        % Loop through trials and get masked value for each
                        % trial:
                        clear meanTrialTimeEffect
                        for t=1:size(relData,3)
                            if trialEEGMeasure ==3
                            tData=resData(:,:,t);
                            else    
                            tData=relData(:,:,t);
                            end
                            % use Anne Collins dot product method:
                            meanTrialTimeEffect(t)=tData(:)'*mask(:);
                        end
                         
                        % Fit model:
                        physData=zscore(meanTrialTimeEffect)'; 
                        if tSeriesRes==false
                            rawX=[ones(size(relBehav.PE)), relBehav.PE, relBehav.PE.*meanCentX(relBehav.cpCond),...
                                relBehav.PE.* meanCentX(physData), relBehav.PE.* meanCentX(physData).*meanCentX(relBehav.cpCond)];
                        else
                            rawX=[ones(size(relBehav.PE)), relBehav.PE, relBehav.PE.*zscore(relBehav.LR_hat),...
                                relBehav.PE.* meanCentX(physData), relBehav.PE.* meanCentX(physData).*meanCentX(relBehav.cpCond)]; 
                        end
                                
                        sel=isfinite(relBehav.UP) & all(isfinite(rawX),2);
                        data.Y=relBehav.UP(sel);
                        data.X=(rawX(sel,:));
                        
                        data.includeUniform=true;
                        data.priorWidth=[2, 2, 1, .1, .1];  % strong regularization on phys terms, weak regularization on everything else.
                        data.priorMean= [0, .5, 0, 0, 0];
                        data.startPoint=[5, 0, 1, 0, 0, 0, unifMix];
                        data.whichParams=logical([1, 1, 1, 1,1, 1, 0]);

                        [coeffByTime(ttt,:), negLogLike]=fitLinearModWCircErrs(data);
                    end
                    
                    allPredLR.minTime=minTime;
                    allPredLR.maxTime=maxTime;
                    allPredLR.tSeries_mainEffect(j,:)=coeffByTime(:,5);
                    allPredLR.tSeries_intEffect(j,:)=coeffByTime(:,6);
                end
                
                clear ROI_params_wModTerms 
                for rr=1:size(relROI.sign);

                    if useAllTrialsForAddingEEG

                       % New strategy for preventing bad data from
                       % impacting model-based analysis --
                       % we will run circular regression over all trials
                       % but set physiological signal equal to zero on trials
                       % where we don't have measurements of it. 

                       % regress model predicted variance out of the
                       % physiological data:
                       
                       physData=zScoreX(meanTrialEffect(:,rr)); % normalize variance
                       physDataFull=zeros(size(gT));            % insert into vector with length equal to total trials
                       physDataFull(gT)=physData;               % where "bad trials" have physiological signal equal to zero. 

                       % pull predicted variance out of direct learning
                       % effect term:
                       pMat=[ones(length(allBehavData.LR_hat),1), allBehavData.LR_hat];
                       [B,BINT,directLearningResidual] = regress(physDataFull,pMat);

                       
                       % pull predicted variance out of conditional learning
                       % effect term:
                       pMat=[ones(length(allBehavData.LR_hat),1), allBehavData.LR_hat];
                       [B,BINT,conditionalLearningResidual] = regress(physDataFull.*meanCentX(allBehavData.cpCond),pMat);
                       
                       
                       
                       rawX=[ones(size(allBehavData.PE)), allBehavData.PE, allBehavData.PE.*zscore(allBehavData.LR_hat),...
                           allBehavData.PE.*zscore(directLearningResidual), allBehavData.PE.*zscore(conditionalLearningResidual) ...
                           ]; % LR_hat is predicted LR from pure behavioral model.
                       
                       sel=isfinite(allBehavData.UP) & all(isfinite(rawX), 2) &~allBehavData.postOdd;
                     
                       data.Y=allBehavData.UP(sel);
                       data.X=(rawX(sel,:));

                       
                       data.includeUniform=true;
                       data.priorWidth=[2, 2, 1, .1, .1];  % strong regularization on phys terms, weak regularization on everything else. 
                       data.priorMean= [0, .5, 0, 0, 0];
                       data.startPoint=[5, 0, 1, 0, 0, 0, unifMix];
                       data.whichParams=logical([1, 1, 1, 1,1, 1, 0]);
                    
                       [ROI_params_wModTerms(rr,:), negLogLike]=fitLinearModWCircErrs(data);
                       
                    else
                        
                        % This was old way -- only uses trials for which we
                        % have EEG data... We should be able to get better
                        % estimates of all terms in model other than the
                        % EEG terms by adding more trials. This should
                        % reduce the uncertainty on the EEG term through
                        % parameter tradeoff. New way also regresses out
                        % covariance between the predicted learning rate
                        % and the EEG signal. 
                        
                        physData=zScoreX(meanTrialEffect(:,rr));
                        
                         if makeEEGHist
                            hist(physData)
                            hist_fn=sprintf('hist_of_eegEffect_%g_%g.eps', j, rr)
                            saveas(gcf, hist_fn)
                            close all
                        end
                        
                        
                        
                        rawX=[ones(size(relBehav.PE)), relBehav.PE, relBehav.PE.*zscore(relBehav.LR_hat),...
                            relBehav.PE.* meanCentX(physData), relBehav.PE.* meanCentX(physData).*meanCentX(relBehav.cpCond) ...
                            ]; % LR_hat is predicted LR from pure behavioral model.
 
                      
                        sel=isfinite(relBehav.UP) & all(isfinite(rawX),2) &~relBehav.postOdd;
                        data.Y=relBehav.UP(sel);
                        data.X=(rawX(sel,:));
                        
                        yHatCPCorr(j,rr)=corr(relBehav.LR_hat(sel&relBehav.cpCond==1), physData(sel&relBehav.cpCond==1));
                        yHatOddCorr (j,rr)=corr(relBehav.LR_hat(sel&relBehav.cpCond==0), physData(sel&relBehav.cpCond==0));
                        
                        
                        data.includeUniform=true;
                        data.priorWidth=[2, 2, 1, .1, .1];  % strong regularization on phys terms, weak regularization on everything else.
                        data.priorMean= [0, .5, 0, 0, 0];
                        data.startPoint=[5, 0, 1, 0, 0, 0, unifMix];
                        data.whichParams=logical([1, 1, 1, 1,1, 1, 0]);
                              
                        [ROI_params_wModTerms(rr,:), negLogLike]=fitLinearModWCircErrs(data);

                    end
                    
 
                    if bootstrapSEM
                        
                        permData=data;
                         permData.priorWidth=[2, 2, 1, 100, 100];
                       
                        for bb=1:numBootstraps
                        toKeep=randi(length(data.Y), length(data.Y), 1);
                        permData.Y=data.Y(toKeep);
                        permData.X=data.X(toKeep,:);
                        [permBoot(bb,:)]=fitLinearModWCircErrs(permData);
                        end
                        ROI_params_wModTerms_bootStd(rr,:)=  nanstd(permBoot); 
                    end
                    
                end


                % here j is indexing the subject, f is indexing the ROI
                % analysis.
                
                if j ==1
                relROI.ROI_params=[];
                relROI.ROI_params_wModTerms=[];    
                relROI.ROI_params_sepConditions=[];  
                end    

                relROI.ROI_params(j,:,:)=ROI_params;
                relROI.ROI_params_wModTerms(j,:,:)=ROI_params_wModTerms;
                relROI.ROI_params_sepConditions(j,:,:)=ROI_params_sepConditions;
                
                if bootstrapSEM
                relROI.ROI_params_wModTerms_bootStd(j,:,:)=ROI_params_wModTerms_bootStd;
                relROI.ROI_params_bootStd(j,:,:)=ROI_params_bootStd;
                relROI.ROI_params_sepConditions_bootStd(j,:,:)=ROI_params_sepConditions_bootStd;
                end
                
                eval(['allROIs.' ROIs_forSingleTrialAnalysis{f} '=relROI;']) ;
                disp('got here')
                
                
            end
        end
        clear data
    end
    
    toc
end

%% And now to see what we've got:

input('Done with analysis, hit enter to continue:')

% select participants that meet criterion to be included in EEG analysis:
goodSub=sub_goodEEGfrac>goodTrialThresh;

% Ditch black in the plotting colors. 
if sum(cbColors(1,:))==0
    cbColors(1,:)=[]; % ditch black...
end

  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Make ERP figures:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('storedERPs')

    makeFig_ERP_figs;
    
    for ch = 1:length(storedERPs)
        erpChan=erpChans{ch};
        
       
        figure
        hold on
        
        plot([-100, 600], [0,0], '--k', 'linewidth', 2)
        
        meanERP=nanmean(storedERPs(ch).cpERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).cpERP(goodSub,:))./sqrt(sum(goodSub))';
        H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(2,:)},false);
        meanERP=nanmean(storedERPs(ch).oddERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).oddERP(goodSub,:))./sqrt(sum(goodSub))';
        H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(3,:)},false);
        meanERP=nanmean(storedERPs(ch).normalERP(goodSub,:));
        semERP=nanstd(storedERPs(ch).normalERP(goodSub,:))./sqrt(sum(goodSub))';
        H1=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(1,:)},false);
        xlim([-100, 600])
        
        
        ylabel(sprintf('Event related potential (%s)', erpChan))
        xlabel('time relative to outcome')
        ff=legend([H1.mainLine,H2.mainLine, H3.mainLine], {'expected', 'change point', 'oddball'})
        set(ff, 'location', 'northwest', 'box', 'off')
        
        fig_fn=fullfile(figDir, ['ERP_fig_', erpChan, '_' date, '.eps']);
        set(gca, 'box','off')
        ylim([-5, 5])
        saveas(gcf, fig_fn,'epsc2')
        close all
        
        % Now, lets look at CP and oddball Difference curves (subtract out
        % response to expected outcomes)
        
        cpDiff=storedERPs(ch).cpERP-storedERPs(ch).normalERP;
        oddDiff=storedERPs(ch).oddERP-storedERPs(ch).normalERP;
        
        figure
        
        hold on
        plot([-100, 600], [0, 0], '--k')
        meanERP=nanmean(cpDiff(goodSub,:));
        semERP=nanstd(cpDiff)./sqrt(sum(goodSub))';
        H2=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(2,:)},false);
        meanERP=nanmean(oddDiff(goodSub,:));
        semERP=nanstd(oddDiff(goodSub,:))./sqrt(sum(goodSub))';
        H3=shadedErrorBar(downSampTimes, meanERP, semERP ,{'color', cbColors(3,:)},false);
        xlim([-100, 600])
        
        
        ylabel(sprintf('Event related difference (%s)', erpChan))
        xlabel('time relative to outcome')
        ff=legend([H2.mainLine, H3.mainLine], { 'change point', 'oddball'})
        set(ff, 'location', 'northwest', 'box', 'off')
        set(gca, 'box','off')
        fig_fn=fullfile(figDir, ['ERD_fig_', erpChan, '_', date, '.eps']);
        saveas(gcf, fig_fn,'epsc2')
        close all
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save regression data if you ran it, otherwise, load it! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if forceRunRegs
    fn=fullfile(scriptDir, 'just_regs_latestRegressionData_altHP.mat');
    save(fn, 'basicCoeffs', 'allModBasedCoeffs',  'allLRCoeffs')
else
    fn=fullfile(scriptDir, 'just_regs_latestRegressionData_altHP.mat');
    a=load(fn);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Do "region" of interest analyses on spatiotemporal clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get labels for times we actually used:
timeLabels=downSampTimes(timesToLookAt);

if runRoiAnalysis
    
    
    if doReviewerRequest
        % Do behavioral tests requested by reviewers (6/25/19):
        [h, p, bint, stats]=ttest(subBehavCoeffs_sepTerms)

        xLabs_revFig={'PE*surprise (changepoint)', 'PE*surprise (oddball)'}
        
        selCoeffs=5:6;
        sepTermCoeffs_toPlot=subBehavCoeffs_sepTerms(:, selCoeffs)
        Scale=max(abs(sepTermCoeffs_toPlot(:)));
        xJit = smartJitter(sepTermCoeffs_toPlot,.03,.1);
        % is there any covariance among parameter estimates?
        [r, p]=corr(sepTermCoeffs_toPlot)
        ll=size(sepTermCoeffs_toPlot, 1);
        plot([0, size(sepTermCoeffs_toPlot, 2)], [0 0], '--k')
        colorInds=[2, 1];
        hold on
        for i = 1:size(sepTermCoeffs_toPlot, 2)
            plot(ones(ll, 1).*colorInds(i)+xJit(:,i) , sepTermCoeffs_toPlot(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(colorInds(i),:), 'markerEdgeColor', 'k', 'lineWidth', 1);
        end
        ylim([-Scale, Scale]);
        xlim([.5, size(sepTermCoeffs_toPlot, 2)+.5]);
        ylabel('Coefficient')
        set(gca, 'xtick', 1:length(xLabs_revFig)*2, 'xticklabel', xLabs_revFig, 'box', 'off')
        
        fn=fullfile(figDir, ['behavFigForReviewrs',  '.eps']);
        saveas(gcf, fn, 'epsc2')
        
        
        mean(sepTermCoeffs_toPlot)
        std(sepTermCoeffs_toPlot)./sqrt(length((sepTermCoeffs_toPlot)))
        
    end
    
    
    %load sepCoeffsWorkspace.mat
    
    
    % Create an image of the template:
    exTempBin=(allROIs.basicCoeffs_cpPlusOdd.maps(:,:,allROIs.basicCoeffs_cpPlusOdd.peakTime==390))
    
    
    %     MAKE AN EXAMPLE TEMPLATE:
    %     threshT=allROIs.basicCoeffs_cpPlusOdd.fullTMap.*exTempBin
    %     imagesc(threshT, [0, 8])
    %     saveas(gcf, 'exTemplate.eps', 'epsc2')
    
    %  keyboard
    for f = 1   % f = 2 gets CP-oddball figures...
        
        disp(['getting ' ROIs_forSingleTrialAnalysis{f} ])
        
        eval(['relROI=allROIs.' ROIs_forSingleTrialAnalysis{f} ';']) ;
        
        if isfield(relROI, 'ROI_params')
            
            % Hand select ROIs. Ugh. 
            selROIs=[1]; % if f = 2
             
             
            % 1) direct effect on learning, 2) interaction learning*condition
            
            
            % ROI number 6 seems to be specific to change-points!!!
            %selROIs=[1];
            allRelCoeffs=[relROI.ROI_params(:,selROIs,5), relROI.ROI_params(:,selROIs,6)]
            allRelCoeffs=allRelCoeffs(goodSub,:);
            
            xLabs=(repmat(relROI.peakTime(selROIs), 1, 2))
            
            
            relROI.peakTime(selROIs)
            Scale=max(abs(allRelCoeffs(:)));
            xJit = smartJitter(allRelCoeffs,.03,.1);
            % is there any covariance among parameter estimates?
            [r, p]=corr(allRelCoeffs)
            ll=size(allRelCoeffs, 1);
            plot([0, size(allRelCoeffs, 2)], [0 0], '--k')
            nROI=length(selROIs);
            hold on
            for i = 1:size(allRelCoeffs, 2)
                if i <=nROI
                    plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(i,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                else
                    plot(ones(ll, 1).*i+xJit(:,i) , allRelCoeffs(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(i-nROI,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                end
            end
            ylim([-Scale, Scale]);
            xlim([.5, size(allRelCoeffs, 2)+.5]);
            ylabel('Coefficient')
            set(gca, 'xtick', 1:length(xLabs)*2, 'xticklabel', xLabs, 'box', 'off')
            
            fn=fullfile(figDir, ['ROI_learning_figure_', ROIs_forSingleTrialAnalysis{f}, '.eps']);
            saveas(gcf, fn, 'epsc2')
            close all
         
            % REMOVED Y-HAT stuff from here.

            % GET BEHAVIORAL REGRESSION COEFFICIENTS OF INTEREST: 
            surpriseEffectOnLR=subBehavCoeffs(:,5);
            contextualSurpriseEffectOnLR=subBehavCoeffs(:,6);
            
               
            make_p300_learningFig_r1;
            %make_p300_resLearningFig;
            %makeIndDiffFig;
            
            
            makeSepCondFig
            
            
            % resave ROI data that now includes single trial coefficients.
            %fn=fullfile(scriptDir, 'ROI_data.mat')
            
            fn=fullfile(scriptDir, 'ROI_data_7-24-19.mat')
            save(fn, 'allROIs');
        end

    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test whether ODD+CPP signal predicts learning beyond predictions made by behavioral model  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TEST physiological LR predictions over time:
if exist('allPredLR')
    
    relRange=[300, 700];
    numPerms=50000;

    meanTime=(allPredLR.maxTime+allPredLR.minTime)./2;
    isDat=isfinite(sum(allPredLR.tSeries_intEffect));

    % PLOT DIRECT EFFECT ON LEARNING:
    % TEST physiological LR predictions over time:
    % Convolve with gaussian (width = 2*interval = 16 ms)
    B=normpdf([-6:6], 0, 2 )
    B=B./sum(B);
    for sub=1:size(allPredLR.tSeries_intEffect, 1)
        smMainEffect(sub,:)= conv(allPredLR.tSeries_mainEffect(sub,:), B, 'same');
    end
    meanBeta1=nanmean(smMainEffect(goodSub,:));
    semBeta1=nanstd(smMainEffect(goodSub,:))./sqrt(sum(goodSub));

    hold on
    plot([meanTime(1), meanTime(end)], [0, 0], '--k')
    H1=shadedErrorBar(meanTime(isDat),meanBeta1(isDat),semBeta1(isDat),{'color', cbColors(3,:)},false);
    ylim([-.06, .06]);
    set(gca, 'box', 'off');
    xlim([0, 1000]);
  

    [h, p, int, stats]=ttest(smMainEffect(goodSub,:));
    meanTime=(allPredLR.maxTime+allPredLR.minTime)./2;
    isDat=isfinite(stats.tstat);
    plot([min(meanTime), max(meanTime)], [0, 0], '--k');
    plot([relRange' relRange']', [-4, 4; -4, 4]', '--k');
    xlim([0, 1000]);
 
    % Do significance testing:
    maskedDat=smMainEffect(goodSub,meanTime>=relRange(1)& meanTime<=relRange(2));
    [h, p, int, stats]=ttest(maskedDat);
    totCount=0;
    clustMass=zeros(size(maskedDat,2),1);
    for i = 1:size(maskedDat,2)
        if h(i)==1
            totCount=totCount+1;
        elseif i>1&& h(i-1)==1
            clustMass((i-totCount):(i-1))=abs(sum(stats.tstat((i-totCount):(i-1)))) ;
            totCount=0;
        end
    end
    realMass=clustMass;
    
    maxClustMass=nan(numPerms,1);
    for j = 1:numPerms
        % flip random subject coefficients:
        permMaskedDat= maskedDat;
        doFlip=logical(binornd(1,.5,size(maskedDat,1),1));
        
        permMaskedDat(doFlip,:)=  permMaskedDat(doFlip,:).*-1;
        [h, p, int, stats]=ttest(permMaskedDat);
        totCount=0;
        clustMass=zeros(size(permMaskedDat,2),1);
        for i = 1:size(permMaskedDat,2)
            if h(i)==1
                totCount=totCount+1;
            elseif i>1&& h(i-1)==1
                clustMass((i-totCount):(i-1))=abs(sum(stats.tstat((i-totCount):(i-1)))) ;
                totCount=0;
            end
        end
        maxClustMass(j)=max(clustMass);
    end
    nanmean(max(realMass)>maxClustMass) 
    
    maskedTime=meanTime(meanTime>=relRange(1)& meanTime<=relRange(2));
    maskedBeta1=meanBeta1(meanTime>=relRange(1)& meanTime<=relRange(2));
    clustCorrSig1=realMass> prctile(maxClustMass, 95);
    plot(maskedTime(clustCorrSig1), maskedBeta1(clustCorrSig1), '.r');
    disp('Main effect stats:')
    disp(sprintf('Cluster Mass = %g \n Cluster corrected p = %g', max(realMass), 1-nanmean(max(realMass)>maxClustMass)));
    
    
    % PLOT CONDITIONAL LEARNING EFFECT:
    

    % Convolve with gaussian (width = 2*interval = 16 ms)
    B=normpdf([-6:6], 0, 2 );
    B=B./sum(B);
    for sub=1:size(allPredLR.tSeries_intEffect, 1)
      smIntEffect(sub,:)= conv(allPredLR.tSeries_intEffect(sub,:), B, 'same');
    end

        
    
    % plot interaction:

    [h, p, int, stats]=ttest(smIntEffect(goodSub,:));
    meanBeta2=nanmean(smIntEffect(goodSub,:));
    semBeta2=nanstd(smIntEffect(goodSub,:))./sqrt(sum(goodSub));
    hold on
    plot([meanTime(1), meanTime(end)], [0, 0], '--k');
    H2=shadedErrorBar(meanTime(isDat),meanBeta2(isDat),semBeta2(isDat),{'color', cbColors(4,:)},false);
    ylim([-.06, .06]);
    set(gca, 'box', 'off');
    xlim([0, 1000]);

%     hold on
%     plot(meanTime(isDat), stats.tstat(isDat))
%     ylim([-4, 4])
%     plot([min(meanTime), max(meanTime)], [0, 0], '--k')
    ylabel('Residual learning coefficient');
    xlabel('Time');
    %title('Interaction')
    plot([relRange' relRange']', [-4, 4; -4, 4]', '--k') ;
    xlim([0, 1000]);
    ff=legend([H1.mainLine, H2.mainLine],  'Direct', 'Conditional');
    set(ff, 'box', 'off');

    
    
    % Do significance testing:
    maskedDat=smIntEffect(goodSub,meanTime>=relRange(1)& meanTime<=relRange(2));
    [h, p, int, stats]=ttest(maskedDat);
    totCount=0;
    clustMass=zeros(size(maskedDat,2),1);
    for i = 1:size(maskedDat,2)
        if h(i)==1
            totCount=totCount+1;
        elseif i>1&& h(i-1)==1
            clustMass((i-totCount):(i-1))=abs(sum(stats.tstat((i-totCount):(i-1)))) ;
            totCount=0;
        end
    end
    realMass=clustMass;
    
    maxClustMass=nan(numPerms,1);
    for j = 1:numPerms
        % flip random subject coefficients:
        permMaskedDat= maskedDat;
        doFlip=logical(binornd(1,.5,size(maskedDat,1),1));
        
        permMaskedDat(doFlip,:)=  permMaskedDat(doFlip,:).*-1;
        [h, p, int, stats]=ttest(permMaskedDat);
        totCount=0;
        clustMass=zeros(size(permMaskedDat,2),1);
        for i = 1:size(permMaskedDat,2)
            if h(i)==1
                totCount=totCount+1;
            elseif i>1&& h(i-1)==1
                clustMass((i-totCount):(i-1))=abs(sum(stats.tstat((i-totCount):(i-1)))) ;
                totCount=0;
            end
        end
        maxClustMass(j)=max(clustMass);
    end
    nanmean(max(realMass)>maxClustMass);
       
    maskedTime=meanTime(meanTime>=relRange(1)& meanTime<=relRange(2));
    maskedBeta2=meanBeta2(meanTime>=relRange(1)& meanTime<=relRange(2));
    clustCorrSig2=realMass> prctile(maxClustMass, 95);
    
    disp('Interaction effect stats:')
    disp(sprintf('Cluster Mass = %g \n Cluster corrected p = %g', max(realMass), 1-nanmean(max(realMass)>maxClustMass)));

    

    if tSeriesRes==true
        fn='trialByTrialEffectsAcrossTimeWindows_res.eps';
    else
        fn='trialByTrialEffectsAcrossTimeWindows.eps';
    end
    
    saveas(gcf, fn, 'epsc2')
    close all


    % get a t-map of the CP+Oddball effect for schematic:
    make_p300_resLearningFig_d2;
    

    
end


%% Cluster-based multiple comparisons corrected hypothesis testing
% The idea here is to identify "clusters" of activation in the space of electrodes
% and time.


% GET CONNECTION MATRIX specifying which electrodes are connected to which:
if exist('connectionMat.mat') & ~recreateConnectionMat
    load connectionMat.mat % if we've already got one made, load it.
else
    % otherwise, create one from scratch:
    
    % Then get the XYZ coordinates for each channel:
    allLocas=[[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]]' ;
    % Loop through the channels and get the distance between that channel and
    % all other channels.
    for i = 1:length(allLocas)
        relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
    end
    
    % Set a threshold on distance... and mark channels that fall within that
    % threshold:
    connectionMat=relDist<connectThresh;
    connectionMat=connectionMat-eye(length(connectionMat));
    
    % CHECK OUT CONNECTIONS:
    %imagesc(connectionMat);
    %close all
    
    % check average number of connections:
    %numCons=nanmean(sum(connectionMat));
    
    % Look at connectivity weights across scalp, alla AC:
    % figure;
    % for e = 1:64
    % subplot(8,8,e)
    % topoplot(connectionMat(e,:),EEG.chanlocs)
    % end
    % saveas(gcf, 'connectionMatFig.eps', 'epsc2')
    % close all
    save connectionMat.mat connectionMat
end



% Create a sorting variable for displaying regression coefficients over
% channels and time:
allLocas=[[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]]' ;
for i = 1:length(allLocas)
    relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
end
sortChannelIndex=strmatch( sortChannel, {EEG.chanlocs.labels});
[dist, displayIndex]=sort(relDist(sortChannelIndex,:), 'descend');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Do significance testing?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRN changed this march 22 so it only uses "good subjects"

if ~exist('allLRCoeffs')
    allLRCoeffs=a.allLRCoeffs;
    allModBasedCoeffs=a.allModBasedCoeffs;
    basicCoeffs=a.basicCoeffs;
end


    
%fn=fullfile(scriptDir, 'ROI_data.mat')
fn=fullfile(scriptDir, 'ROI_data_7-24-19.mat')
if doSigTesting
    % GET ALL ROIs
    clear allROIs
    for i = 1: length(regMatsForROI)
        eval(['regDat=' regMatsForROI{i} '.data;']);
        eval(['regLabels=' regMatsForROI{i} '.labels;']);
        
        % if we are looking at basic regression coefficients, create
        % contrasts for oddball + CP and oddball - CP:
        if strcmp('basicCoeffs', regMatsForROI{i})
 
            cpRegNum=strmatch('changePoint' , regLabels) 
            oddballRegNum=strmatch('oddball' , regLabels) 
            
            % Create CP + Oddball contrast:
            regDat(:,:,end+1,:)= regDat(:,:,cpRegNum,:)+regDat(:,:,oddballRegNum,:); 
            regLabels{end+1}='cpPlusOdd';
            
            
            % Create CP - Oddball contrast:
            regDat(:,:,end+1,:)= regDat(:,:,cpRegNum,:)-regDat(:,:,oddballRegNum,:); 
            regLabels{end+1}='cpMinusOdd';
            
          end
        
        % Get 8th ROI
        
        for j = 1:size(regDat, 3)
            
            EEG_dat=permute(regDat(:,:,j,goodSub), [4, 1, 2, 3]);
            % get clusters, cluster sizes, cluster masses for positive clusters (ie
            % p<.05 in a one tailed positive test):
            posClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'right');

            
            aaa=figure
            hold on
            imagesc(timeLabels, 1:64, posClusterInfo.tMap(displayIndex,:))
            %imagesc(timeLabels, 1:64, posClusterInfo.tMap(:,:))
            plot([0, 0], [1, 64], '-k')
            plot([500, 500], [1, 64], '-r')
            set(gca, 'box', 'off')
            
            % get clusters, cluster sizes, cluster masses for negative clusters (ie
            % p<.05 in a one tailed negative test):
            negClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'left');
            
            if i ==2&j==2;
                
                % MRN 3-6-19 work from here!!!!
                
                % CP - oddball:
                cp_odd_contrast=permute(regDat(:,:,3,goodSub), [4, 1, 2, 3])-permute(regDat(:,:,2,goodSub), [4, 1, 2, 3]);
            
                %[H,P,CI,STATS]=ttest(cp_odd_contrast);
                
                posDifference=getEEG_clusterSize(cp_odd_contrast, connectionMat, thresh, 'right');
                negDifference=getEEG_clusterSize(cp_odd_contrast, connectionMat, thresh, 'left');
 
                %t_map=posDifference.tMap;
                
                figure
                hold on
                imagesc(timeLabels, 1:64, posDifference.tMap)
                imagesc(timeLabels, 1:64, posClusterInfo.tMap(:,:))
                plot([0, 0], [1, 64], '-k')
                plot([500, 500], [1, 64], '-r')
                set(gca, 'box', 'off')
                ylabel('channel')
                xlabel('time relative to outcome')
                colorbar
                close all
                
            end
                
                
            if byHand
                figure
                hold on
                imagesc(timeLabels, 1:64, posClusterInfo.tMap(displayIndex,:))
                %imagesc(timeLabels, 1:64, posClusterInfo.tMap(:,:))
                plot([0, 0], [1, 64], '-k')
                plot([500, 500], [1, 64], '-r')
                set(gca, 'box', 'off')
                ylabel('channel')
                xlabel('time relative to outcome')
                colorbar
                close all
                
                [h]=viewEEG_timeseries(posClusterInfo.tMap(displayIndex,:), timeLabels, EEG.chanlocs, .00000001, 10, [0, 500, 1500, 2000], false, displayIndex)
                close all
                
            end
            
            % OK... now we have our statistics of interest... we just need a null
            % distribution to compare them to. Lets run through a loop that:
            % 1) flips the signs of the data for each subject according to a fair coin
            % toss.
            % 2) find the maximum statistic values that you see in the entire dataset
            % from these "permuted" datasets.
            
            
            clusterInfo=posClusterInfo;
            allTs=0:25:max(timeLabels);
            
            for k=1:numPerm
                % create a subject length array containing randomly assigned -1's and
                % 1s:
                permArray=ones(size(EEG_dat, 1), 1);
                permArray(logical(binornd(1, .5, size(EEG_dat, 1), 1)))=-1;
                
                % multiply each subject timeseries by the -1 or 1 assigned randomly on
                % this trial
                sz=size(EEG_dat);
                permMat=repmat(permArray, [1 sz(2:end)]);
                
                % get cluster statistics for permuted dataset:
                permClusterInfo=getEEG_clusterSize(EEG_dat.*permMat, connectionMat, thresh, 'right');
                
                % store the maximum of the statistics, for use in null distribution:
                maxSize(k)=(max(permClusterInfo.clustSizeMap(:)));
                maxWt(k)=(max(permClusterInfo.clustWtMap(:)));
            end
            
            % For a two tailed test, find the the minimum cluster statistics necessary
            % to beat 97.5% of the null distribution
            % Based on "mass"
            % I like this statistic better!!!
            massTh = prctile(maxWt,  97.5 );
            gPos=unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap	>massTh));
            gNeg=unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap	>massTh));
            threshMap=zeros(size(negClusterInfo.tMap));
            threshMap(posClusterInfo.clustWtMap	>massTh)=posClusterInfo.tMap(posClusterInfo.clustWtMap	>massTh);
            threshMap(negClusterInfo.clustWtMap	>massTh)=negClusterInfo.tMap(negClusterInfo.clustWtMap	>massTh);
            nn=[regMatsForROI{i}, '_', regLabels{j}];
            
            
            
            if byHand
                figure
                hold on
                lims=max(abs(threshMap(:)));
                imagesc(timeLabels, 1:64, threshMap(displayIndex,:), [-lims, lims])
                colorbar
                set(gca, 'box', 'off')
                title(nn);
                xlabel('time (ms)')
                ylabel('channel (Fz top)')
                saveas(gcf, fullfile(figDir, [nn 'Fig.eps']), 'epsc2')
                close all
                
                [h]=viewEEG_timeseries(threshMap, timeLabels, EEG.chanlocs, .00000001, 10, [0, 500, 1500, 2000], true, displayIndex);
                close all
                
            end
            
            
            % STORE ROIs
            % we'll want:
            % 1) map --
            % sign(direction) of effect
            % t map
            
            clear ROIs
            ROIs.sign=[ones(size(gPos)); -ones(size(gNeg))];
            k =1
            while k <=length(gPos)
                ROIs.maps(:,:,k)=posClusterInfo.ID_map==gPos(k);
                k=k+1
            end
            
            while k <= length(gPos) +length(gNeg)
                ROIs.maps(:,:,k)=negClusterInfo.ID_map==gNeg(k-length(gPos));
                k=k+1
            end
            
            for k = 1:length(ROIs.sign)
                clustTs=ROIs.maps(:,:,k).*abs(posClusterInfo.tMap);
                [ROIs.peakChannel(k),J] = find(clustTs==max(clustTs(:)));
                ROIs.peakTime(k)=timeLabels(J);
                
                
                clustIsSig(:,k)=sum(ROIs.maps(:,:,k))>1;
                
            end
            
            ROIs.fullTMap=posClusterInfo.tMap;
            eval(['allROIs.', nn, '=ROIs']);
            

            % MAKE basic plots:           
            subplot(2, 1,1)
            hold on
            for k = 1:length(ROIs.sign)
                plot(timeLabels', clustIsSig(:,k)*100+k, 'o', 'markerFaceColor', cbColors(k+1,:), 'markerEdgeColor', 'none', 'markerSize', 8)
            end
            ylim([100 100+length(ROIs.sign)+1])
            xlim([min(timeLabels), max(timeLabels)])
            set(gca, 'box', 'off')
            colorbar

            subplot(2, 1,2)
            cLim=[8]
            hold on
            imagesc(timeLabels, 1:length(EEG.chanlocs), posClusterInfo.tMap)
            set(gca, 'clim', [-cLim, cLim])
            colorbar
            xlim([min(timeLabels), max(timeLabels)])
            ylim([0, 65])
            set(gca, 'box', 'off', 'ytickLabels', '')
            
            
            
            
            
            if isfield(ROIs, 'peakTime')
                plot(ROIs.peakTime,  ROIs.peakChannel, 'ok', 'markerSize', 6, 'markerFaceColor', 'none', 'markerEdgeColor', 'w', 'lineWidth', 1)
            end
            
            eventTimes=[0, 500, 1500, 2000]
            for event=1:length(eventTimes)
                plot([eventTimes(event), eventTimes(event)], [1, 64], '-k')
            end
            
            
            fn=[nn, 'just_time.eps'];
            saveas(gcf, fullfile(figDir,fn), 'epsc2')
            close all
            

            if isfield(ROIs, 'peakTime')
                
                % Lets make plots to show NO difference effect at key
                % timepoints:
                if i ==2 & j ==9
                    ROIs.peakTime=[ROIs.peakTime, 494, 670];
                end
                
               
                if i ==2 & j ==8
                     ROIs.peakTime=[ROIs.peakTime];
                end
                
                
                for k = 1:length(ROIs.peakTime)
                    
                    % which channel is max:
                    
%                     peakPos=find(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k))== max(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k))))                   
%                     peakNeg=find(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k))== min(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k))))                   
% 
%                     EEG.chanlocs(peakPos)
%                    disp(sprintf('At timepoint %s the peak positive effect is on channel %s and peak negative on channel %s', num2str(ROIs.peakTime(k)), EEG.chanlocs(peakPos).labels, EEG.chanlocs(peakNeg).labels));
          
                    figure
                    topoplot(posClusterInfo.tMap(:,timeLabels==ROIs.peakTime(k)),EEG.chanlocs );
                    title([nn, ' ' num2str(ROIs.peakTime(k))]);
                    fn=[nn, ' ' num2str(ROIs.peakTime(k)), '.eps']
                    set(gca, 'clim', [-5, 5])
                    colorbar    
                    saveas(gcf, fullfile(figDir, fn), 'epsc2')     
                    close all
                end
                
            end
            
            
           % do same loop at 150 msec, 
           rawList=[150, 350, 502 654]
           
           
                for k = 1:length(rawList)
                    figure
                    topoplot(posClusterInfo.tMap(:,timeLabels==rawList(k)),EEG.chanlocs );
                    title([nn, ' ' num2str(rawList(k))]);
                    fn=[nn, '_' num2str(rawList(k)), '_', date, '.eps']
                    set(gca, 'clim', [-5, 5])
                    colorbar  
                    saveas(gcf, fullfile(figDir, fn), 'epsc2')     
                    close all
                end

            
            % MAKE plot of individual subject ROI values:
            
            
            % get subject ROI coefficients:
            
            if isfield(ROIs, 'maps')
 
                subCoefMat=nan(size(EEG_dat, 1), size(ROIs.maps,3));
                for k = 1:size(EEG_dat, 1);
                    subDat=squeeze(EEG_dat(k,:,:));
                    for z = 1:size(ROIs.maps,3)
                        subCoefMat(k,z)=nanmean(subDat(ROIs.maps(:,:,z)));
                    end
                end
                
                ROIs.subCoeffs=subCoefMat;
                
                Scale=max(abs(subCoefMat(:)));
                xJit = smartJitter(subCoefMat,.1,.1);
                % is there any covariance among parameter estimates?
                [r, p]=corr(subBehavCoeffs)
                
                ll=size(subCoefMat, 1);
                plot([0, size(subCoefMat, 2)], [0 0], '--k')
                
                
                nPlot=min([size(subCoefMat, 2), 8]);
                
                
                hold on
                plot([0, nPlot+.5], [0, 0], '--k')
                
                for k = 1:nPlot
                    plot(ones(ll, 1).*k+xJit(:,k) , subCoefMat(:,k), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(9-k,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
                end
                ylim([-Scale, Scale]);
                xlim([.5, size(subCoefMat, 2)+.5]);
                ylabel('Coefficient')
                set(gca, 'xticklabel', {}, 'box', 'off')
                
                fn=fullfile(figDir, [nn, '_indSubCoeffs.eps'])
                saveas(gcf, fn, 'epsc2')
                close all
            end
            
        end
    end
    
    fn=fullfile(scriptDir, 'ROI_data_7-24-19.mat')
    save(fn, 'allROIs');
    
else
    load(fn);
end



% 6-20-19:
% 1: Get separate measures of CP and ODDball coeffs from ROIs
% 2: Make "funnel plot" on number of included trials


% Original submission:
% earlyP300_ind=find(allROIs.basicCoeffs_cpPlusOdd.peakTime==494);
% lateP300_ind=find(allROIs.basicCoeffs_cpPlusOdd.peakTime==670);

% New clusters after including all participants:
earlyP300_ind=find(allROIs.basicCoeffs_cpPlusOdd.peakTime==390);
%lateP300_ind=find(allROIs.basicCoeffs_cpPlusOdd.peakTime==678);


for i = 1:size(a.basicCoeffs.data, 4)
    % Get CP and ODD coefficients for each subject for all time/electrode
    % pairs
    subCpData=a.basicCoeffs.data(:,:,3,i);
    subOddData=a.basicCoeffs.data(:,:,2,i);

    meanCoeff_earlyCp(i)  = nanmean(subCpData(allROIs.basicCoeffs_cpPlusOdd.maps(:,:,earlyP300_ind)));
    meanCoeff_earlyOdd(i) = nanmean(subOddData(allROIs.basicCoeffs_cpPlusOdd.maps(:,:,earlyP300_ind)));
%     meanCoeff_lateCp(i)   = nanmean(subCpData(allROIs.basicCoeffs_cpPlusOdd.maps(:,:,lateP300_ind)));
%     meanCoeff_lateOdd(i)   = nanmean(subOddData(allROIs.basicCoeffs_cpPlusOdd.maps(:,:,lateP300_ind)));

end

% Do stats for each individual contrast/cluster:
[h1, p1, s1, STATS1]=ttest(meanCoeff_earlyCp(goodSub))
[h2, p2, s2, STATS2]=ttest(meanCoeff_earlyOdd(goodSub))
% [h3, p3, s3, STATS3]=ttest(meanCoeff_lateCp(goodSub))
% [h4, p4, s4, STATS4]=ttest(meanCoeff_lateOdd(goodSub))




subDataToPlot=[meanCoeff_earlyCp(goodSub)', meanCoeff_earlyOdd(goodSub)',]

Scale=max(abs(subDataToPlot(:)));
xJit = smartJitter(subDataToPlot,.1,.05);
% is there any covariance among parameter estimates?
[r, p]=corr(subBehavCoeffs)

ll=size(subDataToPlot, 1);
plot([0, size(subDataToPlot, 2)], [0 0], '--k')


nPlot=min([size(subDataToPlot, 2), 8]);



FS=14
hold on
plot([0, nPlot+.5], [0, 0], '--k')

cpOddLabels={'Changepoint', 'Oddball'}
clear H
for k = 1:nPlot
    if mod(k, 2)==0;
        pointColor=cbColors(2,:);
    else
         pointColor=cbColors(1,:);
         
    end
        
    H(k)=plot(ones(ll, 1).*k+xJit(:,k) , subDataToPlot(:,k), 'o', 'markerSize', 10, 'markerFaceColor', pointColor, 'markerEdgeColor', 'k', 'lineWidth', 1);
end
ylim([-Scale, Scale]);
xlim([.5, size(subDataToPlot, 2)+.5]);
ylabel('Coefficient')

set(gca, 'xtick', [1:4], 'xticklabel', cpOddLabels, 'box', 'off', 'fontSize', FS)
%ff=legend(H(1:2), 'Changepoint coefficient', 'Oddball Coefficient')
%set(ff, 'location', 'southwest', 'fontSize', FS, 'box', 'off')
fn=fullfile(figDir, [nn, '_indSub_CP_ODD_P300_COEFFS.eps'])
saveas(gcf, fn, 'epsc2')
close all












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     MAKE FIGURES FOR SFN!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if makeSFN_figs
    
    
    exPE=-pi:.1:pi
    simUncOdd=    nanmean(allBehavData.modRU(~allBehavData.cpCond));
    simUncCP=    nanmean(allBehavData.modRU(allBehavData.cpCond))
    
    for i = 1:length(exPE)
        
        [modSurp, modRU, errBased_LR_oddball(i), errBased_UP]=getTrialVarsFromPEs_cannon...
            (gaussStd(1), exPE(i), simHaz,  true, false,...
            simUncOdd, 0, 1, 1, gaussDriftStd(1), true, 2*pi);
        
        
        [modSurp, modRU, errBased_LR_CP(i), errBased_UP]=getTrialVarsFromPEs_cannon...
            (gaussStd(1), exPE(i), simHaz,  true, false,...
            simUncCP, 0, 1, 1, gaussDriftStd(1), false, 2*pi);
        
        
    end
    
    
    hold on
    plot([-pi, pi], [-pi, pi], '--k')
    plot([-pi, pi], [0, 0], '--k')
    for i = .1:.1:.9
        plot([-pi, pi], i.*[-pi, pi]+ (1-i).* [0, 0], 'color', [.5 .5 .5])   ;
    end
    plot(exPE, exPE.*errBased_LR_oddball, 'color', cbColors(3,:))
    plot(exPE, exPE.*errBased_LR_CP, 'color', cbColors(2,:))
    ylabel('Update')
    xlabel('Prediction error')
    xlim([-pi, pi])
    ylim([-pi, pi])
    set(gca, 'box', 'off')
    fn='LR_schematic_Fig.eps';
    saveas(gcf, fullfile(figDir, fn), 'epsc2')
    close all
    
    
    
    
    % Schematic to explain basic behavioral model:
    
    % load figureWorkspace_10-7-15.mat
    
    for i= 1:max(allSubBehavData.subNum)
        
        %keyboard
        
        subSel=allSubBehavData.subNum==i;
        hold on
        
        sel=subSel&~strcmp(allSubBehavData.cond, 'main')
        plot([-pi, pi], [-pi, pi], '--k')
        plot([-pi, pi], [0, 0], '--k')

        
        % MRN changed to plot points in random order so that you can see
        % both blue and orange ones. 
        
        
        
        
        hold on
        yDat=allSubBehavData.UP(subSel);
        xDat=allSubBehavData.PE(subSel);
        colorDat=strcmp(allSubBehavData.cond(subSel), 'main');
        Order=randperm(length(yDat));
        for j = 1:length(Order)
            if abs(xDat(Order(j)))<.5
                ms=4;
            elseif abs(xDat(Order(j)))<1
                ms=6;
            elseif abs(xDat(Order(j)))<1.5    
                 ms=8;
            else
                ms=10;
            end
            
            if colorDat(Order(j))==1
                
                plot(xDat(Order(j)), yDat(Order(j)), 'o', 'markerFaceColor', cbColors(2,:), 'markerSize', ms, 'markerEdgeColor', 'none', 'lineWidth', 1);
                
            else
            plot(xDat(Order(j)), yDat(Order(j)), 'o', 'markerFaceColor', cbColors(3,:), 'markerSize', ms, 'markerEdgeColor', 'none', 'lineWidth', 1);
            end
        end

        xlabel('Prediction error')
        ylabel('Update')
        
        ylim([-pi, pi])
        xlim([-pi, pi])
        set(gca, 'box', 'off')
        fn=[num2str(i) '_exBehavAnalysisFig.eps'];
        saveas(gcf, fullfile(figDir, fn), 'epsc2')
        
        
        
        
        % Zoom in on [-1, 1]
        ylim([-.3, .3])
        xlim([-.3, .3])
        set(gca, 'xticklabels', '', 'yticklabels', '')
        xlabel('')
        ylabel('')
        fn=[num2str(i) '_exBehavAnalysisFig_inset.eps'];
        saveas(gcf, fullfile(figDir, fn), 'epsc2')
        
        close all
        
        
        
        subplot(2, 1, 1)
        
        sel=subSel&~strcmp(allSubBehavData.cond, 'main')&allSubBehavData.block==4;
        
        hold on
        plot(allSubBehavData.distMean(sel), '--', 'color', cbColors(1,:))
        plot(allSubBehavData.outcome(sel), 'o', 'markerSize', 8, 'markerFaceColor', cbColors(4,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
        plot(allSubBehavData.pred(sel), '-', 'color', cbColors(7,:))
        ylabel('Angle (deg)')
        xlabel('Trial')
        ff=legend('Cannon', 'Cannon ball', 'Shield position')
        set(ff, 'location', 'northeast', 'box', 'off')
        set(gca, 'box', 'off')
        
        ylim([0, 450])
        
        subplot(2, 1, 2)
        
        hold on
        plot(allSubBehavData.modSurp(sel), '-', 'color', cbColors(2,:))
        plot(allSubBehavData.modRU(sel), '-', 'color', cbColors(3,:))
        xlabel('Trial')
        gg=legend('Oddball probability', 'Uncertainty')
        set(gg, 'location', 'northeast', 'box', 'off')
        set(gca, 'box', 'off')
        ylim([0, 1.2])
        
        fn=fullfile(figDir, [num2str(i), '_oddballExRun.eps']);
        saveas(gcf, fn, 'epsc2')
        close all
        
        
        
        subplot(2, 1, 1)
        
        sel=subSel&strcmp(allSubBehavData.cond, 'main')&allSubBehavData.block==4;
        
        hold on
        plot(allSubBehavData.distMean(sel), '--', 'color', cbColors(1,:))
        plot(allSubBehavData.outcome(sel), 'o', 'markerSize', 8, 'markerFaceColor', cbColors(4,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
        plot(allSubBehavData.pred(sel), '-', 'color', cbColors(7,:))
        ylabel('Angle (deg)')
        xlabel('Trial')
        ff=legend('Cannon', 'Cannon ball', 'Shield position')
        set(ff, 'location', 'northeast', 'box', 'off')
        set(gca, 'box', 'off')
        ylim([0, 450])
        
        
        
        
        
        subplot(2, 1, 2)
        
        hold on
        plot(allSubBehavData.modSurp(sel), '-', 'color', cbColors(2,:))
        plot(allSubBehavData.modRU(sel), '-', 'color', cbColors(3,:))
        xlabel('Trial')
        gg=legend('Change-point probability', 'Uncertainty')
        set(gg, 'location', 'northeast', 'box', 'off')
        set(gca, 'box', 'off')
        ylim([0, 1.2])
        
        
        fn=fullfile(figDir, [num2str(i), '_changepointExRun.eps']);
        saveas(gcf, fn, 'epsc2')
        
        close all
        
    end
    
    
    
    behavLabels;
    gs=subBehavCoeffs(:,1)<10;
    goodSubBehavCoeffs=subBehavCoeffs(gs,:);
    goodSubBehavCoeffs_anyUp=subBehavCoeffs_anyUp(gs,:);
    
    % STATS FOR MANUSCRIPT REPORTING!
    [H,P,CI,STATS]= ttest(subBehavCoeffs)
    meanBeta=nanmean(subBehavCoeffs);
    semBeta=nanstd(subBehavCoeffs)./sqrt(sum(gs))
    
    for zz = 1:length(meanBeta)
        disp(sprintf('STATS for %s', behavLabels{zz}))
        disp(sprintf('mean/SEM beta = %g/%g, t = %g, dof= %g, p = %g', meanBeta(zz), semBeta(zz), STATS.tstat(zz), STATS.df(zz), P(zz)))
        
    end
    
    
    % BASIC BEHAVIORAL MODEL FIGURE:
    
    varNormCoeffs=goodSubBehavCoeffs./repmat((std(goodSubBehavCoeffs)), size(goodSubBehavCoeffs, 1), 1);
    Scale=max(abs(goodSubBehavCoeffs(:)));
    
    xJit = smartJitter(varNormCoeffs,.06,.6);
    % is there any covariance among parameter estimates?
    [r, p]=corr(goodSubBehavCoeffs)
    
    
    
    
    % Make behavioral figure for paper (MRN 1-20-19)
    figure(1)
    ll=size(varNormCoeffs, 1);
    plot([0, size(varNormCoeffs, 2)], [0 0], '--k')
    hold on
    for i = 2:size(varNormCoeffs, 2)
        plot(ones(ll, 1).*i+xJit(:,i) , goodSubBehavCoeffs(:,i), 'o', 'markerSize', 6, 'markerFaceColor', cbColors(i,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
        signRankP(i)=signrank(goodSubBehavCoeffs(:,i))
    end
    ylim([-Scale, Scale]);
    %ylim([-2, 4.2]);
    
    xlim([1.5, size(varNormCoeffs, 2)+.5]);
    set(gcf, 'Position', [10,10, 520, 235])
    
    
    ylabel('Coefficient')
    ylim([-2, 2])
    set(gca, 'xtick', 1:size(varNormCoeffs, 2), 'xticklabel', {}, 'box', 'off')
    saveas(gcf, 'basicBehaviorFig_reSize_7-28-19.eps', 'epsc2')
    
    
    close all
    
    behavLabels
    
    cbColors(end+1,:)=[.1, .1, .1]
    % BASIC BEHAVIORAL MODEL FIGURE... but excluding non-update trials:
    
    Scale=max(abs(goodSubBehavCoeffs_anyUp(:)));
    xJit = smartJitter(varNormCoeffs,.1,.8);
    % is there any covariance among parameter estimates?
    [r, p]=corr(goodSubBehavCoeffs_anyUp)
    
    ll=size(varNormCoeffs, 1);
    plot([0, size(varNormCoeffs, 2)], [0 0], '--k')
    
    hold on    
    behavLabels
    
    cbColors(end+1,:)=[.1, .1, .1]
    % BASIC BEHAVIORAL MODEL FIGURE... but excluding non-update trials:
    
    Scale=max(abs(goodSubBehavCoeffs_anyUp(:)));
    xJit = smartJitter(varNormCoeffs,.1,.8);
    % is there any covariance among parameter estimates?
    [r, p]=corr(goodSubBehavCoeffs_anyUp)
    
    ll=size(varNormCoeffs, 1);
    plot([0, size(varNormCoeffs, 2)], [0 0], '--k')
    
    hold on
    for i = 1:size(varNormCoeffs, 2)
        plot(ones(ll, 1).*i+xJit(:,i) , goodSubBehavCoeffs_anyUp(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
    end
    ylim([-Scale, Scale]);
    xlim([1.5, size(varNormCoeffs, 2)+.5]);
    ylabel('Coefficient')
    set(gca, 'xtick', 1:size(varNormCoeffs, 2), 'xticklabel', behavLabels, 'box', 'off')
    saveas(gcf, 'basicBehaviorFig_anyUp.eps', 'epsc2')

    for i = 1:size(varNormCoeffs, 2)
        plot(ones(ll, 1).*i+xJit(:,i) , goodSubBehavCoeffs_anyUp(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
    end
    ylim([-Scale, Scale]);
    xlim([1.5, size(varNormCoeffs, 2)+.5]);
    ylabel('Coefficient')
    set(gca, 'xtick', 1:size(varNormCoeffs, 2), 'xticklabel', behavLabels, 'box', 'off')
    saveas(gcf, 'basicBehaviorFig_anyUp.eps', 'epsc2')
    close all
    
    
    %     binX=[meanCentX(blockCond), meanCentX(modRU), ...
    %       meanCentX(modSurp), meanCentX(modSurp).* meanCentX(blockCond), ...
    %         meanCentX(allBehavData.hit)]
    %
    
    
    %         binX=[meanCentX(blockCond), meanCentX(modRU), ...
    %         meanCentX(modSurp), meanCentX(modSurp).* meanCentX(blockCond), ...
    %         meanCentX(allBehavData.hit)]
    %
    
    anyUpLabels={'int', 'context', 'RU', 'surprise', 'surprise*context', 'hit'}
    
    Scale=max(abs(anyUpBehavCoeffs(:)));
    xJit = smartJitter(anyUpBehavCoeffs,.12,.6);
    % is there any covariance among parameter estimates?
    [r, p]=corr(goodSubBehavCoeffs_anyUp)
    
    ll=size(anyUpBehavCoeffs, 1);
    plot([0, size(anyUpBehavCoeffs, 2)], [0 0], '--k')
    hold on
    for i = 1:size(anyUpBehavCoeffs, 2)
        plot(ones(ll, 1).*i+xJit(:,i) , anyUpBehavCoeffs(:,i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(i,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
    end
    ylim([-10, 10]);
    xlim([.5, size(anyUpBehavCoeffs, 2)+.5]);
    ylabel('Coefficient')
    set(gca, 'xticklabel', anyUpLabels, 'box', 'off')
    saveas(gcf, 'basicBehaviorFig_upOrNot.eps', 'epsc2')
    close all
    
    %% WHOLE BRAIN CORRECTED HYPOTHESIS TESTING:
    
    % right now this is setup to run analyses using a single coefficient map
    % (EEG_dat):
    
end





