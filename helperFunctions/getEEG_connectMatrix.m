% Get channel matrix based on 

% First load some EEG data, ie:
% EEG=pop_loadset(fn)

% Then get the XYZ coordinates for each channel:
allLocas=[[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]]' ;
% Loop through the channels and get the distance between that channel and
% all other channels.
for i = 1:length(allLocas)
    relDist(:,i)=sqrt(sum((allLocas-repmat(allLocas(i,:), length(allLocas), 1)).^2, 2));
end

% Set a threshold on distance... and mark channels that fall within that
% threshold:
thresh=40;
connectionMat=relDist<thresh

connectionMat=connectionMat-eye(length(connectionMat));

% visualize connectivity matrix:
imagesc(connectionMat)
    
% check to make sure that on average, things are connected to ~4 other nodes:    
nanmean(sum(connectionMat))

% Look at connectivity weights across scalp, alla AC:
figure;
for e = 1:64
subplot(8,8,e)
topoplot(connectionMat(e,:),EEG.chanlocs)
end
saveas(gcf, 'connectionMatFig.eps', 'epsc2')
close all
% save connectionMat
save connectionMat.mat connectionMat
close all

%% Cluster-based multiple comparisons corrected hypothesis testing 
% The idea here is to identify "clusters" of activation in the space of electrodes
% and time. 


% load the connection matrix if you didn't just generate it:
% load connectionMat.mat

% load some data... ERP is a dataset that Anne gave me to test this code.
% Anne's data had subject specific t-stats for error versus correct. 
% Standard fMRI-style approach would be to use regression coefficients
% directly, rather than t-stats (coefficient/standard error of
% coefficient). I'm not really going to think too hard about the
% implications of using one or the other right now... 
EEG_dat=ERP(sujets,:,:,3);  % EEG_dat = subjects X channels X timepoints.

EEG_dat=EEG_dat(:,:,751:1500); % chop dataset down so that 
timesCut=times(751:1500)

thresh=.025; % cluster forming threshold. any ttest across subjects getting 
%               a p value lower than this will be included as cluster (or piece of a cluster) 
             % NOTE: the threshold is set for one-tailed tests in each
             % direction, so the "true" threshold is actually thresh*2;
             
             
% get clusters, cluster sizes, cluster masses for positive clusters (ie
% p<.05 in a one tailed positive test):
posClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'right');

posClusterInfo_clean=getEEG_clusterSize_clean(EEG_dat, connectionMat, thresh, 'right');

imagesc(posClusterInfo.ID_map)
figure
imagesc(posClusterInfo_clean.ID_map)

figure
plot(posClusterInfo.ID_map(:), posClusterInfo_clean.ID_map(:), '.')

figure
imagesc(posClusterInfo_clean.ID_map==2)

overlay=double(posClusterInfo_clean.ID_map==6);
overlay(posClusterInfo_clean.ID_map==7)=7;

imagesc(overlay, [0 10])


figure
imagesc(posClusterInfo.ID_map==5)

unique(posClusterInfo_clean.ID_map(posClusterInfo.ID_map==10))


% get clusters, cluster sizes, cluster masses for negative clusters (ie
% p<.05 in a one tailed negative test):
negClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'left');


% OK... now we have our statistics of interest... we just need a null
% distribution to compare them to. Lets run through a loop that:
% 1) flips the signs of the data for each subject according to a fair coin
% toss.
% 2) find the maximum statistic values that you see in the entire dataset
% from these "permuted" datasets.

numPerm=200;
for i=1:numPerm
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
    maxSize(i)=(max(permClusterInfo.clustSizeMap(:)));
    maxWt(i)=(max(permClusterInfo.clustWtMap(:)));
end

% For a two tailed test, find the the minimum cluster statistics necessary
% to beat 97.5% of the null distribution
sizeTh = prctile(maxSize,97.5) 
massTh = prctile(maxWt,97.5) 


% Look for clusters that beat the null:

% Based on "mass"  
% I like this statistic better!!!

gPos=unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap	>massTh));
gNeg=unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap	>massTh));

threshMap=zeros(size(negClusterInfo.tMap));
threshMap(posClusterInfo.clustWtMap	>massTh)=posClusterInfo.tMap(posClusterInfo.clustWtMap	>massTh);
threshMap(negClusterInfo.clustWtMap	>massTh)=negClusterInfo.tMap(negClusterInfo.clustWtMap	>massTh);


% look at thresholded maps across time:
figure
imagesc(threshMap, [-5 5])

figure % look at all positive clusters
imagesc(posClusterInfo.ID_map)

figure % and just suprathreshold clusters
sigClust=posClusterInfo.ID_map;
sigClust(posClusterInfo.clustWtMap	<massTh)=-1;
imagesc(sigClust)
close all



% We can also test significance based on "size". 
% I'm less into this stat, as it doesn't give you any credit for being
% "really signifant"... but we can look at it anyway:

sgPos=unique(posClusterInfo.ID_map(posClusterInfo.clustSizeMap	>sizeTh));
sgNeg=unique(negClusterInfo.ID_map(negClusterInfo.clustSizeMap	>sizeTh));

sthreshMap=zeros(size(negClusterInfo.tMap));
sthreshMap(posClusterInfo.clustSizeMap	>sizeTh)=posClusterInfo.tMap(posClusterInfo.clustSizeMap	>sizeTh);
sthreshMap(negClusterInfo.clustSizeMap	>sizeTh)=negClusterInfo.tMap(negClusterInfo.clustSizeMap	>sizeTh);




% Its sort of nice to look at the clusters in topoplots... which requires
% chanLocs... i was loading a separate EEGlab dataset just to get this, but
% i'm sure you guys have something like this lying around.
% load EEG_demoFile.mat
a=figure
timtopo(threshMap, EEG.chanlocs, 'title', 'mass thresh')
b=figure
timtopo(sthreshMap, EEG.chanlocs, 'title', 'size thresh')
c=figure
timtopo(posClusterInfo.tMap, EEG.chanlocs, 'title', 'no thresh')

% you should be able to see the clusters in the topo maps... for anything that i've
% run the clusters look exactly the same using a mass threshold as they do
% using the size threshold.



