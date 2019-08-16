function clusterInfo=getEEG_clusterSize(EEG_dat, connectMat, thresh, tail);

% The goal of this function is to compute the size of thresholded
% clusters in EEG data. This size can be tested against a permutation
% distribution to determine whole brain corrected significance.


% Inputs:
% EEG_dat: subject number X electrode ID X  timepoint (summary statistic from each
% subject should be one column).


% connectMat: N*N matrix of connections between electrodes
% thresh = cluster forming threshold; 
% tail = tail of t-test used to form clusters:
%        'right' = look for clusters of positive effects
%       'left' = look for clusters of negative effects
%       'both' = look for clusters that contain positive or negative
%       effects
%


% Output:

% clusterInfo: Structure containing the following fields:
% clustSizeMap: Electrode X Timepoints map of clusters, each labeled with
% the number of pairs in that cluster

% ID_map: Electrode X Timepoints map of clusters, each labeled with a
% unique identifier. Unfortunately, do to my new and improved algorith...
% the IDs are not necesarily contigous... so some IDs might be skipped...
% but the non-zero IDs will always be unique.

% clustWtMap: Electrode X Timepoints map of clusters, each labeled with
% the "cluster mass." Cluster mass =  mean|(t-stat)|*pairs in cluster.

% tMap: T-statistic from t-test for each electrode*timepoint pair. 



% Algorithm: This is a bit clugey... I was a bit wrongheaded in my original
% approach... then I added some things that have made it faster and (as far
% as I can tell) stil accurate. here is the post-hoc description of whats happening:

% 0) Identify electrode-time pairs that meet the cluster forming threshold
% 1) find the indices of the first such pair. 
% 2) list the neighboring electrodes to that pair.
% 3) look through that list:
    % A) if you find a electrode-time that is significant, add it to the
    % cluster (ie incremement cluster size by 1, set cluster ID
    % accordingly)
    % B) if you find an electrode-time that is already labeled with an ID
    % (ie it is part of another cluster) combine the clusters (sum the
    % cluster sizes and choose the lower cluster ID). 
    % C) if the list is empty, step forward 1 timepoint and create a new
    % list containing all electrodes that were IDed in the previous
    % timestep. If nothing was ID-ed for this cluster on the previous
    % timestep then stop moving forward in time. And instead start walking
    % backward.
    % D) If you have been walking backward and you get to the time of the
    % original significant pair (ie from step 1) then stop moving
    % backward... break from the loop... and go back to step 1 again. 
    
    % I find this set of operations annoyingly complicated... but it seems
    % much faster than the original way i was doing it. 
    



%% HERE we go (The code):
[tmpSigMap, ~, ~, d]=ttest(EEG_dat, [], thresh, tail);
tmpSigMap=squeeze(tmpSigMap);
tmpWt=squeeze(abs(d.tstat));
tMap=squeeze(d.tstat);



ptsRemain=sum(tmpSigMap(:)); % check to see that there is something significant..
ID=0; % initialize ID

ID_map=zeros(size(tmpSigMap));
clustSizeMap=zeros(size(tmpSigMap));
clustWtMap=zeros(size(tmpSigMap));


usedArray=false(size(ID_map,1), 1);

while ptsRemain>0
    ID=max(ID_map(:))+1;
    [I,J]=find(tmpSigMap, 1); % ok, found first significant point:\
    ID_map(I,J)=ID;           % set map to ID
    clustSize=1;              % set cluster size to 1;
    tmpSigMap(I,J)=false;     % don't reuse this point;
    
    listNeighbors=find(connectMat(I,:)); % initial list of neighbors
    
    % empty list of points we've checked out.
    usedArray(listNeighbors)=true; % set the initial set of neighbors to used.
    stillLooking=true;  % keep looking till you are done.
    fwdDone =false; % we haven't run forward through all timepoints
    
    
    firstSig=J;
    
    
    while stillLooking
        if isempty(listNeighbors)&&~fwdDone
            % if we're in the forward search, we want to look for
            % significant points at the next timepoint that match the
            % electroJ=5des that were in this cluster.
            listNeighbors=find(ID_map(:,J)==ID);
            J=J+1;
            usedArray=false(size(ID_map,1), 1);
            usedArray(listNeighbors)=true;
        elseif isempty(listNeighbors)&&fwdDone
            % if we've finished the forward sweep, lets work our way
            % backwards:
            listNeighbors=find(ID_map(:,J)==ID);
            J=J-1;
            usedArray=false(size(ID_map,1), 1);
            usedArray(listNeighbors)=true;
        end
        
        % now that we have a list of neighbors, lets look through them:
        
        if ~isempty(listNeighbors)&ID_map(listNeighbors(1),J)~=0&ID_map(listNeighbors(1),J)~=ID
            % if we've located a cluster... but that cluster has a
            % different ID than we are currently using, lets combine the
            % clusters:
            nID=ID_map(listNeighbors(1),J);
            clustSize=clustSize+clustSizeMap(listNeighbors(1),J);
            ID_map(ID_map==ID)=nID;
            ID=nID;
            
            %disp('combining clusters');
            
        elseif isempty(listNeighbors)&~fwdDone;
            % if we ran past the data, lets put it in reverse:
            fwdDone=true;
        elseif ~isempty(listNeighbors)&tmpSigMap(listNeighbors(1),J)
            I=listNeighbors(1);
            ID_map(I,J)=ID;           % set map to ID
            clustSize=clustSize+1;    % increment cluster size;
            tmpSigMap(I,J)=false;     % don't reuse this point;
            
            % if we haven't already used it, list the new neighbor in
            % our list.
            newNeighbors=find(connectMat(I,:)'&~usedArray&(ID_map(I,J)~=ID));
            listNeighbors=[listNeighbors, newNeighbors];
            usedArray(listNeighbors)=true;
        end
        
        if length(listNeighbors)>1
            listNeighbors=listNeighbors(2:end);
        else % We've dealt with this timepoint
            listNeighbors=[];
            if J==size(ID_map, 2)&&fwdDone==false;
                fwdDone=true;
            elseif J<=firstSig&&fwdDone==true
                stillLooking=false;
            end
            
        end
        
    end
    clustSizeMap(ID_map==ID)=clustSize;
    clustWtMap(ID_map==ID)=clustSize.*nanmean(tmpWt(ID_map==ID));

    ptsRemain=sum(tmpSigMap(:)); % check to see that there is something significant..
    
end


clusterInfo.clustSizeMap=clustSizeMap;
clusterInfo.ID_map=ID_map;
clusterInfo.clustWtMap=clustWtMap;
clusterInfo.tMap=tMap;
