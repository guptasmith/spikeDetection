function spikeTimes = detectSpikes(dataSet, newparams)
% [spikeTimes] = detectSpikes(dataSet, newparams)
% Finds the positions of spikes in a recording with multiple trials
% Output is the set of spike indices in a cell data structure
% Second Derivative threshold : a threshold value for second derviate is set
%                               and peak-points crossing this threshold is
%                               used for further clustering
% Derivative clustering : Peak points are clustered on the basis of their
%                         positions in space defined by maximum of signal,
%                         the first derivative and the second derivative
%                         for each peak
% If you need to change some parameters, provide them in a structure newparams
% Default parameters:
params.Fs = 20000; % The sampling frequency
params.lowPassFreq = 500; % frequencies above this are removed when calculating the derivatives
params.spikeSdThresh = 1.5;
params.numTrials = 7;
params.onSet = 2;    % Used for threshold calculation
params.showDuration = 10;
params.clustSep = .7; % Parameter that governs how separated the clusters should be
params.maxclust = 2*ones(1,size(dataSet,1)); %maxclust parameter for cluster function. Put expected number of spike types.
params.spikeDebug = 0;
params.ssr = 0;

% update the parameters as supplied (only the supplied parameters are changed)
if nargin > 1
    params = updateParams(params,newparams);
end

%% second derivative thresholding
numRows = size(dataSet,1);
featureSetAll = cell(numRows,1);
spikeIndexDerAll = cell(1,numRows);

for dataRowInd=1:numRows
    
    data=dataSet(dataRowInd,:);
    dataLP = filterLowPass(data,params.lowPassFreq,struct('Fs',params.Fs));
    dif=diff(dataLP);
    difdif=diff(diff(dataLP));
    difNeg = (dif<=0);
    difBoundary = diff(difNeg);
    difBoundaryDownIndices = find(difBoundary==-1);  %% indices of all troughs
    difBoundaryUpIndices = find(difBoundary==1);  %%indices of all crests
    
    % detecting spikes in the signal
    sdWin = min(params.Fs*params.onSet,length(difdif));
    sd = std(difdif(1:sdWin));
    thresh=-params.spikeSdThresh*sd;
    significantData=difdif<thresh;  % mark points (as 1) that are below the threshold
    boundaries = diff(significantData); % the down crossing point (at the threshold) of the secondary derivative will be 1 and the upward crossing point will be -1
    boundaryIndices = find(boundaries == 1); % finds the indices of the downward crossing points only.
    boundaryUpIndices = find(boundaries == -1); % finds the indices of the upward crossing points only.
    
    N = length(boundaryIndices);% number of spikes
    NUp = length(boundaryUpIndices);
    
    if N > 0
        if NUp > 0
            % if the first spike happens to be cut at the start time, remove it
            if(boundaryUpIndices(1) < boundaryIndices(1))
                boundaryUpIndices = boundaryUpIndices(2:end);
            end
            % then, if the last spike is cut at the ending time, remove it
            if(length(boundaryUpIndices) < length(boundaryIndices))
                boundaryIndices = boundaryIndices(1:end-1);
                N=N-1;
            end
        end
    end
    
    
    % detect the spike positions
    spikeIndexDer = [];
    if N > 0
        if NUp > 0
            for i = 1:N
                lastTrough = difBoundaryDownIndices(find(boundaryIndices(i) - difBoundaryDownIndices > 0,1,'last'));
                nextTrough = difBoundaryDownIndices(find(difBoundaryDownIndices - boundaryIndices(i) > 0,1));
                [val,index] = max(dataLP(lastTrough:nextTrough));
                index = index + lastTrough - 1;
                if ismember(index-1,difBoundaryUpIndices) & ~ismember(index,difBoundaryDownIndices)
                    spikeIndexDer = [spikeIndexDer, index];
                end
            end
        end
    end
    
    spikeIndexDer = unique(spikeIndexDer);
    if ~isempty(spikeIndexDer)
        spikeIndexDer([1 end]) = []; %% removing first and last identified spikes
        % as some recording has a noise peaks at
        % these location with high derivative values
    end
    
    %% Derivative clustering
    if params.ssr==1
        featureSet = zeros(length(spikeIndexDer),3);
    else
        featureSet = zeros(length(spikeIndexDer),3);
    end
    for i=1:length(spikeIndexDer)
        lastTrough = difBoundaryDownIndices(find(difBoundaryDownIndices ...
            < spikeIndexDer(i), 1, 'last'));
        nextTrough = difBoundaryDownIndices(find(difBoundaryDownIndices ...
            > spikeIndexDer(i), 1));
        if isempty(nextTrough)
            nextTrough = length(difdif);
        end
        
        
        if params.ssr==1
            difVal = abs(min(dif(lastTrough:nextTrough)))+max(dif(lastTrough:nextTrough));
            difdifVal = abs(min(difdif(lastTrough:nextTrough)))+max(difdif(lastTrough:nextTrough));
            halfMax = (dataLP(spikeIndexDer(i))+ max(dataLP(lastTrough),dataLP(nextTrough)))/2;
            fwhm = find(dataLP(spikeIndexDer(i):nextTrough)<=halfMax,1) -...
                find(dataLP(lastTrough:spikeIndexDer(i))>=halfMax,1) + spikeIndexDer(i) - lastTrough;
            heightSpike = min(dataLP(spikeIndexDer(i))-dataLP(lastTrough),dataLP(spikeIndexDer(i))-dataLP(nextTrough));
            featureSet(i,:) = [heightSpike,difVal,difdifVal];
        else
            difdifMin = min(difdif(lastTrough:nextTrough));
            heightSpike = min(dataLP(spikeIndexDer(i))-dataLP(lastTrough),dataLP(spikeIndexDer(i))-dataLP(nextTrough));
            featureSet(i,:) = [dataLP(spikeIndexDer(i)),heightSpike,abs(difdifMin)];
        end
        
    end
    featureSetAll{dataRowInd} = featureSet;
    spikeIndexDerAll{dataRowInd} = spikeIndexDer;
end

[m,std_sett] = mapstd(cell2mat(featureSetAll)');
if isempty(m)
    spikeTimes = cell(1,numRows);
    return
else
    m=m';
end
[val,minDifPointInd] = min(vecnorm(m'-min(m,[],1)'));
z = linkage(m,'average','cityblock');
t = cluster(z,'criterion','distance','maxclust',min(params.maxclust));
clust1Label = t(minDifPointInd);
clust1Set = m(t==clust1Label,:);
clust2Set = setdiff(m,clust1Set,'rows');

distMatrixBwClust  = pdist2(clust1Set,clust2Set);
[minDistBwClust, minDistInd] = min(distMatrixBwClust,[],'all','linear');
[i,j] = ind2sub(size(distMatrixBwClust),minDistInd);
clust1Dist = vecnorm(clust1Set' - clust1Set(i,:)');
clust1Dist(clust1Dist==0)=Inf;
clust1MinDist = min(clust1Dist);
clust2Dist = vecnorm(clust2Set' - clust2Set(j,:)');
clust2Dist(clust2Dist==0)=Inf;
clust2MinDist = min(clust2Dist);
if clust2MinDist==Inf
    clust2MinDist=0;
end


if (minDistBwClust<params.clustSep*(clust1MinDist+clust2MinDist) && size(dataSet,1)>1 ) || (length(unique(params.maxclust))>1)
    if ~all(unique(params.maxclust)==2)
        if params.spikeDebug
            fprintf('Trials are considered separately, using the provided vector values for maxclust.\n');
        end
        for reiterInd=1:size(dataSet,1)
            newparams = params;
            newparams.maxclust = params.maxclust(reiterInd);
            spikeTimes(reiterInd) =  detectSpikes(dataSet(reiterInd,:),newparams);
        end
    else
        if params.spikeDebug
            fprintf(['Trials are considered separately. Same value of 2 used for maxclust for each trial. ',...
                'In case of errors, try using a vector of value for maxclust.\n']);
        end
        for reiterInd=1:size(dataSet,1)
            newparams = params;
            newparams.maxclust = 2;
            spikeTimes(reiterInd) =  detectSpikes(dataSet(reiterInd,:),newparams);
        end
    end
else
    if params.spikeDebug && size(dataSet,1)>1
        if max(params.maxclust)==2
            fprintf(['Trials are consisdered together with default value of 2 for maxclust. ',...
                'In case of errors, try using a different single value for maxclust.\n']);
        else
            fprintf(['Trials are considered together with value of %d for maxclust. ',...
                'In case of errors try splitting the trials by providing vector of values for maxlcust.\n'],params.maxclust(1));
        end
    end
    spikeTimes = cell(1,numRows);
    start_index=0;
    end_index=0;
    for dataRowInd=1:numRows
        spikeIndexDerRow = spikeIndexDerAll{dataRowInd};
        start_index = end_index + 1;
        end_index = end_index + length(spikeIndexDerRow);
        spikeIndex = spikeIndexDerRow(ismember(m(start_index:end_index,:),clust2Set,'rows'));
        spikeTimes{dataRowInd} = spikeIndex/params.Fs;
        
        %         %debug plots
        if params.spikeDebug
            figure
            data = dataSet(dataRowInd,:);
            plot(0:1/params.Fs:size(data,2)/params.Fs-1/params.Fs,data);
            hold on
            scatter(spikeIndex/params.Fs, data(spikeIndex))
            xlabel('Time (ms)');
            ylabel('Membrane Potential (mV)');
            set(gca,'xtick',0:params.showDuration,'xticklabel',0:params.showDuration);
            figure
            hold on
            grid on
            if params.ssr==1
                scatter3(clust1Set(:,1),clust1Set(:,2),clust1Set(:,3),'r')
                scatter3(clust2Set(:,1),clust2Set(:,2),clust2Set(:,3),'b')
                xlabel('Height');
                ylabel('First Derivative');
                zlabel('Second Derivative');
            else
                scatter3(clust1Set(:,1),clust1Set(:,2),clust1Set(:,3),'r')
                scatter3(clust2Set(:,1),clust2Set(:,2),clust2Set(:,3),'b')
                xlabel('Absolute');
                ylabel('Height');
                zlabel('Second Derivative');
            end
        end
    end
end
end
