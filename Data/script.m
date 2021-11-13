data=[];
currPath = mfilename('fullpath');
filepath=  [currPath(1:end-7) filesep 'Recordings'];
file_list = { 'SG_B606/SG_B606_0014.mat','PS_B403/PS_B403_0004.mat','PS_B411/PS_B411_0016.mat',...
    'PS_B186/PS_B186_0015.mat','SG_B561/SG_B561_0034.mat','PS_B181/PS_B181_0012.mat',...
    };
for i=1:length(file_list)
    filename = file_list{i};
    load([filepath filesep filename])
    params.lowPassFreq = 500;
    params.spikeSdThresh = 1.5;
    params.clustSep = 0.7;
    spikeTimes = detectSpikes(data,params);
end