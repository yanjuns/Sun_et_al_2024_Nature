function generate_spatial_ratemap(datapath, smoothocc)
% This is a general function of calculating the spatial firing ratemap from
% neuron and behav data.
% yanjuns@stanford.edu
if ~exist('smoothocc','var')||isempty(smoothocc)
    smoothocc = true;
end
if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end

load ('neuronIndividualsf.mat','neuronIndividualsf');
load ('behavIndividualsf.mat','behavIndividualsf');
load ('thresh.mat','thresh');
meanFRindividuals = cell(1,length(neuronIndividualsf));
firingrateAll = cell(1,length(neuronIndividualsf));
countAll = cell(1,length(neuronIndividualsf));
countTime = cell(1,length(neuronIndividualsf));
smoccu = cell(1,length(neuronIndividualsf)); %smoothed countTime
ratemap = cell(1,length(neuronIndividualsf));
% Calculate spatial map and mean firing rate for each individual sessions
for k = 1:length(neuronIndividualsf)
    [firingrateAll{k},countAll{k},countTime{k},smoccu{k}] = calculate_firing_ratemap(neuronIndividualsf{k},behavIndividualsf{k},thresh,2,smoothocc);
    for m = 1:size(firingrateAll{1,k},2)
        meanFRindividuals{1,k}(m,1)= sum(sum(countAll{1,k}{1,m}))/sum(sum(countTime{1,k}));
    end
end
% Calculate smoothed firing rate maps
for ii = 1:length(firingrateAll)
    fr = firingrateAll{ii};
    % to remove potential inf in the unsmoothed ratemap
    idx = cellfun(@isinf, fr,'UniformOutput',false);
    for jj = 1:length(fr)
        fr{jj}(idx{jj}) = 0;
    end
    % smoothed ratemap
    ratemap{ii} = cellfun(@(x) filter2DMatrices(x,1), fr, 'uni', 0);
end

save('firingrateAll.mat', 'firingrateAll', 'countAll', 'countTime', 'smoccu', 'ratemap', 'meanFRindividuals', '-v7.3');

end

