load ('neuronIndividualsf.mat');
load ('behavIndividualsf.mat');
threshscale = [1.1:0.1:2.5];
cell2plot = 415;
smoothocc = true;
for n = 1:length(threshscale)
    load ('thresh.mat');
    thresh = thresh* threshscale(n);
    firingrateAll = cell(1,length(neuronIndividualsf));
    countAll = cell(1,length(neuronIndividualsf));
    countTime = cell(1,length(neuronIndividualsf));
    smoccu = cell(1,length(neuronIndividualsf)); %smoothed countTime
    ratemap = cell(1,length(neuronIndividualsf));
    % Calculate spatial map and mean firing rate for each individual sessions
    for k = 1:length(neuronIndividualsf)
        [firingrateAll{k},countAll{k},countTime{k},smoccu{k}] = calculate_firing_ratemap(neuronIndividualsf{k},behavIndividualsf{k},thresh,2,smoothocc);
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
    plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,ratemap,cell2plot,thresh,'S','auto'); %1:size(neuronIndividualsf{1}.C,1);
    text(80,20,['threshscale = ',num2str(threshscale(n))])
end