function [PlaceCells,Tinfo,infoScoreThresh] = permuting_spikes(neuron,behav,thresh,smoothocc,deltaTall,occThresh,nboot,binsize,Stype,Method)
% This function aims to identify place cells by shuffling neuron
% spikes over a course of random spike shift. Get a distribution of shuffled
% information scores and make 95% of the distribution as threshold for
% place cells.
% Input argument
% filepath: path of filtered neuron var for each individual sessions
% occThresh: minimum time animal spend in each bin. Usually 0.1-0.2s
% nboot: times of shuffling
% binsize: 1 by default
% Stype: shuffle type, based on either info/sec or info/spike, or both
% Method: individual threshold for each cell vs. universal threshold for
% all cells
% Note: Pay attention if calculating by using neuron.trace or neuron.S in line 26 and line 55
% Developed by Suoqin Jin, modified by Yanjun Sun 10/19/18, 5/1/19, 9/9/19
if ~exist('smoothocc','var') || isempty(smoothocc)
    smoothocc = false;
end
if ~exist('occThresh','var') || isempty(occThresh)
    occThresh = 0.1;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 1000;
end
if ~exist('deltaTall','var') || isempty(deltaTall)
    deltaTall = randi([round(0.05*length(neuron.time)),length(neuron.time)-round(0.05*length(neuron.time))],nboot,1); %msec
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end
if ~exist('Stype','var') || isempty(Stype)
    Stype = 'both';
end
if ~exist('Method','var') || isempty(Method)
    Method = 'individual';
end
% calculate regular spatial information
[firingrateAll,countAll,countTime,smoccu] = calculate_firing_ratemap(neuron,behav,thresh,binsize,smoothocc);
infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = zeros(length(firingrateAll),1);
meanfr = [];
for jj = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{1,jj}))/sum(sum(countTime));
    if smoothocc
        [infoPerSecondnull(jj), infoPerSpikenull(jj)] = Doug_spatialInfo(firingrateAll{jj},MeanFiringRateAll,smoccu,occThresh);
    else
        [infoPerSecondnull(jj), infoPerSpikenull(jj)] = Doug_spatialInfo(firingrateAll{jj},MeanFiringRateAll,countTime,occThresh);
    end
    meanfr(jj,1) = MeanFiringRateAll;
end
% calculate peak firing rate (optional)
ratemap = cell(size(firingrateAll,1), size(firingrateAll,2));
for ii = 1:length(firingrateAll)
    ratemap{ii} = filter2DMatrices(firingrateAll{ii}, 1);
end
peakfr = NaN(length(ratemap),1);
for m = 1:length(ratemap)
    peakfr(m,1) = max(max(ratemap{m}));
end

% calculate shuffled spatial information
infoPerSecondboot = zeros(length(firingrateAll),nboot);
infoPerSpikeboot = zeros(length(firingrateAll),nboot);
for nE = 1:nboot
    deltaT = deltaTall(nE);
    neuronboot = neuron;
    timeboot = (1:length(neuronboot.time)) + deltaT;
    idx = timeboot > length(neuronboot.time);
    timeboot(idx) = timeboot(idx) - length(neuronboot.time);
    [~,index] = sort(timeboot);
    neuronboot.S = neuronboot.S(:,index);
    neuronboot.trace = neuronboot.trace(:,index);
    [firingrateAll,countAll,countTime,smoccu] = calculate_firing_ratemap(neuronboot,behav,thresh,binsize,smoothocc);
    for jj = 1:length(firingrateAll)
        MeanFiringRateAll= sum(sum(countAll{1,jj}))/sum(sum(countTime));
        if smoothocc
            [infoPerSecondboot(jj,nE), infoPerSpikeboot(jj,nE)] = Doug_spatialInfo(firingrateAll{jj},MeanFiringRateAll,smoccu,occThresh);
        else
            [infoPerSecondboot(jj,nE), infoPerSpikeboot(jj,nE)] = Doug_spatialInfo(firingrateAll{jj},MeanFiringRateAll,countTime,occThresh);
        end
    end
end

switch Stype
    case {'info/sec', 'Info/Sec', 'bit/sec', 'Bit/Sec'}
        infoScore = infoPerSecondboot;
        if strcmp(Method, 'individual')
            infoScoreThresh = NaN(size(thresh,1),1);
            for ii = 1:size(infoScore,1)
            infoScoreThresh(ii,1) = quantile(infoScore(ii,:),0.95);
            end
            PlaceCells = find(infoPerSecondnull > infoScoreThresh);
            allneuron = [1:length(infoPerSecondnull)]';
            Lia = ismember(allneuron,PlaceCells);
            Tinfo = table(allneuron,infoPerSecondnull,Lia,'VariableNames',{'neuron','infoPerSecond','PlaceCell'});
            Tinfo = sortrows(Tinfo,{'infoPerSecond'},{'descend'});
        else
            infoScoreThresh = quantile(infoScore(:),0.95);
            figure
            histogram(infoScore(:),'Normalization','probability');
            hold on
            ylim = get(gca,'ylim');
            line([infoScoreThresh infoScoreThresh],get(gca,'ylim'),'LineStyle','--','Color','r','LineWidth',1)
            xlim([min(infoScore(:)),max(infoScore(:))])
            text(infoScoreThresh+0.05,range(ylim)/2,['Thresh = ',num2str(infoScoreThresh,'%.2f')],'Color','red','FontSize',8)

            Tinfo = table([1:length(infoPerSecondnull)]',infoPerSecondnull,'VariableNames',{'neuron','infoPerSecond'});
            Tinfo = sortrows(Tinfo,{'infoPerSecond'},{'descend'});
            PlaceCells = Tinfo.neuron(Tinfo.infoPerSecond > infoScoreThresh);
        end

    case {'info/spike','Info/Spike','info/event','Info/Event','bit/spike','Bit/Spike'};
        infoScore = infoPerSpikeboot;
        if strcmp(Method, 'individual')
            infoScoreThresh = NaN(size(thresh,1),1);
            for ii = 1:size(infoScore,1);
            infoScoreThresh(ii,1) = quantile(infoScore(ii,:),0.95);
            end
            PlaceCells = find(infoPerSpikenull > infoScoreThresh);
            allneuron = [1:length(infoPerSpikenull)]';
            Lia = ismember(allneuron,PlaceCells);
            Tinfo = table(allneuron,infoPerSpikenull,Lia,'VariableNames',{'neuron','infoPerSpike','PlaceCell'});
            Tinfo = sortrows(Tinfo,{'infoPerSpike'},{'descend'});
        else
        infoScoreThresh = quantile(infoScore(:),0.95);
        figure
        histogram(infoScore(:),'Normalization','probability');
        hold on
        ylim = get(gca,'ylim');
        line([infoScoreThresh infoScoreThresh],get(gca,'ylim'),'LineStyle','--','Color','r','LineWidth',1)
        xlim([min(infoScore(:)),max(infoScore(:))])
        text(infoScoreThresh+0.05,range(ylim)/2,['Thresh = ',num2str(infoScoreThresh,'%.2f')],'Color','red','FontSize',8)

        Tinfo = table([1:length(infoPerSpikenull)]',infoPerSpikenull,'VariableNames',{'neuron','infoPerSpike'});
        Tinfo = sortrows(Tinfo,{'infoPerSpike'},{'descend'});
        PlaceCells = Tinfo.neuron(Tinfo.infoPerSpike > infoScoreThresh);
        end
        
    case {'both','Both'};
        infoScore1 = infoPerSecondboot;
        infoScoreThresh1 = NaN(size(thresh,1),1);
        for ii = 1:size(infoScore1,1);
        infoScoreThresh1(ii,1) = quantile(infoScore1(ii,:),0.95);
        end
        PlaceCells1 = find(infoPerSecondnull > infoScoreThresh1);
        allneuron1 = [1:length(infoPerSecondnull)]';
        Lia1 = ismember(allneuron1,PlaceCells1);
%         Tinfo1 = table(allneuron1,infoPerSecondnull,Lia1,'VariableNames',{'neuron','infoPerSecond','PlaceCell'});
%         Tinfo1 = sortrows(Tinfo1,{'infoPerSecond'},{'descend'});
        
        infoScore2 = infoPerSpikeboot;
        infoScoreThresh2 = NaN(size(thresh,1),1);
        for ii = 1:size(infoScore2,1);
        infoScoreThresh2(ii,1) = quantile(infoScore2(ii,:),0.95);
        end
        PlaceCells2 = find(infoPerSpikenull > infoScoreThresh2);
        allneuron2 = [1:length(infoPerSpikenull)]';
        Lia2 = ismember(allneuron2,PlaceCells2);
%         Tinfo2 = table(allneuron2,infoPerSpikenull,Lia2,'VariableNames',{'neuron','infoPerSpike','PlaceCell'});
%         Tinfo2 = sortrows(Tinfo2,{'infoPerSpike'},{'descend'});

        %output all the information into a table
        Tinfo = table(allneuron2,infoPerSecondnull,Lia1,infoPerSpikenull,Lia2,meanfr,peakfr,...
            'VariableNames',{'neuron','bitpsec','bitpsec_pcidx','bitpspike','bitpspike_pcidx','meanfr','peakfr'});
        PlaceCells = {}; infoScoreThresh = {};
        PlaceCells{1} = PlaceCells1; PlaceCells{2} = PlaceCells2;
        infoScoreThresh{1} = infoScoreThresh1; infoScoreThresh{2} = infoScoreThresh2;
end
end
