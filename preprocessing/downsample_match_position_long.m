function [neuronIndivMch,behavIndivMch] = downsample_match_position_long...
    (neuronIndivLR,behavIndivLR,binsize)
% This function aims to matche two sessions based on behav position
% based on Hardcastle et al (2017) and Kiah Hardcastle's code
% Yanjun Sun, Stanford University, 9/11/2019
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end
%% get parameters and initial the matrices
X = cellfun(@(x) x.position(:,2),behavIndivLR,'uni',0);
Y = cellfun(@(x) x.position(:,1),behavIndivLR,'uni',0);
xmax = max(max(cell2mat(cellfun(@(x) max(x.position(:,2)),behavIndivLR,'uni',0))));
ymax = max(max(cell2mat(cellfun(@(x) max(x.position(:,1)),behavIndivLR,'uni',0))));
%initialize matrices
xAxis = 0:binsize:ceil(xmax+binsize);
yAxis = 0:binsize:ceil(ymax+binsize);
% [timeframe,~,~] = cellfun(@(x,y) histcounts2(x, y, xAxis, yAxis), X,Y, 'uni',0);
%% Determine number of observations in each bin in session1
dsidx = [];
for i = 1:numel(xAxis)-1
    for j = 1:numel(yAxis)-1
        ind = cell(1,length(neuronIndivLR));
        for n = 1:length(neuronIndivLR)
            ind{n} = find(X{n} >= xAxis(i) & X{n} < xAxis(i+1) & Y{n} >= yAxis(j) & Y{n} < yAxis(j+1));
        end
        % Number of observations to keep in each map
        numToKeep = min(cellfun(@numel,ind));
        dsidx = [dsidx; cell2mat(cellfun(@(x)...
            datasample(x, numToKeep, 'replace',false), ind,'uni',0))];
    end
end
%% result output
neuronIndivMch = cell(1,length(neuronIndivLR));
behavIndivMch = cell(1,length(neuronIndivLR));
for ii = 1:length(neuronIndivLR)
    dsidxea = sort(dsidx(:,ii));
    neuron1_mch = struct; behav1_mch = struct;
    %output behavior
    behav1_mch.position = behavIndivLR{ii}.position(dsidxea,:);
    behav1_mch.time = behavIndivLR{ii}.time(dsidxea,:);
    behav1_mch.speed = behavIndivLR{ii}.speed(:,dsidxea);
    %output neuron
    behavtime1 = behavIndivLR{ii}.time;
    idx1 = knnsearch(behavtime1, neuronIndivLR{ii}.time);
    neurontime1 = behavtime1(idx1);
    idxn1 = ismember(neurontime1,behav1_mch.time);
    if isempty(neuronIndivLR{ii}.pos)
        neuronIndivLR{ii}.pos = interp1(behavIndivLR{ii}.time, behavIndivLR{ii}.position, neuronIndivLR{ii}.time);
    end
    neuron1_mch.pos = neuronIndivLR{ii}.pos(idxn1,:);
    neuron1_mch.time = neuronIndivLR{ii}.time(idxn1,:);
    neuron1_mch.S = neuronIndivLR{ii}.S(:,idxn1);
    neuron1_mch.trace = neuronIndivLR{ii}.trace(:,idxn1);
    %final
    neuronIndivMch{ii} = neuron1_mch;
    behavIndivMch{ii} = behav1_mch;
end

end