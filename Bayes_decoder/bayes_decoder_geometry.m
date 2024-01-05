function decodeResults = bayes_decoder_geometry(datapath,condition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('condition','var') || isempty(condition)
    condition = 'ctr';
end
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd;
end
cd(datapath);
% set parameters
session_name = {'circle','triangle','square','hex'}; %expA
areasize = 4;
frame_avg = 6;
numFolds = 10;
% load the data
load('neuronIndividualsf.mat','neuronIndividualsf')
% load('behavIndividualsf.mat','behavIndividualsf')
load('thresh.mat','thresh')
load('env_geometry.mat','S')
% prepare for data
data = struct;
data.neuron = []; data.position = [];
data.neuron_ds = []; data.position_ds = [];
for ii = 1:length(session_name)
    idx = ismember(S, session_name{ii});
    idx = find(idx == 1);
    idx = idx(1);
    neurodata = neuronIndividualsf{idx}.S;
    behavdata = neuronIndividualsf{idx}.pos;
    ctr = [(max(behavdata(:,1)) - min(behavdata(:,1)))/2 + min(behavdata(:,1)),...
        (max(behavdata(:,2)) - min(behavdata(:,2)))/2 + min(behavdata(:,2))]; %binsize = 2
    if strcmp(condition, 'cor') || strcmp(condition, 'corner')
        if ii == 4
            ctr(:,1) = ctr(:,1) +10; ctr(:,2) = ctr(:,2) +15;
        else
            ctr(:,1) = ctr(:,1) +14; ctr(:,2) = ctr(:,2) +14;
%             ctr(:,1) = ctr(:,1) -14; ctr(:,2) = ctr(:,2) +14;
        end
    end
    roi = find(behavdata(:,1) > ctr(1)-areasize & behavdata(:,1) < ctr(1) +areasize &...
        behavdata(:,2) > ctr(2)-areasize & behavdata(:,2) < ctr(2)+areasize);
    neuron = neurodata(:,roi);
    neuron = neuron > repmat(thresh,1,size(neuron,2));
    neuron = double(neuron');
    position = repmat(ii, size(neurodata(:, roi),2), 1);
    ds = 1:frame_avg:length(position);
    position_ds = position(ds,:);
    fun = @(block_struct) sum(block_struct.data);
    neuron_ds = blockproc(neuron, [frame_avg, size(neuron,2)], fun);
    data.neuron = [data.neuron;neuron];
    data.neuron_ds = [data.neuron_ds;neuron_ds];
    data.position = [data.position;position];
    data.position_ds = [data.position_ds;position_ds];
end

decodeResults = cell(numFolds,1);
for k = 1:numFolds
    sections = numFolds*5;
    % divide the data up into 5*num_folds pieces
    edges = round(linspace(1,numel(data.position_ds)+1,sections+1));
    % get test data from edges - each test data chunk comes from entire session
    test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1];
    train_ind = setdiff(1:numel(data.position_ds),test_ind);

    % training dataset
    neuron_train = data.neuron_ds(train_ind,:);
    behav_train = data.position_ds(train_ind,:);

    % prepare for test data and decoding results
    testdata = struct;
    testdata.neuron_ds = data.neuron_ds(test_ind,:);
    testdata.position_ds = data.position_ds(test_ind,:);

    %% train the decoder and make predictions
    Mdl = fitcnb(neuron_train, behav_train', 'DistributionNames', 'mn');
    [label,Posterior,Cost] = predict(Mdl,testdata.neuron_ds);
    testdata.decodePosition = label;
    testdata.posterior = Posterior;
    testdata.accuracy = numel(find(testdata.position_ds == label))/numel(label);
    decodeResults{k} = testdata; 
end

end