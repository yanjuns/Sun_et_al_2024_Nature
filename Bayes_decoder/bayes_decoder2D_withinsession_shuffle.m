function decodeShuffle = bayes_decoder2D_withinsession_shuffle(neuronIndividualsf, behavIndividualsf, thresh, placecell, nboot, numFolds, frame_avg, binsize, binWidth)
% This function takes calcium imaging data to feed in a naive bayes
% classifier in order to use population neuronal activities to predict
% animals location.
% Input args:
% ts: session number of training data set, if more than 2 elements, will
% combine all the sessions into training set
% tts: session number of testing data set
% frame_avg: number of frames for temporal gathering of the data
% tportion: training portion if training and testing sessions are the same
% if twostep is true, the code will perform the 2 step Bayesian decoding by
% adding a continuity constraint of animals behavior into account.
% K is a scaling factor of sigma on transition matrix, usually 2 to 3 based
% on my testing situation.
% Output:
% traindata: shows the data used for training
% testdata: include the testdata and all the decoded results
% Ref: Kechen Zhang et al., J Neurophysi, 1998
% Yanjun Sun, Stanford University, 11/11/2019
% updated with 2 step Bayesian option, 4/13/2020
% updated within session condition with cross valiation 12/24/22

if ~exist('binWidth', 'var') || isempty(binWidth)
    binWidth = 1.6; %cm
end
if ~exist('binsize', 'var') || isempty(binsize)
    binsize = 2;
end
if ~exist('frame_avg', 'var') || isempty(frame_avg)
    frame_avg = 12;
end
if ~exist('numFolds', 'var') || isempty(numFolds)
    numFolds = 10;
end
if ~exist('nboot', 'var') || isempty(nboot)
    nboot = 100;
end
if ~exist('placecell', 'var') || isempty(placecell)
    placecell = repmat({[1:length(thresh)]'},1,length(neuronIndividualsf));
end

decodeShuffle = cell(numFolds, length(neuronIndividualsf));
for n = 1:length(neuronIndividualsf)
    %% prepare for training data
    %preprocess the dataset to get binarized neuronal spikes
    data = struct;
    data.neuron = neuronIndividualsf{1,n}.S > repmat(thresh,1,size(neuronIndividualsf{1,n}.S,2));
    data.neuron = double(data.neuron(placecell{n},:)');
    if isempty(neuronIndividualsf{1,n}.pos)
        data.position = interp1(behavIndividualsf{1,n}.time,...
            behavIndividualsf{1,n}.position, neuronIndividualsf{1,n}.time);
    else
        data.position = neuronIndividualsf{1,n}.pos;
    end
    if any(isnan(data.position(:,1)))
        idx = find(isnan(data.position(:,1)));
        didx = diff(idx);
        if ~isempty(didx) && any(didx ~= 1)
            ind = find(didx ~= 1);
            data.position(idx(1:ind),:) = repmat(data.position(idx(ind)+1,:),size(idx(1:ind),1),1);
            data.position(idx(ind+1:end),:) = repmat(data.position(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
        else
            if idx(1) == 1
                data.position(idx,:) = repmat(data.position(idx(end)+1,:),size(idx,1),1);
            else
                data.position(idx,:) = repmat(data.position(idx(1)-1,:),size(idx,1),1);
            end
        end
    end
    
    %bin and vectorize 2d position
    xAxis = 0:binsize:ceil(max(data.position(:,1)+binsize)); %binedges of x position
    yAxis = 0:binsize:ceil(max(data.position(:,2)+binsize)); %binedges of y position
    [N,~,~,binY,binX] = histcounts2(data.position(:,2), data.position(:,1), yAxis, xAxis);
    data.positionBin = [binX, binY];
    position_vec = 1:(size(N,1)*size(N,2));
    position_mat = reshape(position_vec, [size(N,1),size(N,2)]);
    data.positionMat = position_mat;
    data.positionVec = [];
    for ii = 1:length(data.positionBin)
        data.positionVec(ii,1) = position_mat(data.positionBin(ii,2), data.positionBin(ii,1));
    end
    
    %downsample behavior and average neuron data using a time window
    ds = 1:frame_avg:length(data.positionVec);
    data.positionVec_ds = data.positionVec(ds,:);
    data.positionBin_ds = data.positionBin(ds,:);
    fun = @(block_struct) sum(block_struct.data);
    data.neuron_ds = blockproc(data.neuron, [frame_avg, size(data.neuron,2)], fun);
    
    %% train and predct with a cross validation
    for k = 1:numFolds
        sections = numFolds*5;
        % divide the data up into 5*num_folds pieces
        edges = round(linspace(1,numel(data.positionVec_ds)+1,sections+1));
        % get test data from edges - each test data chunk comes from entire session
        test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
            edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
            edges(k+4*numFolds):edges(k+4*numFolds+1)-1];
        train_ind = setdiff(1:numel(data.positionVec_ds),test_ind);
        
        % training dataset
        neuron_train = data.neuron_ds(train_ind,:);
        behav_train = data.positionVec_ds(train_ind,:);
        
        % prepare for test data and decoding results
        testdata = struct;
        testdata.neuron_ds = data.neuron_ds(test_ind,:);
        testdata.positionVec_ds = data.positionVec_ds(test_ind,:);
        testdata.positionBin_ds = data.positionBin_ds(test_ind,:);
        
        %% train the decoder and make predictions
        Mdl = fitcnb(neuron_train, behav_train', 'DistributionNames', 'mn');
        
        %% Shuffle testdata.neuron_ds
        neuron_ds = testdata.neuron_ds';
        positionMat = data.positionMat;
        positionBin_ds = testdata.positionBin_ds;
        datalength = size(neuron_ds,2);
        deltaT = randi([round(0.05*datalength),datalength-round(0.05*datalength)],nboot,1);
        dEshuffle = [];
        parfor iter = 1:nboot
            dT = deltaT(iter);
            neuron_dsboot = [1:datalength] + dT;
            idx = neuron_dsboot > datalength;
            neuron_dsboot(idx) = neuron_dsboot(idx) - datalength;
            [~,index] = sort(neuron_dsboot); % sort shuffled spikes
            neuron_shuffle = neuron_ds(:,index);
            neuron_shuffle = neuron_shuffle';
            % make predictions using shuffled data
            [label,~,~] = predict(Mdl,neuron_shuffle);
            %decode postion in 2D
            decodePosition = [];
            for ii = 1:length(label)
                [decodeposY, decodeposX] = find(positionMat == label(ii));
                decodePosition = [decodePosition; [decodeposX, decodeposY]];
            end
            % caculate decoding error
            decode_error = [];
            for ii = 1:length(decodePosition)
                decode_error(ii,:) = norm(positionBin_ds(ii,:) - decodePosition(ii,:));
            end
            decode_error = decode_error * binWidth;
            dEshuffle = [dEshuffle;median(decode_error)];
        end
        
        %% Output final data
        decodeShuffle{k,n} = dEshuffle;
    end
end

