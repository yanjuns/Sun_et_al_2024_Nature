function decodeResults = bayes_decoder2D_withinsession(neuronIndividualsf, behavIndividualsf, thresh, placecell, KO, KOcell, numFolds, twostep, K, frame_avg, binsize, binWidth)
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
if ~exist('K', 'var') || isempty(K)
    K = 2.5;
end
if ~exist('twostep', 'var') || isempty(twostep)
    twostep = true;
end
if ~exist('numFolds', 'var') || isempty(numFolds)
    numFolds = 10;
end
if ~exist('KO', 'var') || isempty(KO)
    KO = false;
end
if ~exist('placecell', 'var') || isempty(placecell)
    placecell = repmat({[1:length(thresh)]'},1,length(neuronIndividualsf));
end

decodeResults = cell(numFolds, length(neuronIndividualsf));
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
        if KO
            testdata.neuron_ds = data.neuron_ds(test_ind,:);
            if iscell(KOcell)
                KOcell_ea = KOcell{n};
                testdata.neuron_ds(:,KOcell_ea) = 0;
            else
                testdata.neuron_ds(:,KOcell) = 0;
            end
        else
            testdata.neuron_ds = data.neuron_ds(test_ind,:);
        end
        testdata.positionVec_ds = data.positionVec_ds(test_ind,:);
        testdata.positionBin_ds = data.positionBin_ds(test_ind,:);
        
        %% train the decoder and make predictions
        Mdl = fitcnb(neuron_train, behav_train', 'DistributionNames', 'mn');
        [label,Posterior,Cost] = predict(Mdl,testdata.neuron_ds);
        %decode postion in 2D
        decodePosition = [];
        for ii = 1:length(label)
            [decodeposY, decodeposX] = find(data.positionMat == label(ii));
            decodePosition = [decodePosition; [decodeposX, decodeposY]];
        end
        testdata.decodePosition = decodePosition;
        testdata.decodePositionVec = label;
        testdata.posterior = Posterior;
        
        %% perform 2 step Bayesian decoding by adding the constraint of animals location probability
        tau = frame_avg/15;
        if twostep
            %calculate the instantaneous speed after downsampling
            positionBin_ds = [testdata.positionBin_ds; data.positionBin(end,:)];
            testdata.instspeed = [];
            for ii = 1:length(testdata.positionBin_ds)
                disttravel = norm(positionBin_ds(ii+1,:) - positionBin_ds(ii,:));
                instspeed = disttravel/tau;
                testdata.instspeed = [testdata.instspeed;instspeed];
            end
            minspeed = min(testdata.instspeed(testdata.instspeed > 0));
            testdata.instspeed(testdata.instspeed == 0) = minspeed;
            %build the transition matrix based on previous decoded position
            transitionMat_raw = [];
            %   K is a scale factor of sigma, see Kechen Zhang et al., J Neurophysi, 1998, equation (44)
            for m = 1:length(testdata.positionBin_ds)
                for ii = 1:size(data.positionMat,1)
                    for jj = 1: size(data.positionMat,2)
                        p = exp(-norm([jj,ii] - testdata.decodePosition(m,:)).^2/(2*(K*testdata.instspeed(m)).^2));
                        transitionMat_raw(ii,jj,m) = p;
                    end
                end
            end
            % testdata.transitionMat_raw = transitionMat_raw;
            %normalize the transition matrix
            %     transitionMat = zeros(size(transitionMat_raw,1),size(transitionMat_raw,2),size(transitionMat_raw,3));
            %     for ii = 1:size(transitionMat_raw,3)
            %         p = transitionMat_raw(:,:,ii);
            %         p = p./sum(sum(p));
            %         transitionMat(:,:,ii) = p;
            %     end
            %     testdata.transitionMat = transitionMat;
            %vectorize the transition matrix
            transitionVec = [];
            for ii = 1:size(transitionMat_raw,3)
                p = squeeze(transitionMat_raw(:,:,ii));
                transitionVec(ii,:) = p(:);
            end
            transitionVec = transitionVec(:,Mdl.ClassNames);
            transitionVec = transitionVec ./ sum(transitionVec,2);% normalize each vector
            firstframe = ones(1, size(transitionVec,2));% first decoding position is independent on the previous one
            transitionVec = [firstframe; transitionVec];% adding ones to ahead of all the datapoint
            transitionVec = transitionVec(1:end-1,:);% remove the last datapoint
            % testdata.transitionVec = transitionVec;
            Posterior2step = testdata.posterior .* transitionVec;
            % Posterior2step = log(testdata.posterior) + log(transitionVec);
            testdata.posterior2step = Posterior2step;
            % find the decoded position
            Label2step = [];
            for ii = 1:size(Posterior2step,1)
                [~,idxx] = max(Posterior2step(ii,:));
                Label2step(ii,:) = Mdl.ClassNames(idxx);
            end
            testdata.decodePositionVec2step = Label2step;
            decodePosition2step = [];
            for ii = 1:length(label)
                [decodeposY, decodeposX] = find(data.positionMat == Label2step(ii));
                decodePosition2step = [decodePosition2step; [decodeposX, decodeposY]];
            end
            testdata.decodePosition2step = decodePosition2step;
            %quick check the improvment of 2step bayes on eliminating the erratic jumping
            %     figure;
            %     subplot(2,1,1)
            %     plot(testdata.positionVec_ds)
            %     hold on
            %     plot(testdata.decodePositionVec)
            %     subplot(2,1,2)
            %     plot(testdata.positionVec_ds)
            %     hold on
            %     plot(testdata.decodePositionVec2step)
        end
        
        %% caculate decoding error
        decode_error = [];
        for ii = 1:length(decodePosition)
            decode_error(ii,:) = norm(testdata.positionBin_ds(ii,:) - decodePosition(ii,:));
        end
        decode_error = decode_error * binWidth;
        testdata.decodeError = decode_error;
        if twostep
            decode_error2step = [];
            for ii = 1:length(decodePosition)
                decode_error2step(ii,:) = norm(testdata.positionBin_ds(ii,:) - testdata.decodePosition2step(ii,:));
            end
            decode_error2step = decode_error2step * binWidth;
            testdata.decodeError2step = decode_error2step;
        end
        % decode error in 2d position map
        decode_error2 = NaN(size(data.positionMat,1),size(data.positionMat,2),length(testdata.decodeError));
        for ii = 1:length(testdata.decodeError)
            decode_error2(testdata.positionBin_ds(ii,2),testdata.positionBin_ds(ii,1),ii) = testdata.decodeError(ii);
        end
        testdata.decodeError2 = decode_error2;
        if twostep
            decode_error2_2step = NaN(size(data.positionMat,1),size(data.positionMat,2),length(testdata.decodeError2step));
            for ii = 1:length(testdata.decodeError2step)
                decode_error2_2step(testdata.positionBin_ds(ii,2),testdata.positionBin_ds(ii,1),ii) = testdata.decodeError2step(ii);
            end
            testdata.decodeError2_2step = decode_error2_2step;
        end
        
        %% Output final data
        decodeResults{k,n} = testdata;
    end
end

