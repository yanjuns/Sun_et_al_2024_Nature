function [traindata, testdata] = bayes_decoder2D_xsession(neuronIndividualsf, behavIndividualsf, thresh, ts, tts, twostep, K, placecell, KO, KOcell, frame_avg, tportion, binsize, binWidth, tau)
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
% Yanjun Sun, Stanford University, 11/11/2019
% updated with 2 step Bayesian option, 4/13/2020
% Ref: Kechen Zhang et al., J Neurophysi, 1998

if ~exist('binWidth', 'var') || isempty(binWidth)
    binWidth = 1.7; %cm
end
if ~exist('binsize', 'var') || isempty(binsize)
    binsize = 2;
end
% if ~exist('filepath', 'var') || isempty(filepath)
%     filepath = pwd;
% end
if ~exist('tportion', 'var') || isempty(tportion) %portion of data take to train Bayes in the training set
    tportion = 1;
end
if ~exist('frame_avg', 'var') || isempty(frame_avg)
    frame_avg = 12;
end
if ~exist('placecell', 'var') || isempty(placecell)
    placecell = (1:length(thresh))';
end
if ~exist('KO', 'var') || isempty(KO)
    KO = false;
end
%% preprocess training dataset to get binarized neuronal spikes
traindata = struct;
if isscalar(ts)
    traindata.neuron = neuronIndividualsf{1,ts}.S > repmat(thresh,1,size(neuronIndividualsf{1,ts}.S,2));
    traindata.neuron = double(traindata.neuron(placecell,:)');
    if isempty(neuronIndividualsf{1,ts}.pos)
        traindata.position = interp1(behavIndividualsf{1,ts}.time,...
            behavIndividualsf{1,ts}.position, neuronIndividualsf{1,ts}.time);
    else
        traindata.position = neuronIndividualsf{1,ts}.pos;
    end
    if isnan(traindata.position(end,:))
        traindata.position(end,:) = behavIndividualsf{1,ts}.position(end,:);
    elseif isnan(traindata.position(1,:))
        traindata.position(1,:) = behavIndividualsf{1,ts}.position(1,:);
    end
else %if ts is more than one session, combine data from multiple sessions
    tneuron = {}; tposition = {};
    for ii = 1:length(ts)
        trainneuron = neuronIndividualsf{1,ts(ii)}.S > repmat(thresh,1,size(neuronIndividualsf{1,ts(ii)}.S,2));
        tneuron{ii} = double(trainneuron(placecell,:)');
        if isempty(neuronIndividualsf{1,ts(ii)}.pos)
            tposition{ii} = interp1(behavIndividualsf{1,ts(ii)}.time,...
                behavIndividualsf{1,ts(ii)}.position, neuronIndividualsf{1,ts(ii)}.time);
        else
            tposition{ii} = neuronIndividualsf{1,ts(ii)}.pos;
        end
        if isnan(tposition{ii}(end,:))
            tposition{ii}(end,:) = behavIndividualsf{1,ts(ii)}.position(end,:);
        elseif isnan(tposition{ii}(1,:))
            tposition{ii}(1,:) = behavIndividualsf{1,ts(ii)}.position(1,:);
        end
    end
    traindata.neuron = vertcat(tneuron{:});
    traindata.position = vertcat(tposition{:});
end
%% bin and vectorize 2d position
xAxis = 0:binsize:ceil(max(traindata.position(:,1)+binsize)); %binedges of x position
yAxis = 0:binsize:ceil(max(traindata.position(:,2)+binsize)); %binedges of y position
[N,~,~,binY,binX] = histcounts2(traindata.position(:,2), traindata.position(:,1), yAxis, xAxis);
traindata.positionBin = [binX, binY];
position_vec = 1:(size(N,1)*size(N,2));
position_mat = reshape(position_vec, [size(N,1),size(N,2)]);
traindata.positionMat = position_mat;
traindata.positionVec = [];
for ii = 1:length(traindata.positionBin)
    traindata.positionVec(ii,1) = position_mat(traindata.positionBin(ii,2), traindata.positionBin(ii,1));
end
%% downsample behavior and average neuron data using a time window
ds = 1:frame_avg:length(traindata.positionVec);
traindata.positionVec_ds = traindata.positionVec(ds,:);
traindata.positionBin_ds = traindata.positionBin(ds,:);
fun = @(block_struct) sum(block_struct.data);
traindata.neuron_ds = blockproc(traindata.neuron, [frame_avg, size(traindata.neuron,2)], fun);
%% prepare for test data
testdata = struct;
if KO
    testdata.neuron = neuronIndividualsf{1,tts}.S > repmat(thresh,1,size(neuronIndividualsf{1,tts}.S,2));
    testdata.neuron = double(testdata.neuron(placecell,:)');
    testdata.neuron(:,KOcell) = 0;
else
    testdata.neuron = neuronIndividualsf{1,tts}.S > repmat(thresh,1,size(neuronIndividualsf{1,tts}.S,2));
    testdata.neuron = double(testdata.neuron(placecell,:)');
end
if isempty(neuronIndividualsf{1,tts}.pos)
    testdata.position = interp1(behavIndividualsf{1,tts}.time,...
        behavIndividualsf{1,tts}.position, neuronIndividualsf{1,tts}.time);
else
    testdata.position = neuronIndividualsf{1,tts}.pos;
end
if isnan(testdata.position(end,:))
    testdata.position(end,:) = behavIndividualsf{1,tts}.position(end,:);
elseif isnan(testdata.position(1,:))
    testdata.position(1,:) = behavIndividualsf{1,tts}.position(1,:);
end
% testdata.speed = interp1(behavIndividualsf{1,tts}.time,...
%     behavIndividualsf{1,tts}.speed, neuronIndividualsf{1,tts}.time);

%bin the test data
% binsize = 2;
xAxis = 0:binsize:ceil(max(traindata.position(:,1)+binsize)); %binedges of x position
yAxis = 0:binsize:ceil(max(traindata.position(:,2)+binsize)); %binedges of y position
if max(testdata.position(:,1)) > max(xAxis)
    idxO = find(testdata.position(:,1) > max(xAxis));
    testdata.position(idxO,1) = max(xAxis);
end
if max(testdata.position(:,2)) > max(yAxis)
    idxO = find(testdata.position(:,2) > max(yAxis));
    testdata.position(idxO,2) = max(yAxis);
end
[N,~,~,binY,binX] = histcounts2(testdata.position(:,2), testdata.position(:,1), yAxis, xAxis);
testdata.positionBin = [binX, binY];
% position_vec = 1:(size(N,1)*size(N,2));
% position_mat = reshape(position_vec, [size(N,1),size(N,2)]);
testdata.positionMat = position_mat;
testdata.positionVec = [];
for ii = 1:length(testdata.positionBin)
    testdata.positionVec(ii,1) = position_mat(testdata.positionBin(ii,2), testdata.positionBin(ii,1));
end
ds2 = 1:frame_avg:length(testdata.positionVec);
testdata.positionVec_ds = testdata.positionVec(ds2,:);
testdata.positionBin_ds = testdata.positionBin(ds2,:);
% testdata.instspeed = testdata.speed(ds2,:);
fun = @(block_struct) sum(block_struct.data);
testdata.neuron_ds = blockproc(testdata.neuron, [frame_avg, size(testdata.neuron,2)], fun);

%% train the decoder and make predictions
if tportion < 1
    traindata.subset = round(size(traindata.neuron_ds,1)*tportion);
    Mdl = fitcnb(traindata.neuron_ds(1:traindata.subset,:), traindata.positionVec_ds(1:traindata.subset,:)', 'DistributionNames', 'mn');
%     [label,Posterior,Cost] = predict(Mdl,traindata.neuron_ds);
else
    Mdl = fitcnb(traindata.neuron_ds, traindata.positionVec_ds', 'DistributionNames', 'mn');
end
[label,Posterior,Cost] = predict(Mdl,testdata.neuron_ds);
%decode postion in 2D
decodePosition = [];
for ii = 1:length(label)
    [decodeposY, decodeposX] = find(testdata.positionMat == label(ii));
    decodePosition = [decodePosition; [decodeposX, decodeposY]];
end
testdata.decodePosition = decodePosition;
testdata.decodePositionVec = label;
testdata.posterior = Posterior;
testdata.Mdl = Mdl;
%% perform 2 step Bayesian decoding by adding the constraint of animals location probability
if ~exist('tau', 'var') || isempty(tau)
    tau = frame_avg/15;
end
if ~exist('K', 'var') || isempty(K)
    K = 3;
end
if twostep
    %calculate the instantaneous speed after downsampling
    positionBin_ds = [testdata.positionBin_ds;testdata.positionBin(end,:)];
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
%     K = 1; % a scale factor of sigma, see Kechen Zhang et al., J Neurophysi, 1998, equation (44)
    for n = 1:length(testdata.positionBin_ds)
        for ii = 1:size(testdata.positionMat,1)
            for jj = 1: size(testdata.positionMat,2)
                p = exp(-norm([jj,ii] - testdata.decodePosition(n,:)).^2/(2*(K*testdata.instspeed(n)).^2));
                transitionMat_raw(ii,jj,n) = p;
            end
        end
    end
    testdata.transitionMat_raw = transitionMat_raw;
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
    testdata.transitionVec = transitionVec;
    Posterior2step = testdata.posterior .* transitionVec;
%     Posterior2step = log(testdata.posterior) + log(transitionVec);
    testdata.posterior2step = Posterior2step;
    % find the decoded position
    Label2step = [];
    for ii = 1:size(Posterior2step,1)
        [~,idx] = max(Posterior2step(ii,:));
        Label2step(ii,:) = Mdl.ClassNames(idx);
    end
    testdata.decodePositionVec2step = Label2step;
    decodePosition2step = [];
    for ii = 1:length(label)
        [decodeposY, decodeposX] = find(testdata.positionMat == Label2step(ii));
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
decode_error2 = NaN(size(testdata.positionMat,1),size(testdata.positionMat,2),length(testdata.decodeError));
for ii = 1:length(testdata.decodeError)
    decode_error2(testdata.positionBin_ds(ii,2),testdata.positionBin_ds(ii,1),ii) = testdata.decodeError(ii);
end
testdata.decodeError2 = decode_error2;
if twostep
    decode_error2_2step = NaN(size(testdata.positionMat,1),size(testdata.positionMat,2),length(testdata.decodeError2step));
    for ii = 1:length(testdata.decodeError2step)
        decode_error2_2step(testdata.positionBin_ds(ii,2),testdata.positionBin_ds(ii,1),ii) = testdata.decodeError2step(ii);
    end
    testdata.decodeError2_2step = decode_error2_2step;
end

end

