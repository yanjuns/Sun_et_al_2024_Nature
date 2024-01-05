function neuronIndividuals = splitMultiConditionNeuronData(neuron0, dataName, datasetNames, miniscope_pos)
% this function aims to split combine processed neuron data from CNMF-E based
% on frame numbers of each sessions.
% Input vars
% neuron0: CNMF-E generated and num2read updated neuron var
% dataName: file name of the neuron var
% datasetNames: session names or sessions.mat
neuronIndividuals = cell(1,length(neuron0.num2read)-1);
for i = 2:length(neuron0.num2read)
    neuron0 = importdata(dataName); %Change the Neuron var name here
    if i == 2
        start = 1;
    else
        start = sum(neuron0.num2read(2:i-1))+1;
    end
    neuronIndividuals{i-1} = neuron0;
    neuronIndividuals{i-1}.C = neuron0.C(:,start:sum(neuron0.num2read(2:i)));
    neuronIndividuals{i-1}.C_raw = neuron0.C_raw(:,start:sum(neuron0.num2read(2:i)));
    neuronIndividuals{i-1}.S = neuron0.S(:,start:sum(neuron0.num2read(2:i)));
    neuronIndividuals{i-1}.trace = neuron0.trace(:,start:sum(neuron0.num2read(2:i)));
    neuronIndividuals{i-1}.trace_raw = neuron0.trace_raw(:,start:sum(neuron0.num2read(2:i)));
    neuronIndividuals{i-1}.num2read = neuron0.num2read(i);
end
% save neuronIndividuals.mat neuronIndividuals
save ('neuronIndividuals.mat', 'neuronIndividuals', '-v7.3', '-nocompression'); % for saving large vars
% add the time information
for i = 1:length(neuronIndividuals)
    fid=fopen(['timestamp_' datasetNames{i} '.dat'],'r');
    timedata = textscan(fid, '%d%d%d%d', 'HeaderLines', 1);
    if miniscope_pos == 1
        t = timedata{3}(timedata{1} == 1);t(1) = 0;   %% make sure the miniscope is for USB port 1 %%behav cam port 0
    elseif miniscope_pos == 0
        t = timedata{3}(timedata{1} == 0);t(1) = 0;  %% Giocomo lab, miniscope is 0, behav is 1
    end
    idxt = find(diff(t)<=0);
    t(idxt+1) = t(idxt)+1;
    % the first is trainning and the second is testing
    neuronIndividuals{i}.time = double(t);
end
% save neuronIndividuals.mat neuronIndividuals
save ('neuronIndividuals.mat', 'neuronIndividuals', '-v7.3', '-nocompression'); %for saving large vars
%Code for checking time points 
% fid=fopen('Y:\Steve\miniscope_imaging\0216\trial3\timestamp.dat','r');
% timedata = textscan(fid, '%d%d%d%d', 'HeaderLines', 1);
% t = timedata{3}(timedata{1} == 1);t(1) = 0;
% figure;hist(diff(double(t)),100)
% t2 = timedata{3}(timedata{1} == 0);t2(1) = 0;
% figure;hist(diff(double(t2)),100)
