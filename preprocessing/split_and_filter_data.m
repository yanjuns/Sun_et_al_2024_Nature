%% Split and filter neuron and behav data
 % This process contains multiple sections to split combine processed
 % longitudinal neuron data to their originally defined sessions based
 % on their frames numbers; 
 % then filter slow speed data for both neuron and behav vars;
 % finally this process will define a universal amplitude threshold for
 % Calcium signals. 
 % 1/22/2019 by Yanjun Sun
 % 6/7/2022 updated by Yanjun Sun
%% Get frame number and session information
filepath = pwd;
d = dir(filepath);
n = {d.name}';
load('sessions.mat')
splitname = strsplit(sessions{1},'_');
mouse = splitname{1};
%% load the varibale neuron, add frames information for each session.
nid = contains(n,'Neuron.mat');
filename = cell2mat(n(nid));
load (fullfile(filepath, filename)); 
neuron.num2read = [sum(frames),frames];
save(fullfile(filepath, filename), 'neuron', '-v7.3', '-nocompression'); %Save the neuron variable update again in order to run

%% Input sessions information and split neuron data. 
neuronIndividuals = splitMultiConditionNeuronData(neuron, filename, sessions, 0); %last arg for miniscope_pos

%% 2nd step: scale behavior positions and filter out data associated with low running speed
scale_behav(filepath)
preprocessingData_batch(filepath, filename) %Suoqin's version better

%% 3rd step: get a universal neuron threshold
findneuron = dir('*_Neuron_filter.mat');
load (findneuron.name); % And load behav if under different names
thresh = determiningFiringEventThresh(neuron,'std','S'); %determine the neuron firing threshold, 3SD
save thresh.mat thresh;

% datapath = pwd;
% [noisethresh,noise] = estimate_noisethresh(datapath, 'C');

%% draw spatial footprint of neurons
% findfile = dir('*_Neuron.mat');
% load(findfile.name);
% [footprint]=spatial_footprint_calculation(neuron,0.4);
% figure;
% imagesc(footprint)
% axis image
% colormap pink