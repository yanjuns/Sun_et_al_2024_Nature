% load data
%%% temporal dynamical analysis of Ca2+ imaging data
% addpath(genpath('E:\Google Drive synchronization\projects\Project2_ Xu\Miniscope_Matlab_code'))
% filefolder = 'H:\Yanjun_and_Suoqin\Miniscope data analysis_2018\M3244F\0217ctr_0221CNO_0224postCtrl_align';
filefolder = pwd;
%% load the individual neuron data and extract the time information in each condition
load(fullfile(filefolder,'neuronIndividuals.mat'))
mice = 'M3914F';
% sessions = {'0716C','0716S','0717C','0717S','0718C','0718S','0719C','0719S','0720C','0720S','0721C','0721S',...
%     '0722C','0722S','0723C','0723S','0724C','0724S','0725C','0725S','0726C','0726S'};
sessions = {'0815C','0815S','0817C','0817S','0819C','0819S','0820C','0820S','0822C','0822S','0823C','0823S',...
    '0824C','0824S','0826C','0826S','0827C','0827S','0828C','0828S','0830C','0830S'};
behavIndividuals = cell(1,length(neuronIndividuals));
for i = 1:length(neuronIndividuals)
    behav0 = importdata(fullfile(filefolder,[mice,'_',sessions{i},'_Behav.mat']));
%     behav0 = importdata(fullfile(filefolder,[sessions{i},'_Behav.mat']));
    t = find(diff(behav0.time)<=0);
    while ~isempty(t)
        behav0.time(t+1) = behav0.time(t)+1;
        t = find(diff(behav0.time)<=0);
    end
    dx = [0; diff(behav0.position(:,1))];
    dy = [0; diff(behav0.position(:,2))];
    behav0.speed = sqrt((dx).^2+(dy).^2)/behav0.dt;
    ceil(1/behav0.dt/5)
    behav0.speed = smoothts(behav0.speed','b',ceil(1/behav0.dt/5));
    behavIndividuals{i} = behav0;
end
%% select data for filtering
neuron_selected = 2; %change this number for every session
behav0 = behavIndividuals{neuron_selected};
neuron0 = neuronIndividuals{neuron_selected};
% check whether the neuron data has been downsampled
downsampling = length(neuron0.time)/size(neuron0.trace,2)
if downsampling ~= 1
    neuron0.time = double(neuron0.time);
    neuron0.time = neuron0.time(1:downsampling:end);
end
%Calculate indices which satisfy our velocity criterion. Exclude outliers
thresh_speed_min = 2; thresh_speed_trans = 0;
outliers_removal = false;factor_upper = 1.5;factor_lower = 1.05; % for '0217Ctrl','0221CNO',
% outliers_removal = true;factor_upper = 1.5;factor_lower = 1.25; % for '0224PCtrl'
behav0.filterPara.thresh_speed_min = thresh_speed_min;behav0.filterPara.thresh_speed_trans = thresh_speed_trans;
behav0.filterPara.factor_upper = factor_upper;behav0.filterPara.factor_lower = factor_lower;
frames_pass = behavDataFiltering(behav0,behav0.speed,thresh_speed_min,thresh_speed_trans,outliers_removal,factor_upper,factor_lower);
% correct the behav data and plot the data
behavf = behavPlot(behav0,behav0.speed,frames_pass);
% correct the neuron data and plot the data
[neuronf,frame_pass_neuron] = neuronDataFiltering(neuron0,frames_pass,behav0); % % Note: this code will change the variables neuron0 and neuronIndividuals. Thus, one can only run once.

%% save data
% behavIndividuals{neuron_selected} = behav;
% neuronIndividuals{neuron_selected} = neuron0;
behavIndividualsf{neuron_selected} = behavf;
neuronIndividualsf{neuron_selected} = neuronf;

% save the filtered data after running all the sessions
save Neuronf.mat neuronf
save Behavf.mat behavf

clear neuronf & behavf

%% save after the done all the sessions
save ('neuronIndividualsf.mat', 'neuronIndividualsf', '-v7.3', '-nocompression');
save behavIndividualsf.mat behavIndividualsf;

%% reconstruct the combined neuron variable after filtering
% load the combined neuron data file
load(fullfile(filefolder,'M3924F_0723-0803CS_Neuron2.mat'))

neuron.C = [];neuron.C_raw = [];neuron.C_df = [];neuron.S = [];neuron.trace = [];neuron.trace_raw = [];neuron.time = [];neuron.num2read = [];
for i = 1:length(neuronIndividualsf)
    neuron.C = [neuron.C,neuronIndividualsf{i}.C];
    neuron.C_raw = [neuron.C_raw,neuronIndividualsf{i}.C_raw];
    neuron.C_df = [neuron.C_df,neuronIndividualsf{i}.C_df];
    neuron.S = [neuron.S,neuronIndividualsf{i}.S];
    neuron.trace = [neuron.trace,neuronIndividualsf{i}.trace];
    neuron.trace_raw = [neuron.trace_raw,neuronIndividualsf{i}.trace_raw];
    neuron.time = [neuron.time;neuronIndividualsf{i}.time];
    neuron.num2read = [neuron.num2read,neuronIndividualsf{i}.num2read];
end
neuron.num2read = [sum(neuron.num2read),neuron.num2read];
% save 'M3924F_0716-0726CS_Neuron2.mat' neuron
save ('M3924F_0716-0726CS_filter_Neuron2.mat', 'neuron', '-v7.3', '-nocompression');
