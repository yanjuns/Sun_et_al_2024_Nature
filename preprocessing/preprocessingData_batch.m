function preprocessingData_batch(filefolder, filename)
% load data
%%% temporal dynamical analysis of Ca2+ imaging data
% addpath(genpath('E:\Google Drive synchronization\projects\Project2_ Xu\Miniscope_Matlab_code'))
% filefolder = 'H:\Yanjun_and_Suoqin\Miniscope data analysis_2018\M3244F\0217ctr_0221CNO_0224postCtrl_align';
% filefolder = pwd;
%% load the individual neuron data and extract the time information in each condition
load(fullfile(filefolder,'neuronIndividuals.mat'));
load(fullfile(filefolder,'behavIndividuals.mat'));
load(fullfile(filefolder,'sessions.mat'));

%% select data for filtering
for neuron_selected = 1:length(neuronIndividuals)
    behav0 = behavIndividuals{neuron_selected};
    neuron0 = neuronIndividuals{neuron_selected};
    % check whether the neuron data has been downsampled
    downsampling = round(length(neuron0.time)/size(neuron0.trace,2));
    if downsampling ~= 1
        neuron0.time = double(neuron0.time);
        neuron0.time = neuron0.time(1:downsampling:end);
    end
    %Calculate indices which satisfy our velocity criterion. Exclude outliers
    thresh_speed_min = 2; thresh_speed_trans = 0;
    outliers_removal = false;
    factor_upper = 1.5;factor_lower = 1.05; % for '0217Ctrl','0221CNO',
    behav0.filterPara.thresh_speed_min = thresh_speed_min;
    behav0.filterPara.thresh_speed_trans = thresh_speed_trans;
    behav0.filterPara.factor_upper = factor_upper;
    behav0.filterPara.factor_lower = factor_lower;
    % frames_pass = behavDataFiltering(behav0,behav0.speed,thresh_speed_min,thresh_speed_trans,outliers_removal,factor_upper,factor_lower);
    frames_pass = speed_filtering(behav0.speed,thresh_speed_min);
    
    % correct the behav data and plot the data
    behavf = behavPlot(behav0,behav0.speed,frames_pass);
    % correct the neuron data and plot the data
    [neuronf,frame_pass_neuron] = neuronDataFiltering(neuron0,frames_pass,behav0); % % Note: this code will change the variables neuron0 and neuronIndividuals. Thus, one can only run once.
    
    %% save data into each individual folders
    behavIndividualsf{neuron_selected} = behavf;
    neuronIndividualsf{neuron_selected} = neuronf;
    % save the filtered data after running all the sessions
    dirname = num2str(sessions{neuron_selected});
    if ~exist(dirname,'dir')
        mkdir (dirname);
    end
    savepath = [filefolder,'\',num2str(sessions{neuron_selected})];
    save (fullfile(savepath, 'Neuronf.mat'), 'neuronf', '-v7.3');
    save (fullfile(savepath, 'Behavf.mat'), 'behavf', '-v7.3');
    savefig (1, fullfile(savepath, 'filterBehav.fig'));
    close all;
    clear neuronf & behavf;
    
end
%% add interpolated position to neuronIndividualsf and save variables
for jj = 1:length(neuronIndividualsf)
    neuronIndividualsf{jj}.pos = interp1(behavIndividualsf{jj}.time,...
        behavIndividualsf{jj}.position,neuronIndividualsf{jj}.time);
end
save('neuronIndividualsf.mat', 'neuronIndividualsf', '-v7.3');
save('behavIndividualsf.mat', 'behavIndividualsf', '-v7.3');

%% reconstruct the combined neuron variable after filtering
% load the combined neuron data file
load(fullfile(filefolder,filename)); %change to correct file name
neuron.C = [];neuron.C_raw = [];neuron.C_df = [];neuron.S = [];neuron.trace = [];neuron.trace_raw = [];neuron.time = [];neuron.num2read = [];
for i = 1:length(neuronIndividualsf);
    neuron.C = [neuron.C,neuronIndividualsf{i}.C];
    neuron.C_raw = [neuron.C_raw,neuronIndividualsf{i}.C_raw];
    %     neuron.C_df = [neuron.C_df,neuronIndividualsf{i}.C_df];
    neuron.S = [neuron.S,neuronIndividualsf{i}.S];
    neuron.trace = [neuron.trace,neuronIndividualsf{i}.trace];
    neuron.trace_raw = [neuron.trace_raw,neuronIndividualsf{i}.trace_raw];
    neuron.time = [neuron.time;neuronIndividualsf{i}.time];
    neuron.num2read = [neuron.num2read,neuronIndividualsf{i}.num2read];
end
neuron.num2read = [sum(neuron.num2read),neuron.num2read];
% save 'M3924F_0716-0726CS_Neuron2.mat' neuron
splitname = strsplit(filename,'.');
savename = splitname{1};
save (fullfile(filefolder, [savename '_filter.mat']), 'neuron', '-v7.3'); %change to correct file name
