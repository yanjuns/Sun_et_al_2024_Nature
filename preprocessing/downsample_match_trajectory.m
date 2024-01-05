function [neuronc, behavc, neuronb, behavb] = downsample_match_trajectory(neuron1,behav1,neuron2,behav2)
% This function aims to match two sessions behavior and neuron activities
% by tightly match the animals running trajectory. The idea is to simulate
% the running trajectory of less data session on the session with
% more data. The core function is knnsearch to find the nearest point in 2d
% position. The resulted running trajectory should be matched at each time
% point.
% Yanjun Sun, Stanford University, 8/1/2020
if isempty(neuron1.pos)
    neuron1.pos = interp1(behav1.time,behav1.position,neuron1.time);
    neuron2.pos = interp1(behav2.time,behav2.position,neuron2.time);
end

if length(neuron1.pos) > length(neuron2.pos)
    neurona = neuron1; neuronb = neuron2;
    behava = behav1; behavb = behav2;
else
    neurona = neuron2; neuronb = neuron1;
    behava = behav2; behavb = behav1;
end

IDX = knnsearch(neurona.pos,neuronb.pos);
neuronc = struct;
neuronc.pos = neurona.pos(IDX,:);
neuronc.S = neurona.S(:, IDX);
neuronc.trace = neurona.trace(:, IDX); % neuron.trace may not be good after this reorganization
neuronc.time = neurona.time(IDX);
[~,I] = sort(neuronc.time); % sort variables according to the time
neuronc.pos = neuronc.pos(I,:);
neuronc.S = neuronc.S(:, I);
neuronc.trace = neuronc.trace(:, I); 
neuronc.time = neuronc.time(I);

idx = knnsearch(behava.position,behavb.position);
behavc = struct;
behavc.position = behava.position(idx,:);
behavc.speed = behava.speed(idx);
behavc.time = behava.time(idx);
[~,I] = sort(behavc.time);
behavc.position = behavc.position(I,:);
behavc.speed = behavc.speed(:, I);
behavc.time = behavc.time(I);

%% check the matching results
% [firingrateAll,countAll,countTime] = calculate_subset_ratemap(neuronb,behavb,thresh);
% figure;
% ratemap = filter2DMatrices(firingrateAll{29},1);
% pcolor(ratemap);
% shading flat
% axis image
% colormap jet
% caxis([0 2])
end

