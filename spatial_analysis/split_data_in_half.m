function [neuron1,behav1,neuron2,behav2] = split_data_in_half(neuron,behav)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes he
if ~isfield(neuron,'pos') || isempty(neuron.pos)
    neuron.pos = interp1(behav.time,behav.position,neuron.time);
end

dl_neuron = round(length(neuron.time)/2);
dl_behav = round(length(behav.time)/2);

neuron1 = struct;
neuron1.S = neuron.S(:, 1:dl_neuron);
neuron1.time = neuron.time(1:dl_neuron);
neuron1.pos = neuron.pos(1:dl_neuron,:);

neuron2 = struct;
neuron2.S = neuron.S(:, dl_neuron+1:end);
neuron2.time = neuron.time(dl_neuron+1:end);
neuron2.pos = neuron.pos(dl_neuron+1:end,:);

behav1 = struct;
behav1.time = behav.time(1:dl_behav);
behav1.position = behav.position(1:dl_behav,:);
if isfield(behav,'hd') && ~isempty(behav.hd)
    behav1.hd = behav.hd(1:dl_behav);
end

behav2 = struct;
behav2.time = behav.time(dl_behav+1:end);
behav2.position = behav.position(dl_behav+1:end,:);
if isfield(behav,'hd') && ~isempty(behav.hd)
    behav2.hd = behav.hd(dl_behav+1:end);
end

end

