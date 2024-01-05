function [neuron1,behav1,neuron2,behav2] = split_data_interleaved(neuron,behav,sections)
%This function aims to split data into two in an interleaved manner. It
%split the data into equally spaced chunks and combine the even vs. odd
%chunks in the data. 
%yanjuns@stanford.edu
%8/30/2023

if ~exist('sections','var') || isempty(sections)
    sections = 20;
end
if ~isfield(neuron,'pos') || isempty(neuron.pos)
    neuron.pos = interp1(behav.time,behav.position,neuron.time);
end

neuron_edges = round(linspace(1,size(neuron.S,2)+1,sections+1));
behav_edges = round(linspace(1,length(behav.time)+1,sections+1));

neuron1 = struct;
neuron1.S = [];
neuron1.time = [];
neuron1.pos = [];
for ii = 1:2:length(neuron_edges)-1
    neuron1.S = [neuron1.S,neuron.S(:,[neuron_edges(ii):neuron_edges(ii+1)-1])];
    neuron1.time = [neuron1.time; neuron.time([neuron_edges(ii):neuron_edges(ii+1)-1])];
    neuron1.pos = [neuron1.pos; neuron.pos([neuron_edges(ii):neuron_edges(ii+1)-1],:)];
end
neuron2 = struct;
neuron2.S = [];
neuron2.time = [];
neuron2.pos = [];
for ii = 2:2:length(neuron_edges)
    neuron2.S = [neuron2.S,neuron.S(:,[neuron_edges(ii):neuron_edges(ii+1)-1])];
    neuron2.time = [neuron2.time; neuron.time([neuron_edges(ii):neuron_edges(ii+1)-1])];
    neuron2.pos = [neuron2.pos; neuron.pos([neuron_edges(ii):neuron_edges(ii+1)-1],:)];
end


behav1 = struct;
behav1.time = [];
behav1.position = [];
behav1.hd = [];
for ii = 1:2:length(behav_edges)-1
    behav1.time = [behav1.time; behav.time([behav_edges(ii):behav_edges(ii+1)-1])];
    behav1.position = [behav1.position; behav.position([behav_edges(ii):behav_edges(ii+1)-1],:)];
    if isfield(behav,'hd') && ~isempty(behav.hd)
        behav1.hd = [behav1.hd; behav.hd([behav_edges(ii):behav_edges(ii+1)-1],:)];
    end
end
behav2 = struct;
behav2.time = [];
behav2.position = [];
behav2.hd = [];
for ii = 2:2:length(behav_edges)
    behav2.time = [behav2.time; behav.time([behav_edges(ii):behav_edges(ii+1)-1])];
    behav2.position = [behav2.position; behav.position([behav_edges(ii):behav_edges(ii+1)-1],:)];
    if isfield(behav,'hd') && ~isempty(behav.hd)
        behav2.hd = [behav2.hd; behav.hd([behav_edges(ii):behav_edges(ii+1)-1],:)];
    end
end



end

