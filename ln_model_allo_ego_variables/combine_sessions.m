function [neuron,behav] = combine_sessions(neuronIndividualsf,behavIndividualsf,row,col)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
neuron = struct;
behav = struct;

neuron.S = [];
neuron.time = [];
behav.position = [];
behav.positionblue = [];
behav.time = [];
behav.speed = [];
ss = length(col);
for ii = 1:ss
    neuron.S = [neuron.S,neuronIndividualsf{row(ii),col(ii)}.S];
    if isempty(neuron.time)
        neuron.time = [neuron.time;neuronIndividualsf{row(ii),col(ii)}.time];
    else
        neuron.time = [neuron.time;(neuronIndividualsf{row(ii),col(ii)}.time + max(neuron.time))];
    end
    behav.position = [behav.position;behavIndividualsf{row(ii),col(ii)}.position];
    behav.positionblue = [behav.positionblue;behavIndividualsf{row(ii),col(ii)}.positionblue];
    behav.speed = [behav.speed,behavIndividualsf{row(ii),col(ii)}.speed];
    if isempty(behav.time)
        behav.time = [behav.time;behavIndividualsf{row(ii),col(ii)}.time];
    else
        behav.time = [behav.time;(behavIndividualsf{row(ii),col(ii)}.time + max(behav.time))];
    end
end

