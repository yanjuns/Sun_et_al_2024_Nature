function [noisethresh,noise] = estimate_noisethresh(datapath, datasource, tscale)
% This function aims to determine the noise threshold of calcium signals in
% neuron.C by estimating the noise level using neuron.C_raw - neuron.C.
% Then use 2 x SD of noise as threshold.

if ~exist('tscale','var') || isempty(tscale)
    tscale = 2;
end
if ~exist('datasource','var') || isempty(datasource)
    datasource = 'trace';
end
cd(datapath);
% find the concatenated and speed filtered neuron variable
% findfile = dir('*_filter.mat');
findfile = dir('*_Neuron.mat');
load(findfile.name);
%estinate the noise and noise threshold
if strcmp(datasource, 'trace')
    noise = neuron.trace_raw - neuron.trace;
else
    noise = neuron.C_raw - neuron.C;
end
noisethresh = tscale * std(noise,[],2);
save('noisethresh.mat','noisethresh','noise')

end