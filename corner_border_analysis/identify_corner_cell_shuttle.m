function C = identify_corner_cell_shuttle(datapath)
%This function aims to identify corner cells for each session in each mouse
%yanjuns@stanford.edu, 10/12/22
%C is the final structure that contains all the information

if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end

cd(datapath)
%determine corner cell threshold by shuffling.
load ('NBIndivLR.mat','neuronIndivLR','behavIndivLR','env_coor','ratemap');
load ('thresh.mat','thresh');
neuronIndividualsf = neuronIndivLR;
behavIndividualsf = behavIndivLR;
%compute corner metrics including corner scores
cmetrics = cell(1,length(ratemap));
for n = 1:length(ratemap)
    cmetrics{n} = compute_corner_score(ratemap{n},env_coor{n});
end

cscore_shuffle = cell(1, length(neuronIndividualsf));
parfor ii = 1:length(neuronIndividualsf)
    neuron = neuronIndividualsf{1,ii};
    behav = behavIndividualsf{1,ii};
    if isempty(env_coor{ii})
        cscore_shuffle_ea = [];
    else
        cscore_shuffle_ea = permuting_spikes_corner_shuffle(neuron,behav,thresh,env_coor{ii})
    end
    cscore_shuffle{ii} = cscore_shuffle_ea;
end
%method 1: compare cornerscore with cscore_thresh to determine corner cells
cscore_thresh = cellfun(@(x) quantile(x,0.95,2), cscore_shuffle, 'uni', 0);
cornerscore = cellfun(@(x) x.cornerscore, cmetrics, 'uni', 0);
cornercell = cellfun(@(x,y) find(x > y), cornerscore, cscore_thresh, 'uni', 0);
%method 2: use a fixed threshold to determine corner cells
% cornercell = cellfun(@(x) find(x.cscore > 0.5), cmetrics, 'uni', 0);

%Adding another constraint: pair-wise distances between any fields should
%be greater than the corner-center distance/2
%to compute pair-wise field distances and distance threshold
d2 = cell(1,length(cmetrics));
dist_thresh = cell(1,length(cmetrics));
d_env = cellfun(@pdist, env_coor, 'UniformOutput', false);
d_env_square = cellfun(@squareform, d_env, 'UniformOutput', false);
for n = 1:length(cmetrics)
    if isempty(cmetrics{n}.field_coor)
        d2{n} = [];
        dist_thresh{n} = [];
    else
        d2{n} = cellfun(@(x) pdist(x), cmetrics{n}.field_coor, 'uni', 0);
        dist_thresh{n} = mean(d_env_square{n}(end,1:end-1))/2;
    end
end
%find cells that their paired-wise field distances are all greater than the
%certain threshold
cell_largedist = cell(1, length(cmetrics));
for n = 1:length(d2)
    d2_ea = d2{n};
    if isempty(d2_ea)
        cell_largedist{n} = [];
    else
        d2_thresh = dist_thresh{n};
        idx = cellfun(@(x) all(x > d2_thresh), d2_ea);
        cell_largedist{n} = find(idx == 1);
    end
end
cornercellf = cellfun(@(x,y) intersect(x,y), cornercell, cell_largedist, 'uni',0);

%% output vars into a structure
C = struct;
C.cmetrics = cmetrics;
C.cscore = cornerscore;
C.cscore_shuffle = cscore_shuffle;
C.cscore_thresh = cscore_thresh;
C.cornercell_raw = cornercell;
C.cell_largedist = cell_largedist;
C.cornercell = cornercellf;

end

