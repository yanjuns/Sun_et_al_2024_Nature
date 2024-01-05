function C = identify_corner_cell(mask, datapath, threshlevel, env_coor)
%This function aims to identify corner cells for each session in each mouse
%yanjuns@stanford.edu, 10/12/22
%C is the final structure that contains all the information
if ~exist('threshlevel','var')||isempty(threshlevel)
    threshlevel = 0.3;
end
if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)
%load data
if ~exist('env_coor','var')||isempty(env_coor)
    load ('env_geometry.mat','env_coor');
end
load ('neuronIndividualsf.mat','neuronIndividualsf');
load ('behavIndividualsf.mat','behavIndividualsf');
load ('thresh.mat','thresh');
load ('firingrateAll','ratemap')
%if only run specific sessions
% ss = 6;
% neuronIndividualsf = neuronIndividualsf(ss);
% behavIndividualsf = behavIndividualsf(ss);
% ratemap = ratemap(ss);
% env_coor = env_coor(ss);

if ~exist('mask','var')||isempty(mask)
    mask = repmat({false},1,length(ratemap));
end

%compute corner metrics including corner scores
cmetrics = cell(1,length(ratemap));
for n = 1:length(ratemap)
    cmetrics{n} = compute_corner_score(ratemap{n},env_coor{n},mask{n},true,threshlevel);
end

cscore_shuffle = cell(1, length(neuronIndividualsf));
parfor ii = 1:length(neuronIndividualsf)
    neuron = neuronIndividualsf{1,ii};
    behav = behavIndividualsf{1,ii};
    if isempty(env_coor{ii})
        cscore_shuffle_ea = [];
    else
        cscore_shuffle_ea = permuting_spikes_corner_shuffle(neuron,behav,thresh,env_coor{ii},mask{ii},threshlevel)
    end
    cscore_shuffle{ii} = cscore_shuffle_ea;
end
%method 1: compare cornerscore with cscore_thresh to determine corner cells
cscore_thresh = cellfun(@(x) quantile(x,0.95,2), cscore_shuffle, 'uni', 0);
cornerscore = cellfun(@(x) x.cornerscore, cmetrics, 'uni', 0);
cornercell = cellfun(@(x,y) find(x > y), cornerscore, cscore_thresh, 'uni', 0);
%method 2: use a fixed threshold to determine corner cells
% cornercell = cellfun(@(x) find(x.cscore > 0.5), cmetrics, 'uni', 0);

%% Adding field distance constraint
%pair-wise distances between any major fields should be greater than 
%the corner-center distance/2. 
%Note: major fields are fields with largest n corner scores, which n is the
%number of corners.
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
        % to only consider major fields, not extra fields
        ncorners = size(d_env_square{n},1)-1;
        [~,idx] = cellfun(@(x) maxk(x,ncorners), cmetrics{n}.cornerscore_ea, 'uni', 0);
        fieldLocation = cellfun(@(x,y) x(y,:), cmetrics{n}.field_coor, idx, 'uni', 0);
        d2{n} = cellfun(@(x) pdist(x), fieldLocation, 'uni', 0);
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

%% plot corner score distribution for illustration purpose
% msize = 25;
% x  = zeros(msize);
% corners = [1,1;1,msize;msize,msize;msize,1];
% ctr = [round(msize/2),round(msize/2)];
% cscore = NaN(msize);
% 
% for ii = 1:size(x,1)
%     for jj = 1:size(x,2)
%         d = pdist2([ii,jj], corners);
%         dmin = min(d);
%         dc = pdist2([ii,jj], ctr);
%         cscore(ii,jj) = (dc-dmin)/(dmin+dc);
%     end
% end
% 
% figure
% imagesc(cscore);
% axis image
% colormap summer
end

