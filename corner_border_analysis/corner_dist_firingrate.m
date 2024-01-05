function [dist2corner,dist2cornerb] = corner_dist_firingrate(datapath,binsize)
% this function aims to find the firing rate of each bin in the
% ratemap and the distance of that bin to the nearest corner. This prepares
% the quantification of linear distance to corners vs. the firing rate for
% each neuron.
% dist2cornerb is the binned version of dist2corner
% Yanjun Sun, yanjuns@stanford.edu, 11/22/2022

if ~exist('binsize','var')||isempty(binsize)
    binsize = 1;
end
if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end

cd(datapath);
load('firingrateAll.mat','ratemap')
load('env_geometry.mat','env_coor')
% %expC concave
% load('env_geometry.mat','ccav_coor')
% env_coor = ccav_coor;

dist2corner = struct;
dist2corner.fr = cell(1,length(ratemap));
dist2corner.d = cell(1,length(ratemap));
for n = 1:length(ratemap)
    ratemap_ss = ratemap{n};
    env_coor_ea = env_coor{n}(1:end-1,:); %last row is centroid
    %flip the xy coordinates to match the position on the matrix
    env_coor_ea = [env_coor_ea(:,2),env_coor_ea(:,1)]; 
    d_ea = [];
    for ii = 1:size(ratemap_ss{1},1)
        for jj = 1:size(ratemap_ss{1},2)
            if ~isnan(ratemap_ss{1}(ii,jj))
                d_ea = [d_ea;pdist2([ii,jj],env_coor_ea)];
            end
        end
    end
    dist2corner_fr = cell(1,length(ratemap_ss));
    for k = 1:length(ratemap_ss)
        ratemap_ea = ratemap_ss{k};
        fr_ea = [];
        for ii = 1:size(ratemap_ea,1)
            for jj = 1:size(ratemap_ea,2)
                if ~isnan(ratemap_ea(ii,jj))
                    fr_ea = [fr_ea;ratemap_ea(ii,jj)];
                end
            end
        end
        dist2corner_fr{k} = fr_ea;
    end
    dist2corner.d{n} = d_ea;
    dist2corner.fr{n} = dist2corner_fr;
end

% bin the x data (distance) into 1cm bins and average the y (fr) data
% within each bin. Easy for making plots.
dist2cornerb = struct;
dist2cornerb.fr = cell(1,length(ratemap));
dist2cornerb.d = cell(1,length(ratemap));
for ii = 1:length(dist2corner.fr)
    dist2corner_frea = dist2corner.fr{ii};
    dist2corner_dea = dist2corner.d{ii};
    frb = [];
    x = min(dist2corner_dea,[],2);
    for jj = 1:length(dist2corner_frea)
        y = dist2corner_frea{1,jj};
        [xbin,ybin] = bin_x_y(x,y,binsize);
        frb = [frb;ybin];
    end
    db = [1:length(xbin)]*binsize;
    dist2cornerb.d{ii} = db;
    dist2cornerb.fr{ii} = frb;
end

end

