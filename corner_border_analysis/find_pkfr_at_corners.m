function pkfr_corners = find_pkfr_at_corners(ratemap, convsize, datapath, condition)
%This function aims to find the peak fr at each phyical corner for each
%cell. 
%yanjuns@stanford.edu, 12/14/22
if ~exist('condition','var')||isempty(condition)
    condition = 'normal';
end
if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)
if ~exist('convsize','var')||isempty(convsize)
    convsize = 12;
end
if ~exist('ratemap','var')||isempty(ratemap)
    load('firingrateAll.mat', 'ratemap')
end

if strcmp(condition, 'shuttle')
    load('NBindivLR.mat','env_coori')
elseif strcmp(condition, 'obj')
    load('env_geometry.mat','obj_coori','S')
    env_coori = obj_coori;
else
    load('env_geometry.mat','env_coori','S')
end
envc = cellfun(@(x) x(1:end-1,:), env_coori, 'uni', 0);
pkfr_corners = cell(1,length(envc)); %pkfr at each corner
for n = 1:length(envc)
    envc_ea = envc{1,n};
    pkfr_ss = zeros(size(envc_ea,1),length(ratemap{1,n}));
    for jj = 1:length(ratemap{1,n})
    ratemap_ea = ratemap{1,n}{1,jj};
    v = ones(convsize);
    pkfr_ea = [];
    for ii = 1:size(envc_ea,1)
        M = zeros(size(ratemap_ea,1),size(ratemap_ea,2));
        M(envc_ea(ii,2),envc_ea(ii,1)) = 1;
        corner_regidx = logical(conv2(M,v,'same'));
        corner_region = ratemap_ea(corner_regidx);
        pkfr_ea = [pkfr_ea;max(corner_region)];
%         %to plot each extracted region
%         row = sum(any(corner_regidx > 0, 2));
%         col = sum(any(corner_regidx > 0, 1));
%         corner_region = reshape(corner_region, [row,col]);
%         figure
%         imagesc(corner_region)
%         axis image
    end
    pkfr_ss(:,jj) = pkfr_ea;
    end
    pkfr_corners{n} = pkfr_ss;
end

end

