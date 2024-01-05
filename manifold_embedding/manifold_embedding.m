function [reduction_all, umap_template_all, clusterIds_all] = manifold_embedding(neuronIndividualsf,behavIndividualsf,thresh,env_coor,S,ss2plot)
% This function performs manifold embedding for neural data. It runs PCA
% for neural data first, and takes the first 10 dimensions of PCA to run a
% 3D embedding using UMAP. 
% This function needs UMAP package for matlab, which can be downloaded with
% the following link
% https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap

% load data
reduction_all = {};
umap_template_all = {};
clusterIds_all = {};

ax = figure;
set(ax, 'Position', [0, 200, 1200, 800]);
hold on;
for n = 1:length(ss2plot)
% sessions to analyze
ss = ss2plot(n);
ssname = S{ss};
corners = env_coor{ss};
numCorner = size(corners,1)-1;
%identify each corner
binsize = 2;
if numCorner <= 4
    areasize = 5;
else
    areasize = 4;
end
if isempty(neuronIndividualsf{ss}.pos)
    neuronIndividualsf{ss}.pos = interp1(behavIndividualsf{ss}.time,behavIndividualsf{ss}.position,neuronIndividualsf{ss}.time);
end
position = neuronIndividualsf{ss}.pos/binsize;
if any(isnan(position(:,1)))
    idx = find(isnan(position(:,1)));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        position(idx(1:ind),:) = repmat(position(idx(ind)+1,:),size(idx(1:ind),1),1);
        position(idx(ind+1:end),:) = repmat(position(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            position(idx,:) = repmat(position(idx(end)+1,:),size(idx,1),1);
        else
            position(idx,:) = repmat(position(idx(1)-1,:),size(idx,1),1);
        end
    end
end

%% coloring the position by distance to corners.
dist2corner = pdist2(position,corners(1:end-1,:));
[D,I] = min(dist2corner,[],2);
Dcorner.dist = cell(1,numCorner);
Dcorner.idx = cell(1,numCorner);
Dcorner.colors = cell(1,numCorner);
for ii = 1:numCorner
    idx1 = find(I == ii);
    idx2 = find(D < areasize);
    idx = intersect(idx1,idx2);
    dist = D(idx);
    Dcorner.idx{ii} = idx;
    Dcorner.dist{ii} = dist;
end
colors_ori = {};
colors_ori{1} = [188 190 192]/255;%grey
colors_ori{2} = [0 0.4470 0.7410]; %blue
colors_ori{3} = [0.8500 0.3250 0.0980]; %orange
colors_ori{4} = [0.9290 0.6940 0.1250]; %yellow
colors_ori{5} = [0.4940 0.1840 0.5560]; %purple
colors_ori{6} = [0.4660 0.6740 0.1880]; %green
colors_ori{7} = [0.3010 0.7450 0.9330]; %light blue

k = 1;
for ii = 2:length(Dcorner.dist)+1
    clength = length(Dcorner.dist{k});
    colors_grd = [linspace(colors_ori{ii}(1),colors_ori{1}(1),clength)',...
        linspace(colors_ori{ii}(2),colors_ori{1}(2),clength)', ...
        linspace(colors_ori{ii}(3),colors_ori{1}(3),clength)'];
    [~, ind] = sort(Dcorner.dist{k});
    Dcorner.colors{k} = colors_grd(ind, :);
    k = k+1;
end

% plot regions
subplot(round(length(ss2plot)/2),4,2*n-1)
plot(position(:,1),position(:,2),'Color',[188 190 192]/255)
axis image
hold on
for ii = 1:numCorner
    plot(position(Dcorner.idx{ii},1),position(Dcorner.idx{ii},2));
end
title(ssname)
%% perform PCA first to reduce noise in the neural data
data_raw = neuronIndividualsf{ss}.S;
data = data_raw > repmat(thresh,1,size(data_raw,2));
data = double(data');
for ii = 1:size(data,2)
    data_ea = data(:,ii);
    data_filt = imgaussfilt(data_ea,5);
    data_z = zscore(data_filt);
    data(:,ii) = data_z;
end
[pcaIndiv.coeff,pcaIndiv.score,pcaIndiv.latent,...
    pcaIndiv.tsquared,pcaIndiv.explained,pcaIndiv.mu] = pca(data);

%% Perform UMAP and plot
nDim = 10;
[reduction, umap_template, clusterIds]=run_umap(pcaIndiv.score(:,1:nDim), 'min_dist', 0.1, 'n_neighbors', 100, 'n_components', 3, 'verbose', 'none');

reduction_bsl = reduction;
for ii = 1:length(Dcorner.idx)
    reduction_bsl(Dcorner.idx{ii},:) = NaN;
end

subplot(round(length(ss2plot)/2),4,2*n)
hold on
scatter3(reduction_bsl(:,1),reduction_bsl(:,2),reduction_bsl(:,3), 100, [188 190 192]/255, '.')
axis square
xlabel('UMAP1')
ylabel('UMAP2')
zlabel('UMAP3')
scatter3(reduction(Dcorner.idx{1},1),reduction(Dcorner.idx{1},2),reduction(Dcorner.idx{1},3),200, Dcorner.colors{1}, '.')
scatter3(reduction(Dcorner.idx{2},1),reduction(Dcorner.idx{2},2),reduction(Dcorner.idx{2},3),200, Dcorner.colors{2}, '.')
scatter3(reduction(Dcorner.idx{3},1),reduction(Dcorner.idx{3},2),reduction(Dcorner.idx{3},3),200, Dcorner.colors{3}, '.')
if numCorner == 4
    scatter3(reduction(Dcorner.idx{4},1),reduction(Dcorner.idx{4},2),reduction(Dcorner.idx{4},3),200, Dcorner.colors{4}, '.')
end
if numCorner == 6
    scatter3(reduction(Dcorner.idx{4},1),reduction(Dcorner.idx{4},2),reduction(Dcorner.idx{4},3),200, Dcorner.colors{4}, '.')
    scatter3(reduction(Dcorner.idx{5},1),reduction(Dcorner.idx{5},2),reduction(Dcorner.idx{5},3),200, Dcorner.colors{5}, '.')
    scatter3(reduction(Dcorner.idx{6},1),reduction(Dcorner.idx{6},2),reduction(Dcorner.idx{6},3),200, Dcorner.colors{6}, '.')
end
title(ssname)

% %plot individually
% figure
% hold on
% scatter3(reduction(:,1),reduction(:,2),reduction(:,3), 100, '.','MarkerEdgeColor',[188 190 192]/255)
% axis square
% xlabel('UMAP1')
% ylabel('UMAP2')
% zlabel('UMAP3')
% scatter3(reduction(Dcorner.idx{1},1),reduction(Dcorner.idx{1},2),reduction(Dcorner.idx{1},3),200, '.','MarkerEdgeColor','#0072BD')
% scatter3(reduction(Dcorner.idx{2},1),reduction(Dcorner.idx{2},2),reduction(Dcorner.idx{2},3),200, '.','MarkerEdgeColor','#D95319')
% scatter3(reduction(Dcorner.idx{3},1),reduction(Dcorner.idx{3},2),reduction(Dcorner.idx{3},3),200, '.','MarkerEdgeColor','#EDB120')
% if numCorner == 4
%     scatter3(reduction(Dcorner.idx{4},1),reduction(Dcorner.idx{4},2),reduction(Dcorner.idx{4},3),200, '.','MarkerEdgeColor','#7E2F8E')
% end
% if numCorner == 6
%     scatter3(reduction(Dcorner.idx{4},1),reduction(Dcorner.idx{4},2),reduction(Dcorner.idx{4},3),200, '.','MarkerEdgeColor','#77AC30')
%     scatter3(reduction(Dcorner.idx{5},1),reduction(Dcorner.idx{5},2),reduction(Dcorner.idx{5},3),200, '.','MarkerEdgeColor','#4DBEEE')
%     scatter3(reduction(Dcorner.idx{6},1),reduction(Dcorner.idx{6},2),reduction(Dcorner.idx{6},3),200, '.','MarkerEdgeColor','#A2142F')
% end

%% output data
reduction_all{n} = reduction;
umap_template_all{n} = umap_template;
clusterIds_all{n} = clusterIds;
end

end