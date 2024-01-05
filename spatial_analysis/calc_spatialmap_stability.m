function Rxy = calc_spatialmap_stability(neuron,behav,thresh,binsize)
%calculation within session rate map stability (between 1st and 2nd half)
if ~exist('binsize','var') || isempty(binsize)
    binsize = 2;
end
[neuron1,behav1,neuron2,behav2] = split_data_interleaved(neuron,behav,30);
% [neuron1,behav1,neuron2,behav2] = split_data_in_half(neuron,behav);
%get a map axis
xAxis = 0:binsize:ceil(max(behav.position(:,1)+binsize)); %binedges of x position
yAxis = 0:binsize:ceil(max(behav.position(:,2)+binsize)); %binedges of y position
%get ratemap for 1st and 2nd half of the data
[ratemap1,~,~] = calculate_subset_ratemap(neuron1,behav1,thresh,xAxis,yAxis);
[ratemap2,~,~] = calculate_subset_ratemap(neuron2,behav2,thresh,xAxis,yAxis);
%smooth ratemap
ratemapS1 = cellfun(@(x) filter2DMatrices(x,1), ratemap1, 'uni',0);
ratemapS2 = cellfun(@(x) filter2DMatrices(x,1), ratemap2, 'uni',0);
%get rid of NaNs in the rate map
% idx1 = cellfun(@isnan,ratemapS1,'uni',0);
% for ii = 1:length(ratemapS1)
%     ratemapS1{ii}(idx1{ii}) = 0;
% end
% idx2 = cellfun(@isnan,ratemapS2,'uni',0);
% for ii = 1:length(ratemapS2)
%     ratemapS2{ii}(idx2{ii}) = 0;
% end
%compute the correlation
Rxy = cellfun(@(x,y) corr(x(:),y(:),'Rows','complete'), ratemapS1, ratemapS2);
Rxy = Rxy';
end
