function wfd=weighted_firing_distance(fields,binsize)
map=fields{1};
for i=2:length(fields)
    map=map+fields{i};
end
map=map/nansum(map(:));
[ly,lx]=size(map);
my=(1:ly)'*ones(1,lx);
mx=ones(ly,1)*(1:lx);
distance_matrix=min(min(my,mx),min(flipud(my),fliplr(mx)))*binsize;
wfd=nansum(nansum(map.*distance_matrix))/nansum(nansum(map.*ones(ly,lx)));
