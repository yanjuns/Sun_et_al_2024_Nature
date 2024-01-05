function R = calc_spatialmap_stabilityx(map1,map2)
%This function computes 2d correlations for spatial maps from different
%sessions. map1 and map2 are two cell arrays that each contains all the
%ratemaps from the session. 

%remove any rows and cols that are all NaNs
map1 = cellfun(@(x) x(~all(isnan(x),2),~all(isnan(x),1)), map1, 'uni',0);
map2 = cellfun(@(x) x(~all(isnan(x),2),~all(isnan(x),1)), map2, 'uni',0);
map1 = cellfun(@convertNaN2zero, map1, 'uni',0);
map2 = cellfun(@convertNaN2zero, map2, 'uni',0);

%find the max number of row and col
maxrow = max(size(map1{1},1),size(map1{2},1));
maxcol = max(size(map1{1},2),size(map1{2},2));

%scale the maps to the max size
mapA = cellfun(@(x) imresize(x,[maxrow,maxcol],'method','nearest'), map1, 'uni', 0);
mapB = cellfun(@(x) imresize(x,[maxrow,maxcol],'method','nearest'), map2, 'uni', 0);

%compute the correlation
R = cellfun(@(x,y) corr2(x,y), mapA, mapB)';

end

