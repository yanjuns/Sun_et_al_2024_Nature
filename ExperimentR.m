%% Experiment R
% Experiment C aims to identify neurons that might encode convex corners
load('F:\analysis_folders.mat','expR')
datapath = expR;
%% Identify corner cells
%determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('env_geometry.mat','S')
    mask = cell(1,length(S));
    for jj = 1:length(S)
        if strcmp(S{jj}, 'largeXm1d0') || strcmp(S{jj}, 'largeXm1d90')
            mask{jj} = true;
        else
            mask{jj} = false;
        end
    end
    %determine corner cell for each session by spike shuffling.
    %NOTE, for expR, C was using a 0.3/0.35 threshold for identifying
    %corner cells. C2 was using a 0.4 threshold for identifying corner
    %cells. C2 works slightly better becuase the data in the large
    %environment is bit noiser. 
    %in this experiment, C2 is equal to CR1;
    C2 = identify_corner_cell(mask,[],0.4);
    save('corner_metrics.mat','C2','-append')
end
%% Incorporate within sessions stability for corner cells
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat')
    load('behavIndividualsf.mat')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndividualsf, behavIndividualsf, 'uni',0);
    save('spatial_metrics.mat', 'map_stb')
end
%% Revised definition of corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'C2')
    C = C2;
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C', '-append')    
end

%% OBJECTS
%% largesqobj1
%% Convex objects in large square
%% mannually identify the location of objects in the largeRTcvx session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [~,Nobj,obj_coor,obj_coori] = identify_obj_geometry;
    save('env_geometry.mat','Nobj','obj_coor','obj_coori','-append')
end
% Note: there are 8 sampling points in total. The first 3 are the convex
% corners of the triangle objects, the next 3 are the linear wall of the
% triangle objects, and the last 2 are from the cylindar. 
%% Quantify the peak firing rate around objects ***
for n = 1:length(datapath)
    cd(datapath{n})
    %find the pkfr at corners for corner cells
    pkfr_obj = find_pkfr_at_corners([],3,[],'obj');
    %simulate a ratemap with constant firing
    load('neuronSim.mat','ratemapSim')
    %find the pkfr at corners of the simulated neuron
    pkfrsim_obj = find_pkfr_at_corners(ratemapSim,3,[],'obj');
    save('corner_metrics.mat','pkfr_obj','pkfrsim_obj','-append')
end
%% whether cvex corner cells respond to object curvature ***
load('F:\analysis_folders.mat','expR')
datapath = expR([1,3:7,9:10]); %***
session_name = {'largesqobj1'}; %expC
bsl = {'largeXm1d90'};
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','pkfr_obj','pkfrsim_obj')
    load('corner_metricsR1.mat','C')
    idx = contains(S, session_name);
    pkfr_obj = pkfr_obj(idx);
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell(idxbsl);
    cornercell = cc{end};
    allcell = 1:length(pkfr_obj{1});
    ind = ~ismember(allcell, cornercell);
    ncc = allcell(ind);
    %corner cell
    pkfrcc_atobj = cellfun(@(x) x(:,cornercell), pkfr_obj, 'uni', false);
    %corrected using simulated cell
    pkfrsim_atobj = pkfrsim_obj(idx);
    pkfrccObj = cellfun(@(x,y) x./y, pkfrcc_atobj, pkfrsim_atobj, 'uni', false);
    %non-corner cell
    pkfrncc_atobj = cellfun(@(x) x(:,ncc), pkfr_obj, 'uni', false);
    pkfrnccObj = cellfun(@(x,y) x./y, pkfrncc_atobj, pkfrsim_atobj, 'uni', false);
    save('corner_metricsR1.mat','pkfrccObj','pkfrnccObj','-append')
end

%gather data from all mice
pkfrObj.cc = [];
pkfrObj.ncc = [];
numcc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    % to exclude mice that have less than 3 convex corner cells
    bsl = {'largeXm1d90'};
    load('env_geometry.mat','S')
    load('corner_metricsR1.mat','C')
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell{idxbsl};
    numcc = [numcc;numel(cc)];
    load('corner_metricsR1.mat','pkfrccObj','pkfrnccObj')
    pkfrcc = cellfun(@(x) mean(x,2), pkfrccObj, 'uni', 0);
    pkfrObj.cc = [pkfrObj.cc, cell2mat(pkfrcc)];
    pkfrncc = cellfun(@(x) mean(x,2), pkfrnccObj, 'uni', 0);
    pkfrObj.ncc = [pkfrObj.ncc, cell2mat(pkfrncc)];
end
% Note: the first two corner locations in the oval are high curvatures. 
pkfrObj.cchicurve = nanmean(pkfrObj.cc(1:3,:))';
pkfrObj.cclocurve = nanmean(pkfrObj.cc(4:6,:))';
pkfrObj.ccobj = nanmean(pkfrObj.cc(7:8,:))';
pkfrObj.ncchicurve = nanmean(pkfrObj.ncc(1:3,:))';
pkfrObj.ncclocurve = nanmean(pkfrObj.ncc(4:6,:))';
pkfrObj.nccobj = nanmean(pkfrObj.ncc(7:8,:))';
pkfrObj.numcc = numcc;
save('F:\Results_experimentR\pkfrObj_R1.mat','pkfrObj')
data.objtri = [pkfrObj.cchicurve - pkfrObj.cclocurve,...
    pkfrObj.ncchicurve - pkfrObj.ncclocurve];
data.objcyl = [pkfrObj.ccobj - pkfrObj.cclocurve,...
    pkfrObj.nccobj - pkfrObj.ncclocurve];
save('F:\Results_experimentR\pkfrObj_R1.mat','data','-append')

figure
subplot(2,1,1)
plot(data.objtri')
subplot(2,1,2)
plot(data.objcyl')


