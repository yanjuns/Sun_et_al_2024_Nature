%% Experiment K
% Experiment C aims to identify neurons that might encode convex corners
load('F:\analysis_folders.mat','expK')
datapath = expK;

%% Identify corner cells
%determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('env_geometry.mat','S')
    mask = cell(1,length(S));
    for jj = 1:length(S)
        if strcmp(S{jj}, 'largeXbsl')
            mask{jj} = true;
        else
            mask{jj} = false;
        end
    end
    %determine corner cell for each session by spike shuffling.
    %NOTE, for expK, C was using a 0.3 threshold for identifying
    %corner cells.
    CR1 = identify_corner_cell(mask);
    save('corner_metrics.mat','CR1','-append')
end
%% Incorporate within sessions stability for corner cells
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat')
    load('behavIndividualsf.mat')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndividualsf, behavIndividualsf, 'uni',0);
    save('spatial_metrics.mat', 'map_stb','-append')
end
%% Revised definition of corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'CR1')
    C = CR1;
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C')    
end
%% Quantify the peak firing rate at corners for all sessions
for n = 1:length(datapath)
    cd(datapath{n})
    %find the pkfr at corners for corner cells
    pkfr_corners = find_pkfr_at_corners([],8);
    %simulate a ratemap with constant firing
    [neuronSim,ratemapSim] = simulate_ratemap;
    %find the pkfr at corners of the simulated neuron
    pkfrsim_corners = find_pkfr_at_corners(ratemapSim,8);
    save('corner_metrics.mat','pkfr_corners','pkfrsim_corners','-append')
    save('neuronSim.mat','neuronSim','ratemapSim','-v7.3')
end

%% OVAL
%% Oval experiment
%% whether ccav corner cells respond to curvature in the oval ***
session_name = {'oval','ovald90'}; %expC
bsl = {'squarebsl'};
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
    load('corner_metricsR1.mat','C')
    load('border_metrics.mat','bcell')
    idx = contains(S, session_name);
    pkfr_corners = pkfr_corners(idx);
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell(idxbsl);
    cornercell = cc{end};
    allcell = 1:length(pkfr_corners{1});
    ind = ~ismember(allcell, cornercell);
    ncc = allcell(ind);
    bc = bcell(idxbsl);
    bvc = bc{end};
%     idxbvc = ~ismember(bc, cornercell);
%     bvc = bc(idxbvc);
%     idxcc = ~ismember(cornercell,bc);
%     cornercell = cornercell(idxcc);
    
    %corner cell
    pkfrcc_atcorner = cellfun(@(x) x(:,cornercell), pkfr_corners, 'uni', false);
    %corrected using simulated cell
    pkfrsim_atcorner = pkfrsim_corners(idx);
    pkfrccOval = cellfun(@(x,y) x./y, pkfrcc_atcorner, pkfrsim_atcorner, 'uni', false);
    %non-corner cell
    pkfrncc_atcorner = cellfun(@(x) x(:,ncc), pkfr_corners, 'uni', false);
    pkfrnccOval = cellfun(@(x,y) x./y, pkfrncc_atcorner, pkfrsim_atcorner, 'uni', false);
    %boundary vector cell (2b continued)
    pkfrbvc_atcorner = cellfun(@(x) x(:,bvc), pkfr_corners, 'uni', false);
    pkfrbvcOval = cellfun(@(x,y) x./y, pkfrbvc_atcorner, pkfrsim_atcorner, 'uni', false);
    save('corner_metricsR1.mat','pkfrccOval','pkfrnccOval','pkfrbvcOval','-append')
end

pkfrOval.cc = [];
pkfrOval.ncc = [];
pkfrOval.bvc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','pkfrccOval','pkfrnccOval','pkfrbvcOval')
    pkfrcc = cellfun(@(x) mean(x,2), pkfrccOval, 'uni', 0);
    pkfrOval.cc = [pkfrOval.cc, mean(cell2mat(pkfrcc),2)];
    
    pkfrncc = cellfun(@(x) mean(x,2), pkfrnccOval, 'uni', 0);
    pkfrOval.ncc = [pkfrOval.ncc, mean(cell2mat(pkfrncc),2)];
    
    pkfrbvc = cellfun(@(x) mean(x,2), pkfrbvcOval, 'uni', 0);
    pkfrOval.bvc = [pkfrOval.bvc, mean(cell2mat(pkfrbvc),2)];
end
% Note: the first two corner locations in the oval are high curvatures. 
pkfrOval.cchicurve = mean(pkfrOval.cc(1:2,:))';
pkfrOval.cclocurve = mean(pkfrOval.cc(3:4,:))';
pkfrOval.ncchicurve = mean(pkfrOval.ncc(1:2,:))';
pkfrOval.ncclocurve = mean(pkfrOval.ncc(3:4,:))';
pkfrOval.bvchicurve = mean(pkfrOval.bvc(1:2,:))';
pkfrOval.bvclocurve = mean(pkfrOval.bvc(3:4,:))';
save('F:\Results_experimentK\pkfrOval_R1.mat','pkfrOval')

[h,p] = signrank(pkfrOval.cchicurve,pkfrOval.cclocurve)
[h,p] = signrank(pkfrOval.ncchicurve,pkfrOval.ncclocurve)
[h,p] = signrank(pkfrOval.bvchicurve,pkfrOval.bvclocurve)
data = [pkfrOval.cchicurve - pkfrOval.cclocurve, pkfrOval.ncchicurve - pkfrOval.ncclocurve];
figure
plot(data')
%use GraphPad to plot and perform statistics

%% re-run oval data using raw fr data
session_name = {'oval','ovald90'}; %expC
bsl = {'squarebsl'};
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
    load('corner_metricsR1.mat','C')
    idx = contains(S, session_name);
    pkfr_corners = pkfr_corners(idx);
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell(idxbsl);
    cornercell = cc{end};
    allcell = 1:length(pkfr_corners{1});
    ind = ~ismember(allcell, cornercell);
    ncc = allcell(ind);
    %corner cell
    pkfrccOval_raw = cellfun(@(x) x(:,cornercell), pkfr_corners, 'uni', false);
    %non-corner cell
    pkfrnccOval_raw  = cellfun(@(x) x(:,ncc), pkfr_corners, 'uni', false);
    save('corner_metricsR1.mat','pkfrccOval_raw','pkfrnccOval_raw','-append')
end

pkfrOval.cc = [];
pkfrOval.ncc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','pkfrccOval_raw','pkfrnccOval_raw')
    pkfrcc = cellfun(@(x) mean(x,2), pkfrccOval_raw, 'uni', 0);
    pkfrOval.cc = [pkfrOval.cc, mean(cell2mat(pkfrcc),2)];
    
    pkfrncc = cellfun(@(x) mean(x,2), pkfrnccOval_raw, 'uni', 0);
    pkfrOval.ncc = [pkfrOval.ncc, mean(cell2mat(pkfrncc),2)];
end
% Note: the first two corner locations in the oval are high curvatures. 
pkfrOval.cchicurve = mean(pkfrOval.cc(1:2,:))';
pkfrOval.cclocurve = mean(pkfrOval.cc(3:4,:))';
pkfrOval.ncchicurve = mean(pkfrOval.ncc(1:2,:))';
pkfrOval.ncclocurve = mean(pkfrOval.ncc(3:4,:))';
pkfrOval_raw = pkfrOval;
save('F:\Results_experimentK\pkfrOval_R1.mat','pkfrOval_raw','-append')

data = [pkfrOval_raw.cchicurve - pkfrOval_raw.cclocurve, pkfrOval_raw.ncchicurve - pkfrOval_raw.ncclocurve];
figure
plot(data')
%use GraphPad to plot and perform statistics

%% OBJECTS
%% LargeRTcvx
%% Convex objects in large rectangle
%% mannually identify the location of objects in the largeRTcvx session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [~,Nobj,obj_coor,obj_coori] = identify_obj_geometry;
    save('env_geometry.mat','Nobj','obj_coor','obj_coori','-append')
end
% Note: there are 8 sampling points for the large obj and 4 sampling points
% for the small obj. 
%% Quantify the peak firing rate around objects for all sessions
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
%% Revised definition of corner cells ***
%Note: the orignal C was correct by using normal mask size, CR1 used
%slightly larger mask size than C. 
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'C')
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    Cvx = C;
    save('corner_metricsR1.mat', 'Cvx', '-append')    
end
%% whether cvex corner cells respond to object curvature ***
% Note: the first 8 points are for the large obj. 
session_name = {'largeRTcvx'}; %expC
bsl = {'largeXbsl'};
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','pkfr_obj','pkfrsim_obj')
    load('corner_metricsR1.mat','Cvx')
    C = Cvx;
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
for n = 1:length(datapath)
    cd(datapath{n})
    % to exclude mice that have less than 3 convex corner cells
    bsl = {'largeXbsl'};
    load('env_geometry.mat','S')
    load('corner_metrics.mat','C')
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell{idxbsl};
    if numel(cc) >=3
        load('corner_metricsR1.mat','pkfrccObj','pkfrnccObj')
        pkfrcc = cellfun(@(x) mean(x,2), pkfrccObj, 'uni', 0);
        pkfrObj.cc = [pkfrObj.cc, mean(cell2mat(pkfrcc),2)];
        pkfrncc = cellfun(@(x) mean(x,2), pkfrnccObj, 'uni', 0);
        pkfrObj.ncc = [pkfrObj.ncc, mean(cell2mat(pkfrncc),2)];
    end
end
% Note: the first two corner locations in the oval are high curvatures. 
pkfrObj.cclocurve = mean(pkfrObj.cc(1:8,:))';
pkfrObj.cchicurve = mean(pkfrObj.cc(9:12,:))';
pkfrObj.ncclocurve = mean(pkfrObj.ncc(1:8,:))';
pkfrObj.ncchicurve = mean(pkfrObj.ncc(9:12,:))';
save('F:\Results_experimentK\pkfrObj_R1.mat','pkfrObj')

[h,p] = signrank(pkfrObj.cchicurve,pkfrObj.cclocurve)
[h,p] = signrank(pkfrObj.ncchicurve,pkfrObj.ncclocurve)
[h,p] = signrank(pkfrObj.bvchicurve,pkfrObj.bvclocurve)
data = [pkfrObj.cchicurve - pkfrObj.cclocurve, pkfrObj.ncchicurve - pkfrObj.ncclocurve];
figure
plot(data')

%% re-run analysis using raw fr data
% Note: the first 8 points are for the large obj. 
session_name = {'largeRTcvx'}; %expK
bsl = {'largeXbsl'};
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','pkfr_obj','pkfrsim_obj')
    load('corner_metricsR1.mat','Cvx')
    C = Cvx;
    idx = contains(S, session_name);
    pkfr_obj = pkfr_obj(idx);
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell(idxbsl);
    cornercell = cc{end};
    allcell = 1:length(pkfr_obj{1});
    ind = ~ismember(allcell, cornercell);
    ncc = allcell(ind);
    %corner cell
    pkfrccObj_raw = cellfun(@(x) x(:,cornercell), pkfr_obj, 'uni', false);
    %non-corner cell
    pkfrnccObj_raw = cellfun(@(x) x(:,ncc), pkfr_obj, 'uni', false);
    save('corner_metricsR1.mat','pkfrccObj_raw','pkfrnccObj_raw','-append')
end

%gather data from all mice
pkfrObj.cc = [];
pkfrObj.ncc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    % to exclude mice that have less than 3 convex corner cells
    bsl = {'largeXbsl'};
    load('env_geometry.mat','S')
    load('corner_metrics.mat','C')
    idxbsl = strcmp(S, bsl);
    cc = C.cornercell{idxbsl};
    if numel(cc) >=3
        load('corner_metricsR1.mat','pkfrccObj_raw','pkfrnccObj_raw')
        pkfrcc = cellfun(@(x) mean(x,2), pkfrccObj_raw, 'uni', 0);
        pkfrObj.cc = [pkfrObj.cc, mean(cell2mat(pkfrcc),2)];
        pkfrncc = cellfun(@(x) mean(x,2), pkfrnccObj_raw, 'uni', 0);
        pkfrObj.ncc = [pkfrObj.ncc, mean(cell2mat(pkfrncc),2)];
    end
end
% Note: the first two corner locations in the oval are high curvatures. 
pkfrObj.cclocurve = mean(pkfrObj.cc(1:8,:))';
pkfrObj.cchicurve = mean(pkfrObj.cc(9:12,:))';
pkfrObj.ncclocurve = mean(pkfrObj.ncc(1:8,:))';
pkfrObj.ncchicurve = mean(pkfrObj.ncc(9:12,:))';
pkfrObj_raw = pkfrObj; 
save('F:\Results_experimentK\pkfrObj_R1.mat','pkfrObj_raw','-append')

data = [pkfrObj_raw.cchicurve - pkfrObj_raw.cclocurve, pkfrObj_raw.ncchicurve - pkfrObj_raw.ncclocurve];
figure
plot(data')
