%% Experiment F
% manipulations of visual appearance of the corner. 
load('F:\analysis_folders.mat','expF')
datapath = expF;

%% SHUTTLE BOX EXPERIMENT
% proportion of corner cells
%prepare for the data
%split the data into left and right compartments
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat');
    load('behavIndividualsf.mat');
    load('thresh.mat');
    load('sessions.mat', 'sessions')
    S = cellfun(@(x) strsplit(x,'_'), sessions, 'uni', 0);
    S = cellfun(@(x) x{end}, S, 'uni', 0);
    idx = strcmp(S, 'shuttle') | strcmp(S, 'shuttleb');
    neuronShuttle = neuronIndividualsf(idx);
    behavShuttle = behavIndividualsf(idx);
    neuronIndivLR = cell(1,length(neuronShuttle)*2);
    behavIndivLR = cell(1,length(neuronShuttle)*2);
    midline = NaN(1, length(neuronShuttle));
    for ii = 1:length(neuronShuttle)
        [neuronm, neuronR, behavm, behavR, autocenter] = split_neuron_behav_LR(neuronShuttle{ii}, behavShuttle{ii});
        jj = 2*ii-1;
        neuronIndivLR{1,jj} = neuronm; neuronIndivLR{1,jj+1} = neuronR;
        behavIndivLR{1,jj} = behavm; behavIndivLR{1,jj+1} = behavR;
        midline(1,ii) = autocenter;
    end
    close all
    %get the raw ratemap and smoothed ratemap
    meanFRindividuals = cell(1,length(neuronIndivLR));
    firingrateAll = cell(1,length(neuronIndivLR));
    countAll = cell(1,length(neuronIndivLR));
    countTime = cell(1,length(neuronIndivLR));
    smoccu = cell(1,length(neuronIndivLR)); %smoothed countTime
    ratemap = cell(1,length(neuronIndivLR));
    % Calculate spatial map and mean firing rate for each individual sessions
    for k = 1:length(neuronIndivLR)
        [firingrateAll{k},countAll{k},countTime{k},smoccu{k}] = ...
            calculate_firing_ratemap(neuronIndivLR{k},behavIndivLR{k},thresh,2,true);
        for m = 1:size(firingrateAll{1,k},2)
            meanFRindividuals{1,k}(m,1)= sum(sum(countAll{1,k}{1,m}))/sum(sum(countTime{1,k}));
        end
    end
    % Calculate smoothed firing rate maps
    for ii = 1:length(firingrateAll)
        fr = firingrateAll{ii};
        % to remove potential inf in the unsmoothed ratemap
        idx = cellfun(@isinf, fr,'UniformOutput',false);
        for jj = 1:length(fr)
            fr{jj}(idx{jj}) = 0;
        end
        % smoothed ratemap
        ratemap{ii} = cellfun(@(x) filter2DMatrices(x,1), fr, 'uni', 0);
    end
    save('NBindivLR.mat', 'neuronIndivLR', 'behavIndivLR', 'firingrateAll',...
        'countAll', 'countTime', 'smoccu', 'ratemap', 'meanFRindividuals', '-v7.3');
    clearvars -except datapath n
end
%% mannually identify the location of corners for each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [N,env_coor,env_coori] = identify_env_geometry_shuttle;
    save('NBindivLR.mat','N','env_coor','env_coori','-append')
end
%% determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    CR1 = identify_corner_cell_shuttle;
    save('NBindivLR.mat','CR1','-append')
end
%% within sessions stability for cells
for k = 1:length(datapath)
    cd(datapath{k});
    load('NBindivLR.mat','neuronIndivLR','behavIndivLR')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndivLR, behavIndivLR, 'uni',0);
    save('NBindivLR.mat', 'map_stb','-append')
end
%% Revised definition of corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat', 'C', 'map_stb')
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('NBindivLRR1.mat', 'C')    
end
%% for plot
propCC_shuttle = [];
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('NBindivLRR1.mat','C')
    load('thresh.mat','thresh')
    data_raw = cellfun(@numel, C.cornercell);
    data_raw = reshape(data_raw, [2,2])';
    %find sessions that contain less than 5 corner cells
    idx = any(data_raw < 3, 2); 
    data = cellfun(@numel, C.cornercell)./numel(thresh);
    data = reshape(data, [2,2])';
    %filter out sessions that contain less than 5 corner cells
    data(idx,:) = NaN;
    propCC_shuttle = [propCC_shuttle;nanmean(data)];
end
save('F:\Results_experimentF\prop_cornercell_shuttle_R1.mat','propCC_shuttle')

%% peak fr between the grey and black compartments
% xsession corner cells in shuttle box
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('NBindivLRR1.mat','C')
    cornercellx = cell(1,4);
    cornercellx{1} = intersect(C.cornercell{1},C.cornercell{2});
    cornercellx{2} = intersect(C.cornercell{1},C.cornercell{2});
    cornercellx{3} = intersect(C.cornercell{3},C.cornercell{4});
    cornercellx{4} = intersect(C.cornercell{3},C.cornercell{4});
    save('NBindivLRR1.mat','cornercellx','-append')
end

%find peak fr at each corner for each neuron
for n = 1:length(datapath)
    cd(datapath{n})
    %find the pkfr at corners for corner cells
    load('NBindivLR.mat','ratemap')
    pkfr_corners = find_pkfr_at_corners(ratemap,[],'shuttle');
    %simulate a ratemap with constant firing
    [neuronSim,ratemapSim] = simulate_ratemap([],'shuttle');
    %find the pkfr at corners of the simulated neuron
    pkfrsim_corners = find_pkfr_at_corners(ratemapSim,[],'shuttle');
    save('NBindivLR.mat','pkfr_corners','pkfrsim_corners','neuronSim','ratemapSim','-append')
end

% gather data for ploting
% plot corner cells that defined in both compartments
pkfrcc_atC = []; pkfrcc_atCs = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat','pkfr_corners','pkfrsim_corners')
    load('NBindivLRR1.mat','cornercellx')
    cornercell = cornercellx;
    %corner cell
    pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
    pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
    pkfrcc_atcorner = [pkfrcc_atcorner{:}];
    pkfrcc_atC = [pkfrcc_atC;mean(pkfrcc_atcorner)]; 
    %correction using simulated cell
    pkfrsim_atcorner = pkfrsim_corners;
    pkfrsim_atcorner = [pkfrsim_atcorner{:}];
    pkfrcc_atCs = [pkfrcc_atCs;mean(pkfrcc_atcorner./pkfrsim_atcorner)];
end
data_xctxcc = [nanmean(pkfrcc_atCs(:,[1,3]),2),nanmean(pkfrcc_atCs(:,[2,4]),2)];
save('F:\Results_experimentF\pkfr_shuttle_R1.mat','pkfrcc_atC','pkfrcc_atCs','data_xctxcc')

% plot corner cells that defined in grey compartments
pkfrcc_atC = []; pkfrcc_atCs = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat','pkfr_corners','pkfrsim_corners')
    load('NBindivLRR1.mat','C')
    cornercell = C.cornercell;
    cornercell{2} = cornercell{1};
    cornercell{4} = cornercell{3};
    %corner cell
    pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
    pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
    pkfrcc_atcorner = [pkfrcc_atcorner{:}];
    pkfrcc_atC = [pkfrcc_atC;mean(pkfrcc_atcorner)]; 
    %correction using simulated cell
    pkfrsim_atcorner = pkfrsim_corners;
    pkfrsim_atcorner = [pkfrsim_atcorner{:}];
    pkfrcc_atCs = [pkfrcc_atCs;mean(pkfrcc_atcorner./pkfrsim_atcorner)];
end
data_greycc = [nanmean(pkfrcc_atCs(:,[1,3]),2),nanmean(pkfrcc_atCs(:,[2,4]),2)];
save('F:\Results_experimentF\pkfr_shuttle_R1.mat','data_greycc','-append')

% plot corner cells that defined in grey compartments
pkfrcc_atC = []; pkfrcc_atCs = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat','pkfr_corners','pkfrsim_corners')
    load('NBindivLRR1.mat','C')
    cornercell = C.cornercell;
    cornercell{1} = cornercell{2};
    cornercell{3} = cornercell{4};
    %corner cell
    pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
    pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
    pkfrcc_atcorner = [pkfrcc_atcorner{:}];
    pkfrcc_atC = [pkfrcc_atC;mean(pkfrcc_atcorner)]; 
    %correction using simulated cell
    pkfrsim_atcorner = pkfrsim_corners;
    pkfrsim_atcorner = [pkfrsim_atcorner{:}];
    pkfrcc_atCs = [pkfrcc_atCs;mean(pkfrcc_atcorner./pkfrsim_atcorner)];
end
data_blackcc = [nanmean(pkfrcc_atCs(:,[1,3]),2),nanmean(pkfrcc_atCs(:,[2,4]),2)];
save('F:\Results_experimentF\pkfr_shuttle_R1.mat','data_blackcc','-append')

%% representative rate maps for shuttle box
% % plot representative ratemaps
% % good examples:
% % M4101Fa2, cell 15
% % M4113Fa, cell 303
% % M4115, cell 47
% load('NBindivLR.mat','ratemap')
% cell2plot = 15;
% figure
% subplot(2,2,1)
% pcolor(ratemap{1}{cell2plot})
% shading flat
% axis image
% subplot(2,2,2)
% pcolor(ratemap{2}{cell2plot})
% shading flat
% axis image
% subplot(2,2,3)
% pcolor(ratemap{3}{cell2plot})
% shading flat
% axis image
% subplot(2,2,4)
% pcolor(ratemap{4}{cell2plot})
% shading flat
% axis image
% colormap jet

%% TRIANGLE
%% MODIFIED TRIANGLE EXPERIMENT
% load('F:\analysis_folders.mat','expD')
% datapath = expD;
% 
% %Percentage of corner cells for the triangle and modified triangle session
% session_name = {'rightTri','rightTrim1'}; %expD
% prop_cornercell = cell(length(datapath),length(session_name));
% for n = 1:length(datapath)
%     cd(datapath{n})
%     load('corner_metrics.mat','C')
%     load('env_geometry.mat','S') %load all session name
%     load('thresh.mat','thresh')
%     numNeurons = length(thresh);
%     pcc = cell(1,length(session_name));
%     for ii = 1:length(session_name)
%         idx = strcmp(S, session_name{ii});
%         C_ea = C.cornercell(idx);
%         C_ea = C_ea(end);
%         C_ea = cellfun(@numel, C_ea);
%         pcc{ii} = C_ea(:)/numNeurons;
%     end
%     prop_cornercell(n,:) = pcc;
% end
% propCC_rightTri = cell2mat(prop_cornercell);
% save('F:\Results_experimentF\prop_cornercell_shuttle.mat','propCC_rightTri','-append')
% 
% % peak firing rate of corner cell in each corner
% session_name = {'rightTri','rightTrim1'}; %expD
% pkfrcc_atCnc = cell(1,length(session_name));
% pkfrcc_atCs = cell(1,length(session_name));
% pkfrcc_atC = cell(1,length(session_name));
% for k = 1:length(session_name)
%     pkfrcc_atC = []; pkfrcc_atCornc = []; pkfrcc_atCs = [];
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('env_geometry.mat','S')
%         load('corner_metrics.mat','pkfr_corners','C','cornercellx','pkfrsim_corners')
%         idx1 = strcmp(S, session_name{k});
%         temp = cellfun(@numel, C.cornercell);
%         idx2 = temp >= 5;
%         idx = idx1 & idx2;
%         idxs = 1:length(idx);
%         idx = idxs(idx);%select the last session
%         idx = idx(end);
%         pkfr_corners = pkfr_corners(idx);
%         cornercell = C.cornercell(idx);
%         allcell = 1:length(pkfr_corners{1,1});
%         ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
%         ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
%         %corner cell
%         pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
%         pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
%         pkfrcc_atcorner = [pkfrcc_atcorner{:}];
%         pkfrcc_atC = [pkfrcc_atC,pkfrcc_atcorner];
%         %correction using non-corner cell
%         pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
%         pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
%         pkfrncc_atcorner = [pkfrncc_atcorner{:}];
%         pkfrcc_atCornc = [pkfrcc_atCornc,pkfrcc_atcorner./pkfrncc_atcorner];
%         %correction using simulated cell
%         pkfrsim_atcorner = pkfrsim_corners(idx);
%         pkfrsim_atcorner = [pkfrsim_atcorner{:}];
%         pkfrcc_atCs = [pkfrcc_atCs,pkfrcc_atcorner./pkfrsim_atcorner];
%     end
%     pkfrcc_atCnc{k} = pkfrcc_atCornc';
%     pkfrcc_atCs{k} = pkfrcc_atCs';
%     pkfrcc_atC{k} = pkfrcc_atC';
% end
% save('F:\Results_experimentF\pkfr_rightTrim1.mat','pkfrcc_atC','pkfrcc_atCs','pkfrcc_atCnc')

%% Experiment G, DARK RECORDING AND WHISKER TRIMMING
%% Experiment G, dark recording and whisker trimming
%% Experiment G
% dark recording and whisker trimming
load('F:\analysis_folders.mat','expG')
datapath = expG;

%% determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    CR1 = identify_corner_cell;
    save('corner_metrics.mat','CR1','-append')
end
%% within sessions stability for cells
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

%% Percentage of corner cells for each session
prop_cornercell = [];
prop_placecell = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('spatial_metrics.mat','placecell')
    load('thresh.mat','thresh')
    if contains(datapath{n},'M4113F') || contains(datapath{n},'M4131F')
        idx = [1,2,3,4];
    else
        idx = [1,2,3,5];
    end
    cc_all = C.cornercell(idx);
    pc_all = placecell(idx);
    numNeurons = length(thresh);
    cc = cellfun(@numel, cc_all)./numNeurons;
    prop_cornercell = [prop_cornercell;cc];
    pc = cellfun(@numel, pc_all)./numNeurons;
    prop_placecell = [prop_placecell;pc];
end
%note for saved data: data were organized as baseline1, baseline2, dark,
%whisker. 
save('F:\Results_experimentF\prop_cell_dark_whisker_R1.mat','prop_cornercell','prop_placecell')
%% compare the firing rate of corner cells across baseline,dark and whisker sessions
%find peak fr at each corner for each neuron
for n = 1:length(datapath)
    cd(datapath{n})
    %find the pkfr at corners for corner cells
    pkfr_corners = find_pkfr_at_corners;
    %simulate a ratemap with constant firing
    [neuronSim,ratemapSim] = simulate_ratemap;
    %find the pkfr at corners of the simulated neuron
    pkfrsim_corners = find_pkfr_at_corners(ratemapSim);
    save('corner_metrics.mat','pkfr_corners','pkfrsim_corners','-append')
    save('neuronSim.mat','neuronSim','ratemapSim','-v7.3')
end
%gather data for plotting
session_name = {'square','squareNoLight','squareNoSense'}; %expA
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCors = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
        load('corner_metricsR1.mat','C')
        if contains(datapath{n},'M4113F') || contains(datapath{n},'M4131F')
            idx_all = [2,3,4];
        else
            idx_all = [2,3,5];
        end
        idx = idx_all(k);
        pkfr_corners = pkfr_corners(idx);
        cornercell = C.cornercell(idx);
        allcell = 1:length(pkfr_corners{1,1});
        ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
        ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        %corrected using non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,pkfrcc_atcorner./pkfrncc_atcorner];
        %corrected using simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,pkfrcc_atcorner./pkfrsim_atcorner];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
end
save('F:\Results_experimentF\pkfrcc_atCorner_R1.mat','pkfrcc_atCnc','pkfrcc_atCs')
data = cellfun(@(x) mean(x,2), pkfrcc_atCs, 'uni', 0);

%% Percentage of boundary cells for each session
prop_bcell = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('border_metrics.mat','bcell_raw')
    bcell = bcell_raw;
    load('thresh.mat','thresh')
    if contains(datapath{n},'M4113F') || contains(datapath{n},'M4131F')
        idx = [1,2,3,4];
    else
        idx = [1,2,3,5];
    end
    bc_all = bcell(idx);
    numNeurons = length(thresh);
    bc = cellfun(@numel, bc_all)./numNeurons;
    prop_bcell = [prop_bcell;bc];
end
%note for saved data: data were organized as baseline1, baseline2, dark,
%whisker. 
save('F:\Results_experimentF\prop_cell_dark_whisker.mat','prop_cornercell','prop_placecell')

%% OCCUPANCY MATCHING
%% get place cell, cornercell, and boundary cell with occupancy matching
%% placecellm
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat');
    if contains(datapath{k},'M4113F') || contains(datapath{k},'M4131F')
        idx = [2,3,4];
    else
        idx = [2,3,5];
    end
    neuronIndividualsf = neuronIndividualsf(idx);
    behavIndividualsf = behavIndividualsf(idx);
    niter = 24;
    nss = length(neuronIndividualsf);
    metricsm = cell(niter,length(neuronIndividualsf));
    metricsm_avg = cell(1,length(neuronIndividualsf));
    parfor n = 1:niter
        [neuronIndivm, behavIndivm] = match_neuron_behav(neuronIndividualsf, behavIndividualsf, 'longitudinal_full');
        for ii = 1:nss
            neuronm = neuronIndivm{1,ii};
            behavm = behavIndivm{1,ii};
            [~,TinfoL,~] = permuting_spikes(neuronm,behavm,thresh,true);
            metricsm{n,ii} = table2array(TinfoL); %defined by bit/sec
        end
    end
    %average the value for each iteration
    for jj = 1:size(metricsm,2)
        metrics_each = squeeze(metricsm(:,jj));
        metrics_each = cat(3, metrics_each{:});
        metricsm_avg{1,jj} = array2table(nanmean(metrics_each,3),...
            'VariableNames',{'neuron','bitpsec','bitpsec_pcidx','bitpspike','bitpspike_pcidx','meanfr','peakfr'});
    end
    %get place cells
    placecellm = cellfun(@(x) find(x.bitpspike_pcidx >= 0.3 &...
        x.meanfr > quantile(x.meanfr,0.05)), metricsm_avg, 'uni',0);
    %save data
    save('spatial_metrics.mat', 'metricsm', 'metricsm_avg','placecellm', '-append');
end
%results were stored on Prism for plot
prop_placecellm = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat','placecellm')
    load('thresh.mat','thresh')
    pc_all = placecellm;
    numNeurons = length(thresh);
    pc = cellfun(@numel, pc_all)./numNeurons;
    prop_placecellm = [prop_placecellm;pc];
end

%% cornercellm
for m = 1:length(datapath)
    cd(datapath{m});
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat');
    load('env_geometry.mat', 'env_coor');
    if contains(datapath{m},'M4113F') || contains(datapath{m},'M4131F')
        idx = [2,3,4];
    else
        idx = [2,3,5];
    end
    neuronIndividualsf = neuronIndividualsf(idx);
    behavIndividualsf = behavIndividualsf(idx);
    niter = 24;
    nss = length(neuronIndividualsf);
    Cm = cell(niter,1);
    parfor n = 1:niter
        [neuronIndivm, behavIndivm] = match_neuron_behav(neuronIndividualsf, behavIndividualsf, 'longitudinal_full');
        firingrateAll = cell(1,nss);
        ratemap = cell(1,nss);
        % Calculate spatial map and mean firing rate for each individual sessions
        for k = 1:nss
            [firingrateAll{k},~,~,~] = calculate_firing_ratemap(neuronIndivm{k},behavIndivm{k},thresh,2,true);
        end
        % Calculate smoothed firing rate maps
        for ii = 1:nss
            fr = firingrateAll{ii};
            % to remove potential inf in the unsmoothed ratemap
            idx = cellfun(@isinf, fr,'UniformOutput',false);
            numcell = length(fr);
            for jj = 1:numcell
                fr{jj}(idx{jj}) = 0;
            end
            % smoothed ratemap
            ratemap{ii} = cellfun(@(x) filter2DMatrices(x,1), fr, 'uni', 0);
        end
        Cm{n} = identify_corner_cell_match(neuronIndivm, behavIndivm, thresh, env_coor, ratemap);
    end
    save('corner_metrics.mat', 'Cm', '-append');
end
%results were stored on Prism for plot
ncc = [];
for m = 1:length(datapath)
    cd(datapath{m});
    load('corner_metrics.mat', 'Cm');
    load('thresh.mat','thresh');
    numN = length(thresh);
    ncc_ea = [];
    for ii = 1:length(Cm)
        C_ea = Cm{ii};
        ncc_ea = [ncc_ea;cellfun(@numel, C_ea.cornercell)];
    end
    ncc = [ncc;mean(ncc_ea)./numN];
end

