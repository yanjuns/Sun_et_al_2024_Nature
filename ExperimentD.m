%% Experiment D
% This code is for analyzing corner coding in the subiculum for animals
% running in a right triangle and trapezoid. 
load('F:\analysis_folders.mat','expD')
datapath = expD;
%% Corner cell with different threshold
% get corner cells at different threshold level
threshlevel = [0.1:0.1:0.6];
for n = 1:length(datapath)
    Cthresh = cell(1,length(threshlevel));
    cd(datapath{n})
    for ii = 1:length(threshlevel)
        %determine corner cell for each session by spike shuffling.
        Cthresh{ii} = identify_corner_cell([], [], threshlevel(ii));
    end
    save('corner_metrics.mat','Cthresh','-append')
end
%% Stability
%map_stb: using interleaved method
%map_stability: using split to half method
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat')
    load('behavIndividualsf.mat')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndividualsf, behavIndividualsf, 'uni',0);
    if ~exist('spatial_metrics.mat','file')
        save('spatial_metrics.mat','map_stb')
    else
        save('spatial_metrics.mat','map_stb','-append')
    end
end
%% Revised definition of corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'Cthresh')
    C = Cthresh{4};
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C')    
end

%% Percentage of corner cells for each session
session_name = {'rightTri','trapezoid','rightTrim1'}; %expD
prop_cornercell = cell(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    load('thresh.mat','thresh')
    numNeurons = length(thresh);
    pcc = cell(1,length(session_name));
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        C_ea = C.cornercell(idx);
        C_ea = cellfun(@numel, C_ea);
        pcc{ii} = C_ea(:)./numNeurons;
    end
    prop_cornercell(n,:) = pcc;
end
prop_cornercell_ea = cellfun(@mean, prop_cornercell);
save('F:\Results_experimentD\prop_cornercell_R1.mat','prop_cornercell','prop_cornercell_ea')
%plot using Graphpad by adding the data from the circle sessions

%% Corner cell firing rate of each corner in the rightTri and trapezoid environments
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
session_name = {'rightTri','trapezoid'}; %expD
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
pkfrcc_atC = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCornc = []; pkfrcc_atCors = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
        load('corner_metricsR1.mat','C')
        idx = strcmp(S, session_name{k});
        temp = cellfun(@numel, C.cornercell);
%         idx2 = temp >= 3;
%         idx = idx1 & idx2;
        pkfr_corners = pkfr_corners(idx);
        cornercell = C.cornercell(idx);
        allcell = 1:length(pkfr_corners{1,1});
        ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
        ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,nanmean(pkfrcc_atcorner,2)];
        %correction using non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCornc = [pkfrcc_atCornc,nanmean(pkfrcc_atcorner./pkfrncc_atcorner,2)];
        %correction using simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,nanmean(pkfrcc_atcorner./pkfrsim_atcorner,2)];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCornc';
    pkfrcc_atCs{k} = pkfrcc_atCors';
    pkfrcc_atC{k} = pkfrcc_atCor';
end
save('F:\Results_experimentD\pkfrcc_atCorner_R1mouse0.4.mat','pkfrcc_atCnc','pkfrcc_atCs','pkfrcc_atC')
%plot using GraphPad with some rearrangement of the data

%% DETERMINE RANDOM LEVEL OF CORNER CELLS (simulated corners or CA1 data)
%% Using randomly assigned corners
load('F:\analysis_folders.mat','expD')
datapath = expD;
session_name = {'rightTri','trapezoid','rightTrim1'}; %expD
%mannually assign the location of corners for each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [~,~,env_coorS,env_coorSi] = identify_env_geometry;
    save('env_geometry.mat','env_coorS','env_coorSi','-append')
end
%determin corner cells in each session using simulated corner locations
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('env_geometry.mat','env_coorS')
    Csim = identify_corner_cell([], [], [], env_coorS);
    save('corner_metrics.mat','Csim','-append')
end
% proportion of corner cells
propCC_randC = cell(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metrics.mat','Csim')
    C = Csim;
    load('spatial_metrics.mat', 'map_stb')
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    load('env_geometry.mat','S') %load all session name
    load('thresh.mat','thresh')
    numNeurons = length(thresh);
    pcc = cell(1,length(session_name));
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        C_ea = C.cornercell(idx);
        C_ea = cellfun(@numel, C_ea);
        pcc{ii} = C_ea(:)/numNeurons;
    end
    propCC_randC(n,:) = pcc;
end
propCC_randC_ea = cellfun(@mean, propCC_randC);
save('F:\Results_experimentD\prop_cornercell_R1random.mat','propCC_randC','propCC_randC_ea')

%% ARCHIVE
%% archived code
%% Corner cell firing rate of each corner in the trapezoid environment (old version)
% load('F:\analysis_folders.mat','expD')
% datapath = expD;
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('corner_metrics.mat', 'C')
%     load('env_geometry.mat', 'env_coor', 'S')
%     load('firingrateAll.mat', 'ratemap')
%     idx = strcmp(S, 'trapezoid');
%     corners = env_coor(idx);
%     cornercell = C.cornercell(idx);
%     ratemap = ratemap(idx);
%     ratemapc = cellfun(@(x,y) x(y), ratemap, cornercell, 'uni',0);
%     trapzPeaks = {};
%     for n = 1:length(ratemapc)
%         ratemapss = ratemapc{n};
%         c1 = [ceil(corners{n}(1,2)),ceil(corners{n}(1,1))];
%         c2 = [floor(corners{n}(2,2)),ceil(corners{n}(2,1))];
%         c3 = [floor(corners{n}(3,2)),floor(corners{n}(3,1))];
%         c4 = [ceil(corners{n}(4,2)),floor(corners{n}(4,1))];
%         c1area = [c1(1):c1(1)+4;c1(2):c1(2)+4];
%         c2area = [c2(1)-4:c2(1);c2(2):c2(2)+4];
%         c3area = [c3(1)-4:c3(1);c3(2)-4:c3(2)];
%         c4area = [c4(1):c4(1)+4;c4(2)-2:c4(2)+2];
%         cornerpeaks = [];
%         for ii = 1:length(ratemapss)
%             rtmap = ratemapss{ii};
%             c1peak = max(max(rtmap(c1area(1,:),c1area(2,:))));
%             c2peak = max(max(rtmap(c2area(1,:),c2area(2,:))));
%             c3peak = max(max(rtmap(c3area(1,:),c3area(2,:))));
%             c4peak = max(max(rtmap(c4area(1,:),c4area(2,:))));
%             cpeak = [c1peak,c2peak,c3peak,c4peak];
%             cornerpeaks = [cornerpeaks;cpeak];
%         end
%         trapzPeaks{n} = cornerpeaks;
%     end
%     save('corner_metrics.mat','trapzPeaks','-append')
% end
% 
% peaks_tpz = [];
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('corner_metrics.mat', 'trapzPeaks')
%     peaksz = cellfun(@mean, trapzPeaks, 'uni', 0);
%     peaksz = vertcat(peaksz{:});
%     peaks_tpz = [peaks_tpz;peaksz];
% end
% peaks_tpz(:,2) = mean(peaks_tpz(:,1:2),2);
% peaks_tpz = peaks_tpz(:,2:end);
% peaks_tpz(:,1:2) = [peaks_tpz(:,2),peaks_tpz(:,1)];
%% determine threshold 0.1-0.6
% % arrange the results
% session_name = {'rightTri','trapezoid'}; %expA
% propCC = cell(1,length(session_name));
% for ii = 1:length(session_name)
%     propCC_ea = [];
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('env_geometry.mat','S')
%         load('corner_metrics.mat','Cthresh')
%         idx = strcmp(S, session_name{ii});
%         numCell = length(Cthresh{1}.cscore{1});
%         numCC = [];
%         for jj = 1:length(Cthresh)
%             Cea = Cthresh{jj};
%             Cea_ss = Cea.cornercell(idx);
%             Cea_ss = cellfun(@numel, Cea_ss)/numCell;
%             numCC = [numCC,Cea_ss'];
%         end
%         propCC_ea = [propCC_ea;numCC];
%     end
%     propCC{ii} = propCC_ea;
% end
% % compute the the firing rate as a function of distance to corners for each
% % threshold
% % frCC is number of threshold by type of sessions. 
% threshlevel = [0.1:0.1:0.6];
% frCC = cell(length(threshlevel),length(session_name));
% for ii = 1:length(threshlevel)
%     frcc = {};
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('corner_metrics.mat','Cthresh')
%         d2cornerq = corner_dist_firingrate_quantify(session_name,[],Cthresh{ii});
%         frcc = [frcc;d2cornerq.frcc];
%     end
%     dsize = cellfun(@(x) size(x,2), frcc);
%     msize = min(dsize);
%     msize = repmat(msize, length(dsize), 1);
%     msize = num2cell(msize);
%     frcc = cellfun(@(x,y) x(:,1:y), frcc, msize, 'uni', false); 
%     % allocate each session's data into corresponding cells
%     frcc1 = frcc(:,1);
%     frCC{ii,1} = vertcat(frcc1{:});
%     frcc2 = frcc(:,2);
%     frCC{ii,2} = vertcat(frcc2{:});
% end
% data = propCC; 
% figure
% % plot prop of corner cells for each threshold
% subplot(3,2,1)
% stdshade(data{1},0.3, [65 68 135]/255)
% xlim([0.5,6.5])
% % ylim([0.02, 0.12])
% axis square
% ylabel('prop of neurons')
% xlabel('threshold (x 0.1)')
% title('triangle')
% subplot(3,2,2)
% stdshade(data{2},0.3, [43 120 142]/255)
% xlim([0.5,6.5])
% % ylim([0.02, 0.12])
% ylabel('prop of neurons')
% xlabel('threshold (x 0.1)')
% axis square
% title('square')
% % plot firing rate for each threshold
% subplot(3,2,3)
% hold on
% for ii = 1:size(frCC,1)
%     plot([1:length(mean(frCC{ii,1}))]*1.55, smoothdata(nanmean(frCC{ii,1}), 2,'movmean',3)')
% end
% axis square
% xlim([0,20])
% % ylim([0.05,0.5])
% xlabel('distance to corners (cm)')
% ylabel('spike rate (Hz)')
% subplot(3,2,4)
% hold on
% for ii = 1:size(frCC,1)
%     plot([1:length(mean(frCC{ii,2}))]*1.55, smoothdata(nanmean(frCC{ii,2}), 2,'movmean',3)')
% end
% axis square
% xlim([0,20])
% % ylim([0.05,0.5])
% xlabel('distance to corners (cm)')
% ylabel('spike rate (Hz)')
% % plot the difference between max and min firing rate 
% subplot(3,2,5)
% frCC_diff = cellfun(@(x) max(nanmean(x)) - min(nanmean(x)), frCC);
% plot(frCC_diff(:,1)', 'Color',[65 68 135]/255, 'LineWidth', 2)
% axis square
% xlim([0.5, 6.5])
% % ylim([0.25, 0.45])
% xlabel('threshold (x 0.1)')
% ylabel('spike rate difference (Hz)')
% subplot(3,2,6)
% plot(frCC_diff(:,2)', 'Color',[43 120 142]/255, 'LineWidth', 2)
% axis square
% xlim([0.5, 6.5])
% % ylim([0.25, 0.45])
% xlabel('threshold (x 0.1)')
% ylabel('spike rate difference (Hz)')
