%% Experiment E
% This code is for analyzing corner coding in the subiculum for inserted
% corner experiments (icm)
load('F:\analysis_folders.mat','expE')
datapath = expE;
%% Re-define corner cells for revision
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    CR1 = identify_corner_cell([],[],0.4);
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
    save('spatial_metrics.mat', 'map_stb')
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

%% compute the peak fr at the inserted corner across different sessions for corner cells
%simulate neuron to correct the firing rate
%note: M4114F was used simulate_ratemap_xsession, as the firing rate
%difference between icm0 and icm60 was greater than 10%
for n = 1:length(datapath)
    cd(datapath{n})
    %simulate a ratemap with constant firing
    [neuronSim,ratemapSim] = simulate_ratemap;
    [neuronSimx,ratemapSimx] = simulate_ratemap_xsession;
    save('neuronSim.mat','neuronSim','ratemapSim','neuronSimx','ratemapSimx','-v7.3')
end
%for corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('neuronSim.mat','ratemapSim')
    load('env_geometry.mat','env_coori','S')
    load('firingrateAll.mat','ratemap')
    if strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4131Fe')
        load('corner_metrics.mat','CR1')
        C = CR1;
    else
        load('corner_metricsR1.mat','C')
    end
    idx = strcmp(S,'icm0');
    corners_icm0 = env_coori(idx);
    icm_loc = corners_icm0{1}(end-1,:);
    pkfr_icm = cell(1,length(ratemap));
    pkfr_icm_Sim = cell(1,length(ratemap));
    for ii = 1:length(ratemap)
        ratemap_ea = ratemap{ii};
        ratemapSim_ea = ratemapSim{ii};
        pkfr_icm{ii} = cellfun(@(x) max(max(x(icm_loc(2)-1:icm_loc(2)+1,...
            icm_loc(1)-1:icm_loc(1)+1))), ratemap_ea);
        pkfr_icm_Sim{ii} = cellfun(@(x) max(max(x(icm_loc(2)-1:icm_loc(2)+1,...
            icm_loc(1)-1:icm_loc(1)+1))), ratemapSim_ea);
    end
    %for corner cells
    if strcmp(datapath{n}, 'F:\subiculum_mice_May2022\M4104e')
        cornercell = union(C.cornercell{1},C.cornercell{2});
    else
        idx_cc = strcmp(S,'icmb');
        cornercell = cell2mat(C.cornercell(idx_cc));
    end
    pkfrcc_icm = cellfun(@(x) x(cornercell), pkfr_icm, 'uni', 0);
    pkfrcc_icmS = cellfun(@(x,y) x./y, pkfrcc_icm, pkfr_icm_Sim, 'uni', 0);
    save('corner_metricsR1.mat','pkfrcc_icm','pkfrcc_icmS','-append')
end
% for non-corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','env_coori','S')
    load('firingrateAll.mat','ratemap')
    load('corner_metricsR1.mat','C')
    load('neuronSim.mat','ratemapSimx')
    idx = strcmp(S,'icm0');
    corners_icm0 = env_coori(idx);
    icm_loc = corners_icm0{1}(end-1,:);
    pkfr_icm = cell(1,length(ratemap));
    pkfr_icm_Sim = cell(1,length(ratemap));
    for ii = 1:length(ratemap)
        ratemap_ea = ratemap{ii};
        ratemapSim_ea = ratemapSimx{ii};
        pkfr_icm{ii} = cellfun(@(x) max(max(x(icm_loc(2)-1:icm_loc(2)+1,...
            icm_loc(1)-1:icm_loc(1)+1))), ratemap_ea);
        pkfr_icm_Sim{ii} = cellfun(@(x) max(max(x(icm_loc(2)-1:icm_loc(2)+1,...
            icm_loc(1)-1:icm_loc(1)+1))), ratemapSim_ea);
    end
    idx_cc = strcmp(S,'icmb');
    cornercell = cell2mat(C.cornercell(idx_cc));
    allcell = 1:length(ratemap{1});
    idx2 = ~ismember(allcell, cornercell);
    noncc = allcell(idx2);
    pkfrncc_icm = cellfun(@(x) x(noncc), pkfr_icm, 'uni', 0);
    pkfrncc_icmS = cellfun(@(x,y) x./y, pkfrncc_icm, pkfr_icm_Sim, 'uni', 0);
    save('corner_metricsR1.mat','pkfrncc_icm','pkfrncc_icmS','-append')
end

% Organize the data for plot
pkfricmS = [];
pkfricmS_ncc = [];
pkfricm = [];
pkfricm_ncc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','pkfrcc_icmS','pkfrncc_icmS','pkfrcc_icm','pkfrncc_icm')
    %corrected firing rate
    pkfrS_avg = cellfun(@mean, pkfrcc_icmS);
    pkfrnccS_avg = cellfun(@mean, pkfrncc_icmS);
    if length(pkfrS_avg) > 5
        if strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4113Fe')||...
                strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4114Fe')
            pkfrS_avg = pkfrS_avg([2,3,4,5,6]);
            pkfrnccS_avg = pkfrnccS_avg([2,3,4,5,6]);
        else
            pkfrS_avg = pkfrS_avg([2,1,4,5,6]);
            pkfrnccS_avg = pkfrnccS_avg([2,1,4,5,6]);
        end
    end
    pkfricmS = [pkfricmS;pkfrS_avg];    
    pkfricmS_ncc = [pkfricmS_ncc;pkfrnccS_avg];
    %uncorrected firing rate
    pkfr_avg = cellfun(@mean, pkfrcc_icm);
    pkfrncc_avg = cellfun(@mean, pkfrncc_icm);
    if length(pkfr_avg) > 5
        pkfr_avg = pkfr_avg([2,1,4,5,6]);
        pkfrncc_avg = pkfrncc_avg([2,1,4,5,6]);
    end
    pkfricm = [pkfricm;pkfr_avg];
    pkfricm_ncc = [pkfricm_ncc;pkfrncc_avg];
end
save('F:\Results_experimentE\pkfr_icm_R1.mat','pkfricmS','pkfricmS_ncc','pkfricm','pkfricm_ncc')

%% compute the peak fr at the regular environmental corners across different sessions for corner cells
load('F:\analysis_folders.mat','expE')
datapath = expE;
%find peak fr at each corner for each neuron
for n = 1:length(datapath)
    cd(datapath{n})
    %find the pkfr at corners for corner cells
    pkfr_corners = find_pkfr_at_corners;
    %find the pkfr at corners of the simulated neuron
    load('neuronSim.mat','ratemapSimx')
    pkfrsim_corners = find_pkfr_at_corners(ratemapSimx);
    save('corner_metricsR1.mat','pkfr_corners','pkfrsim_corners','-append')
end

%gather data for plot
pkfrcc.srdC = [];
pkfrcc.ctrC = [];%this number is not accurate, do not use for analysis. 
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','env_coori','S')
    if strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4131Fe')
        load('corner_metrics.mat','CR1')
        C = CR1;
    else
        load('corner_metricsR1.mat','C')
    end
    if strcmp(datapath{n}, 'F:\subiculum_mice_May2022\M4104e')
        cornercell = union(C.cornercell{1},C.cornercell{2});
    else
        idx = strcmp(S,'icmb');
        cornercell = C.cornercell{idx};
    end
    load('corner_metricsR1.mat','pkfr_corners','pkfrsim_corners')
    if length(pkfr_corners) > 5
        if strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4113Fe')||...
                strcmp(datapath{n}, 'F:\subiculum_mice_Aug2022\M4114Fe')
            pkfr_corners = pkfr_corners([2,3,4,5,6]);
            pkfrsim_corners = pkfrsim_corners([2,3,4,5,6]);
        else
            pkfr_corners = pkfr_corners([2,1,4,5,6]);
            pkfrsim_corners = pkfrsim_corners([2,1,4,5,6]);
        end
    end
    pkfr_atCors = cellfun(@(x,y) x./y, pkfr_corners, pkfrsim_corners, 'uni', 0);
    pkfrcc_atCors = cellfun(@(x) mean(x(:,cornercell),2), pkfr_atCors, 'uni', 0);
    pkfrcc_atCors{1,1}(5,1) = NaN;
    pkfrcc_atSrdC = cellfun(@(x) mean(x(1:4)), pkfrcc_atCors);
    pkfrcc_atCtrC = cellfun(@(x) x(5,1), pkfrcc_atCors);
    pkfrcc.srdC = [pkfrcc.srdC;pkfrcc_atSrdC];
    pkfrcc.ctrC = [pkfrcc.ctrC;pkfrcc_atCtrC];
end
save('F:\Results_experimentE\pkfr_icm_R1.mat','pkfrcc','-append')

%% Plot the whole data
%plot
load('F:\Results_experimentE\pkfr_icm_R1.mat')
data = pkfrcc.srdC(:,2:end);
figure
subplot(2,3,1)
plot([1.15:1:4.15], data',':o','Color',[128 130 133]/255,'markersize',6,'linewidth',0.1)
hold on
err = std(data)/sqrt(size(data,1));
errorbar(mean(data), err,'Marker','.','markersize',30,'Color','k','linewidth',1)
xlim([0.5,4.5])
ylim([0,3.5])
axis square
xticklabels({'0 cm','1.5 cm','3.0 cm','6.0 cm'})
ylabel('corrected firng rate (fold)')
subplot(2,3,2)
plot([1.15:1:5.15], pkfricmS',':o','markersize',6,'linewidth',0.1)
hold on
err = nanstd(pkfricmS)/sqrt(size(pkfricmS,1));
errorbar(nanmean(pkfricmS), err,'Marker','.','markersize',30,'Color','k','linewidth',1)
xlim([0.5,5.5])
ylim([0,3.5])
axis square
xticklabels({'baseline','0 cm','1.5 cm','3.0 cm','6.0 cm'})
ylabel('corrected firng rate (fold)')
subplot(2,3,3)
plot([1.15:1:5.15], pkfricmS_ncc',':o','markersize',6,'linewidth',0.1)
hold on
err = nanstd(pkfricmS_ncc)/sqrt(size(pkfricmS_ncc,1));
errorbar(nanmean(pkfricmS_ncc), err,'Marker','.','markersize',30,'Color','k','linewidth',1)
xlim([0.5,5.5])
ylim([0,3.5])
axis square
xticklabels({'baseline','0 cm','1.5 cm','3.0 cm','6.0 cm'})
ylabel('corrected firng rate (fold)')
subplot(2,3,5)
plot([1.15:1:5.15], pkfricm',':o','Color',[128 130 133]/255,'markersize',6,'linewidth',0.1)
hold on
err = std(pkfricm)/sqrt(size(pkfricm,1));
errorbar(nanmean(pkfricm), err,'Marker','.','markersize',30,'Color','k','linewidth',1)
xlim([0.5,5.5])
ylim([0,1.5])
axis square
xticklabels({'baseline','0 cm','1.5 cm','3.0 cm','6.0 cm'})
ylabel('firng rate (Hz)')
subplot(2,3,6)
plot([1.15:1:5.15], pkfricm_ncc',':o','Color',[128 130 133]/255,'markersize',6,'linewidth',0.1)
hold on
err = std(pkfricm_ncc)/sqrt(size(pkfricm_ncc,1));
errorbar(nanmean(pkfricm_ncc), err,'Marker','.','markersize',30,'Color','k','linewidth',1)
xlim([0.5,5.5])
ylim([0,1.5])
axis square
xticklabels({'baseline','0 cm','1.5 cm','3.0 cm','6.0 cm'})
ylabel('firng rate (Hz)')

%% CONVEX
%% CONVEX CORNER ANALYSIS
%% Convex corner anlaysis for icm experiment 
%% problem: too few neurons were identified
% load('F:\analysis_folders.mat','expE')
% datapath = expE;
% % identify convex corners from icm0
% for ii = 1:length(datapath)
%     cd(datapath{ii})
%     [cvex_coor,cvex_coori] = identify_env_geometry_icmcvx;
%     save('env_geometry.mat','cvex_coor','cvex_coori','-append')
% end
% 
% %determine corner cells in each session
% for ii = 1:length(datapath)
%     cd(datapath{ii})
%     load('env_geometry.mat','S')
%     mask = cell(1,length(S));
%     for jj = 1:length(S)
%         mask{jj} = true;
%     end
%     %determine corner cell for each session by spike shuffling.
%     %NOTE, for expC, C was using a 0.3/0.35 threshold for identifying
%     %corner cells. C2 was using a 0.4 threshold for identifying corner
%     %cells. C2 works slightly better becuase the data in the large
%     %environment is bit noiser.
%     CV = identify_corner_cell(mask);
%     save('corner_metrics.mat','CV','-append')
% end



