%% Nature revision
%% 1. changes made for the definition of corner cells
% To make penalty linear, the extra field that near the corner will be less
% panlized than the extra field that near the center. To achieve this, the
% panelty score that range from -1 to 1 was subctracted by 1, which bring
% it to -2 to 0, and the absolute value was taken for the panelty. 
% Related functions:
% compute_corner_score.m

% Added a variable 'threshlevel' to determine the best threshold

% For field distance, instead of using all pairs of field's distance, I
% used the distance for only the major fields. 

%% CHANGE THE THRESHOLD FROM 0.1 - 0.6
%% Conclusion: 0.3-0.4 is the best range of threshold
%% Change detection threshold for experiment A
load('F:\analysis_folders.mat','expA')
datapath = expA;
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
% arrange the results
session_name = {'triangle','square','hex','circle'}; %expA
propCC = cell(1,length(session_name));
for ii = 1:length(session_name)
    propCC_ea = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','Cthresh')
        idx = strcmp(S, session_name{ii});
        numCell = length(Cthresh{1}.cscore{1});
        numCC = [];
        for jj = 1:length(Cthresh)
            Cea = Cthresh{jj};
            Cea_ss = Cea.cornercell(idx);
            Cea_ss = cellfun(@numel, Cea_ss)/numCell;
            numCC = [numCC,Cea_ss'];
        end
        propCC_ea = [propCC_ea;numCC];
    end
    propCC{ii} = propCC_ea;
end
save('F:\Results_revision\corner_metrics_thresh.mat','propCC','-v7.3')

% compute the the firing rate as a function of distance to corners for each
% threshold
session_name = {'triangle','square','hex'}; %expA
% frCC is number of threshold by type of sessions. 
frCC = cell(length(threshlevel),length(session_name));
for ii = 1:length(threshlevel)
    frcc = {};
    for n = 1:length(datapath)
        cd(datapath{n})
        load('corner_metrics.mat','Cthresh')
        d2cornerq = corner_dist_firingrate_quantify(session_name,[],Cthresh{ii});
        frcc = [frcc;d2cornerq.frcc];
    end
    dsize = cellfun(@(x) size(x,2), frcc);
    msize = min(dsize);
    msize = repmat(msize, length(dsize), 1);
    msize = num2cell(msize);
    frcc = cellfun(@(x,y) x(:,1:y), frcc, msize, 'uni', false); 
    % allocate each session's data into corresponding cells
    frcc1 = frcc(:,1);
    frCC{ii,1} = vertcat(frcc1{:});
    frcc2 = frcc(:,2);
    frCC{ii,2} = vertcat(frcc2{:});
    frcc3 = frcc(:,3);
    frCC{ii,3} = vertcat(frcc3{:});
end
save('F:\Results_revision\corner_metrics_thresh.mat','frCC','-append')

data = propCC; 
figure
% plot prop of corner cells for each threshold
subplot(3,4,1)
stdshade(data{1},0.3, [65 68 135]/255)
xlim([0.5,6.5])
ylim([0.02, 0.12])
axis square
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
title('triangle')
subplot(3,4,2)
stdshade(data{2},0.3, [43 120 142]/255)
xlim([0.5,6.5])
ylim([0.02, 0.12])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('square')
subplot(3,4,3)
stdshade(data{3},0.3, [34 168 132]/255)
xlim([0.5,6.5])
ylim([0.02, 0.12])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('hexagon')
subplot(3,4,4)
stdshade(data{4},0.3, [68 1 84]/255)
xlim([0.5,6.5])
ylim([0, 0.01])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('circle')
% plot firing rate for each threshold
subplot(3,4,5)
hold on
for ii = 1:size(frCC,1)
    plot([1:length(mean(frCC{ii,1}))]*1.55, smoothdata(mean(frCC{ii,1}), 2,'movmean',3)')
end
axis square
xlim([0,20])
ylim([0.05,0.5])
xlabel('distance to corners (cm)')
ylabel('spike rate (Hz)')
subplot(3,4,6)
hold on
for ii = 1:size(frCC,1)
    plot([1:length(mean(frCC{ii,2}))]*1.55, smoothdata(mean(frCC{ii,2}), 2,'movmean',3)')
end
axis square
xlim([0,20])
ylim([0.05,0.5])
xlabel('distance to corners (cm)')
ylabel('spike rate (Hz)')
subplot(3,4,7)
hold on
for ii = 1:size(frCC,1)
    plot([1:length(mean(frCC{ii,3}))]*1.55, smoothdata(mean(frCC{ii,3}), 2,'movmean',3)')
end
axis square
xlim([0,20])
ylim([0.05,0.5])
xlabel('distance to corners (cm)')
ylabel('spike rate (Hz)')
% plot the difference between max and min firing rate 
subplot(3,4,9)
frCC_diff = cellfun(@(x) max(mean(x)) - min(mean(x)), frCC);
plot(frCC_diff(:,1)', 'Color',[65 68 135]/255, 'LineWidth', 2)
axis square
xlim([0.5, 6.5])
ylim([0.25, 0.45])
xlabel('threshold (x 0.1)')
ylabel('spike rate difference (Hz)')
subplot(3,4,10)
plot(frCC_diff(:,2)', 'Color',[43 120 142]/255, 'LineWidth', 2)
axis square
xlim([0.5, 6.5])
ylim([0.25, 0.45])
xlabel('threshold (x 0.1)')
ylabel('spike rate difference (Hz)')
subplot(3,4,11)
plot(frCC_diff(:,3)', 'Color',[34 168 132]/255, 'LineWidth', 2)
axis square
xlim([0.5, 6.5])
ylim([0.25, 0.45])
xlabel('threshold (x 0.1)')
ylabel('spike rate difference (Hz)')

%% Change detection threshold for experiment C
load('F:\analysis_folders.mat','expC')
datapath = expC;
threshlevel = [0.1:0.1:0.6];
for n = 1:length(datapath)
    Cthresh = cell(1,length(threshlevel));
    cd(datapath{n})
    load('env_geometry.mat','S')
    mask = cell(1,length(S));
    for jj = 1:length(S)
        if strcmp(S{jj}, 'largesq') || strcmp(S{jj}, 'largeRT')
            mask{jj} = false;
        else
            mask{jj} = true;
        end
    end
    %determine corner cell for each session by spike shuffling.
    %NOTE, for expC, C was using a 0.3/0.35 threshold for identifying
    %corner cells. C2 was using a 0.4 threshold for identifying corner
    %cells. C2 works slightly better becuase the data in the large
    %environment is bit noiser.
    for ii = 1:length(threshlevel)
        Cthresh{ii} = identify_corner_cell(mask, [], threshlevel(ii));
    end
    save('corner_metrics.mat','Cthresh','-append')
end

% plot the results
session_name = {'largesq','largeRT','largeX','largeXm1'}; %expA
propCC = cell(1,length(session_name));
for ii = 1:length(session_name)
    propCC_ea = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','Cthresh')
        idx = strcmp(S, session_name{ii});
        numCell = length(Cthresh{1}.cscore{1});
        numCC = [];
        for jj = 1:length(Cthresh)
            Cea = Cthresh{jj};
            Cea_ss = Cea.cornercell(idx);
            Cea_ss = cellfun(@numel, Cea_ss)/numCell;
            numCC = [numCC,Cea_ss'];
        end
        propCC_ea = [propCC_ea;numCC];
    end
    propCC{ii} = propCC_ea;
end
data = propCC; 
figure
subplot(1,4,1)
stdshade(data{1},0.3, [65 68 135]/255)
xlim([0.5,6.5])
% ylim([0.02, 0.12])
axis square
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
title('square')
subplot(1,4,2)
stdshade(data{2},0.3, [43 120 142]/255)
xlim([0.5,6.5])
% ylim([0.02, 0.12])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('rectangle')
subplot(1,4,3)
stdshade(data{3},0.3, [34 168 132]/255)
xlim([0.5,6.5])
% ylim([0.02, 0.12])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('convex')
subplot(1,4,4)
stdshade(data{4},0.3, [68 1 84]/255)
xlim([0.5,6.5])
% ylim([0, 0.01])
ylabel('prop of neurons')
xlabel('threshold (x 0.1)')
axis square
title('convex_m1')

%% DISTRIBUTION OF CORNER SCORES AND SHUFFLE
%% Plot shuffled distribution of corner score and compare to the actual corner score
load('F:\analysis_folders.mat','expA')
datapath = expA;
session_name = {'triangle','square','hex'}; %expA
allcscore = cell(length(datapath), length(session_name));
allcscore_shuffle = cell(length(datapath), length(session_name)); 
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat', 'Cthresh')
    C = Cthresh{4};
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        cornercell = C.cornercell(idx);
        cscore = C.cscore(idx);
        cscore_shuffle = C.cscore_shuffle(idx);
%         CCscore_ea = cellfun(@(x,y) x(y,:), cscore, cornercell, 'UniformOutput',false);
%         CCscore_shuffle_ea = cellfun(@(x,y) x(y,:), cscore_shuffle, cornercell, 'UniformOutput',false);
        CCscore_ea = cscore;
        CCscore_shuffle_ea = cscore_shuffle;
        allcscore{n,ii} = vertcat(CCscore_ea{:});
        allcscore_shuffle{n,ii} = vertcat(CCscore_shuffle_ea{:});
    end
end
% for all cells
% reformat the data
Cscore_all = struct;
score_triangle = allcscore(:,1);
Cscore_all.triangle = vertcat(score_triangle{:});
score_square = allcscore(:,2);
Cscore_all.square= vertcat(score_square{:});
score_hex = allcscore(:,3);
Cscore_all.hex= vertcat(score_hex{:});
% reformat shuffled data
score_triangle = allcscore_shuffle(:,1);
Cscore_all.triangle_shuffle = vertcat(score_triangle{:});
score_square = allcscore_shuffle(:,2);
Cscore_all.square_shuffle= vertcat(score_square{:});
score_hex = allcscore_shuffle(:,3);
Cscore_all.hex_shuffle= vertcat(score_hex{:});
save('F:\Results_revision\corner_metrics_randomlevel.mat','Cscore_all','-append')
% plot (use the plot togather with shuffle comparison in the next section)

% % for corner cell only
% % reformat the data
% Cscore_cc = struct;
% score_triangle = allcscore(:,1);
% Cscore_cc.triangle = vertcat(score_triangle{:});
% score_square = allcscore(:,2);
% Cscore_cc.square= vertcat(score_square{:});
% score_hex = allcscore(:,3);
% Cscore_cc.hex= vertcat(score_hex{:});
% % reformat shuffled data
% score_triangle = allcscore_shuffle(:,1);
% Cscore_cc.triangle_shuffle = vertcat(score_triangle{:});
% score_square = allcscore_shuffle(:,2);
% Cscore_cc.square_shuffle= vertcat(score_square{:});
% score_hex = allcscore_shuffle(:,3);
% Cscore_cc.hex_shuffle= vertcat(score_hex{:});
% save('F:\Results_revision\corner_metrics_randomlevel.mat','Cscore_cc','-append')
% % plot
% thresh_triangle = quantile(Cscore_cc.triangle_shuffle(:), 0.95);
% thresh_square = quantile(Cscore_cc.square_shuffle(:), 0.95);
% thresh_hex = quantile(Cscore_cc.hex_shuffle(:), 0.95);
% figure
% subplot(2,3,1)
% histogram(Cscore_cc.triangle, [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_triangle,thresh_triangle],[0,80],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square
% subplot(2,3,4)
% histogram(Cscore_cc.triangle_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_triangle,thresh_triangle],[0,50000],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square
% subplot(2,3,2)
% histogram(Cscore_cc.square, [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_square,thresh_square],[0,80],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square
% subplot(2,3,5)
% histogram(Cscore_cc.square_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_square,thresh_square],[0,50000],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square
% subplot(2,3,3)
% histogram(Cscore_cc.hex, [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_hex,thresh_hex],[0,80],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square
% subplot(2,3,6)
% histogram(Cscore_cc.hex_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
% line([thresh_hex,thresh_hex],[0,50000],'linewidth',1.5,'color','r')
% xlim([-1.5,1.5])
% xlabel('corner score')
% ylabel('counts')
% axis square

%% Compare the shuffled distribution with and without penalty enabled
load('F:\analysis_folders.mat','expA')
datapath = expA;
% get corner cells at different threshold level
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    %shuffle with penalty
    Cshuffwp = identify_corner_cell;
    save('corner_metrics.mat','Cshuffwp','-append')
end
session_name = {'triangle','square','hex'}; %expA
allcscore = cell(length(datapath), length(session_name));
allcscore_shuffle = cell(length(datapath), length(session_name)); 
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat', 'Cshuffwp')
    C = Cshuffwp;
    cscorep = cell(1,length(C.cmetrics));
    for jj = 1:length(C.cmetrics)
        cscorep{jj} = cellfun(@sum, C.cmetrics{jj}.cornerscorep);
    end
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        cornercell = C.cornercell(idx);
        %cornerscore without penalize
        cscore = cellfun(@(x,y) x+y', C.cscore(idx),cscorep(idx), 'uni', 0);
        cscore_shuffle = C.cscore_shuffle(idx);
        allcscore{n,ii} = vertcat(cscore{:});
        allcscore_shuffle{n,ii} = vertcat(cscore_shuffle{:});
    end
end
% rearrange the data
Cscore_shuffwp = struct;
score_triangle = allcscore(:,1);
Cscore_shuffwp.triangle = vertcat(score_triangle{:});
score_square = allcscore(:,2);
Cscore_shuffwp.square= vertcat(score_square{:});
score_hex = allcscore(:,3);
Cscore_shuffwp.hex= vertcat(score_hex{:});
% reformat shuffled data
score_triangle = allcscore_shuffle(:,1);
Cscore_shuffwp.triangle_shuffle = vertcat(score_triangle{:});
score_square = allcscore_shuffle(:,2);
Cscore_shuffwp.square_shuffle= vertcat(score_square{:});
score_hex = allcscore_shuffle(:,3);
Cscore_shuffwp.hex_shuffle= vertcat(score_hex{:});
save('F:\Results_revision\corner_metrics_randomlevel.mat','Cscore_shuffwp','-append')

% Final plot
load('F:\Results_revision\corner_metrics_randomlevel.mat', 'Cscore_all', 'Cscore_shuffwp')
load('F:\Results_revision\corner_metrics_CA1.mat', 'cscore_CA1')
thresh_triangle = quantile(Cscore_all.triangle_shuffle(:), 0.95);
thresh_square = quantile(Cscore_all.square_shuffle(:), 0.95);
thresh_hex = quantile(Cscore_all.hex_shuffle(:), 0.95);
figure
subplot(5,3,1)
histogram(Cscore_shuffwp.triangle, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_triangle,thresh_triangle],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/o penalty)')
ylabel('counts')
axis square
title('triangle')
subplot(5,3,2)
histogram(Cscore_shuffwp.square, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_square,thresh_square],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/o penalty)')
ylabel('counts')
axis square
title('square')
subplot(5,3,3)
histogram(Cscore_shuffwp.hex, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_hex,thresh_hex],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/o penalty)')
ylabel('counts')
axis square
title('hexagon')
subplot(5,3,4)
histogram(Cscore_all.triangle, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_triangle,thresh_triangle],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/ penalty)')
ylabel('counts')
axis square
title('triangle')
subplot(5,3,7)
histogram(Cscore_all.triangle_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_triangle,thresh_triangle],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/o penalty)')
ylabel('counts')
axis square
subplot(5,3,5)
histogram(Cscore_all.square, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_square,thresh_square],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/ penalty)')
ylabel('counts')
axis square
title('square')
subplot(5,3,8)
histogram(Cscore_all.square_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_square,thresh_square],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/o penalty)')
ylabel('counts')
axis square
subplot(5,3,6)
histogram(Cscore_all.hex, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_hex,thresh_hex],[0,550],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,550])
xlabel('corner score for all the cell(w/ penalty)')
ylabel('counts')
axis square
title('hexagon')
subplot(5,3,9)
histogram(Cscore_all.hex_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_hex,thresh_hex],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/o penalty)')
ylabel('counts')
axis square
subplot(5,3,10)
histogram(Cscore_shuffwp.triangle_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_triangle,thresh_triangle],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/ penalty)')
ylabel('counts')
axis square
subplot(5,3,11)
histogram(Cscore_shuffwp.square_shuffle, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_square,thresh_square],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/ penalty)')
ylabel('counts')
axis square
subplot(5,3,12)
histogram(Cscore_shuffwp.hex_shuffle(:), [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_hex,thresh_hex],[0,800000],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
xlabel('shuffled corner score (w/ penalty)')
ylabel('counts')
axis square
subplot(5,3,14)
histogram(cscore_CA1, [-1.5:0.05:1.5], 'FaceColor', 'k')
line([thresh_square,thresh_square],[0,900],'linewidth',1.5,'color','r')
xlim([-1.5,1.5])
ylim([0,900])
xlabel('corner score, CA1')
ylabel('counts')
axis square

%% CORNER CELL STABILITY
%% Determine an appropriate stability by shuffling
load('F:\analysis_folders.mat','expA')
datapath = expA;
%determin corner cells in each session using simulated corner locations
%map_stb: using interleaved method
%map_stability: using split to half method
for n = 1:length(datapath)
    cd(datapath{n})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat','thresh')
    map_stb_shuffle = cell(1,length(neuronIndividualsf));
    parfor ii = 1:length(neuronIndividualsf)
        neuron = neuronIndividualsf{1,ii};
        behav = behavIndividualsf{1,ii};
        Rxy_shuffle = calc_spatialmap_stability_shuffle(neuron,behav,thresh);
        map_stb_shuffle{1,ii} = Rxy_shuffle; %defined by bit/spike
    end
    save('spatial_metrics.mat','map_stb_shuffle','-append')
end
session_name = {'hex'}; %expA
stability_thresh = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat','map_stb_shuffle')
    load('corner_metrics.mat','Cthresh');
    load('env_geometry.mat','S');
    idx = ismember(S,session_name);
    C = Cthresh{3};
    CC_stability_shuffle = cellfun(@(x,y) x(y,:), map_stb_shuffle(idx), C.cornercell(idx), 'uni', 0);
    CC_stability_shuffle = mean(cellfun(@(x) nanmean(quantile(x',0.95)), CC_stability_shuffle));
%     CC_stability_shuffle = mean(cellfun(@(x) nanmean(quantile(x',0.95)), map_stability_shuffle));
    stability_thresh = [stability_thresh;CC_stability_shuffle];
end
% the stability threshold for shuffling from hex for corner cells is 0.4. 

%% DETERMINE RANDOM LEVEL OF CORNER CELLS (simulated corners or CA1 data)
%% Using randomly assigned corners
%mannually assign the location of corners for each session
load('F:\analysis_folders.mat','expA')
datapath = expA;
for ii = 1:length(datapath)
    cd(datapath{ii})
    [~,~,env_coorS,env_coorSi] = identify_env_geometry;
    save('env_geometry.mat','env_coorS','env_coorSi','-append')
end
%determin corner cells in each session using simulated corner locations
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('env_geometry.mat','env_coorS')
    mask = repmat({true},1,11);
    Csim_msk = identify_corner_cell([], [], [], env_coorS);
    save('corner_metrics.mat','Csim_msk','-append')
end

% proportion of corner cells
session_name = {'circle','triangle','square','hex'}; %expA
propCC_randC = cell(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metrics.mat','Csim_msk')
    C = Csim_msk;
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
save('F:\Results_revision\corner_metrics_randomlevel.mat','propCC_randC','propCC_randC_ea','-append')

%% DETERMINE THE PROPORTION OF CORNER CELLS WERE FILTERED BY THE PAIR-WISE DISTANCE
%% How effective of the field distance constraint for corner cells ***
load('F:\analysis_folders.mat','expA')
datapath = expA;
session_name = {'triangle','square','hex'}; %expA
propfilt = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metricsR1.mat','C')
    idx = ismember(S, session_name);
    bef = cellfun(@numel, C.cornercell_raw(idx));
    aft = cellfun(@numel, C.cornercell(idx));
    propfilt = [propfilt; (bef-aft)/bef];
end
% wrote in response to reviewer letter

%% CA1 CA1 CA1
%% CA1 CA1
%% Using CA1 data
load('F:\CA1_mice\analyze_info.mat', 'CA1data')
datapath_ori = CA1data;
for n = 1:length(datapath_ori)
    cd(datapath_ori{n})
    load('neuronIndividualsf.mat');
    load('behavIndividualsf.mat');
    load('thresh.mat','thresh');
    load('property.mat','mouse');
    load('sessions.mat','sessions');
%     if strcmp(mouse, 'M4056')
%         idx = [5,7,9];
%     elseif strcmp(mouse, 'M4071')
%         idx = [4,6,8];
%     else
%         idx = [3,5,7];
%     end
    if strcmp(mouse, 'M3974F')
        idx = [6,8,10];
    else
        idx = [5,7,9];
    end
    neuronIndividualsf = neuronIndividualsf(idx);
    behavIndividualsf = behavIndividualsf(idx);
    sessions = sessions(idx);
    sessions = cellfun(@(x) [mouse,'_',x,'_square'], sessions, 'uni', 0);
    % correct behav position
    for ii = 1:length(behavIndividualsf)
        behav0 = behavIndividualsf{ii};
        edgecorrect = min(behav0.position);
        behav0.position = (behav0.position - edgecorrect)*1.25;
        behavIndividualsf{ii} = behav0;
    end
    %change to analysis direction
    cd('F:\CA1_mice')
    if ~exist(mouse, 'dir')
        mkdir(mouse)
    end
    cd(['F:\CA1_mice\',mouse])
    save('neuronIndividualsf.mat','neuronIndividualsf','-v7.3');
    save('behavIndividualsf.mat','behavIndividualsf','-v7.3');
    save('thresh.mat','thresh');
    save('sessions.mat','sessions');
    %generate ratemap
    generate_spatial_ratemap([], false)
end

%generate corner coordinates
load('F:\CA1_mice\analyze_info.mat','expCA1')
datapath = expCA1;
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('firingrateAll.mat','ratemap')
    ratemap_size = cellfun(@(x) x{1}, ratemap, 'uni', 0);
    env_coor = cellfun(@(x) [1,1;size(x,2)-1,1;size(x,2)-1,size(x,1)-1;1,size(x,1)-1;size(x,1)/2,size(x,2)/2],...
        ratemap_size, 'uni', 0);
    save('env_geometry.mat','env_coor')
end

% get the number of corner cells
load('F:\CA1_mice\analyze_info.mat','expCA1')
datapath = expCA1;
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    C = identify_corner_cell;
    save('corner_metrics.mat','C','-v7.3')
end

% plot corner cells for a given session
load('firingrateAll.mat','ratemap')
load('corner_metricsR1.mat','C')
session2plot = 3;
plot_ratemap(ratemap{session2plot},C.cornercell{session2plot})
% plot corner score histogram of corner cells

%% Incorporate within sessions stability for corner cells
load('F:\analysis_folders.mat', 'expCA1')
datapath = expCA1;
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat')
    load('behavIndividualsf.mat')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndividualsf, behavIndividualsf, 'uni',0);
    save('spatial_metrics.mat', 'map_stb','-append')
end
%% Final Proportion of corner cells for CA1 data
% final number of corner cells for CA1 data
load('F:\analysis_folders.mat', 'expCA1')
datapath = expCA1;
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'C')
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C', '-v7.3')    
end

% distribution of corner cell stability
CC_stabilityCA1 = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'C')
    CCS = cellfun(@(x,y) x(y), map_stb, C.cornercell, 'uni', 0);
    CC_stabilityCA1 = [CC_stabilityCA1; CCS{end}];
end
save('F:\Results_revision\corner_metrics_CA1.mat','CC_stabilityCA1','-append')
%plot
load('F:\Results_experimentA\stabilitycc_R1.mat','CC_stability')
load('F:\Results_revision\corner_metrics_CA1.mat','CC_stabilityCA1')
figure
subplot(1,2,1)
histogram(CC_stability, [-0.5:0.05:1])
line([0.3,0.3],[0,120],'linewidth',1.5,'color','r')
xlabel('corner cell stability')
ylabel('counts')
axis square
subplot(1,2,2)
histogram(CC_stabilityCA1, [-0.5:0.05:1])
line([0.3,0.3],[0,10],'linewidth',1.5,'color','r')
xlabel('corner cell stability, CA1')
ylabel('counts')
axis square

% organize the data for final corner cells
cscore_CA1 = [];
numCC_CA1 = [];
total_CA1 = [];
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('corner_metricsR1.mat','C');
    load('thresh.mat','thresh');
    cscore_CA1 = [cscore_CA1; C.cscore{end}];
    numCC_CA1 = [numCC_CA1;cellfun(@numel, C.cornercell(end))];
    total_CA1 = [total_CA1;length(thresh)];
end
propCC_CA1 = numCC_CA1./total_CA1;
save('F:\Results_revision\corner_metrics_CA1.mat', 'propCC_CA1', 'cscore_CA1', 'numCC_CA1', 'total_CA1','-append')
%plot in GraphPad






