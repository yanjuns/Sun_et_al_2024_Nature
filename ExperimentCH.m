%% Experiment C
% Experiment C aims to identify neurons that might encode convex corners
load('F:\analysis_folders.mat','expC')
datapath = expC;

%% Identify corner cells
%mannually identify the location of corners for each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [S,N,ccav_coor,ccav_coori] = identify_env_geometry;
    save('env_geometry.mat','ccav_coor','ccav_coori','-append')
end

%determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
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
    C = identify_corner_cell(mask,[],0.4);
    save('corner_metrics.mat','C','-append')
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
    load('corner_metrics.mat', 'Cthresh')
    C = Cthresh{4};
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C')    
end
%% determine corner cells using multiple sessions
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    session_name = {'largeX','largeXm1'}; %expA
    ccellx = cell(1,length(session_name));
    %for experiment C
    for jj = 1:length(session_name)
        idx = strcmp(S, session_name{jj});
        C_square = C.cornercell(idx);
        ccellx{jj} = unique(cat(1, C_square{:}));
    end
    cornercellx = 1:median(cellfun(@numel,C.cscore));
    for jj = 1:length(ccellx)
        cornercellx = intersect(cornercellx,ccellx{jj});
    end
    save('corner_metricsR1.mat','cornercellx','-append')
end
%% Percentage of corner cells for each session
session_name = {'largesq','largeRT','largeX','largeXm1','largeXm2','largeXm3'}; %expC
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
        pcc{ii} = C_ea(:)/numNeurons;
    end
    prop_cornercell(n,:) = pcc;
end
prop_cornercell_ea = cellfun(@mean, prop_cornercell);
save('F:\Results_experimentC\prop_cornercell_R1.mat','prop_cornercell_ea')

%% Percentage of overlap between concave and convex corner cells
session_name = {'largesq','largeRT','largeX','largeXm1'};
prop_overlap = NaN(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    idxb = strcmp(S,'largesq');
    bsl = C.cornercell(idxb);
    bsl = bsl{1};
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        C_ea = C.cornercell(idx);
        C_ea = C_ea{:};
        idxc = ismember(C_ea, bsl);
        p_overlap = sum(idxc)/numel(bsl);
        prop_overlap(n,ii) = p_overlap;
    end
end
save('F:\Results_experimentC\prop_overlap_R1.mat','prop_overlap')

% for generating venn plot
venn = struct;
venn.a = []; venn.b = []; venn.c = [];
venn.ab = []; venn.ac = []; venn.bc =[];
venn.abc = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    idxa = strcmp(S,'largesq');
    a = cell2mat(C.cornercell(idxa));
    idxb = strcmp(S,'largeRT');
    b = cell2mat(C.cornercell(idxb));
    idxc = strcmp(S,'largeX');
    c = cell2mat(C.cornercell(idxc));
    venn.a = [venn.a; numel(a)/numel(a)];
    venn.b = [venn.b; numel(b)/numel(a)];
    venn.c = [venn.c; numel(c)/numel(a)];
    
    venn.ab = [venn.ab; numel(intersect(a,b))/numel(a)];
    venn.ac = [venn.ac; numel(intersect(a,c))/numel(a)];
    venn.bc = [venn.bc; numel(intersect(b,c))/numel(a)];
    venn.abc = [venn.abc; numel(intersect(intersect(a,b),intersect(a,c)))/numel(a)];
end
save('F:\Results_experimentC\prop_overlap_R1.mat','venn','-append')

% determine the random level of overlap
rand_overlap = struct;
rand_overlap.squarerect = [];
rand_overlap.squarecvx = [];
rand_overlap.rectcvx = [];
nboot = 1000;
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    load('thresh.mat','thresh')
    numNeuron = length(thresh);
    idxb = strcmp(S,'largesq');
    bsl = C.cornercell(idxb);
    bsl = numel(bsl{1});
    idxc = strcmp(S, 'largeRT');
    cvx = C.cornercell(idxc);
    cvx = numel(cvx{1});
    rand_overlap_ea = NaN(1,nboot);
    for ii = 1:nboot
        p1 = randperm(numNeuron,bsl);
        p2 = randperm(numNeuron,cvx);
        po = intersect(p1,p2);
        rand_overlap_ea(ii) = numel(po)/numel(p1);
    end
    rand_overlap.squarerect = [rand_overlap.squarerect; rand_overlap_ea];
end
save('F:\Results_experimentC\prop_overlap_R1.mat','rand_overlap','-append')
figure
subplot(1,3,1)
histogram(rand_overlap.rectcvx(:),[0:0.05:0.5])
xlim([0,0.4])
ylim([0,9000])
axis square
xlabel('% overlap')
t3 = quantile(rand_overlap.squarerect(:),0.95);
line([t3,t3],[0,9000],'linewidth',1.5,'color','k')
line([0.311,0.311],[0,9000],'linewidth',1.5,'color','r')
subplot(1,3,2)
histogram(rand_overlap.squarecvx(:),[0:0.05:0.5])
hold on
xlim([0,0.4])
ylim([0,9000])
xlabel('% overlap')
axis square
t1 = quantile(rand_overlap.squarecvx(:),0.95);
line([t1,t1],[0,9000],'linewidth',1.5,'color','k')
line([0.0325,0.0325],[0,9000],'linewidth',1.5,'color','r')
subplot(1,3,3)
histogram(rand_overlap.rectcvx(:),[0:0.05:0.5])
xlim([0,0.4])
ylim([0,9000])
axis square
xlabel('% overlap')
t2 = quantile(rand_overlap.rectcvx(:),0.95);
line([t2,t2],[0,9000],'linewidth',1.5,'color','k')
line([0,0],[0,9000],'linewidth',1.5,'color','r')

%% CONCAVE AND CONVEX CORNER CELLS ARE ANATOMICALLY SEPARATED
%% Plot an example in anatomy
load('neuronIndividualsf.mat');
load('corner_metricsR1.mat','C')
neuronSelected1 = C.cornercell{1};
neuronSelected2 = C.cornercell{3};
ncontour = neuronIndividualsf{1}.Coor;
figure
hold on
for ii = 1:length(ncontour)
    plot(ncontour{ii}(1,:),ncontour{ii}(2,:),'Color',[128 130 133]/255,'LineWidth', 0.5)
    set(gca, 'YDir','reverse')
    axis image
end
for jj = 1:length(neuronSelected1)
    plot(ncontour{neuronSelected1(jj)}(1,:),ncontour{neuronSelected1(jj)}(2,:),'Color',[43 120 142]/255, 'LineWidth', 1)
end
for jj = 1:length(neuronSelected2)
    plot(ncontour{neuronSelected2(jj)}(1,:),ncontour{neuronSelected2(jj)}(2,:),'Color',[147 38 103]/255, 'LineWidth', 1)
end
%% intra- and inter-cluster distance of concave and convex cells
% measure pair-wise distance of intra- and inter-cluster of corner cells. 
% 1px = 0.89um (12.5mm Ach lens miniscope)
load('F:\analysis_folders.mat','expC')
datapath = expC;

pwdist = struct;
pwdist.iaccav = []; pwdist.iacvex = []; pwdist.inter = []; 
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat')
    load('corner_metricsR1.mat','C')
    ctr = struct;
    ctr.ccav = neuronIndividualsf{1}.centroid(union(C.cornercell{1},C.cornercell{2}),:);
    ctr.cvex = neuronIndividualsf{1}.centroid(union(C.cornercell{3},C.cornercell{4}),:);
    pwdist.iaccav = [pwdist.iaccav; mean(pdist(ctr.ccav))];
    pwdist.iacvex = [pwdist.iacvex; mean(pdist(ctr.cvex))];
    pwdist.inter = [pwdist.inter; mean(mean(pdist2(ctr.ccav,ctr.cvex)))];
end
save('F:\Results_experimentC\anatdist_ccav_cvex_R1.mat','pwdist')
%plot
x1 = pwdist.iaccav .* 0.89;
x2 = pwdist.iacvex .* 0.89;
y = pwdist.inter .* 0.89;

figure;
plot(x1,y,'.','markersize',40)
hold on
plot(x2,y,'.','markersize',40)
axis square
legend('concave corner cell', 'convex corner cell')
xlim([40,120])
ylim([40,120])
xlabel('intra-group distance (um)')
ylabel('inter-group distance (um)')

%% FIRING RATE
%% compute the the firing rate as a function of distance to corners
load('F:\analysis_folders.mat','expC')
datapath = expC;
session_name = {'largesq','largeRT','largeX','largeXm1','largeXm2','largeXm3'}; %expC
%get the distance to corners for each bin and associated firing rate
for n = 1:length(datapath)
    cd(datapath{n})
    [d2corner,d2cornerb] = corner_dist_firingrate;
    save('corner_metrics.mat','d2corner','d2cornerb','-append')
end
%average data for each type of sessions
for n = 1:length(datapath)
    cd(datapath{n})
    d2cornerq = corner_dist_firingrate_quantify(session_name);
    save('corner_metricsR1.mat','d2cornerq','-append')
end
%% plot for convex corners
%plot dist2corner firing rate for convex corners, ExpC
data = cell(1,length(session_name));
for ii = 1:length(session_name)
    temp = struct; temp.x=[];temp.ycc=[];temp.yncc=[];temp.yccx=[];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('corner_metricsR1.mat','d2cornerq')
        temp.x = padconcatenation(temp.x, nanmean(d2cornerq.d{ii},1),1);
        temp.ycc = padconcatenation(temp.ycc, nanmean(d2cornerq.frcc{ii},1),1);
        temp.yncc = padconcatenation(temp.yncc, nanmean(d2cornerq.frncc{ii},1),1);
        temp.yccx = padconcatenation(temp.yccx, nanmean(d2cornerq.frccx{ii},1),1);
    end
    data{ii} = temp;
end
[h,p] = signrank(mean(data{1}.ycc(:,1:3),2),mean(data{1}.yncc(:,1:3),2))
[h,p] = signrank(nanmean(data{1}.ycc(:,end-3:end-1),2),nanmean(data{1}.yncc(:,end-3:end-1),2))
[h,p] = signrank(mean(data{2}.ycc(:,1:3),2),mean(data{2}.yncc(:,1:3),2))
[h,p] = signrank(nanmean(data{2}.ycc(:,end-3:end-1),2),nanmean(data{2}.yncc(:,end-3:end-1),2))

%plot
figure
subplot(3,2,1)
stdshade(smoothdata(data{1}.yccx,2,'movmean',3),0.4, [147 38 103]/255, [1:length(data{1}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{1}.yncc,2,'movmean',3),0.4, [109 110 113]/255, [1:length(data{1}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.8])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex in square')
subplot(3,2,2)
stdshade(smoothdata(data{2}.yccx,2,'movmean',3),0.4, [147 38 103]/255,[1:length(data{2}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{2}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:length(data{2}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.75])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex in rectangle')
subplot(3,2,3)
stdshade(smoothdata(data{3}.ycc,2,'movmean',3),0.4, [147 38 103]/255,[1:length(data{3}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{3}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:length(data{3}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.8])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex')
subplot(3,2,4)
stdshade(smoothdata(data{4}.ycc,2,'movmean',3),0.4, [147 38 103]/255,[1:length(data{4}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{4}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:length(data{4}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.8])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex m1')
subplot(3,2,5)
stdshade(smoothdata(data{5}.ycc,2,'movmean',3),0.4, [147 38 103]/255,[1:length(data{5}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{5}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:length(data{5}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.8])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex m2')
subplot(3,2,6)
stdshade(smoothdata(data{6}.ycc,2,'movmean',3),0.4, [147 38 103]/255,[1:length(data{6}.ycc)]*1.55)
hold on
stdshade(smoothdata(data{6}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:length(data{6}.yncc)]*1.55)
axis square
xlim([0,24])
ylim([0, 0.8])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('convex m3')

%% plot for concave corners
% %plot dist2corner firing rate, ExpC
% data = cell(1,length(session_name));
% for ii = 1:length(session_name)
%     temp = struct; temp.x=[];temp.ycc=[];temp.yncc=[];temp.yccx=[];
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('corner_metrics.mat','d2cornerq_ccav')
%         d2cornerq = d2cornerq_ccav;
%         temp.x = padconcatenation(temp.x, nanmean(d2cornerq.d{ii},1),1);
%         temp.ycc = padconcatenation(temp.ycc, nanmean(d2cornerq.frcc{ii},1),1);
%         temp.yncc = padconcatenation(temp.yncc, nanmean(d2cornerq.frncc{ii},1),1);
%         temp.yccx = padconcatenation(temp.yccx, nanmean(d2cornerq.frccx{ii},1),1);
%     end
%     data{ii} = temp;
% end
% %plot
% figure
% subplot(2,2,1)
% stdshade(smoothdata(data{3}.ycc,2,'movmean',3),0.4, [28 117 188]/255)
% hold on
% stdshade(smoothdata(data{3}.yncc,2,'movmean',3),0.4, [109 110 113]/255)
% axis square
% % xlim([0,15])
% ylim([0, 0.75])
% xlabel('distance to corners')
% ylabel('firing rate (spikes/sec)')
% title('convex')
% subplot(2,2,2)
% stdshade(smoothdata(data{4}.ycc,2,'movmean',3),0.4, [28 117 188]/255)
% hold on
% stdshade(smoothdata(data{4}.yncc,2,'movmean',3),0.4, [109 110 113]/255)
% axis square
% % xlim([0,15])
% ylim([0, 0.75])
% xlabel('distance to corners')
% ylabel('firing rate (spikes/sec)')
% title('convex m1')
% subplot(2,2,3)
% stdshade(smoothdata(data{5}.ycc,2,'movmean',3),0.4, [28 117 188]/255)
% hold on
% stdshade(smoothdata(data{5}.yncc,2,'movmean',3),0.4, [109 110 113]/255)
% axis square
% % xlim([0,15])
% ylim([0, 0.75])
% xlabel('distance to corners')
% ylabel('firing rate (spikes/sec)')
% title('convex m2')
% subplot(2,2,4)
% stdshade(smoothdata(data{6}.ycc,2,'movmean',3),0.4, [28 117 188]/255)
% hold on
% stdshade(smoothdata(data{6}.yncc,2,'movmean',3),0.4, [109 110 113]/255)
% axis square
% % xlim([0,15])
% ylim([0, 0.75])
% xlabel('distance to corners')
% ylabel('firing rate (spikes/sec)')
% title('convex m3')
%% Non-geometric impact on convex corner coding
%% compare largeX and largeXm1 at the corner with different decorations ***
load('F:\analysis_folders.mat','expC')
datapath = expC;
session_name = {'largesq','largeRT','largeX','largeXm1','largeXm2','largeXm3'}; %expC
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCors = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
        load('corner_metricsR1.mat','C')
        idx = strcmp(S, session_name{k});
        idx2 = strcmp(S, 'largeX');
        pkfr_corners = pkfr_corners(idx);
        cornercell = C.cornercell(idx2);
        allcell = 1:length(pkfr_corners{1,1});
        ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
        ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        %         pkfrcc_atcorner = cellfun(@(x,y) x(:,cornercellx), pkfr_corners, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        %non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,pkfrcc_atcorner./pkfrncc_atcorner];
        %simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,pkfrcc_atcorner./pkfrsim_atcorner];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
end
%corner cells were defined in largeX session and compared between largeX
%and largeXm1. 
pkfrcc_atCs_m1 = pkfrcc_atCs;
pkfrcc_atCnc_m1 = pkfrcc_atCnc;
save('F:\Results_experimentC\pkfrcc_atCorner_R1.mat','pkfrcc_atCnc_m1','pkfrcc_atCs_m1','-append')

%% 270 vs. 315 degree
%% see the firing rate difference between the 270 degree and 315 degree corner cell ***
load('F:\analysis_folders.mat','expC')
datapath = expC;
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
% reidentify corner cells in the largeXm3 session using silghtly larger
% mask, as the triangluar geometry of the corner configuration. 
for ii = 1:length(datapath)
    cd(datapath{ii})
    mask = cell(1,1);
    mask{1} = true;
    %change inside the function to specify which session to run
    Cvx = identify_corner_cell(mask,[],0.4);
    save('corner_metrics.mat','Cvx','-append')
end
%gather data for plotting
session_name = {'largeXm3'}; %expC
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCors = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
        if strcmp(datapath{n},'F:\subiculum_mice_May2022\M4107c')
            load('corner_metrics.mat','C2')
            cornercell = C2.cornercell(6);
        else
            load('corner_metrics.mat','Cvx')
            cornercell = Cvx.cornercell; %using C and Cvx are the same
        end
        idx = strcmp(S, session_name{k});
        pkfr_corners = pkfr_corners(idx);
        allcell = 1:length(pkfr_corners{1,1});
        ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
        ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        %non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,pkfrcc_atcorner./pkfrncc_atcorner];
        %simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,pkfrcc_atcorner./pkfrsim_atcorner];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
end
pkfrcc_atCnc_cvx = pkfrcc_atCnc;
pkfrcc_atCs_cvx = pkfrcc_atCs;
save('F:\Results_experimentC\pkfrcc_atCorner_R1.mat','pkfrcc_atCnc_cvx','pkfrcc_atCs_cvx','-append')
%plot in GraphPad

%% ALLOCENTRIC VS. EGOCENTRIC
%% EGOCENTRIC CHARACTERIZATION WITH CORNER BEARING
%% Using tuning curve method to identify egocentric tuning
load('F:\analysis_folders.mat','expC')
datapath = expC;
%% get egocentric corner directions. 
for n = 1:length(datapath)
    cd(datapath{n})
    get_behav_directions
end
%% plot egocentric bearing of neurons
load('neuronIndividualsf.mat');
load('thresh.mat');
load('corner_metrics.mat','C2')
load('lnp_whole_fminunccdist.mat', 'lnp_model')
cell2plot = find(lnp_model{1} == 8);
plot_eb_map_longitudinal(neuronIndividualsf,thresh,cell2plot)

%print eb map to eps file
filepath = pwd;
neuronSelected = 145;
fig = openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected), '_ebmap.fig']));
print (fig, '-painters', '-depsc', fullfile(filepath,'RatemapFigures',['Cell',num2str(neuronSelected),'_ebmap.eps']));

%% Using LN model to assess how much egocentric tuning contribute to corner coding
load('F:\analysis_folders.mat','expC')
datapath = expC;
%run LN model
edit run_me
% % OPTIONAL, change the metric for defining the ln model
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminuncpbin4_2.mat','fit_results')
    lnp_model_llh = cell(size(fit_results,1),size(fit_results,2));
    for x = 1:size(fit_results,1)
        for y = 1:size(fit_results,2)
            lnp_model_ea = NaN(length(fit_results{x,y}),1);
            for ii = 1:length(fit_results{x,y})
                [lnp_model_ea(ii),~] = select_best_model(fit_results{x,y}{ii},'llh');
            end
            lnp_model_llh{x,y} = lnp_model_ea;
        end
    end
    save('lnp_whole_fminuncpbin4_2.mat', 'lnp_model_llh','-append')
end

%% gather LN data for plotting
%combine all cells together
ln_ccav = [];
ln_cvex = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminuncpbin4R1.mat','lnp_model')
    load('corner_metricsR1.mat','C')
    idx = C.cornercell(1:4);
    lnp_model = cellfun(@(x,y) x(y), lnp_model, idx, 'UniformOutput',false);
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    ccav = lnp_model{1};
    ln_ccav = [ln_ccav;ccav];
    cvex = lnp_model{3};
    ln_cvex = [ln_cvex;cvex];
end
save('F:\Results_experimentC\LN_resultsR1_cornercell.mat','ln_ccav','ln_cvex')

% for each individual mice
LN = struct;
LNP = {};
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminuncpbin4R1.mat','lnp_model')
    load('corner_metricsR1.mat','C')
    idx = C.cornercell(1:4);
    lnp_model = cellfun(@(x,y) x(y), lnp_model, idx, 'UniformOutput',false);
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    histall = cellfun(@(x) histcounts(x,[-0.5:1:15.5])./length(x),...
        lnp_model,'uni',0);
    for ii = 1:size(histall,2)
        LNP{ii}(n,:) = histall{ii};
    end
end
LN.ccav = LNP{1};
LN.cvex = LNP{3};
save('F:\Results_experimentC\LN_resultsR1_cornercell.mat','LN','-append')

%plot
edge = [-0.5:1:15.5];
figure
subplot(2,2,1)
histogram(ln_ccav,edge)
xlim([-0.5,15.5])
xticks([0:15])
set(gca,'XTickLabel',{'NaN','PHSE','PHS','PHE','PSE','HSE','PH','PS','PE','HS','HE','SE','P','H','S','E'});
set(gca, 'XDir','reverse')
ylim([0,40])
ylabel('num of total corner cell')
subplot(2,2,2)
histogram(ln_cvex,edge,'FaceColor',[147 38 103]/255)
xlim([-0.5,15.5])
xticks([0:15])
set(gca,'XTickLabel',{'NaN','PHSE','PHS','PHE','PSE','HSE','PH','PS','PE','HS','HE','SE','P','H','S','E'});
set(gca, 'XDir','reverse')
ylim([0,40])
ylabel('num of total corner cell')
subplot(2,2,3)
position_O = 0:1:15;
boxplot(LN.ccav,'colors',[28 117 188]/255,'positions',position_O,'Symbol','+');
ylim([-0.1 1])
set(gca,'XTickLabel',{'NaN','PHSE','PHS','PHE','PSE','HSE','PH','PS','PE','HS','HE','SE','P','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')
subplot(2,2,4)
boxplot(LN.cvex,'colors',[147 38 103]/255,'positions',position_O,'Symbol','+');
ylim([-0.1 1])
set(gca,'XTickLabel',{'NaN','PHSE','PHS','PHE','PSE','HSE','PH','PS','PE','HS','HE','SE','P','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')

%% BOUNDARY CELL ANALYSIS
%% Border Cell
%% Identify border cells
load('F:\analysis_folders.mat','expC')
datapath = expC;
for n = 1:length(datapath)
    cd(datapath{n})
    load('firingrateAll.mat','ratemap')
    bscore = {}; bfields = {};
    for ii = 1:length(ratemap)
        ratemapS = ratemap{1,ii};
        %remove rows that are all NaNs
        ratemapS = cellfun(@(x) x(~all(isnan(x),2),:), ratemapS, 'uni', 0);
        %remove cols that are all NaNs
        ratemapS = cellfun(@(x) x(:,~all(isnan(x),1)), ratemapS, 'uni', 0);
        %trim the edge
        ratemapS = cellfun(@(x) x(2:end-1,2:end-1), ratemapS, 'uni', 0);
        [bs,fields] = cellfun(@borderScoreCalculation, ratemapS);
        bscore{ii} = bs';
        bfields{ii} = fields;
    end
    save('border_metrics.mat','bscore','bfields','-v7.3')
end

%% Defining border cells with various threshold
coverage_thresh = [0.5:0.1:0.9];
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat','placecell')
    load('border_metrics.mat','bscore','bfields')
    bc1 = cellfun(@(x) find(x > 0.6), bscore, 'uni',0);
    coverage = cellfun(@(x) {x.coverage}, bfields, 'uni',0);
    idxc = cell(1,length(coverage));
    for jj = 1:length(coverage_thresh)
        for ii = 1:length(coverage)
            idxc{jj,ii} = cellfun(@(x) x.perimeter.max>coverage_thresh(jj), coverage{ii});
        end
    end
    bc2 = cellfun(@(x) find(x == 1), idxc, 'uni', 0);
    bc1 = repmat(bc1, length(coverage_thresh), 1);
    bcell_raw = cellfun(@(x,y) intersect(x,y), bc1, bc2, 'uni', 0);
    placecell = repmat(placecell, length(coverage_thresh), 1);
    bcell_all = cellfun(@(x,y) intersect(x,y), bcell_raw, placecell, 'uni',0);
    save('border_metricsR1.mat','bcell_all')
end

%Add stability criteria to border cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('border_metricsR1.mat', 'bcell_all')
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    stablecell = repmat(stablecell, size(bcell_all, 1),1);
    bcell_stb = cellfun(@(x,y) intersect(x,y), stablecell, bcell_all, 'uni', 0);
    save('border_metricsR1.mat', 'bcell_stb','-append')    
end

%overall proportion of border cellls
prop_bcell = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('border_metricsR1.mat','bcell_stb')
    bcell = bcell_stb;
    load('thresh.mat','thresh')
    numb = cellfun(@numel, bcell(:,[1,2]))./length(thresh);
    prop_bcell = cat(3,prop_bcell,numb);
end
save('F:\Results_experimentC\prop_bcell_R1.mat','prop_bcell')

% proportion of BVCs
load('F:\Results_experimentC\prop_bcell_R1.mat','prop_bcell')
data_sq = squeeze(prop_bcell(:,1,:))';
data_RT = squeeze(prop_bcell(:,2,:))';
figure
subplot(1,2,1)
stdshade(data_sq, 0.3, [43 57 144]/255)
axis square
ylim([0,0.15])
xticks([1:5])
xticklabels({'0.5', '0.6', '0.7', '0.8', '0.9'})
xlabel('boundary coverage threshold')
ylabel('prop of total neuron')
title('BVCs in square')
subplot(1,2,2)
stdshade(data_RT, 0.3, [43 57 144]/255)
axis square
ylim([0,0.15])
xticks([1:5])
xticklabels({'0.5', '0.6', '0.7', '0.8', '0.9'})
xlabel('boundary coverage threshold')
ylabel('prop of total neuron')
title('BVCs in rectangle')

% plot the proportion of BVCs in GraphPad Prism. 

%% percentage overlap between border and corner cells ***
load('F:\analysis_folders.mat','expC')
datapath = expC;

overlap = struct;
overlap.b = [];
overlap.ccav = [];
overlap.cvex = [];
bc_overlap = cell(1,5);
for n = 1:length(bc_overlap)
    bc_overlap{n} = overlap;
end
for n = 1:length(datapath)
    cd(datapath{n})
    load('thresh.mat','thresh')
    load('border_metricsR1.mat','bcell_stb')
    load('corner_metricsR1.mat','C')
    for ii = 1:size(bcell_stb,1)
        b = bcell_stb{ii,1};
        bcell = bcell_stb(ii,:);
        bc_overlap{ii}.b = [bc_overlap{ii}.b; numel(b)/numel(thresh)];
        bc = cellfun(@(x,y) intersect(x,y), bcell, C.cornercell, 'uni',0);
        bc_overlap{ii}.ccav = [bc_overlap{ii}.ccav;numel(bc{1})/numel(b)];
        bc_overlap{ii}.cvex = [bc_overlap{ii}.cvex;numel(bc{3})/numel(b)];
    end
end
save('F:\Results_experimentC\bc_overlap_R1.mat','bc_overlap')
% concave cases
% determine the random level of overlap
rand_overlap = [];
nboot = 1000;
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('border_metricsR1.mat','bcell_stb') %load all session name
    load('thresh.mat','thresh')
    numNeuron = length(thresh);
    cc = numel(C.cornercell{1});
    rand_overlap_ea = NaN(nboot,size(bcell_stb,1));
    for jj = 1:size(bcell_stb,1)
        bcell = numel(bcell_stb{jj,1});
        for ii = 1:nboot
            p1 = randperm(numNeuron,cc);
            p2 = randperm(numNeuron,bcell);
            po = intersect(p1,p2);
            rand_overlap_ea(ii,jj) = numel(po)/numel(p2);
        end
    end
    rand_overlap = cat(3, rand_overlap, rand_overlap_ea);
end
% reshape the data
bc_overlap_rand = [];
for ii = 1:size(rand_overlap,2)
    data = squeeze(rand_overlap(:,ii,:));
    bc_overlap_rand = cat(3, bc_overlap_rand, data);
end
save('F:\Results_experimentC\bc_overlap_R1.mat','bc_overlap_rand','-append')
% convex cases
% determine the random level of overlap
rand_overlap_cvx = [];
nboot = 1000;
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('border_metricsR1.mat','bcell_stb') %load all session name
    load('thresh.mat','thresh')
    numNeuron = length(thresh);
    cc = numel(C.cornercell{3});
    rand_overlap_ea = NaN(nboot,size(bcell_stb,1));
    for jj = 1:size(bcell_stb,1)
        bcell = numel(bcell_stb{jj,1});
        for ii = 1:nboot
            p1 = randperm(numNeuron,cc);
            p2 = randperm(numNeuron,bcell);
            po = intersect(p1,p2);
            rand_overlap_ea(ii,jj) = numel(po)/numel(p2);
        end
    end
    rand_overlap_cvx = cat(3, rand_overlap_cvx, rand_overlap_ea);
end
% reshape the data
bc_overlap_randcvx = [];
for ii = 1:size(rand_overlap_cvx,2)
    data = squeeze(rand_overlap_cvx(:,ii,:));
    bc_overlap_randcvx = cat(3, bc_overlap_randcvx, data);
end
save('F:\Results_experimentC\bc_overlap_R1.mat','bc_overlap_randcvx','-append')
%plot
%plot
load('F:\Results_experimentC\bc_overlap_R1.mat','bc_overlap','bc_overlap_rand','bc_overlap_randcvx')
data = []; datacvx = []; data_rand = [];data_randcvx = [];
for ii = 1:length(bc_overlap)
    data = [data,bc_overlap{ii}.ccav];
    temp1 = bc_overlap_rand(:,:,ii);
    data_rand = [data_rand, quantile(temp1(:), 0.95)];
    datacvx = [datacvx, bc_overlap{ii}.cvex];
    temp2 = bc_overlap_randcvx(:,:,ii);
    data_randcvx = [data_randcvx, quantile(temp2(:), 0.95)];
end
figure
subplot(1,2,1)
stdshade(data, 0.3, [43 57 144]/255)
hold on
plot(data_rand, 'r')
axis square
ylim([0,0.15])
xticks([1:5])
xticklabels({'0.5', '0.6', '0.7', '0.8', '0.9'})
xlabel('boundary coverage threshold')
ylabel('prop of overlap')
title('BVC & cav overlap')
subplot(1,2,2)
stdshade(datacvx, 0.3, [43 57 144]/255)
hold on
plot(data_randcvx, 'r')
axis square
ylim([0,0.15])
xticks([1:5])
xticklabels({'0.5', '0.6', '0.7', '0.8', '0.9'})
xlabel('boundary coverage threshold')
ylabel('prop of overlap')
title('BVC & cvx overlap')

%venn diagram
venn = struct;
venn.a = []; venn.b = []; venn.c = [];
venn.ab = []; venn.ac = []; venn.bc =[];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('border_metricsR1.mat','bcell_stb')
    a = bcell_stb{3,1};
    b = C.cornercell{1};
    c = C.cornercell{3};
    venn.a = [venn.a; numel(a)/numel(a)];
    venn.b = [venn.b; numel(b)/numel(a)];
    venn.c = [venn.c; numel(c)/numel(a)];
    venn.ab = [venn.ab; numel(intersect(a,b))/numel(a)];
    venn.ac = [venn.ac; numel(intersect(a,c))/numel(a)];
    venn.bc = [venn.bc; numel(intersect(b,c))/numel(a)];
end
save('F:\Results_experimentC\bc_overlap_R1.mat','venn','-append')

%% plot anatomical locations of border and corner cells
load('neuronIndividualsf.mat','neuronIndividualsf')
load('corner_metricsR1.mat','C')
load('border_metricsR1.mat','bcell_stb')
neuronSelected1 = bcell_stb{3,1};
neuronSelected2 = C.cornercell{1};
neuronSelected3 = [C.cornercell{3};C.cornercell{4}];
neuron_overlap = intersect(neuronSelected1,neuronSelected2);
ncontour = neuronIndividualsf{1}.Coor;
figure
hold on
for ii = 1:length(ncontour)
    plot(ncontour{ii}(1,:),ncontour{ii}(2,:),'Color',[128,130,133]/255,'LineWidth', 0.5)
    set(gca, 'YDir','reverse')
    axis image
end
for jj = 1:length(neuronSelected1)
    plot(ncontour{neuronSelected1(jj)}(1,:),ncontour{neuronSelected1(jj)}(2,:),'Color',	'#4DBEEE', 'LineWidth', 1.2)
end
for jj = 1:length(neuronSelected2)
    plot(ncontour{neuronSelected2(jj)}(1,:),ncontour{neuronSelected2(jj)}(2,:),'Color',[43 120 142]/255, 'LineWidth', 1.2)
end
for jj = 1:length(neuronSelected3)
    plot(ncontour{neuronSelected3(jj)}(1,:),ncontour{neuronSelected3(jj)}(2,:),'Color',[147 38 103]/255, 'LineWidth', 1.2)
end
for jj = 1:length(neuron_overlap)
    plot(ncontour{neuron_overlap(jj)}(1,:),ncontour{neuron_overlap(jj)}(2,:), 'Color','y', 'LineWidth', 1.2)
end
savefig(gcf,'bc_anatomy.fig')

%% intra- and inter-cluster distance of boundary and corner cells
% measure pair-wise distance of intra- and inter-cluster of corner cells. 
% 1px = 0.89um (12.5mm Ach lens miniscope)
load('F:\analysis_folders.mat','expC')
datapath = expC;

pwdist = struct;
pwdist.iab = []; pwdist.iac = []; pwdist.inter = []; 
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat')
    load('corner_metricsR1.mat','C')
    load('border_metricsR1.mat','bcell_stb')
    bcell = bcell_stb{3,1};
    ctr = struct;
    ctr.b = neuronIndividualsf{1}.centroid(bcell,:);
    ctr.c = neuronIndividualsf{1}.centroid(union(C.cornercell{1},C.cornercell{3}),:);
    pwdist.iab = [pwdist.iab; mean(pdist(ctr.b))];
    pwdist.iac = [pwdist.iac; mean(pdist(ctr.c))];
    pwdist.inter = [pwdist.inter; mean(mean(pdist2(ctr.b,ctr.c)))];
end
save('F:\Results_experimentC\anatdist_border_corner_R1.mat','pwdist')
%plot
x1 = pwdist.iab .* 0.89;
x2 = pwdist.iac .* 0.89;
y = pwdist.inter .* 0.89;

figure;
plot(x1,y,'.','markersize',40)
hold on
plot(x2,y,'.','markersize',40)
axis square
legend('concave corner cell', 'convex corner cell')
xlim([50,120])
ylim([50,120])
xlabel('intra-group distance (um)')
ylabel('inter-group distance (um)')

%% example of summed rate map for border cells
load('firingrateAll.mat','ratemap')
load('corner_metrics.mat','C2')
load('border_metrics.mat','bcell')
ss = 1;
rtmap = ratemap{ss};
bc = bcell{ss};
cc = C2.cornercell{ss};
bccc = union(bc,cc);
sum_ratemap = [];
for ii = 1:length(bccc)
    sum_ratemap = cat(3, sum_ratemap, rtmap{bccc(ii)});
end

figure
imagesc(mean(sum_ratemap,3))
axis image

%% EXTRA PLOT
%% Extra plot for supplementary figures
%% plot cornerscore shuffling for supplementary figure
load ('F:\subiculum_mice_Aug2022\M4115c\corner_metrics.mat','C2')
figure
histogram(C2.cscore_shuffle{1, 2}(63,:),[-0.6:0.05:0.9])
hold on
line([C2.cscore_thresh{2}(63),C2.cscore_thresh{2}(63)],[0,140],'linewidth',1.5,'color',[128 130 133]/255)
line([C2.cscore{2}(63),C2.cscore{2}(63)],[0,140],'linewidth',1.5,'color',[0 148 68]/255)


%% EXPERIMENT H
%% Experiment H
%% Experiment H for introducing a 135 degree convex corners to the environment.
% Experiment H aims to look at the firing rate of convex corner cells in 90
% degree vs. 135 degree convex corners
load('F:\analysis_folders.mat','expH')
datapath = expH;

%% Identify corner cells ***
%mannually identify the location of corners for each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [S,N,env_coor,env_coori] = identify_env_geometry;
    save('env_geometry.mat','S','N','env_coor','env_coori')
end

%determine corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('env_geometry.mat','S')
    mask = cell(1,length(S));
    for jj = 1:length(S)
        mask{jj} = true;
    end
    %determine corner cell for each session by spike shuffling.
    %NOTE, for expC, C was using a 0.3/0.35 threshold for identifying
    %corner cells. C2 was using a 0.4 threshold for identifying corner
    %cells. C2 works slightly better becuase the data in the large
    %environment is bit noiser. 
    %C3 is to use the mask ahead of identification of field peaks, with 0.4
    %threshold and also tight location of corners. 
%     C3 = identify_corner_cell(mask);
%     save('corner_metrics.mat','C3','-append')
    CR1 = identify_corner_cell(mask,[],0.4);
    save('corner_metrics.mat','CR1','-append')
end

%% see the firing rate difference between the 270 degree and 225 degree corner cell ***
load('F:\analysis_folders.mat','expH')
datapath = expH;
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
session_name = {'convex2','convex2d90'}; %expC
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCors = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        if strcmp(datapath{n},'F:\subiculum_mice_Aug2022\M4131Fh')
            load('corner_metrics.mat','C2','pkfr_corners','pkfrsim_corners')
            C = C2;
        else
            load('corner_metrics.mat','CR1','pkfr_corners','pkfrsim_corners')
            C = CR1;
        end
        idx = strcmp(S, session_name{k});
        pkfr_corners = pkfr_corners(idx);
        cornercell = C.cornercell(idx);
        allcell = 1:length(pkfr_corners{1,1});
        ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
        ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        %         pkfrcc_atcorner = cellfun(@(x,y) x(:,cornercellx), pkfr_corners, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        %non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,pkfrcc_atcorner./pkfrncc_atcorner];
        %simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,pkfrcc_atcorner./pkfrsim_atcorner];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
end
data = nanmean(cat(3, pkfrcc_atCs{:}),3);
save('F:\Results_experimentC\pkfrcc_atCorner_expH_R1.mat','pkfrcc_atCnc','pkfrcc_atCs','data')
%plot using graphpad prism

%% Reviewer 1, corner score mask explanation
% plot corner score distribution for illustration purpose
msize = 25;
x  = zeros(msize);
corners = [5,5;5,msize-4;msize-4,msize-4;msize-4,5];
ctr = [round(msize/2),round(msize/2)];
cscore = NaN(msize);

for ii = 1:size(x,1)
    for jj = 1:size(x,2)
        d = pdist2([ii,jj], corners);
        dmin = min(d);
        dc = pdist2([ii,jj], ctr);
        cscore(ii,jj) = (dc-dmin)/(dmin+dc);
    end
end

figure
imagesc(cscore);
axis image
colormap winter

cscore(1:4,1:4) = NaN;
cscore(1:4,22:25) = NaN;
cscore(22:25,1:4) = NaN;
cscore(22:25,22:25) = NaN;

x = diag(cscore);
figure
plot(x(1:13))
axis square

%% Reviewer 2, replot allocentric bearing
load('neuronIndividualsf.mat')
figure
for ii = 1:length(neuronIndividualsf)
    neuron = neuronIndividualsf{ii};
    subplot(3,2,ii)
    scatter(neuron.pos(:,1),neuron.pos(:,2),2,neuron.ab,'filled')
    colormap hsv
    axis image
end
