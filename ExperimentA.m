%% Experiment A
% This code is for analyzing corner coding in the subiculum for animals
% running in the circle, square, triangle, and hexagon environments. 
load('F:\analysis_folders.mat','expA')
datapath = expA;
%% Distribution of corner cell stability
session_name = {'triangle','square','hex'}; %expA
CC_stability = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metrics.mat', 'Cthresh')
    C = Cthresh{3};
    load('spatial_metrics.mat', 'map_stb')
    load('env_geometry.mat','S') %load all session name
    idx = ismember(S, session_name);
    cornercell = C.cornercell(idx);
    map_stb = map_stb(idx);
    CCS = cellfun(@(x,y) x(y), map_stb, cornercell, 'uni',0);
    CC_stability = [CC_stability;vertcat(CCS{:})];
end
save('F:\Results_experimentA\stabilitycc_R1.mat','CC_stability','-append')
%% Revised definition of corner cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'map_stb')
    load('corner_metrics.mat', 'Cthresh')
    C = Cthresh{3};
    stablecell = cellfun(@(x) find(x > 0.3), map_stb, 'uni', 0);
    C.cornercell = cellfun(@(x,y) intersect(x,y), stablecell, C.cornercell, 'uni', 0);
    save('corner_metricsR1.mat', 'C', '-append')    
end
%% Percentage of corner cells for each session
session_name = {'circle','triangle','square','hex'}; %expA
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
save('F:\Results_experimentA\prop_cornercell_R1.mat','prop_cornercell','prop_cornercell_ea','-append')
%% compute the the firing rate as a function of distance to corners
session_name = {'triangle','square','hex'}; %expA
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
%plot dist2corner firing rate, ExpA
data = cell(1,length(session_name));
for ii = 1:length(session_name)
    temp = struct; temp.x=[];temp.ycc=[];temp.yncc=[];temp.yccx=[];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('corner_metricsR1.mat','d2cornerq')
        temp.x = padconcatenation(temp.x, nanmean(d2cornerq.d{ii},1),1);
        temp.ycc = padconcatenation(temp.ycc, d2cornerq.frcc{ii},1);
        temp.yncc = padconcatenation(temp.yncc, d2cornerq.frncc{ii},1);
        temp.yccx = padconcatenation(temp.yccx, d2cornerq.frccx{ii},1);
    end
    data{ii} = temp;
end
save('cornercell_linearFr_R1.mat','data')
%plot
rho_all = cellfun(@(x) corr(x.ycc',x.yncc','Rows','complete'), data, 'uni', 0);
rho = cellfun(@diag, rho_all,'uni', 0);

figure
subplot(1,3,1)
stdshade(smoothdata(data{1}.ycc,2,'movmean',3),0.4, [43 120 142]/255, [1:size(data{1}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{1}.yncc,2,'movmean',3),0.4, [109 110 113]/255, [1:size(data{1}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('triangle')
subplot(1,3,2)
stdshade(smoothdata(data{2}.ycc,2,'movmean',3),0.4, [43 120 142]/255,[1:size(data{2}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{2}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:size(data{2}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('square')
subplot(1,3,3)
stdshade(smoothdata(data{3}.ycc,2,'movmean',3),0.4, [43 120 142]/255,[1:size(data{3}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{3}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:size(data{3}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('hexagon')
%% compare the overall firing rate of corner cells for each corner
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
session_name = {'triangle','square','hex'}; %expA
pkfrcc_atCnc = cell(1,length(session_name));
pkfrcc_atCs = cell(1,length(session_name));
pkfrcc_atC = cell(1,length(session_name));
pkfrcc_atsimC = cell(1,length(session_name));
for k = 1:length(session_name)
    pkfrcc_atCor = []; pkfrcc_atCors = []; pkfrcc_atCo = []; pkfrcc_atsimCor = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
        load('corner_metricsR1.mat','C')
        idx = strcmp(S, session_name{k});
%         temp = cellfun(@numel, C.cornercell);
%         idx2 = temp >= 5;
%         idx = idx1 & idx2;
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
        pkfrcc_atCo = [pkfrcc_atCo,mean(pkfrcc_atcorner,2)];
        %non-corner cell
        pkfrncc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, ncc, 'uni', false);
        pkfrncc_atcorner = cellfun(@(x) mean(x,2), pkfrncc_atcorner, 'uni', 0);
        pkfrncc_atcorner = [pkfrncc_atcorner{:}];
        pkfrcc_atCor = [pkfrcc_atCor,mean(pkfrcc_atcorner./pkfrncc_atcorner,2)];
        %simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_atCors = [pkfrcc_atCors,mean(pkfrcc_atcorner./pkfrsim_atcorner,2)];
        pkfrcc_atsimCor = [pkfrcc_atsimCor,mean(pkfrsim_atcorner,2)];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
    pkfrcc_atC{k} = pkfrcc_atCo';
    pkfrcc_atsimC{k} = pkfrcc_atsimCor';
end
save('F:\Results_experimentA\pkfrcc_atCornerR1.mat',...
    'pkfrcc_atCnc','pkfrcc_atCs','pkfrcc_atC','pkfrcc_atsimC')

%plot the traj and ratemap for the simulated neuron data, if needed
load('behavIndividualsf.mat')
load('neuronSim.mat')
plot_rate_map_longitudinal(neuronSim,behavIndividualsf,ratemapSim,1,0,'S','auto'); %1:size(neuronIndividualsf{1}.C,1);
%plot examples for figure
cd('F:\subiculum_mice_Aug2022\M4115a')
load('neuronSim.mat','ratemapSim')
load('firingrateAll.mat','ratemap')
figure
subplot(2,2,1)
pcolor(ratemap{4}{47})
axis image
shading flat
title('raw ratemap')
subplot(2,2,2)
pcolor(ratemapSim{4}{1})
shading flat
axis image
title('simulated ratemap')
subplot(2,2,3)
pcolor(ratemap{4}{47}./ratemapSim{4}{1})
axis image
shading flat
title('corrected ratemap')
%% compare the peak firing rate variations of corner cells at each corner (i.e. max pkfr/min pkfr)
% session_name = {'triangle','square','hex'}; %expA
% pkfrcc_vary = cell(1,length(session_name));
% for k = 1:length(session_name)
%     pkfrcc_atCors = [];
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('env_geometry.mat','S')
%         load('corner_metrics.mat','pkfr_corners','C','pkfrsim_corners')
%         idx1 = strcmp(S, session_name{k});
%         temp = cellfun(@numel, C.cornercell);
%         idx2 = temp >= 5;
%         idx = idx1 & idx2;
%         pkfr_corners = pkfr_corners(idx);
%         cornercell = C.cornercell(idx);
%         %corner cell
%         pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
%         pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
%         pkfrcc_atcorner = [pkfrcc_atcorner{:}];
%         %simulated cell
%         pkfrsim_atcorner = pkfrsim_corners(idx);
%         pkfrsim_atcorner = [pkfrsim_atcorner{:}];
%         pkfrcc_corrected = pkfrcc_atcorner./pkfrsim_atcorner;
%         pkfrcc_atCors = [pkfrcc_atCors,max(pkfrcc_corrected)./min(pkfrcc_corrected)];
%     end
%     pkfrcc_vary{k} = pkfrcc_atCors';
% end
% save('F:\Results_experimentA\pkfrcc_atCorner.mat','pkfrcc_vary','-append')

%% XSESSION CORNER CELLS
%% xsession corner cell analysis
%% Percentage of xsession corner cells
prop_cornercellx = NaN(length(datapath),1);
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','cornercellx')
    load('thresh.mat','thresh')
    numNeurons = length(thresh);
    numccx = numel(cornercellx);
    prop_cornercellx(n,1) = numccx/numNeurons;
end
save('F:\Results_experimentA\prop_cornercell_R1.mat','prop_cornercellx','-append')
%% compare the firing rate at corners across different environment (triangle, square, hex)
session_name = {'triangle','square','hex'}; %expA
pkfrccx_cmp = NaN(length(datapath),length(session_name));
pkfrccx_cmpcrt = NaN(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metrics.mat','pkfr_corners','pkfrsim_corners')
    load('corner_metricsR1.mat','cornercellx')
    load('env_geometry.mat', 'S')
    for ii = 1:length(session_name)
        %raw firing rate
        idx = strcmp(S, session_name{ii});
        pkfr = pkfr_corners(idx);
        pkfrccx = mean(cellfun(@(x) mean(mean(x(:,cornercellx))), pkfr));
        pkfrccx_cmp(n,ii) = pkfrccx;
        %corrected firing rate
        pkfrs = pkfrsim_corners(idx);
        pkfr_crt = cellfun(@(x,y) x./y, pkfr, pkfrs, 'uni',0);
        pkfrccx_crt = mean(cellfun(@(x) mean(mean(x(:,cornercellx))), pkfr_crt));
        pkfrccx_cmpcrt(n,ii) = pkfrccx_crt;
    end
end
save('F:\Results_experimentA\pkfrcc_atCornerR1.mat','pkfrccx_cmp','pkfrccx_cmpcrt','-append')
%% compute the STABILITY of xsession cornercells
%within session stability
%map_stb: using interleaved method
%map_stability: using split to half method
for k = 1:length(datapath)
    cd(datapath{k});
    load('neuronIndividualsf.mat')
    load('behavIndividualsf.mat')
    load('thresh.mat')
    map_stb = cellfun(@(x,y) calc_spatialmap_stability(x,y,thresh), ...
        neuronIndividualsf, behavIndividualsf, 'uni',0);
    save('spatial_metrics.mat', 'map_stb','-append')
end
%gather data for plotting
session_name = {'triangle','square','hex'}; %expA
stability = cell(1,length(session_name));
for k = 1:length(session_name)
    stab = [];
    for n = 1:length(datapath)
        cd(datapath{n})
        load('env_geometry.mat','S')
        load('corner_metrics.mat','cornercellx')
        load('spatial_metrics.mat', 'map_stability')
        idx = strcmp(S, session_name{k});
        map_stb = map_stb(idx);
        allcell = 1:length(map_stb{1,1});
        ind = ~ismember(allcell, cornercellx);
        ncc = allcell(ind);
        %corner cell
        stabcc = cellfun(@(x) x(cornercellx), map_stb, 'uni', false);
        stabcc = cellfun(@(x) nanmean(x), stabcc, 'uni', 0);
        stabcc = [stabcc{:}];
        %non-corner cell
        stabnc = cellfun(@(x) x(ncc), map_stb, 'uni', false);
        stabnc = cellfun(@(x) nanmean(x), stabnc, 'uni', 0);
        stabnc = [stabnc{:}];
        stab = [stab;[stabcc',stabnc']];
    end
    stability{k} = stab;
end

%cross session stability
session_name = {'triangle','square','hex'}; %expA
for k = 1:length(datapath)
    cd(datapath{k});
    load('env_geometry.mat','S')
    load('firingrateAll.mat','ratemap')
    map_stabilityx = cell(1,length(session_name));
    for n = 1:length(session_name)
        idx = strcmp(S, session_name{n});
        map_all = ratemap(idx);
        if size(map_all,2) > 1
            map1 = map_all{1};
            map2 = map_all{2};
            map_stabilityx{1,n} = calc_spatialmap_stabilityx(map1,map2);
        else
            map_stabilityx{1,n} = NaN(length(map1),1);
        end
    end
    save('spatial_metrics.mat', 'map_stabilityx', '-append')
end
%gather data for plotting
stabilityx = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','cornercellx')
    load('spatial_metrics.mat', 'map_stabilityx')
    %cornercellx
    stabcc = cellfun(@(x) x(cornercellx), map_stabilityx, 'uni', false);
    stabcc = cellfun(@(x) nanmean(x), stabcc);
    stabilityx = [stabilityx;stabcc];
end
save('F:\Results_experimentA\stabilityccxR1.mat','stabilityx')

%% DECODING
%% neural decoding analysis
%% Within session decoding for each session
for k = 1:length(datapath)
    cd(datapath{k})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat','thresh')
    decodeResults = bayes_decoder2D_withinsession(neuronIndividualsf,behavIndividualsf,thresh);
    decodeShuffle = bayes_decoder2D_withinsession_shuffle(neuronIndividualsf,behavIndividualsf,thresh);
    save('decode_wisession.mat','decodeResults','decodeShuffle','-v7.3')
end

session_name = {'circle','triangle','square','hex'}; %expA
decode_wisession = NaN(length(datapath)*2,length(session_name)*2);
for ii = 1:length(session_name)
    k=1;
    for n = 1:length(datapath)
        cd(datapath{n})
        load('decode_wisession.mat')
        load('env_geometry.mat', 'S')
        idx = strcmp(S, session_name{ii});
        de = mean(cellfun(@(x) median(x.decodeError2step), decodeResults));
        de = de(idx);
        ds = mean(cellfun(@median, decodeShuffle));
        ds = ds(idx);
        idx2 = length(find(idx == 1));
        decode_wisession(k:k+idx2-1,ii*2-1:ii*2) = [de',ds'];
        k = k+idx2;
    end
end
save('F:\Results_experimentA\decode_results.mat','decode_wisession')

%plot an example decoding trace
load('F:\subiculum_mice_Aug2022\M4113Fa\decode_wisession.mat','decodeResults')
x1 = decodeResults{8,2};
figure
plot(x1.positionVec_ds,'k','linewidth',0.65)
hold on
plot(x1.decodePositionVec2step,'r')
xlim([1,110])
xticklabels([0:8:88])
xlabel('time (sec)')
ylabel('vectorized position (bins)')

%plot an example decoding trace in 2D
x1 = decodeResults{8,2};
figure
subplot(2,1,1)
plot(x1.positionBin_ds(:,1), 'Color', 'k','linewidth',0.65)
hold on,
plot(x1.decodePosition2step(:,1),'Color','r')
xlim([1,110])
xticklabels([0:8:88])
ylabel('x position (cm)')
xlabel('time (sec)')
subplot(2,1,2)
plot(x1.positionBin_ds(:,2), 'Color', 'k','linewidth',0.65)
hold on,
plot(x1.decodePosition2step(:,2),'Color','r')
xlim([1,110])
xticklabels([0:8:88])
ylabel('y position (cm)')
xlabel('time (sec)')

%% Within session decoding with corner cell KO
for k = 1:length(datapath)
    cd(datapath{k})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat','thresh')
    load('corner_metricsR1.mat', 'C')
    decodeResultsKO = bayes_decoder2D_withinsession(neuronIndividualsf,behavIndividualsf,thresh,[],true,C.cornercell);
    save('decode_wisessionR1.mat','decodeResultsKO')
end

%% Compare decoding errors at the corners between normal and KO decoding
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('decode_wisession.mat','decodeResults','decodeResultsKO')
%     de_corners = find_decodeError_at_corners(decodeResults);
%     de_cornersKO = find_decodeError_at_corners(decodeResultsKO);
%     save('decode_wisession.mat','de_corners','de_cornersKO','-append')
% end
% 
% session_name = {'triangle','square','hex'}; %expA
% %ccRatio: corner to center ratio of the decoding error
% de_ccRatio = NaN(length(datapath)*2,length(session_name)*2);
% for ii = 1:length(session_name)
%     k=1;
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('decode_wisession.mat','de_corners','de_cornersKO')
%         load('env_geometry.mat', 'S')
%         idx = strcmp(S, session_name{ii});
%         % ccRatio: corner to center ratio
%         de = cellfun(@(x) mean(x(1:end-1)./x(end)),de_corners);
%         deKO = cellfun(@(x) mean(x(1:end-1)./x(end)),de_cornersKO);
%         de = de(idx);
%         deKO = deKO(idx);
%         idx2 = length(find(idx == 1));
%         de_ccRatio(k:k+idx2-1,ii*2-1:ii*2) = [de',deKO'];
%         k = k+idx2;
%     end
% end
% save('F:\Results_experimentA\decode_results.mat','de_ccRatio','-append')
% %plot triangle
% figure
% plot(de_ccRatio(:,1:2)','.-','markersize',30,'color','k')
% xlim([0.5,2.5])
% xticks([1:2])
% xticklabels({'full','ccKO'})
% % ylim([0.7,1.4])
% ylabel('corner to center ratio for decoding errors')

%use the decoding error map full decoder / the error map of KO decoder
for k = 1:length(datapath)
    cd(datapath{k})
    load('decode_wisession.mat','decodeResults')
    load('decode_wisessionR1.mat','decodeResultsKO')
    de_fullvsKO = compare_decodeError_at_corners(decodeResults,decodeResultsKO);
    save('decode_wisessionR1.mat','de_fullvsKO','-append')
end

session_name = {'triangle','square','hex'}; %expA
%ccRatio: corner to center ratio of the decoding error
de_ccCmp = NaN(length(datapath)*2,length(session_name)*2);
for ii = 1:length(session_name)
    k=1;
    for n = 1:length(datapath)
        cd(datapath{n})
        load('decode_wisessionR1.mat','de_fullvsKO')
        load('env_geometry.mat', 'S')
        idx = strcmp(S, session_name{ii});
        % ccRatio: corner to center ratio
        deCor = cellfun(@(x) mean(x(1:end-1)),de_fullvsKO);
        deCtr = cellfun(@(x) x(end),de_fullvsKO);
        deCor = deCor(idx);
        deCtr = deCtr(idx);
        idx2 = length(find(idx == 1));
        de_ccCmp(k:k+idx2-1,ii*2-1:ii*2) = [deCor',deCtr'];
        k = k+idx2;
    end
end
save('F:\Results_experimentA\decode_results_R1.mat','de_ccCmp')
%plot
%plot triangle
figure
subplot(1,3,1)
plot(de_ccCmp(:,1:2)','.-','markersize',25,'color','k')
xlim([0.5,2.5])
xticks([1:2])
xticklabels({'corner','center'})
ylim([0.9,2.1])
ylabel('corner to center ratio for decoding errors')
subplot(1,3,2)
plot(de_ccCmp(:,3:4)','.-','markersize',25,'color','k')
xlim([0.5,2.5])
xticks([1:2])
xticklabels({'corner','center'})
ylim([0.9,2.1])
ylabel('corner to center ratio for decoding errors')
subplot(1,3,3)
plot(de_ccCmp(:,5:6)','.-','markersize',25,'color','k')
xlim([0.5,2.5])
xticks([1:2])
xticklabels({'corner','center'})
ylim([0.9,2.1])
ylabel('corner to center ratio for decoding errors')
[h,p] = signrank(de_ccCmp(:,1),de_ccCmp(:,2))

%plot individual session for example figure
%M4111a ss7 is a good example for triangle
%M4102a ss3 is a good example for square
%M4111a ss2 is a good example for square
%M4116a ss4 is a good example for square
%M4104a ss3 is a good example for square
load('decode_wisession.mat','decodeResults')
load('decode_wisessionR1.mat','decodeResultsKO')
ss = 4;
de2 = []; de2KO = [];
for ii = 1:size(decodeResults,1)
    de2 = cat(3,de2, decodeResults{ii,ss}.decodeError2_2step);
    de2KO = cat(3,de2KO, decodeResultsKO{ii,ss}.decodeError2_2step);
end
de2 = nanmean(de2,3);
de2KO = nanmean(de2KO,3);
de = de2KO./de2;
de = gaussianfilter2D(de, 3, 1.5);
figure
imagesc(de)
axis image
colormap bone

%% DECODING USING ONLY CORNER CELLS
%% Decoding using only corner cells to train the decoder
%% Decoding compare the decoded quadrant vs. the true quadrant
load('F:\analysis_folders.mat','expA')
datapath = expA;
session_name = 'square'; %expA
for k = 1:length(datapath)
    cd(datapath{k})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat','thresh')
    load('corner_metricsR1.mat', 'C')
    load('env_geometry.mat', 'S')
    idx = strcmp(S, session_name);
%     cornercellx = {cornercellx};
%     cornercellx = repmat(cornercellx, 1, length(idx));
    decodeResultsCC = bayes_decoder2D_withinsession(neuronIndividualsf(idx),behavIndividualsf(idx),...
        thresh,C.cornercell(idx),false,[],[],false);
    %generate a matrix that have four quadrants
    load('firingrateAll.mat','countTime')
    load('env_geometry.mat', 'env_coor')
    env_coor = env_coor(idx);
    mquad_raw = cellfun(@(x) NaN(size(x,1),size(x,2)), countTime(idx), 'uni', 0);
    mquad = mquad_raw;
    for n = 1:length(mquad_raw)
        D = [];
        for ii = 1:size(mquad_raw{n},2)
            for jj = 1:size(mquad_raw{n},1)
                D = [D; pdist2([ii,jj],env_coor{n}(1:end-1,:))];
            end
        end
        [~,I] = min(D,[],2);
        mquad{n} = reshape(I, [size(mquad_raw{n},1), size(mquad_raw{n},2)]);
    end
    save(['decode_quad_',session_name,'_R1.mat'],'decodeResultsCC','mquad','-v7.3')
    % compute the corner agreement between true and decoded positions (quadrants)
    mquad = repmat(mquad, size(decodeResultsCC,1), 1);
    truequad = cellfun(@(x,y) diag(y(x.positionBin_ds(:,2),x.positionBin_ds(:,1))), decodeResultsCC, mquad,'uni', 0);
    decodedquad = cellfun(@(x,y) diag(y(x.decodePosition(:,2),x.decodePosition(:,1))), decodeResultsCC, mquad,'uni', 0);
    correctquad = cellfun(@(x,y) x==y, truequad, decodedquad,'uni', 0);
    save(['decode_quad_',session_name,'_R1.mat'],'truequad','decodedquad','correctquad','-append')
end
% shuffle condition
for k = 1:length(datapath)
    cd(datapath{k})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load('behavIndividualsf.mat','behavIndividualsf')
    load('thresh.mat','thresh')
    load('corner_metrics.mat', 'C')
    load('env_geometry.mat', 'S')
    load(['decode_quad_',session_name,'_R1.mat'],'mquad')
    idx = strcmp(S, session_name);
    [decodeShuffle,correctquadShuffle] = bayes_decoder2D_withinsession_shuffle_quad(neuronIndividualsf(idx), behavIndividualsf(idx), ...
        thresh, mquad, C.cornercell(idx));
    save(['decode_quad_',session_name,'_R1.mat'],'decodeShuffle','correctquadShuffle','-append')
end

% combine results from all mice
session_name = 'square'; %expA
Cquad = []; CquadShuffle = [];
for k = 1:length(datapath)
    cd(datapath{k})
    load(['decode_quad_',session_name,'_R1.mat'],'correctquad','correctquadShuffle')
    %     load('corner_metrics.mat', 'cornercellx')
    x = mean(cellfun(@(x) sum(x)/numel(x), correctquad));
    Cquad = [Cquad; mean(x(:))];
    xs = mean(cellfun(@(x) sum(x)/numel(x), correctquadShuffle));
    CquadShuffle = [CquadShuffle; mean(xs(:))];
end
figure
plot([Cquad, CquadShuffle]') % find the data in GraphPad
[p,h] = signrank(Cquad, CquadShuffle)

% plot example trace
cd('F:\subiculum_mice_Aug2022\M4113Fa')
load('decode_quad_square_R1.mat', 'truequad', 'decodedquad','correctquad')
x = cellfun(@(x) sum(x)/numel(x), correctquad);
figure
plot(truequad{9}(10:41),'k','LineWidth',1)
hold on
plot(decodedquad{9}(10:41),'r', 'LineStyle', '--','LineWidth',1)
ylim([0.5,4.5])
xlabel('time (x0.8 s)')
ylabel('quadrant')

%% WALL HEIGHT
%% Wall height manipulations
%% Percentage of corner cells in square and low wall square
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    CR1 = identify_corner_cell;
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
%% Organize the data and plot
load('F:\analysis_folders.mat','expA2')
datapath = expA2;
session_name = {'square','squarelw'}; %expA
prop_cornercell = cell(length(datapath),length(session_name));
prop_placecell = cell(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat','C')
    load('spatial_metrics.mat','placecell')
    load('env_geometry.mat','S') %load all session name
    load('thresh.mat','thresh')
    numNeurons = length(thresh);
    pcc = cell(1,length(session_name));
    ppc = cell(1,length(session_name));
    for ii = 1:length(session_name)
        idx = strcmp(S, session_name{ii});
        C_ea = C.cornercell(idx);
        C_ea = cellfun(@numel, C_ea);
        pcc{ii} = C_ea(:)/numNeurons;
        PC_ea = placecell(idx);
        PC_ea = cellfun(@numel, PC_ea);
        ppc{ii} = PC_ea(:)/numNeurons;
    end
    prop_cornercell(n,:) = pcc;
    prop_placecell(n,:) = ppc;
end
prop_cornercell_ea = cellfun(@mean, prop_cornercell);
prop_placecell_ea = cellfun(@mean, prop_placecell);
save('F:\Results_experimentA\prop_cornercell_lowWall_R1.mat','prop_cornercell',...
    'prop_cornercell_ea','prop_placecell','prop_placecell_ea')
%% compare the peak fr at corners 
s1 = 'square';
s2 = 'squarelw';
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat', 'S')
    idx1 = strcmp(s1,S);
    idx2 = strcmp(s2,S);
    cornercell_sq = C.cornercell(idx1);
    cornercell_lw = C.cornercell(idx2);
    cornercell_sq = union(cornercell_sq{1},cornercell_sq{end});
    cornercell_lw = union(cornercell_lw{1},cornercell_lw{end});
    cornercellxlw = intersect(cornercell_sq, cornercell_lw);
    save('corner_metricsR1.mat','cornercellxlw','-append')
end

%gather data for plotting
session_name = {'square','squarelw'}; %expA
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
        pkfr_corners = pkfr_corners(idx);
        cornercell = C.cornercell(idx);
        %corner cell
        pkfrcc_atcorner = cellfun(@(x,y) x(:,y), pkfr_corners, cornercell, 'uni', false);
        pkfrcc_atcorner = cellfun(@(x) mean(x,2), pkfrcc_atcorner, 'uni', 0);
        pkfrcc_atcorner = [pkfrcc_atcorner{:}];
        %simulated cell
        pkfrsim_atcorner = pkfrsim_corners(idx);
        pkfrsim_atcorner = [pkfrsim_atcorner{:}];
        pkfrcc_corrected = pkfrcc_atcorner./pkfrsim_atcorner;
        pkfrcc_atCors = [pkfrcc_atCors,mean(pkfrcc_corrected,2)];
    end
    pkfrcc_atCnc{k} = pkfrcc_atCor';
    pkfrcc_atCs{k} = pkfrcc_atCors';
end
save('F:\Results_experimentA\pkfrcc_atCorner_lowWall_R1.mat','pkfrcc_atCs')

%% off corner cell from square to square low wall. 
load('F:\analysis_folders.mat','expA2')
datapath = expA2;

ss1 = 'square';
ss2 = 'squarelw';
data_ea = struct;
data_ea.x=[];
data_ea.ycc=[];
data_ea.yncc=[];
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','d2cornerb')
    load('corner_metricsR1.mat','C')
    sslength = 1:length(S);
    idx1 = strcmp(S, ss1);
    idx1 = sslength(idx1);
    idx2 = strcmp(S, ss2);
    idx2 = sslength(idx2);
    cornercell1 = C.cornercell{idx1(end)};
    cornercell2 = C.cornercell{idx2(1)};
    ind = ~ismember(cornercell1,cornercell2);
    cornercelloff = cornercell1(ind);
    allcell = [1:length(d2cornerb.fr{1})]';
    allcc = union(cornercelloff,cornercell2);
    noncornercell = allcell(~ismember(allcell, allcc));
    data_ea.x = padconcatenation(data_ea.x, d2cornerb.d{idx2(1)},1);
    data_ea.ycc = padconcatenation(data_ea.ycc, nanmean(d2cornerb.fr{idx2(1)}(cornercelloff,:),1),1);
    data_ea.yncc = padconcatenation(data_ea.yncc, nanmean(d2cornerb.fr{idx2(1)}(noncornercell,:),1),1);
end
datax{1} = data_ea;
%datax organized as square in squarelw
save('F:\Results_experimentA\cornercelloff_lowWall_R1.mat','datax') 

%plot
figure
stdshade(smoothdata(datax{1}.ycc,2,'movmean',3),0.4, [43 120 142]/255, [1:size(datax{1}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(datax{1}.yncc,2,'movmean',3),0.4, [109 110 113]/255, [1:size(datax{1}.yncc,2)]*1.55)
axis square
xlim([0,24])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('square in squarelw')

[h,p] = signrank(mean(datax{1}.ycc(:,1:3),2),mean(datax{1}.yncc(:,1:3),2))
[h,p] = signrank(mean(datax{1}.ycc(:,end-2:end),2),mean(datax{1}.yncc(:,end-2:end),2))

%% BOUNDARY CELL ANALYSIS
%% Border Cell
%% plot anatomical locations of border and corner cells
load('neuronIndividualsf.mat','neuronIndividualsf')
load('env_geometry.mat','S')
load('corner_metrics.mat','C')
load('border_metrics.mat','bcell')
idx = strcmp(S,'square');
ss = [1:length(S)];
ss = ss(idx);
ss = ss(1);
neuronSelected1 = bcell{ss};
neuronSelected2 = C.cornercell{ss};
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
    plot(ncontour{neuronSelected1(jj)}(1,:),ncontour{neuronSelected1(jj)}(2,:),'Color','#0072BD', 'LineWidth', 1.2)
end
for jj = 1:length(neuronSelected2)
    plot(ncontour{neuronSelected2(jj)}(1,:),ncontour{neuronSelected2(jj)}(2,:),'Color','#D95319', 'LineWidth', 1.2)
end
for jj = 1:length(neuron_overlap)
    plot(ncontour{neuron_overlap(jj)}(1,:),ncontour{neuron_overlap(jj)}(2,:), 'Color','#77AC30', 'LineWidth', 1.2)
end
savefig(gcf,'bc_anatomy.fig')


%% PLOT CALCIUM TRACE
%% Plot Calcium trace for supplementary figure
%% plot calcium trace with deconvolved spikes
load('F:\subiculum_mice_Aug2022\M4115a\neuronIndividuals.mat')
neuron = 47;
spike = neuronIndividuals{1,2}.S(neuron,1001:5000) > 0;
spike = double(spike);
spike(spike==0) = NaN;
neuron2 = 54;
spike2 = neuronIndividuals{1,2}.S(neuron2,5001:9000) > 0;
spike2 = double(spike2);
spike2(spike2==0) = NaN;

idx = find(spike ==1);
figure;
subplot(2,1,1)
plot(neuronIndividuals{1,2}.trace(neuron,1001:5000))
hold on
for ii = 1:length(idx)
    line([idx(ii),idx(ii)],[40,45],'color','r')
end
subplot(2,1,2)
plot(neuronIndividuals{1,2}.trace(neuron2,5001:9000))
hold on
plot(spike2, '.')


%% ON/OFF CORNER CELL
%% on/off corner cell across different environments
%% analyze neurons that are corner cells in one environment, but not the other environment
%% across same condition
load('F:\analysis_folders.mat','expA')
datapath = expA;
% datapath = datapath(2:end);

ss = 'square';
data_ea = struct;
data_ea.x=[];
data_ea.ycc=[];
data_ea.yncc=[];
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','d2cornerb')
    load('corner_metricsR1.mat','C')
    sslength = 1:length(S);
    idx = strcmp(S, ss);
    idx = sslength(idx);
    cornercell1 = C.cornercell{idx(1)};
    cornercell2 = C.cornercell{idx(2)};
    ind = ~ismember(cornercell1,cornercell2);
    cornercelloff = cornercell1(ind);
    allcell = [1:length(d2cornerb.fr{1})]';
    allcc = union(cornercelloff,cornercell2);
    noncornercell = allcell(~ismember(allcell, allcc));
    data_ea.x = padconcatenation(data_ea.x, d2cornerb.d{idx(2)},1);
    data_ea.ycc = padconcatenation(data_ea.ycc, nanmean(d2cornerb.fr{idx(2)}(cornercelloff,:),1),1);
    data_ea.yncc = padconcatenation(data_ea.yncc, nanmean(d2cornerb.fr{idx(2)}(noncornercell,:),1),1);
end
% data = cell(1,3);
data{2} = data_ea;
save('F:\Results_experimentA\cornercelloff_R1.mat','data') %organized as trignale, square, hex

%plot
figure
subplot(1,3,1)
stdshade(smoothdata(data{1}.ycc,2,'movmean',3),0.4, [43 120 142]/255, [1:size(data{1}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{1}.yncc,2,'movmean',3),0.4, [109 110 113]/255, [1:size(data{1}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('triangle')
subplot(1,3,2)
stdshade(smoothdata(data{2}.ycc,2,'movmean',3),0.4, [43 120 142]/255,[1:size(data{2}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{2}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:size(data{2}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('square')
subplot(1,3,3)
stdshade(smoothdata(data{3}.ycc,2,'movmean',3),0.4, [43 120 142]/255,[1:size(data{3}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(data{3}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:size(data{3}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('hexagon')

[h,p] = signrank(mean(data{2}.ycc(:,1:3),2),mean(data{2}.yncc(:,1:3),2))
[h,p] = signrank(mean(data{2}.ycc(:,end-3:end-1),2),mean(data{2}.yncc(:,end-3:end-1),2))

%% off corner cell across different condition
load('F:\analysis_folders.mat','expA')
datapath = expA;

ss1 = 'triangle';
ss2 = 'square';
data_ea = struct;
data_ea.x=[];
data_ea.ycc=[];
data_ea.yncc=[];
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat','S')
    load('corner_metrics.mat','d2cornerb')
    load('corner_metricsR1.mat','C')
    sslength = 1:length(S);
    idx1 = strcmp(S, ss1);
    idx1 = sslength(idx1);
    idx2 = strcmp(S, ss2);
    idx2 = sslength(idx2);
    cornercell1 = C.cornercell{idx1(end)};
    cornercell2 = C.cornercell{idx2(1)};
    ind = ~ismember(cornercell1,cornercell2);
    cornercelloff = cornercell1(ind);
    allcell = [1:length(d2cornerb.fr{1})]';
    allcc = union(cornercelloff,cornercell2);
    noncornercell = allcell(~ismember(allcell, allcc));
    data_ea.x = padconcatenation(data_ea.x, d2cornerb.d{idx2(1)},1);
    data_ea.ycc = padconcatenation(data_ea.ycc, nanmean(d2cornerb.fr{idx2(1)}(cornercelloff,:),1),1);
    data_ea.yncc = padconcatenation(data_ea.yncc, nanmean(d2cornerb.fr{idx2(1)}(noncornercell,:),1),1);
end
% datax = cell(1,4);
datax{1} = data_ea;
%datax organized as trignale in square, trignale in hex, square in triangle,
%square in hex
save('F:\Results_experimentA\cornercelloff_R1.mat','datax','-append') 

[h,p] = signrank(mean(datax{2}.ycc(:,1:3),2),mean(datax{2}.yncc(:,1:3),2))
[h,p] = signrank(mean(datax{2}.ycc(:,end-3:end-1),2),mean(datax{2}.yncc(:,end-3:end-1),2))

%plot
figure
subplot(1,3,1)
stdshade(smoothdata(datax{1}.ycc,2,'movmean',3),0.4, [43 120 142]/255, [1:size(datax{1}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(datax{1}.yncc,2,'movmean',3),0.4, [109 110 113]/255, [1:size(datax{1}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('triangle in square')
subplot(1,3,2)
stdshade(smoothdata(datax{2}.ycc,2,'movmean',3),0.4, [43 120 142]/255,[1:size(datax{2}.ycc,2)]*1.55)
hold on
stdshade(smoothdata(datax{2}.yncc,2,'movmean',3),0.4, [109 110 113]/255,[1:size(datax{2}.yncc,2)]*1.55)
axis square
xlim([0,23])
ylim([0.05,0.55])
xlabel('distance to corners')
ylabel('firing rate (spikes/sec)')
title('triangle in hex')

%% CORNER SCORE METHOD
%% Some considerations of corner score method
%% Compute the percentage of corner cells that were corrected by the penalty
% As per Doug's comment, this code is compute the percentage of
% corner cells that went through the penalty process during the
% identification process. 
load('F:\analysis_folders.mat', 'expA')
datapath = expA;
session_name = {'triangle','square','hex'}; %expA
pct_pntyA = cell(1,length(datapath));
for n = 1:length(datapath)
    cd(datapath{n})
    load('env_geometry.mat', 'S')
    load('corner_metrics.mat','Cthresh')
    C = Cthresh{3};
    idx1 = ismember(S, session_name);
    cornercell = C.cornercell(idx1);
    cmetrics = C.cmetrics(idx1);
    cornerp = cellfun(@(x,y) x.cornerscorep(y), cmetrics, cornercell, 'uni',0);
    pct = cellfun(@(x) sum(~isempty(x))./length(x), cornerp);
    pct_pntyA{n} = pct;
end
save('F:\Results_experimentA\penalty_percent_R1.mat','pct_pntyA','-append')
P = cellfun(@mean, pct_pntyA);
mean(P)
std(P)/length(P)

%% plot corner score distribution using field shuffling for comparison (Doug's suggestion)
msize = 26;
x  = zeros(msize);
corners = [1,1;1,msize;msize,msize;msize,1];
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

% get shuffle using another method
load('corner_metrics.mat','C2')
C = C2;
numcell = length(C.cscore_thresh{1});
ss = 1;
nboot = 1000;
shuffle_score = NaN(numcell,nboot);
numcorner = 4;
for n = 1:numcell
    numfield = size(C.cmetrics{ss}.field_coor{n},1);
    for ii = 1:nboot
        coor = randi(msize,numfield,2);
        nscore = [];
        for jj = 1:size(coor,1)
            nscore = [nscore;cscore(coor(jj,1),coor(jj,2))];
        end
        shuffle_score(n,ii) = sum(nscore)/numcorner;
    end
end
thresh = quantile(shuffle_score, 0.95, 2);
cornercell = find(C.cscore{ss} > thresh);
cornercell = intersect(cornercell, C.cell_largedist{ss});
idx = ~ismember(C.cornercell{ss}, cornercell);
othercell = C.cornercell{ss}(idx);

%% EXTRA PLOT
%% Extra plot for Fig. S1, corner score distribution
%% plot cornerscore shuffling for supplementary figure
load ('F:\subiculum_mice_Aug2022\M4116c\corner_metrics.mat','C2')
C = C2;
session2plot = 1;
cell2plot = 22;
C.cmetrics{session2plot}.cornerscore_ea{cell2plot}
figure
histogram(C.cscore_shuffle{1, session2plot}(cell2plot,:),[-0.5:0.05:0.8])
hold on
line([C.cscore_thresh{session2plot}(cell2plot),C.cscore_thresh{session2plot}(cell2plot)],[0,120],'linewidth',1.5,'color',[128 130 133]/255)
line([C.cscore{session2plot}(cell2plot),C.cscore{session2plot}(cell2plot)],[0,120],'linewidth',1.5,'color',[0 148 68]/255)
axis square

%plot corresponding ratemap
load ('F:\subiculum_mice_Aug2022\M4116c\firingrateAll.mat','ratemap')
figure
imagesc(ratemap{session2plot}{cell2plot})
axis image
colormap jet

%% MANIFOLD
%% Manifold embedding (Revision 1)
%% Use PCA and UMAP to reveal the structure of neural manifold
datapath = pwd;
cd(datapath);
session_name = {'triangle','square','hex'}; %expA
load('neuronIndividualsf.mat')
load('behavIndividualsf.mat')
load('thresh.mat','thresh')
load('env_geometry.mat','env_coor','S')
session_number = 1:length(S);
idx = ismember(S, session_name);
ss2plot = session_number(idx);
[reduction_all, umap_template_all, clusterIds_all] = manifold_embedding(...
    neuronIndividualsf,behavIndividualsf,thresh,env_coor,S,ss2plot);
save('manifold_embedding.mat','reduction_all','umap_template_all','clusterIds_all','-v7.3')

%% plot neural acitivty on umap manifold
% load('env_geometry.mat','S')
% cell2plot = 47;
% ax = figure;
% set(ax, 'Position', [0, 200, 1200, 800]);
% hold on;
% for ii = 1:length(reduction_all)
%     reduction = reduction_all{ii};
%     ss = ss2plot(ii);
%     ssname = S{ss};
%     data_raw = neuronIndividualsf{ss}.S(cell2plot,:);
%     data = data_raw > thresh(cell2plot);
%     subplot(round(length(ss2plot)/2),2,ii)
%     scatter3(reduction(:,1),reduction(:,2),reduction(:,3), '.','MarkerEdgeColor',[188 190 192]/255)
%     hold on
%     axis square
%     xlabel('UMAP1')
%     ylabel('UMAP2')
%     zlabel('UMAP3')
%     scatter3(reduction(data,1),reduction(data,2),reduction(data,3),100, '.','MarkerEdgeColor','r')
%     title(ssname)
% end
% 
% % plot all corner cells in a specific session
% load('corner_metricsR1.mat','C');
% idx = 2;
% reduction = reduction_all{idx};
% ss = ss2plot(idx);
% ssname = S{ss};
% cell2plot = C.cornercell{ss};
% figure
% hold on
% for ii = 1:length(cell2plot)
%     data_raw = neuronIndividualsf{ss}.S(cell2plot(ii),:);
%     data = data_raw > thresh(cell2plot(ii));
%     if ii == 1
%         scatter3(reduction(:,1),reduction(:,2),reduction(:,3), 100, '.','MarkerEdgeColor',[188 190 192]/255)
%     end
%     hold on
%     axis square
%     xlabel('UMAP1')
%     ylabel('UMAP2')
%     zlabel('UMAP3')
%     scatter3(reduction(data,1),reduction(data,2),reduction(data,3),100, '.','MarkerEdgeColor','r')
%     title(ssname)
% end

%% DECODE GEOMETRY  
%% geometry decoding
%% Decoding of environmental geometry
load('F:\analysis_folders.mat','expA')
datapath = expA;
%decode at the center
for n = 1:length(datapath)
    cd(datapath{n})
    decode_geo_ctr = bayes_decoder_geometry(datapath{n});
    save('decode_geometry.mat','decode_geo_ctr','-v7.3')
end
%decode at the corner
for n = 1:length(datapath)
    cd(datapath{n})
    decode_geo_cor = bayes_decoder_geometry(datapath{n},'cor');
    save('decode_geometry.mat','decode_geo_cor','-append')
end
%shuffle condition
for n = 1:length(datapath)
    cd(datapath{n})
    decode_geo_ctr_shuff = bayes_decoder_geometry_shuffle(datapath{n});
    save('decode_geometry.mat','decode_geo_ctr_shuff','-append')
end
%shuffle condition
for n = 1:length(datapath)
    cd(datapath{n})
    decode_geo_cor_shuff = bayes_decoder_geometry_shuffle(datapath{n},'cor');
    save('decode_geometry.mat','decode_geo_cor_shuff','-append')
end

%combine data for quantification
decode_acc = struct;
decode_acc.ctr = [];
decode_acc.ctr_shuff = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('decode_geometry.mat','decode_geo_ctr','decode_geo_ctr_shuff')
    decode_acc_ea = mean(cellfun(@(x) x.accuracy, decode_geo_ctr));
    decode_acc.ctr = [decode_acc.ctr;decode_acc_ea];
    decode_acc_ea_shuff = mean(cellfun(@mean, decode_geo_ctr_shuff));
    decode_acc.ctr_shuff = [decode_acc.ctr_shuff;decode_acc_ea_shuff];
end
decode_acc.cor = [];
decode_acc.cor_shuff = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('decode_geometry.mat','decode_geo_cor','decode_geo_cor_shuff')
    decode_acc_ea = mean(cellfun(@(x) x.accuracy, decode_geo_cor));
    decode_acc.cor = [decode_acc.cor;decode_acc_ea];
    decode_acc_ea_shuff = mean(cellfun(@mean, decode_geo_cor_shuff));
    decode_acc.cor_shuff = [decode_acc.cor_shuff;decode_acc_ea_shuff];
end
save('F:\Results_experimentA\geometry_decoding.mat','decode_acc')

% plot exaple traces
load('F:\subiculum_mice_Aug2022\M4115a\decode_geometry.mat')
figure
subplot(2,1,1)
ss1 = 8;
plot(decode_geo_ctr{ss1}.position_ds)
hold on
plot(decode_geo_ctr{ss1}.decodePosition)
subplot(2,1,2)
ss2 = 1;
plot(decode_geo_cor{ss2}.position_ds([6:20,43:51]))
hold on
plot(decode_geo_cor{ss2}.decodePosition([6:20,43:51]))

%% Reviewer 2
%% Were any cells identified as conjunctive corner and place cells?
load('F:\analysis_folders.mat','expA')
datapath = expA;
%% Distribution of corner cell stability
session_name = {'triangle','square','hex'}; %expA
prop_pccell = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat', 'placecell')
    load('corner_metricsR1.mat', 'C')
    load('env_geometry.mat','S')
    idx = ismember(S, session_name);
    placecells = placecell(idx);
    cornercells = C.cornercell(idx);
    cornerplacecell = cellfun(@(x,y) intersect(x,y), placecells, cornercells, 'uni', 0);
    prop_pccell = [prop_pccell; mean(cellfun(@(x,y) numel(x)/numel(y), cornerplacecell, cornercells))];
end

%% Estimate the image size
imagesize = [];
for n = 1:length(datapath)
    cd(datapath{n})
    findneuron = dir('*_Neuron.mat');
    load (findneuron.name); % And load behav if under different names
    [ysize,xsize] = size(neuron.Cn);
    imagesize = [imagesize;[xsize,ysize]];
end
save('F:\Results_revision\imagesize.mat','imagesize')

%% ARCHIVE
%% archived code
%% compute the max pkfr and min pkfr ratio for all the fields within each corner cell
% load('F:\analysis_folders.mat','expA')
% datapath = expA;
% %find the pkfr for each detected field from a corner cell
% for k = 1:length(datapath)
%     cd(datapath{k})
%     [pkfrcc,pkfrccx] = find_pkfr_eachfield;
%     save('corner_metrics.mat','pkfrcc','pkfrccx','-append')
% end
% session_name = {'triangle','square','hex'}; %expA
% %for each session
% dRate = {};
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('corner_metrics.mat','pkfrcc')
%     load('env_geometry.mat','S')
%     for ii = 1:length(session_name)
%         idx = strcmp(S, session_name{ii});
%         pkfr = pkfrcc(idx);
%         for jj = 1:length(pkfr)
%         dpkfr = cellfun(@(x) max(x)/min(x), pkfr{jj});
%         dRate{(k-1)*2+jj,ii} = dpkfr;
%         end
%     end
% end
% data = cellfun(@mean, dRate);
% %for each mouse
% dRate = cell(length(datapath),length(session_name));
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('corner_metrics.mat','pkfrcc')
%     load('env_geometry.mat','S')
%     for ii = 1:length(session_name)
%         idx = strcmp(S, session_name{ii});
%         pkfr = pkfrcc(idx);
%         pkfr = [pkfr{:}];
%         dpkfr = cellfun(@(x) max(x)/min(x), pkfr);
%         dRate{k,ii} = dpkfr;
%     end
% end
% data = cellfun(@mean, dRate);
% %use GraphPad to plot
%% compare the average pkfr of cornercellx across different environment (no effect)
% session_name = {'circle','triangle','square','hex'}; %expA
% dRate = cell(length(datapath),length(session_name));
% for k = 1:length(datapath)
%     cd(datapath{k})
%     load('corner_metrics.mat','pkfrccx')
%     load('env_geometry.mat','S')
%     for ii = 1:length(session_name)
%         idx = strcmp(S, session_name{ii});
%         pkfr = pkfrccx(idx);
%         pkfr = [pkfr{:}];
%         dpkfr = cellfun(@(x) mean(x), pkfr);
%         dRate{k,ii} = dpkfr;
%     end
% end
% data = cellfun(@mean, dRate);
% %use GraphPad to plot
%% STABILITY (within session and xsession using sessionly defined corner cells)
% within session stability comparison
% session_name = {'triangle','square','hex'}; %expA
% stability = cell(1,length(session_name));
% for k = 1:length(session_name)
%     stab = [];
%     for n = 1:length(datapath)
%         cd(datapath{n})
%         load('env_geometry.mat','S')
%         load('corner_metrics.mat','C')
%         load('spatial_metrics.mat', 'map_stability')
%         idx = strcmp(S, session_name{k});
%         map_stability = map_stability(idx);
%         cornercell = C.cornercell(idx);
%         allcell = 1:length(map_stability{1,1});
%         ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
%         ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
%         %corner cell
%         stabcc = cellfun(@(x,y) x(y), map_stability, cornercell, 'uni', false);
%         stabcc = cellfun(@(x) nanmean(x), stabcc, 'uni', 0);
%         stabcc = [stabcc{:}];
%         %non-corner cell
%         stabnc = cellfun(@(x,y) x(y), map_stability, ncc, 'uni', false);
%         stabnc = cellfun(@(x) nanmean(x), stabnc, 'uni', 0);
%         stabnc = [stabnc{:}];
%         stab = [stab;[stabcc',stabnc']];
%     end
%     stability{k} = stab;
% end
% % xsession stablity comparison
% session_name = {'triangle','square','hex'}; %expA
% stabilityx = [];
% for n = 1:length(datapath)
%     cd(datapath{n})
%     load('env_geometry.mat','S')
%     load('corner_metrics.mat','C')
%     load('spatial_metrics.mat', 'map_stabilityx')
%     cornercell = cell(1,length(session_name));
%     for k = 1:length(session_name)
%         idx = strcmp(S, session_name{k});
%         ccell = C.cornercell(idx);
%         cornercell{k} = unique(vertcat(ccell{:}));
%     end
%     allcell = 1:length(map_stabilityx{1,1});
%     ind = cellfun(@(x) ~ismember(allcell, x), cornercell, 'uni',0);
%     ncc = cellfun(@(x) allcell(x), ind, 'uni', 0);
%     %corner cell
%     stabcc = cellfun(@(x,y) x(y), map_stabilityx, cornercell, 'uni', false);
%     stabcc = cellfun(@(x) nanmean(x), stabcc);
%     %non-corner cell
%     stabnc = cellfun(@(x,y) x(y), map_stabilityx, ncc, 'uni', false);
%     stabnc = cellfun(@(x) nanmean(x), stabnc);
%     stabilityx = [stabilityx;[stabcc,stabnc]];
% end
