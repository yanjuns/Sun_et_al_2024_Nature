%% LPHN2 EXPERIMENTS
%% Lphn2 KO experiments
%% Load data
load('F:\subiculum_Lphn_Jan2023\datapath.mat','datapath','Ctrl','KO')

%% PLACE CELL
%% Regular place cells
%% Identify place cells for overall environment
for n = 1:length(datapath)
    cd(datapath{n})
    load ('neuronIndividualsf.mat');
    load ('behavIndividualsf.mat');
    load ('thresh.mat')
    spatial_metrics = cell(1, length(neuronIndividualsf));
    parfor ii = 1:length(neuronIndividualsf)
        neuron = neuronIndividualsf{1,ii};
        behav = behavIndividualsf{1,ii};
        [pc,Tinfo,~] = permuting_spikes(neuron,behav,thresh,true);
        spatial_metrics{1,ii} = Tinfo; %defined by bit/spike
    end
    idx = cellfun(@(x) find(x.bitpspike_pcidx == 1 &...
        x.meanfr > 0.1), spatial_metrics, 'uni',0);
    placecell = cellfun(@(x,y) x.neuron(y), spatial_metrics, idx, 'uni', 0);
    if ~exist('spatial_metrics.mat','file')
        save('spatial_metrics.mat','placecell','spatial_metrics')
    else
        save('spatial_metrics.mat','placecell','spatial_metrics','-append')
    end
end
%use a different fr threshold to define place cells
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat','spatial_metrics')
    idx = cellfun(@(x) find(x.bitpspike_pcidx == 1 &...
        x.meanfr > 0.1), spatial_metrics, 'uni',0);
    placecell = cellfun(@(x,y) x.neuron(y), spatial_metrics, idx, 'uni', 0);
    save('spatial_metrics.mat','placecell','-append')
end

%% firing rate and proportion of place cells
pc_metrics.Ctrlmeanfr = [];
pc_metrics.Ctrlprop = [];
pc_metrics.Ctrlbpspike = [];
pc_metrics.Ctrlpkfr = [];
for n = 1:length(Ctrl)
    cd(Ctrl{n})
    load('spatial_metrics.mat','spatial_metrics','placecell')
    load('sessions.mat','sessions')
    session_name = {'shuttle'};
    idx = contains(sessions, session_name);
    spatial_metrics = spatial_metrics(idx);
    placecell = placecell(idx);
    %meanfr
    meanfr = cellfun(@(x) nanmean(x.meanfr), spatial_metrics);
    pc_metrics.Ctrlmeanfr = [pc_metrics.Ctrlmeanfr;mean(meanfr,2)];
    %prop of place cells
    pc_metrics.Ctrlprop = [pc_metrics.Ctrlprop;...
        mean(cellfun(@numel, placecell)./size(spatial_metrics{1,1},1))];
    %bitpspike
    bpspike = cellfun(@(x,y) mean(table2array(x(y,'bitpspike'))),...
        spatial_metrics, placecell);
    pc_metrics.Ctrlbpspike = [pc_metrics.Ctrlbpspike;nanmean(bpspike)];
    %peakfr
    pkfr = cellfun(@(x,y) mean(table2array(x(y,'peakfr'))),...
        spatial_metrics, placecell);
    pc_metrics.Ctrlpkfr = [pc_metrics.Ctrlpkfr;nanmean(pkfr)];
end

pc_metrics.KOmeanfr = [];
pc_metrics.KOprop = [];
pc_metrics.KObpspike = [];
pc_metrics.KOpkfr = [];
for n = 1:length(KO)
    cd(KO{n})
    load('spatial_metrics.mat','spatial_metrics','placecell')
    load('sessions.mat','sessions')
    session_name = {'shuttle'};
    idx = contains(sessions, session_name);
    spatial_metrics = spatial_metrics(idx);
    placecell = placecell(idx);
    %meanfr
    meanfr = cellfun(@(x) nanmean(x.meanfr), spatial_metrics);
    pc_metrics.KOmeanfr = [pc_metrics.KOmeanfr;mean(meanfr,2)];
    %prop of place cells
    pc_metrics.KOprop = [pc_metrics.KOprop;...
        mean(cellfun(@numel, placecell)./size(spatial_metrics{1,1},1))];
    %bitpspike
    bpspike = cellfun(@(x,y) mean(table2array(x(y,'bitpspike'))),...
        spatial_metrics, placecell);
    pc_metrics.KObpspike = [pc_metrics.KObpspike;nanmean(bpspike)];
    %peakfr
    pkfr = cellfun(@(x,y) mean(table2array(x(y,'peakfr'))),...
        spatial_metrics, placecell);
    pc_metrics.KOpkfr = [pc_metrics.KOpkfr;nanmean(pkfr)];
end
save('F:\subiculum_Lphn_Jan2023\Results\pc_metrics.mat','pc_metrics')

%% Compare anatomical density of place cell in Ctrl and KO
% calculate anatomical density for each mouse
for n = 1:length(datapath)
    cd(datapath{n})
    if contains(datapath{n},'M4136F')
        session_name = {'square'};
    else
        session_name = {'circle','square','triangle','hex','shuttle'};
    end
    load('sessions.mat','sessions')
    load('spatial_metrics.mat','placecell')
    idx = contains(sessions, session_name);
    ssLength = 1:length(sessions);
    ss = ssLength(idx);
    %load centriod locations of neurons
    findneuron = dir('*_Neuron_filter.mat');
    load (findneuron.name);
    centriod = neuron.centroid;
    %preallocate results
    pA = cell(1,length(ss));
    pAs = cell(1,length(ss));
    for ii = 1:length(ss)
        pc = placecell{ss(ii)};
        [A,As,pA{ii},pAs{ii}] = anat_density_map(centriod,pc,20);
    end
    pAs_avg = mean(cat(3, pAs{:}),3);
    if ~exist('anatomy.mat','file')
        save('anatomy.mat','A','As','pA','pAs','pAs_avg')
    else
        save('anatomy.mat','A','As','pA','pAs','pAs_avg','-append')
    end
end

% get a uniform size to align all the data
msize = [];
for k = 1:length(datapath)
    cd(datapath{k});
    m = matfile('anatomy.mat');
    ms = size(m, 'As');
    msize = [msize;ms];
end
targetSize = max(msize);
% gather data for all mice to plot
anat = struct;
%Ctrlmice
pd1 = []; pd2 = []; % place cell dimension1, place cell dimension2
for k = 1:length(Ctrl)
    cd(Ctrl{k});
    load('anatomy.mat','pAs_avg');
    pd = imresize(pAs_avg, targetSize); %interpolate all matrix to same size
    pd1 = [pd1;sum(pd)]; pd2 = [pd2;sum(pd,2)'];
end
anat.Ctrlpd1 = pd1;
anat.Ctrlpd2 = pd2;

%KOmice
pd1 = []; pd2 = []; % place cell dimension1, place cell dimension2
for k = 1:length(KO)
    cd(KO{k});
    load('anatomy.mat','pAs_avg');
    pd = imresize(pAs_avg, targetSize); %interpolate all matrix to same size
    pd1 = [pd1;sum(pd)]; pd2 = [pd2;sum(pd,2)'];
end
anat.KOpd1 = pd1;
anat.KOpd2 = pd2;

%normalize the data for plot and quantification
anat.Ctrl_px = anat.Ctrlpd1./sum(anat.Ctrlpd1,2);
anat.Ctrl_py = anat.Ctrlpd2./sum(anat.Ctrlpd2,2);
anat.KO_px = anat.KOpd1./sum(anat.KOpd1,2);
anat.KO_py = anat.KOpd2./sum(anat.KOpd2,2);
save('F:\subiculum_Lphn_Jan2023\Results\pc_anatomy.mat','anat')

figure
subplot(1,2,1)
stdshade(anat.Ctrl_px,0.4,[65,64,66]/255,[],3)
hold on
stdshade(anat.KO_px,0.4,[66,10,104]/255,[],3)
xlabel('Proximal to Distal (Bins, 1bin = 18um)')
ylabel('Proportion of total density')
axis square
ylim([0, 0.09])
for ii = 1:size(anat.Ctrl_px,2)
    [h,p] = ttest2(anat.Ctrl_px(:,ii),anat.KO_px(:,ii));
    if p < 0.05 && p >0.01
        text(ii,0.08,'*','FontSize',15)
    elseif p <= 0.01
        text(ii,0.08,'**','FontSize',15)
    end
end
subplot(1,2,2)
stdshade(anat.Ctrl_py,0.4,[65,64,66]/255,[],3)
hold on
stdshade(anat.KO_py,0.4,[66,10,104]/255,[],3)
xlabel('Caudal to Rostral (Bins, 1bin = 18um)')
ylabel('Proportion of total density')
axis square
ylim([0, 0.11])
for ii = 1:size(anat.Ctrl_py,1)
    [h,p] = ttest2(anat.Ctrl_py(:,ii),anat.KO_py(:,ii));
    if p < 0.05 && p >0.01
        text(ii,0.12,'*','FontSize',15)
    elseif p <= 0.01
        text(ii,0.12,'**','FontSize',15)
    end
end

figure
subplot(2,1,1)
imagesc(anat.Ctrl_px)
xlabel('proximal to distal subiculum')
ylabel('Ctrl mice')
caxis([0,0.09])
subplot(2,1,2)
imagesc(anat.KO_px)
xlabel('proximal to distal subiculum')
ylabel('Lphn2 cKO mice')
caxis([0,0.09])

%% plot each example for place cell anatomy
load('neuronIndividualsf.mat','neuronIndividualsf')
load('spatial_metrics.mat','placecell')
load('anatomy.mat')
ss2plot = 9;
neuronSelected = placecell{ss2plot};
ncontour = neuronIndividualsf{1}.Coor;
figure
subplot(2,2,1)
hold on
for ii = 1:length(ncontour)
    plot(ncontour{ii}(1,:),ncontour{ii}(2,:),'Color',[128,130,133]/255,'LineWidth', 0.5)
    set(gca, 'YDir','reverse')
    axis image
end
subplot(2,2,2)
imagesc(As)
axis image
subplot(2,2,3)
hold on
for ii = 1:length(ncontour)
    plot(ncontour{ii}(1,:),ncontour{ii}(2,:),'Color',[128,130,133]/255,'LineWidth', 0.5)
    set(gca, 'YDir','reverse')
    axis image
end
for jj = 1:length(neuronSelected)
    plot(ncontour{neuronSelected(jj)}(1,:),ncontour{neuronSelected(jj)}(2,:),'Color',[218,28,92]/255, 'LineWidth', 0.5)
end
subplot(2,2,4)
imagesc(pAs{ss2plot})
axis image


%% SHUTTLE BOX EXPERIMENT
%% shuttle box experiment
%% Load data
load('F:\subiculum_Lphn_Jan2023\datapath.mat','datapath','Ctrl','KO')
%% Shuttle box: split the shuttle box data into left and right compartments
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

%% Shuttle box: Identify place cells with left and right matched data
for k = 1:length(datapath)
    cd(datapath{k});
    load('NBindivLR.mat', 'neuronIndivLR', 'behavIndivLR')
    load('thresh.mat');
    niter = 24;
    lengN = length(neuronIndivLR);
    metricsm = cell(niter, lengN);
    parfor n = 1:niter
        [neuronIndivLRm, behavIndivLRm] = match_neuron_behav(neuronIndivLR, behavIndivLR, 'LR');
        for ii = 1:lengN
            neuronL = neuronIndivLRm{1,ii};
            behavL = behavIndivLRm{1,ii};
            [~,Tinfo,~] = permuting_spikes(neuronL,behavL,thresh,true);
            metricsm{n,ii} = table2array(Tinfo); %defined by bit/sec
        end
    end
    %average the value for each iteration
    metricsm_avg = cell(1,lengN);
    for jj = 1:size(metricsm,2)
        metrics_each = squeeze(metricsm(:,jj));
        metrics_each = cat(3, metrics_each{:});
        metricsm_avg{1,jj} = array2table(nanmean(metrics_each,3),...
            'VariableNames',{'neuron','bitpsec','bitpsec_pcidx','bitpspike','bitpspike_pcidx','meanfr','peakfr'});
    end
    save('spatial_metrics.mat', 'metricsm', 'metricsm_avg', '-append');
end
%define place cells
for k = 1:length(datapath)
    cd(datapath{k});
    load('spatial_metrics.mat', 'metricsm_avg')
    placecellsLRm = cellfun(@(x) find(x.meanfr >= 0.1 & x.bitpspike_pcidx >= 0.3), metricsm_avg,...
        'Uniformoutput', false);
    save('spatial_metrics.mat', 'placecellsLRm', '-append')
end

%% Shuttle box: Compare anatomical density of place cell
% calculate anatomical density for each mouse
for n = 1:length(datapath)
    cd(datapath{n})
    load('sessions.mat','sessions')
    load('spatial_metrics.mat','placecellsLRm')
    ss = 1:length(placecellsLRm);
    %load centriod locations of neurons
    findneuron = dir('*_Neuron_filter.mat');
    load (findneuron.name);
    centriod = neuron.centroid;
    %preallocate results
    pA = cell(1,length(ss));
    pAs = cell(1,length(ss));
    for ii = 1:length(ss)
        pc = placecellsLRm{ss(ii)};
        [A,As,pA{ii},pAs{ii}] = anat_density_map(centriod,pc,20);
    end
    pAs_avg = mean(cat(3, pAs{:}),3);
    if ~exist('anatomy_shuttle.mat','file')
        save('anatomy_shuttle.mat','A','As','pA','pAs','pAs_avg')
    else
        save('anatomy_shuttle.mat','A','As','pA','pAs','pAs_avg','-append')
    end
end

% get a uniform size to align all the data
msize = [];
for k = 1:length(datapath)
    cd(datapath{k});
    m = matfile('anatomy_shuttle.mat');
    ms = size(m, 'As');
    msize = [msize;ms];
end
targetSize = max(msize);
% gather data for all mice to plot
anat = struct;
%Ctrlmice
pd1 = []; pd2 = []; % place cell dimension1, place cell dimension2
for k = 1:length(Ctrl)
    cd(Ctrl{k});
    load('anatomy_shuttle.mat','pAs');
    pdx = cellfun(@(x) sum(imresize(x, targetSize)), pAs, 'uni', 0);
    pdx = cat(1,pdx{:});
    pdy = cellfun(@(x) sum(imresize(x, targetSize),2), pAs, 'uni', 0);
    pdy = cat(2,pdy{:});
    pd1 = [pd1;pdx]; pd2 = [pd2;pdy'];
end
anat.Ctrlpd1 = pd1;
anat.Ctrlpd2 = pd2;

%KOmice
pd1 = []; pd2 = []; % place cell dimension1, place cell dimension2
for k = 1:length(KO)
    cd(KO{k});
    load('anatomy_shuttle.mat','pAs');
    pdx = cellfun(@(x) sum(imresize(x, targetSize)), pAs, 'uni', 0);
    pdx = cat(1,pdx{:});
    pdy = cellfun(@(x) sum(imresize(x, targetSize),2), pAs, 'uni', 0);
    pdy = cat(2,pdy{:});
    pd1 = [pd1;pdx]; pd2 = [pd2;pdy'];
end
anat.KOpd1 = pd1;
anat.KOpd2 = pd2;

%normalize the data for plot and quantification
anat.Ctrl_px = anat.Ctrlpd1./sum(anat.Ctrlpd1,2);
anat.Ctrl_py = anat.Ctrlpd2./sum(anat.Ctrlpd2,2);
anat.KO_px = anat.KOpd1./sum(anat.KOpd1,2);
anat.KO_py = anat.KOpd2./sum(anat.KOpd2,2);

% %if needs to average each mouse
% fun = @(block_struct) mean(block_struct.data);
% anat.Ctrl_px = blockproc(anat.Ctrl_px, [4, size(anat.Ctrl_px,2)], fun);
% anat.Ctrl_py = blockproc(anat.Ctrl_py, [4, size(anat.Ctrl_py,2)], fun);
% anat.KO_px = blockproc(anat.KO_px, [4, size(anat.KO_px,2)], fun);
% anat.KO_py = blockproc(anat.KO_py, [4, size(anat.KO_py,2)], fun);

figure
subplot(1,2,1)
stdshade(anat.Ctrl_px,0.5,[0 0.4470 0.7410],[],3)
hold on
stdshade(anat.KO_px,0.5,[0.8500 0.3250 0.0980],[],3)
xlabel('Proximal to Distal (Bins, 1bin = 25um)')
ylabel('Proportion of total density')
axis square
for ii = 1:size(anat.Ctrl_px,2)
    [h,p] = ttest2(anat.Ctrl_px(:,ii),anat.KO_px(:,ii));
    if p < 0.05 && p >0.01
        text(ii,0.08,'*','FontSize',15)
    elseif p <= 0.01
        text(ii,0.08,'**','FontSize',15)
    end
end
subplot(1,2,2)
stdshade(anat.Ctrl_py,0.5,[0 0.4470 0.7410],[],3)
hold on
stdshade(anat.KO_py,0.5,[0.8500 0.3250 0.0980],[],3)
xlabel('Caudal to Rostral (Bins, 1bin = 25um)')
ylabel('Proportion of total density')
axis square
for ii = 1:size(anat.Ctrl_py,2)
    [h,p] = ttest2(anat.Ctrl_py(:,ii),anat.KO_py(:,ii));
    if p < 0.05 && p >0.01
        text(ii,0.12,'*','FontSize',15)
    elseif p <= 0.01
        text(ii,0.12,'**','FontSize',15)
    end
end

