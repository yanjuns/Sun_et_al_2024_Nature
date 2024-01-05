%% load datapath
load('F:\analysis_folders.mat','expA')
datapath = expA;
%% Plot spatial rate maps for a initial examination of cells
%% Calculate spatial map for plotting purpose
for ii = 1:length(datapath)
    cd(datapath{ii})
    generate_spatial_ratemap
end
%% Visually inspect individual spatial maps longitudinally
load ('neuronIndividualsf.mat');
load ('behavIndividualsf.mat');
load ('thresh.mat');
load ('firingrateAll.mat','ratemap')
load ('corner_metricsR1.mat','C')
cell2plot = union(C.cornercell{2},C.cornercell{4});
plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,ratemap,cell2plot,thresh,'S','auto'); %1:size(neuronIndividualsf{1}.C,1);
% open a group of plots;
filepath = pwd;
neuronSelected = cell2plot;
for ii = 1:length(neuronSelected)
    if exist(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap_auto.fig']), 'file')
        openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap_auto.fig']));
    end
end
% Print figure into vector eps files
filepath = pwd;
neuronSelected =58;
fig = openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected), '_ratemap_auto.fig']));
print (fig, '-painters', '-depsc', fullfile(filepath,'RatemapFigures',['Cell',num2str(neuronSelected),'_ratemap_auto.eps']));
% plot the ratemap by varying the threshold, as needed
edit plot_rate_map_varythresh
%% Visually inspect rate maps for each session
load('firingrateAll.mat','ratemap')
load('corner_metricsR1.mat','C')
session2plot =1;
cell2plot = C.cornercell{session2plot};
plot_ratemap(ratemap{session2plot},cell2plot, C.cscore{session2plot}(cell2plot))

%% find specific neurons under the imaging window
load('neuronIndividualsf.mat','neuronIndividualsf')
load('corner_metrics.mat','C')
load('border_metrics.mat','bcell')
neuronSelected1 = bcell{1};
neuronSelected2 = [C.cornercell{1}];
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

%% FUNCTIONAL CELL TYPES
%% Identify different functional cell types
%% Place cell, Corner cell, Boundary cell
%% Identify place cells
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
        x.meanfr > quantile(x.meanfr,0.05)), spatial_metrics, 'uni',0);
    placecell = cellfun(@(x,y) x.neuron(y), spatial_metrics, idx, 'uni', 0);
    save('spatial_metrics.mat','placecell','-append')
end
%% Identify corner cells
%mannually identify the location of corners for each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    [S,N,env_coor,env_coori] = identify_env_geometry;
    save('env_geometry.mat','S','N','env_coor','env_coori')
end

%determin corner cells in each session
for ii = 1:length(datapath)
    cd(datapath{ii})
    %determine corner cell for each session by spike shuffling.
    C = identify_corner_cell;
    save('corner_metrics.mat','C','-v7.3')
end

%determin corner cells using multiple sessions
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('corner_metricsR1.mat','C')
    load('env_geometry.mat','S') %load all session name
    session_name = {'triangle','square','hex'}; %expA
%     session_name = {'rightTri','trapezoid'}; %selected sessions using to define corner cells
    ccellx = cell(1,length(session_name));
    %for experiment A
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
load('F:\analysis_folders.mat','expD')
datapath = expD;
% session_name = {'circle','triangle','square','hex'}; %expA
session_name = {'rightTri','trapezoid','rightTrim1'}; %expD
prop_cornercell = cell(length(datapath),length(session_name));
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metrics.mat','C')
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
save('F:\Results_experimentD\prop_cornercell.mat','prop_cornercell','prop_cornercell_ea')

%% Identify border cells
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

%intersect border score with place cells to get border cell
for n = 1:length(datapath)
    cd(datapath{n})
    load('spatial_metrics.mat','placecell')
    load('border_metrics.mat','bscore','bfields')
    bc1 = cellfun(@(x) find(x > 0.6), bscore, 'uni',0);
    coverage = cellfun(@(x) {x.coverage}, bfields, 'uni',0);
    idxc = cell(1,length(coverage));
    for ii = 1:length(coverage)
        idxc{ii} = cellfun(@(x) x.perimeter.max>0.7, coverage{ii});
    end
    bc2 = cellfun(@(x) find(x == 1), idxc, 'uni', 0);
    bcell_raw = cellfun(@(x,y) intersect(x,y), bc1, bc2, 'uni', 0);
    bcell = cellfun(@(x,y) intersect(x,y), bcell_raw, placecell, 'uni',0);
    save('border_metrics.mat','bcell','bcell_raw','-append')
end

%% Archived code
%% Identify corner cells (first version, not perfect)
% load('D:\Miniscope_analysis_subiculum_shuttle\datapath.mat')
% for n = 1:length(datapath)
%     cd(datapath{n})
%     load('firingrateAll', 'firingrateAll')
%     load('PlaceCellsLR.mat','placecellsLRm')
%     ratemap = firingrateAll;
%     cmetrics = cell(size(ratemap,1),size(ratemap,2));
%     for ii = 1:length(ratemap)
%             cmetrics{ii} = cellfun(@identify_corner_coding, ratemap{ii});
%     end
%     coridx = cellfun(@(x) find([x(:).acpeaks] >=7 & [x(:).acpeaks] <=9 & [x(:).atcorner] == 1),...
%         cmetrics,'uni',false);
%     placecell = cellfun(@(x,y) union(x,y), placecellsLRm(1,:),placecellsLRm(2,:),'uni',0);
%     placecell = [placecell;placecell];
%     corcell = cellfun(@(x,y) intersect(x,y), coridx, placecell, 'uni',0);
%     save('corner.mat','cmetrics','corcell','-v7.3')
% end