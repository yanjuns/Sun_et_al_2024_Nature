%% Analysis of Egocentric corner coding in the subiculum using LN model
load('F:\analysis_folders.mat','expC')
datapath = expC;
%% obtain egeocentric wall bearing/distance and center bearing/distance
for n = 1:length(datapath)
    cd(datapath{n})
    get_egobearing_wall_center
end
%% compute the correlation between variables
egocorr = cell(1,4);
for n = 1:length(datapath)
    cd(datapath{n})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    neuronIndiv = neuronIndividualsf(1:4);
    egocorr{1} = [egocorr{1}; cellfun(@(x) corr(sin(x.eb), sin(x.wallb),'rows','complete'), neuronIndiv)];
    egocorr{2} = [egocorr{2}; cellfun(@(x) corr(sin(x.eb), sin(x.ctrb),'rows','complete'), neuronIndiv)];
    egocorr{3} = [egocorr{3}; cellfun(@(x) corr(x.edist, x.walldist,'rows','complete'), neuronIndiv)];
    egocorr{4} = [egocorr{4}; cellfun(@(x) corr(x.edist, x.ctrdist,'rows','complete'), neuronIndiv)];
end
save('F:\Results_revision\egocorr_expC.mat','egocorr')
%% Run LN models
% the updated code is in package called 'ln_model_ego_bearingdist_variables'
% the most up to date version is on Sherlock and also copied to local. 
%% Plot result of egocentric bearing and egocentric distance to nearest corners
%combine all cells together
ln_ccav = [];
ln_cvex = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminunccdist.mat','lnp_model')
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    ccav = lnp_model{1};
    ln_ccav = [ln_ccav;ccav];
    cvex = lnp_model{3};
    ln_cvex = [ln_cvex;cvex];
end
save('F:\Results_experimentC\LN_resultsR1_cdist.mat','ln_ccav','ln_cvex')

% for each individual mice
LN = struct;
LNP = {};
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminunccdist.mat','lnp_model')
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    histall = cellfun(@(x) histcounts(x,[-0.5:1:15.5])./length(x),...
        lnp_model,'uni',0);
    for ii = 1:size(histall,2)
        LNP{ii}(n,:) = histall{ii};
    end
end
LN.ccav = LNP{1};
LN.cvex = LNP{3};
save('F:\Results_experimentC\LN_resultsR1_cdist.mat','LN','-append')

%plot corner distance and bearing
edge = [0.5:1:15.5];
figure
subplot(2,2,1)
histogram(ln_ccav,edge)
xlim([0.5,15.5])
xticks([1:15])
set(gca,'XTickLabel',{'DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
% ylim([0,50])
ylabel('num of total corner cell')
subplot(2,2,2)
histogram(ln_cvex,edge,'FaceColor',[147 38 103]/255)
xlim([0.5,15.5])
xticks([1:15])
set(gca,'XTickLabel',{'DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
% ylim([0,50])
ylabel('num of total corner cell')
subplot(2,2,3)
position_O = 1:1:15;
boxplot(LN.ccav(:,2:end),'colors',[28 117 188]/255,'positions',position_O,'Symbol','+');
ylim([-0.01 0.17])
set(gca,'XTickLabel',{'DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')
subplot(2,2,4)
boxplot(LN.cvex(:,2:end),'colors',[147 38 103]/255,'positions',position_O,'Symbol','+');
ylim([-0.01 0.17])
set(gca,'XTickLabel',{'DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')

%% DEFINE EGOCENTRIC CORNER CELL USING GLM RESULTS
%% compare fitting between egeocentric corner, boundary and center to define egocentric corner cells
load('F:\analysis_folders.mat','expC')
datapath = expC;
for n = 1:length(datapath)
    cd(datapath{n})
    %load LN data with different models
    load('lnp_whole_fminuncpbin4R1.mat','fit_results','lnp_model')
    fit_resultsp = fit_results;
    lnp_modelp = lnp_model;
    load('lnp_whole_fminunccdist.mat','fit_results','lnp_model')
    fit_resultsc = fit_results;
    lnp_modelc = lnp_model;
    load('lnp_whole_fminuncwdist.mat','fit_results','lnp_model')
    fit_resultsw = fit_results;
    lnp_modelw = lnp_model;
    load('lnp_whole_fminuncctrdist.mat','fit_results','lnp_model')
    fit_resultsctr = fit_results;
    lnp_modelctr = lnp_model;
    %define egocentric corner cells
    ecornercell = cell(1,length(lnp_modelc));
    for ss = 1:length(lnp_modelc)
        %select all the cells that has corner bearing modulation
        ccell = find(lnp_modelc{ss} == 1 | lnp_modelc{ss} == 3 |lnp_modelc{ss} == 4 |...
            lnp_modelc{ss} == 5|lnp_modelc{ss} == 8|lnp_modelc{ss} == 10|lnp_modelc{ss} == 11|lnp_modelc{ss} == 15);
        %Make sure the center-bearing cell has a stronger tuning to boundary and
        %environmental center, see LaChance 2019, Science
        ecornercell_ea1 = [];
        for ii = 1:length(ccell)
            if isnan(lnp_modelw{ss}(ccell(ii)))
                ecornercell_ea1 = [ecornercell_ea1;ccell(ii)];
            else
                wfit = fit_resultsw{ss}{ccell(ii)}{lnp_modelw{ss}(ccell(ii))}(:,2);
                cfit = fit_resultsc{ss}{ccell(ii)}{lnp_modelc{ss}(ccell(ii))}(:,2);
                p = signrank(cfit, wfit, "tail","right");
                if p < 0.05
                    ecornercell_ea1 = [ecornercell_ea1;ccell(ii)];
                end
            end
        end
        ecornercell_ea2 = [];
        for ii = 1:length(ecornercell_ea1)
            if isnan(lnp_modelctr{ss}(ecornercell_ea1(ii)))
                ecornercell_ea2 = [ecornercell_ea2;ecornercell_ea1(ii)];
            else
                cfit = fit_resultsc{ss}{ecornercell_ea1(ii)}{lnp_modelc{ss}(ecornercell_ea1(ii))}(:,2);
                ctrfit = fit_resultsctr{ss}{ecornercell_ea1(ii)}{lnp_modelctr{ss}(ecornercell_ea1(ii))}(:,2);
                p = signrank(cfit, ctrfit, "tail","right");
                if p < 0.05
                    ecornercell_ea2 = [ecornercell_ea2;ecornercell_ea1(ii)];
                end
            end
        end
        %Make sure the center-bearing tuning was not an artifact of conjunctive
        %head direction and position (PH) or (PHS), see LaChance 2019, Science
        ecornercell_ea3 = [];
        for ii = 1:length(ecornercell_ea2)
            if isnan(lnp_modelp{ss}(ecornercell_ea2(ii)))
                ecornercell_ea3 = [ecornercell_ea3;ecornercell_ea2(ii)];
            elseif lnp_modelp{ss}(ecornercell_ea2(ii)) == 1 || ...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 3 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 4 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 5 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 8 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 10 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 11 ||...
                    lnp_modelp{ss}(ecornercell_ea2(ii)) == 15
                ecornercell_ea3 = [ecornercell_ea3;ecornercell_ea2(ii)];
            else
                cfit = fit_resultsc{ss}{ecornercell_ea2(ii)}{lnp_modelc{ss}(ecornercell_ea2(ii))}(:,2);
                pfit = fit_resultsp{ss}{ecornercell_ea2(ii)}{lnp_modelp{ss}(ecornercell_ea2(ii))}(:,2);
                p = signrank(cfit, pfit, "tail","right");
                if p < 0.05
                    ecornercell_ea3 = [ecornercell_ea3;ecornercell_ea2(ii)];
                end
            end
        end
        ecornercell{ss} = ecornercell_ea3;
    end
    save("corner_metricsR1.mat",'ecornercell','-append')
end

%% plot specific egocentric corner neurons
load ('neuronIndividualsf.mat');
load ('behavIndividualsf.mat');
load ('thresh.mat');
load ('firingrateAll.mat','ratemap')
load ('corner_metricsR1.mat', 'ecornercell')
% load ('lnp_whole_fminunccdist.mat', 'lnp_model')
% cell2plot = find(lnp_model{1} == 1 | lnp_model{1} == 3 |lnp_model{1} == 4 |...
%     lnp_model{1} == 5|lnp_model{1} == 8|lnp_model{1} == 10|lnp_model{1} == 11|lnp_model{1} == 15);
cell2plot = ecornercell{3};
plot_eb_map_longitudinal(neuronIndividualsf,thresh,cell2plot,'eb')
plot_eb_map_longitudinal(neuronIndividualsf,thresh,cell2plot,'ctrb')
plot_eb_map_longitudinal(neuronIndividualsf,thresh,cell2plot,'wallb')
plot_eb_map_longitudinal(neuronIndividualsf,thresh,cell2plot,'hd')
plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,ratemap,cell2plot,thresh,'S','auto'); %1:size(neuronIndividualsf{1}.C,1);

% Print figure into vector eps files
filepath = pwd;
neuronSelected = 416;
fig = openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected), '_ebmap.fig']));
print (fig, '-vector', '-depsc', fullfile(filepath,'RatemapFigures',['Cell',num2str(neuronSelected),'_hdmap.eps']));

%% percentage of ecornercell
% proportion of total egeocentric bearing modulated cells
prop_ecell = cell(1,2);
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminunccdist.mat','lnp_model')
    load('thresh.mat','thresh')
    ss1 = 1; ss2 = 3;
    cell2plot1 = find(lnp_model{1} == 1 | lnp_model{1} == 3 |lnp_model{1} == 4 |...
        lnp_model{1} == 5|lnp_model{1} == 8|lnp_model{1} == 10|lnp_model{1} == 11|lnp_model{1} == 15);
    cell2plot2 = find(lnp_model{ss2} == 1 | lnp_model{ss2} == 3 |lnp_model{ss2} == 4 |...
        lnp_model{ss2} == 5|lnp_model{ss2} == 8|lnp_model{ss2} == 10|lnp_model{ss2} == 11|lnp_model{ss2} == 15);
    prop_ecell{1} = [prop_ecell{1};length(cell2plot1)/length(thresh)];
    prop_ecell{2} = [prop_ecell{2};length(cell2plot2)/length(thresh)];
end
save('F:\Results_revision\prop_ego_cells.mat','prop_ecell')

prop_ecornercell = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load("corner_metricsR1.mat",'ecornercell')
    load('thresh.mat','thresh')
    prop_ecornercell = [prop_ecornercell; cellfun(@numel, ecornercell)/numel(thresh)];
end
save('F:\Results_revision\prop_ego_cells.mat','prop_ecornercell','-append')

%% overlap between allocentric and egocentric corner cells
load('F:\analysis_folders.mat','expC')
datapath = expC;
ae_overlap = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('corner_metricsR1.mat', 'C', 'ecornercell')
    cornercell = C.cornercell(1:4);
    overlap_ea = cellfun(@(x,y) numel(intersect(x,y)), cornercell, ecornercell);
    ae_overlap = [ae_overlap;overlap_ea];
end
save('F:\Results_revision\allo_ego_overlap_expC.mat','ae_overlap')

%% Experiment D
%% Experiment D
%% Experiment D
load('F:\analysis_folders.mat','expD')
datapath = expD(1:end-1);
for n = 1:length(datapath)
    cd(datapath{n})
    get_behav_directions
end
for n = 1:length(datapath)
    cd(datapath{n})
    get_egobearing_wall_center
end
%% compute the correlation between variables
egocorr = cell(1,4);
for n = 1:length(datapath)
    cd(datapath{n})
    load('neuronIndividualsf.mat','neuronIndividualsf')
    load("env_geometry.mat",'S')
    idx1 = ismember(S, 'rightTri');
    neuronIndiv1 = neuronIndividualsf(idx1);
    idx2 = ismember(S, 'trapezoid');
    neuronIndiv2 = neuronIndividualsf(idx2);
    egocorr{1} = [egocorr{1}; cellfun(@(x) corr(sin(x.eb), sin(x.wallb),'rows','complete'), neuronIndiv1),...
        cellfun(@(x) corr(x.eb, x.wallb,'rows','complete'), neuronIndiv2)];
    egocorr{2} = [egocorr{2}; cellfun(@(x) corr(sin(x.eb), sin(x.ctrb),'rows','complete'), neuronIndiv1),...
        cellfun(@(x) corr(x.eb, x.ctrb,'rows','complete'), neuronIndiv2)];
    egocorr{3} = [egocorr{3}; cellfun(@(x) corr(x.edist, x.walldist,'rows','complete'), neuronIndiv1),...
        cellfun(@(x) corr(x.edist, x.walldist,'rows','complete'), neuronIndiv2)];
    egocorr{4} = [egocorr{4}; cellfun(@(x) corr(x.edist, x.ctrdist,'rows','complete'), neuronIndiv1),...
        cellfun(@(x) corr(x.edist, x.ctrdist,'rows','complete'), neuronIndiv2)];
end
save('F:\Results_revision\egocorr_expD.mat','egocorr')

%% Plot result of egocentric bearing and egocentric distance to nearest corners
%combine all cells together
ln_ccav = [];
ln_cvex = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminuncpbin4R1.mat','lnp_model')
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    ccav = lnp_model{1};
    ln_ccav = [ln_ccav;ccav];
    cvex = lnp_model{2};
    ln_cvex = [ln_cvex;cvex];
end
save('F:\Results_experimentD\LN_resultsR1_cdist.mat','ln_ccav','ln_cvex')

% for each individual mice
LN = struct;
LNP = {};
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_fminuncpbin4R1.mat','lnp_model')
    lnp_model = cellfun(@convertNaN2zero,lnp_model,'uni',0);
    histall = cellfun(@(x) histcounts(x,[-0.5:1:15.5])./length(x),...
        lnp_model,'uni',0);
    for ii = 1:size(histall,2)
        LNP{ii}(n,:) = histall{ii};
    end
end
LN.ccav = LNP{1};
LN.cvex = LNP{2};
save('F:\Results_experimentD\LN_resultsR1_cdist.mat','LN','-append')

%plot corner distance and bearing
edge = [-0.5:1:15.5];
figure
subplot(2,2,1)
histogram(ln_ccav,edge)
xlim([-0.5,15.5])
xticks([0:15])
set(gca,'XTickLabel',{'NaN','DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
% ylim([0,50])
ylabel('num of total corner cell')
subplot(2,2,2)
histogram(ln_cvex,edge,'FaceColor',[147 38 103]/255)
xlim([-0.5,15.5])
xticks([0:15])
set(gca,'XTickLabel',{'NaN','DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
% ylim([0,50])
ylabel('num of total corner cell')
subplot(2,2,3)
position_O = 0:1:15;
boxplot(LN.ccav,'colors',[28 117 188]/255,'positions',position_O,'Symbol','+');
ylim([-0.1 1])
set(gca,'XTickLabel',{'NaN','DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')
subplot(2,2,4)
boxplot(LN.cvex,'colors',[147 38 103]/255,'positions',position_O,'Symbol','+');
ylim([-0.1 1])
set(gca,'XTickLabel',{'NaN','DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'});
set(gca, 'XDir','reverse')
ylabel('prop of corner cell')

%% distribution of egocentric corner cells ***
% plot the group distribution of characterized egocentric corner cells
load('F:\analysis_folders.mat','expD')
datapath = expD(1:end-1);
ecornercell_all = [];
for n = 1:length(datapath)
    cd(datapath{n})
    load("corner_metricsR1.mat",'ecornercell')
    load('lnp_whole_fminunccdist.mat','lnp_model')
    ecornercell_all = [ecornercell_all; lnp_model{1}(ecornercell{1})];
end
load('F:\analysis_folders.mat','expC')
datapath = expC;
for n = 1:length(datapath)
    cd(datapath{n})
    load("corner_metricsR1.mat",'ecornercell')
    load('lnp_whole_fminunccdist.mat','lnp_model')
    ecornercell_all = [ecornercell_all; lnp_model{1}(ecornercell{1});...
        lnp_model{2}(ecornercell{2});lnp_model{4}(ecornercell{4})];
end
group = {'DHSE','DHS','DHE','DSE','HSE','DH','DS','DE','HS','HE','SE','D','H','S','E'};
N = histcounts(ecornercell_all, [0.5:1:15.5]);
data = N/sum(N);
idx = data ~= 0; 
data(data ==0) = [];
figure
pie(data)
labels = group(idx);
legend(labels);

%% ARCHIVE 
%% Archived code
%% Andy's code
% rx = root.x;      % x position in cm or pixels
% ry = root.y;      % y position in cm or pixels
% md = root.md;     % movement (or head) direction, in radians
% ts = root.ts;     % time stamps (seconds)
% spk = root.spike
% 
% root.x = neuron.pos(:,1);
% root.y = neuron.pos(:,2);
% root.md = neuron.hd;
% root.ts = neuron.time;
% root.spike = double(idx');
% 
% out = EgocentricRatemap(root);
% plotEBC(root, out, 1)
