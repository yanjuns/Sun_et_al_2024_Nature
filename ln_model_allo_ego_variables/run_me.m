%% Description of run_me
% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017
% V2(current version): Yanjun Sun, modified for miniscope data, Feb 2022
% with an implementation of fitting using glmnet, which performs faster
% than fminunc. Of note, the two methods tend to give slightly different
% fitting results

% This script is segmented into several parts. First, the data (an
% example cell) is loaded. Then, 15 LN models are fit to the
% cell's spike train. Each model uses information about 
% position, head direction, running speed, theta phase,
% or some combination thereof, to predict a section of the
% spike train. Model fitting and model performance is computed through
% 10-fold cross-validation, and the minimization procedure is carried out
% through fminunc. Next, a forward-search procedure is
% implemented to find the simplest 'best' model describing this spike
% train. Following this, the firing rate tuning curves are computed, and
% these - along with the model-derived response profiles and the model
% performance and results of the selection procedure are plotted.

%% set up global parameters
p = struct;
p.method = 'fminunc'; % glmnet or fminunc
p.dataType = 'whole'; % choose from whole or subset(LR)
p.session_name = {'largesq','largeRT','largeX','largeXm1'}; %expC
p.n2compute = 'cornercell'; % put 'all' if you want to compute all neurons
p.c = false; % if need to combine data from different sessions
% p.c_vector = [1,1,2,2;3,3,4,4]; % specify the sessions need to be combined, e.g [1,1,2,2;3,3,4,4]
% p.c_vector = [1,1,2,2,3,3;4,4,5,5,6,6];
p.pos_bin_size = 2; % cm
p.n_dir_bins = 18;
p.n_speed_bins = 10;
method = p.method;
%% prepare for LNP data
poolobj = parpool;
poolobj = poolobj.NumWorkers;
for n = 1:length(datapath)
    cd(datapath{n})
    load_data;
    %select sessions
    idx = ismember(S, p.session_name);
    neuronIndiv = neuronIndiv(idx);
    behavIndiv = behavIndiv(idx);
    %preallocate results
    lnp_model = cell(size(neuronIndiv,1),size(neuronIndiv,2));
    llh = cell(size(neuronIndiv,1),size(neuronIndiv,2));
    fit_results = cell(size(neuronIndiv,1),size(neuronIndiv,2));
    if strcmp(p.n2compute,'all')
        n2compute = repmat({1:length(thresh)},1,length(p.session_name));
    else
        load('corner_metrics.mat','C2')
        n2compute = C2.cornercell(idx);
    end
    t0 = tic();
    for x = 1:size(neuronIndiv,1)
        for y = 1:size(neuronIndiv,2)
            neuron = neuronIndiv{x,y};
            behav = behavIndiv{x,y};
            cell2use = n2compute{x,y};
            llh_ea = cell(1,length(cell2use));
            lnp_model_ea = NaN(length(cell2use),1);
            testFit = cell(1,length(cell2use));
            %calculate required vars
            [A, modelType,spiketrain_all, dt, xnbins, ynbins, direction, speed, egob]...
                = prep_lnp_data(neuron,behav,thresh,p.pos_bin_size,p.n_dir_bins,p.n_speed_bins);
            %start modeling for each neuron
            parfor ii = 1:length(cell2use)
                jj = cell2use(ii);
                %fit the model
                fprintf('Fitting all linear-nonlinear (LN) models\n')
                [testFit{ii},~,~,~] = fit_all_ln_models(A,...
                    modelType,spiketrain_all(:,jj),dt,[],[],method);
                %find the simplest model that best describes the spike train
                fprintf('Performing forward model selection\n')
                [lnp_model_ea(ii),llh_ea{ii}] = select_best_model(testFit{ii},'correlation');
            end
            fit_results{x,y} = testFit;
            lnp_model{x,y} = lnp_model_ea;
            llh{x,y} = llh_ea;
%             %in case of aborted workers during computation
%             pp = gcp;
%             if pp.NumWorkers < poolobj
%                 delete(gcp('nocreate'));
%                 parpool;
%             end
        end
    end
    t1 = toc(t0)
    %allocate and save data
    save(['lnp_',p.dataType,'_',method,'.mat'], 'lnp_model', 'llh', 'fit_results','p','-v7.3')
end

%% Kiah's code for ploting individual neurons
% Compute the firing-rate tuning curves
load_data;
fprintf('Computing tuning curves\n')
session2plotx = 1;
session2ploty = 1;
cell2plot = 145;
compute_all_tuning_curves
% plot the results
fprintf('Plotting performance and parameters\n')
plot_performance_and_parameters

%% calculate and plot batch results from multiple mice
load('E:\Miniscope_MorphineCPA_LNP\datapath.mat')
LNP = {};
for n = 1:length(datapath)
    cd(datapath{n})
    load('property.mat','drugside','mouse','bt2')
    if any(strcmp(mouse, {'M4051','M4017F','M4024F','M4028','M4029','M922'}))
        load('lnp_LR_fminunc_ea.mat','lnp_model')
        histea = cellfun(@(x) histcounts(x,[0.5:1:7.5])./length(x),...
            lnp_model,'uni',0);
        histall = cell(2,round(size(lnp_model,2)/2));
        for x = 1:size(histea,1)
            k = 1;
            for y = 1:2:size(histea,2)
                histall{x,k} = mean([histea{x,y};histea{x,y+1}]);
                k = k+1;
            end
        end
        histall = histall(:,bt2);
    else
        load('lnp_LR_fminunc.mat','lnp_model')
        lnp_model = lnp_model(:,bt2);
        histall = cellfun(@(x) histcounts(x,[0.5:1:7.5])./length(x),...
            lnp_model,'uni',0);
    end
    if strcmp(drugside, 'R')
        for ii = 1:2
            for jj = 1:2
                LNP{ii,jj}(n,:) = histall{ii,jj};
            end
        end
    else
        histall = flipud(histall);
        for ii = 1:2
            for jj = 1:2
                LNP{ii,jj}(n,:) = histall{ii,jj};
            end
        end
    end
end
LNPf = LNP;
save('LNP_MOCPP2.mat','LNPf')
% LNPg = LNP;
% save('LNP_LR.mat','LNPg','-append')

data = LNPf;
figure
subplot(2,1,1)
position_O = 1:1:7;
boxplot(data{1,1},'colors','k','positions',position_O,'width',0.18,'Whisker',10);
position_S = 0.7:1:6.7;
hold on
boxplot(data{1,2},'colors',[127 63 152]/255,'positions',position_S,'width',0.18,'Whisker',10);
ylim([-0.05 0.55])
set(gca,'XTickLabel',{'PHS','PH','PS','HS','P','H','S'});
set(gca, 'XDir','reverse')
ylabel('proportion of neuron')
title('Preferred context')
subplot(2,1,2)
position_O = 1:1:7;
boxplot(data{2,1},'colors','k','positions',position_O,'width',0.18,'Whisker',10);
hold on
position_S = 0.7:1:6.7;
boxplot(data{2,2},'colors',[127 63 152]/255,'positions',position_S,'width',0.18,'Whisker',10);
ylim([-0.05 0.55])
set(gca,'XTickLabel',{'PHS','PH','PS','HS','P','H','S'});
set(gca, 'XDir','reverse')
ylabel('proportion of neuron')
title('Non-preferred context')

%%%%%%%%% For whole data %%%%%%%%%%
load('E:\Miniscope_MACPP_LNP\datapath.mat','datapath_MA')
datapath = datapath_MA;
LNP = {};
for n = 1:length(datapath)
    cd(datapath{n})
    load('lnp_whole_glmnet.mat')
    histall = cellfun(@(x) histcounts(x,[0.5:1:7.5])./length(x),...
        lnp_model,'uni',0);
        for ii = 1:size(histall,2)
                LNP{ii}(n,:) = histall{ii};
        end
end
LNPf = LNP;
save('LNP_whole.mat','LNPf')
% LNPg = LNP;
% save('LNP_whole.mat','LNPg','-append')

data = LNPf;
figure
position_O = 1:1:7;
boxplot(data{1,1},'colors','k','positions',position_O,'width',0.18,'Whisker',10);
position_S = 0.7:1:6.7;
hold on
boxplot(data{1,2},'colors','r','positions',position_S,'width',0.18,'Whisker',10);
ylim([-0.05 0.55])
set(gca,'XTickLabel',{'PHS','PH','PS','HS','P','H','S'});
set(gca, 'XDir','reverse')
ylabel('proportion of neuron')
