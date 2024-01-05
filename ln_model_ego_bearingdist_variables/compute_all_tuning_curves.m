%% Description
% This will compute the firing rate tuning curves for position, head
% direction, running speed, and theta phase.
% compute tuning curves for position, head direction, speed, and theta phase
%% select data
neuron = neuronIndiv{session2plotx,session2ploty};
behav = behavIndiv{session2plotx,session2ploty};
%% re-run the specific cell
[A, modelType,spiketrain_all, dt, xnbins, ynbins, direction, speed, egob]...
    = prep_lnp_data(neuron,behav,thresh,p.pos_bin_size,p.n_dir_bins,p.n_speed_bins);

spktrain = spiketrain_all(:,cell2plot);
[testFit,~,param,numFolds] = fit_all_ln_models(A,modelType,spktrain,dt,[],[],p.method);
[lnp_model_ea,llh_ea] = select_best_model(testFit,'correlation');

%% calcualte rate maps
[ratemapAll,~,~] = compute_firing_ratemap(neuron,...
    behav,thresh);
ratemap = ratemapAll{cell2plot};
[pos_curve] = filter2DMatrices(ratemap,1);

filter = gaussmf(-4:4,[2 0]);
filter = filter/sum(filter); 
fr = spiketrain_all(:,cell2plot)/dt;
smooth_fr = conv(fr,filter,'same');
[hd_curve] = compute_1d_tuning_curve(direction,smooth_fr,p.n_dir_bins,0,2*pi);
[speed_curve] = compute_1d_tuning_curve(speed,smooth_fr,p.n_speed_bins,0,ceil(max(speed)));
[eb_curve] = compute_1d_tuning_curve(egob,smooth_fr,p.n_dir_bins,0,2*pi);
