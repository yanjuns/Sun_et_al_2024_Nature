function CVerr = fit_model_CV(A,spiketrain)
%% Description
% This code will run a test trial to estimate the best lambda for glmnet to
% fit. 
% run a test trial to get initial lambda for glmnet
options = glmnetSet;
options.alpha = 0.1;
options.nlambda = 50;
% options.maxit = 10000;
% options.standardize = false;
CVerr = cvglmnet(A, spiketrain, 'poisson', options, 'deviance', 3);