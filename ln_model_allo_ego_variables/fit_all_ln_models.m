function [testFit,trainFit,param,numFolds] = fit_all_ln_models(...
    A,modelType,spiketrain,dt,numFolds,numModels,method)
%% Description
% The model: r = exp(W*theta), where r is the predicted # of spikes, W is a
% matrix of one-hot vectors describing variable (P, H, S, or T) values, and
% theta is the learned vector of parameters.
if ~exist('method','var') || isempty(method)
    method = 'glmnet';
end
if ~exist('numModels','var')||isempty(numModels)
    numModels = 15;
end
if ~exist('numFolds','var')||isempty(numFolds)
    numFolds = 10;
end

%% Fit all 15 LN models
testFit = cell(numModels,1);
trainFit = cell(numModels,1);
param = cell(numModels,1);
numPos = size(A{12},2);
numHD = size(A{13},2);
numSpd = size(A{14},2);
numEgob = size(A{15},2);

% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]);
filter = filter/sum(filter); 

% compute the number of folds we would like to do
for n = 1:numModels
    % to find a resonable lambda value to initiate glmnet
    if strcmp(method, 'glmnet')
        A_sub = A{n};
        spiketrain_sub = spiketrain;
        CVerr = fit_model_CV(A_sub,spiketrain_sub);
        lambda = CVerr.lambda_min;
    else
        lambda = 1e-3;
    end
    fprintf('\t- Fitting model %d of %d\n', n, numModels);
    [testFit{n},trainFit{n},param{n}] = fit_model(A{n},dt,spiketrain,filter,modelType{n},numFolds,...
        numPos,numHD,numSpd,numEgob,method,lambda);
end

end
