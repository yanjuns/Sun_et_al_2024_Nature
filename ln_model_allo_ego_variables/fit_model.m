function [testFit,trainFit,param_mean] = fit_model(...
    A,dt,spiketrain,filter,modelType,numFolds,numPos,numHD,numSpd,numEgob,method,lambda)

%% Description
% This code will section the data into 10 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 9 sections, and test the model performance on the
% remaining section. This procedure will be repeated 10 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.
% Yanjun edits on Mar 2022
% Introduced a faster regression using glmnet,but this result may give
% slightly different results compared to fminunc. 
% if faster is true, it will not using glmnet
if ~exist('lambda','var') || isempty(lambda)
    lambda = 1e-3;
end
if ~exist('method','var') || isempty(method)
    method = 'glmnet';
end
%% Initialize matrices and section the data for k-fold cross-validation

[~,numCol] = size(A);
sections = numFolds*5;

% divide the data up into 5*num_folds pieces
edges = round(linspace(1,numel(spiketrain)+1,sections+1));

% initialize matrices
testFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of test data
trainFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of train data
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
for k = 1:numFolds
    fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
    
    % get test data from edges - each test data chunk comes from entire session
    test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1]   ;
    
    test_spikes = spiketrain(test_ind); %test spiking
    smooth_spikes_test = conv(test_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_test = smooth_spikes_test./dt;
    test_A = A(test_ind,:);
    
    % training data
    train_ind = setdiff(1:numel(spiketrain),test_ind);
    train_spikes = spiketrain(train_ind);
    smooth_spikes_train = conv(train_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_train = smooth_spikes_train./dt;
    train_A = A(train_ind,:);
    data{1} = train_A; data{2} = train_spikes;
    switch method
        case('fminunc')
            %method 1, using fminunc
            opts = optimoptions('fminunc','Gradobj','on','Hessian','on','MaxIter',400,'Display','off','UseParallel',true);
            if k == 1
                init_param = 1e-3*randn(numCol, 1);
            else
                init_param = param;
            end
            [param] = fminunc(@(param) ln_poisson_model(param,data,modelType,numPos,numHD,numSpd,numEgob),init_param,opts);
        case('glmnet')
            %method 2, using glmnet
            options = glmnetSet;
            options.alpha = 0.1;
            options.lambda = linspace(lambda,lambda/0.1,12);
%             options.maxit = 10000;
%             CVerr = cvglmnet(data{1}, data{2}, 'poisson', options);
            B = glmnet(data{1}, data{2}, 'poisson', options);
            [~,I] = min(abs(B.lambda -lambda)); %find the closest lambda
            param = B.beta(:,I); %estimated param from cvglmnet
    end
    
    %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_test = exp(test_A * param)/dt;
    smooth_fr_hat_test = conv(fr_hat_test,filter,'same'); %returns vector same size as original
%     switch method
%         case('fminunc')
%             fr_hat_test = exp(test_A * param)/dt;
%             smooth_fr_hat_test = conv(fr_hat_test,filter,'same'); %returns vector same size as original
%         case('glmnet')
%             result = glmnetPredict(B, test_A, lambda,'response');
%             fr_hat_test = result/dt;
%             smooth_fr_hat_test = conv(fr_hat_test,filter,'same');
%     end

    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
    sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % compute correlation
    dataL = length(smooth_fr_test);
    nboot = 500;
    deltaTall = randi([round(0.1*dataL),dataL-round(0.1*dataL)],nboot,1); %msec
    nullcorr = [];
    for nE = 1:nboot
        deltaT = deltaTall(nE);
        timeboot = (1:dataL) + deltaT;
        idx = timeboot > dataL;
        timeboot(idx) = timeboot(idx) - dataL;
        [~,index] = sort(timeboot);
        shuffled_fr = smooth_fr_test(index);
        nullcorr = [nullcorr;corr(shuffled_fr, smooth_fr_test,'type','Pearson')];
    end
    nullcorr95 = quantile(nullcorr,0.95);
    correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,...
        'type','Pearson') - nullcorr95;
    
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
    r = exp(test_A * param);n = test_spikes; meanFR_test = nanmean(test_spikes);
    log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
    log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
    log_llh_test = (-log_llh_test_model + log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test;
    
    % compute MSE
    mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
    
    % fill in all the relevant values for the test fit cases
    testFit(k,:) = [varExplain_test correlation_test log_llh_test mse_test sum(n) numel(test_ind) ];
    
    %% training data
    %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_train = exp(train_A * param)/dt;
    smooth_fr_hat_train = conv(fr_hat_train,filter,'same'); %returns vector same size as original
%     switch method
%         case('fminunc')
%             fr_hat_train = exp(train_A * param)/dt;
%             smooth_fr_hat_train = conv(fr_hat_train,filter,'same'); %returns vector same size as original
%         case('glmnet')
%             result_train = glmnetPredict(B, train_A, lambda,'response');
%             fr_hat_train = result_train/dt;
%             smooth_fr_hat_train = conv(fr_hat_train,filter,'same');
%     end

    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_train-smooth_fr_train).^2);
    sst = sum((smooth_fr_train-mean(smooth_fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % compute correlation
    correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson');
    
    % compute log-likelihood
    r_train = exp(train_A * param); n_train = train_spikes; meanFR_train = nanmean(train_spikes);
    log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train;
    
    % compute MSE
    mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
    
    trainFit(k,:) = [varExplain_train correlation_train log_llh_train mse_train sum(n_train) numel(train_ind)];

    % save the parameters
    paramMat(k,:) = param;

end

param_mean = nanmean(paramMat);

return
