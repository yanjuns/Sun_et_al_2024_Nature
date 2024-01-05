function [neuronIndivLRm, behavIndivLRm] = match_neuron_behav(neuronIndivLR, behavIndivLR, condition)
% This function performs multiple steps for ratemap analysis
% 1. downsample data to match bahavior coverages across left and right cpp
% experment chambers or across different sessions;
% 2. calculate raw firing ratemap by using matched neuron and behav data;
% 3. smooth the raw firing ratemap;
% 4. trim the smoothed raw firing ratemap in order to prepare for
% correlation computation.
% Input args:
%  -neuronIndivLR, splitted neuron data from left and right cpp chambers
%  -behavIndivLR, splitted behav data from left and right cpp chamebrs
% Output:
%  -neuronIndivLRm, splitted neuron data from left and right cpp chambers
%  -behavIndivLRm, splitted behav data from left and right cpp chamebrs
% note: this function contains customized functions naming as follows
%  -downsample_match_position
% Yanjun Sun, Stanford University, 9/16/2019

if ~exist('condition', 'var') || isempty(condition)
    condition = 'xsession';
end

switch condition
    case {'LR', 'left right', 'Left Right', 'left and right', 'Left and Right'}
        % preallocate
        neuronIndivLRm = cell(size(neuronIndivLR,1), size(neuronIndivLR,2));
        behavIndivLRm = cell(size(behavIndivLR,1), size(behavIndivLR,2));
        % downsample the data points to match behaviors of left and right chambers
%         for ii = 1:length(neuronIndivLR)
%             [neuronLm, behavLm, neuronRm, behavRm] = downsample_match_position(neuronIndivLR{1,ii},behavIndivLR{1,ii},neuronIndivLR{2,ii},behavIndivLR{2,ii});
%             neuronIndivLRm{1,ii} = neuronLm; behavIndivLRm{1,ii} = behavLm;
%             neuronIndivLRm{2,ii} = neuronRm; behavIndivLRm{2,ii} = behavRm;
%         end
        % downsample when data were organized in one row as LRLR
            for ii = 1:length(neuronIndivLR)/2
                jj = 2*ii-1;
                [neuronLm, behavLm, neuronRm, behavRm] = downsample_match_position(neuronIndivLR{1,jj},behavIndivLR{1,jj},neuronIndivLR{1,jj+1},behavIndivLR{1,jj+1});
                neuronIndivLRm{1,jj} = neuronLm; behavIndivLRm{1,jj} = behavLm;
                neuronIndivLRm{1,jj+1} = neuronRm; behavIndivLRm{1,jj+1} = behavRm;
            end
    case{'xsession', 'cross-session', 'x-session', 'longitudinal'}
        % downsample the data points to match behaviors of longitudinal recordings
        % data will be stored in a cell array, with (1,2) stores matched data of
        % session 1 from 1 and 2 matching, (2,1) stores matched data of session 2
        % from 1 and 2 matching, so on and so forth. dimension will be how
        % many subsegments of the data are required
        dimension = size(neuronIndivLR,1);
        if dimension > 1
            neuronIndivLRm = cell(length(neuronIndivLR),length(neuronIndivLR), dimension);
            behavIndivLRm = cell(length(behavIndivLR),length(behavIndivLR), dimension);
            for d = 1:dimension
                for ii = 1:length(neuronIndivLR)
                    for jj = ii:length(neuronIndivLR)
                    [neuronL1, behavL1, neuronL2, behavL2] = downsample_match_position(neuronIndivLR{d,ii},behavIndivLR{d,ii},neuronIndivLR{d,jj},behavIndivLR{d,jj});
                    neuronIndivLRm{ii,jj,d} = neuronL1; behavIndivLRm{ii,jj,d} = behavL1;
                    neuronIndivLRm{jj,ii,d} = neuronL2; behavIndivLRm{jj,ii,d} = behavL2;
                    end
                end
            end
        else
            neuronIndivLRm = cell(length(neuronIndivLR),length(neuronIndivLR));
            behavIndivLRm = cell(length(behavIndivLR),length(behavIndivLR));
            for ii = 1:length(neuronIndivLR)
                for jj = ii:length(neuronIndivLR)
                [neuronL1, behavL1, neuronL2, behavL2] = downsample_match_position(neuronIndivLR{1,ii},behavIndivLR{1,ii},neuronIndivLR{1,jj},behavIndivLR{1,jj});
                neuronIndivLRm{ii,jj} = neuronL1; behavIndivLRm{ii,jj} = behavL1;
                neuronIndivLRm{jj,ii} = neuronL2; behavIndivLRm{jj,ii} = behavL2;
                end
            end
        end
        
    case{'longitudinal_full','Longitudinal_full'}
        % match all the longitudinal sessions together from each side
        % of the CPP.
        dimension = size(neuronIndivLR,1);
        if dimension > 1
            [neuronIndivLm, behavIndivLm] = downsample_match_position_long(neuronIndivLR(1,:), behavIndivLR(1,:));
            [neuronIndivRm, behavIndivRm] = downsample_match_position_long(neuronIndivLR(2,:), behavIndivLR(2,:));
            neuronIndivLRm = [neuronIndivLm;neuronIndivRm];
            behavIndivLRm = [behavIndivLm;behavIndivRm];
        else
            [neuronIndivLRm, behavIndivLRm] = downsample_match_position_long(neuronIndivLR, behavIndivLR);
        end

end
end