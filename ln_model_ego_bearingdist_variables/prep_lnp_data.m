function [A,modelType,spiketrain_all,dt]...
    = prep_lnp_data(neuron,behav,thresh,pos_bin_size,n_dir_bins,n_speed_bins,frame_avg,numModels)
if ~exist('numModels','var')||isempty(numModels)
    numModels = 15;
end
if ~exist('frame_avg','var') || isempty(frame_avg)
    frame_avg = 8;
end
dt = frame_avg * 0.0667;

spiketrain_all = double(neuron.S > repmat(thresh,1,size(neuron.S,2)))';
position = interp1(behav.time, behav.position, neuron.time);
if any(isnan(position(:,1)))
    idx = find(isnan(position(:,1)));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        position(idx(1:ind),:) = repmat(position(idx(ind)+1,:),size(idx(1:ind),1),1);
        position(idx(ind+1:end),:) = repmat(position(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            position(idx,:) = repmat(position(idx(end)+1,:),size(idx,1),1);
        else
            position(idx,:) = repmat(position(idx(1)-1,:),size(idx,1),1);
        end
    end
end
positionblue = interp1(behav.time, behav.positionblue, neuron.time);
if any(isnan(positionblue(:,1)))
    idx = find(isnan(positionblue(:,1)));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        positionblue(idx(1:ind),:) = repmat(positionblue(idx(ind)+1,:),size(idx(1:ind),1),1);
        positionblue(idx(ind+1:end),:) = repmat(positionblue(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            positionblue(idx,:) = repmat(positionblue(idx(end)+1,:),size(idx,1),1);
        else
            positionblue(idx,:) = repmat(positionblue(idx(1)-1,:),size(idx,1),1);
        end
    end
end
speed = interp1(behav.time, behav.speed, neuron.time);
if any(isnan(speed))
    idx = find(isnan(speed));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        speed(idx(1:ind),:) = repmat(speed(idx(ind)+1,:),size(idx(1:ind),1),1);
        speed(idx(ind+1:end),:) = repmat(speed(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            speed(idx,:) = repmat(speed(idx(end)+1,:),size(idx,1),1);
        else
            speed(idx,:) = repmat(speed(idx(1)-1,:),size(idx,1),1);
        end
    end
end
eb = neuron.eb;
if any(isnan(eb))
    idx = find(isnan(eb));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        eb(idx(1:ind),:) = repmat(eb(idx(ind)+1,:),size(idx(1:ind),1),1);
        eb(idx(ind+1:end),:) = repmat(eb(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            eb(idx,:) = repmat(eb(idx(end)+1,:),size(idx,1),1);
        else
            eb(idx,:) = repmat(eb(idx(1)-1,:),size(idx,1),1);
        end
    end
end
cornerdist = neuron.edist;
if any(isnan(cornerdist))
    idx = find(isnan(cornerdist));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        cornerdist(idx(1:ind),:) = repmat(cornerdist(idx(ind)+1,:),size(idx(1:ind),1),1);
        cornerdist(idx(ind+1:end),:) = repmat(cornerdist(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            cornerdist(idx,:) = repmat(cornerdist(idx(end)+1,:),size(idx,1),1);
        else
            cornerdist(idx,:) = repmat(cornerdist(idx(1)-1,:),size(idx,1),1);
        end
    end
end
wallb = neuron.wallb;
if any(isnan(wallb))
    idx = find(isnan(wallb));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        wallb(idx(1:ind),:) = repmat(wallb(idx(ind)+1,:),size(idx(1:ind),1),1);
        wallb(idx(ind+1:end),:) = repmat(wallb(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            wallb(idx,:) = repmat(wallb(idx(end)+1,:),size(idx,1),1);
        else
            wallb(idx,:) = repmat(wallb(idx(1)-1,:),size(idx,1),1);
        end
    end
end
walldist = neuron.walldist;
if any(isnan(walldist))
    idx = find(isnan(walldist));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        walldist(idx(1:ind),:) = repmat(walldist(idx(ind)+1,:),size(idx(1:ind),1),1);
        walldist(idx(ind+1:end),:) = repmat(walldist(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            walldist(idx,:) = repmat(walldist(idx(end)+1,:),size(idx,1),1);
        else
            walldist(idx,:) = repmat(walldist(idx(1)-1,:),size(idx,1),1);
        end
    end
end
ctrb = neuron.ctrb;
if any(isnan(ctrb))
    idx = find(isnan(ctrb));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        ctrb(idx(1:ind),:) = repmat(ctrb(idx(ind)+1,:),size(idx(1:ind),1),1);
        ctrb(idx(ind+1:end),:) = repmat(ctrb(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            ctrb(idx,:) = repmat(ctrb(idx(end)+1,:),size(idx,1),1);
        else
            ctrb(idx,:) = repmat(ctrb(idx(1)-1,:),size(idx,1),1);
        end
    end
end
ctrdist = neuron.ctrdist;
if any(isnan(ctrdist))
    idx = find(isnan(ctrdist));
    didx = diff(idx);
    if ~isempty(didx) && any(didx ~= 1)
        ind = find(didx ~= 1);
        ctrdist(idx(1:ind),:) = repmat(ctrdist(idx(ind)+1,:),size(idx(1:ind),1),1);
        ctrdist(idx(ind+1:end),:) = repmat(ctrdist(idx(ind+1)-1,:),size(idx(ind+1:end),1),1);
    else
        if idx(1) == 1
            ctrdist(idx,:) = repmat(ctrdist(idx(end)+1,:),size(idx,1),1);
        else
            ctrdist(idx,:) = repmat(ctrdist(idx(1)-1,:),size(idx,1),1);
        end
    end
end
%position
ds = 1:frame_avg:length(position);
position = position(ds,:);
positionblue = positionblue(ds,:);
%speed
speed = speed(ds,:);
%ego bearing
eb = eb(ds,:);
%corner distance
cornerdist = cornerdist(ds,:);
%ego bearing to walls
wallb = wallb(ds,:);
%wall distance
walldist = walldist(ds,:);
%ego bearing to the center
ctrb = ctrb(ds,:);
%center distance
ctrdist = ctrdist(ds,:);


%neuron
fun = @(block_struct) sum(block_struct.data);
spiketrain_all = blockproc(spiketrain_all, [frame_avg, size(spiketrain_all,2)], fun);

[posgrid,xnbins,ynbins] = pos_map(position,pos_bin_size);
[hdgrid,hdVec,direction] = hd_map(positionblue(:,1),position(:,1),positionblue(:,2),position(:,2),n_dir_bins);
[speedgrid,speedVec] = speed_map(speed,n_speed_bins);
[ebgrid,~,egob] = eb_map(eb,n_dir_bins);
[cornerdistgrid,cdistVec] = corner_distance_map(cornerdist,pos_bin_size);
[wallbgrid,~,~] = eb_map(wallb,n_dir_bins);
[walldistgrid,~] = corner_distance_map(walldist,pos_bin_size);
[ctrbgrid,~,~] = eb_map(ctrb,n_dir_bins);
[ctrdistgrid,~] = corner_distance_map(ctrdist,pos_bin_size);

%% Build all 15 models
A = cell(numModels,1);
modelType = cell(numModels,1);

% % ALL VARIABLES
% A{1} = [ posgrid hdgrid speedgrid ebgrid]; modelType{1} = [1 1 1 1];
% % THREE VARIABLES
% A{2} = [ posgrid hdgrid speedgrid ]; modelType{2} = [1 1 1 0];
% A{3} = [ posgrid hdgrid  ebgrid]; modelType{3} = [1 1 0 1];
% A{4} = [ posgrid  speedgrid ebgrid]; modelType{4} = [1 0 1 1];
% A{5} = [  hdgrid speedgrid ebgrid]; modelType{5} = [0 1 1 1];
% % TWO VARIABLES
% A{6} = [ posgrid hdgrid]; modelType{6} = [1 1 0 0];
% A{7} = [ posgrid  speedgrid ]; modelType{7} = [1 0 1 0];
% A{8} = [ posgrid   ebgrid]; modelType{8} = [1 0 0 1];
% A{9} = [  hdgrid speedgrid ]; modelType{9} = [0 1 1 0];
% A{10} = [  hdgrid  ebgrid]; modelType{10} = [0 1 0 1];
% A{11} = [  speedgrid ebgrid]; modelType{11} = [0 0 1 1];
% % ONE VARIABLE
% A{12} = posgrid; modelType{12} = [1 0 0 0];
% A{13} = hdgrid; modelType{13} = [0 1 0 0];
% A{14} = speedgrid; modelType{14} = [0 0 1 0];
% A{15} = ebgrid; modelType{15} = [0 0 0 1];

% % using corner distance instead of position
% A{1} = [ cornerdistgrid hdgrid speedgrid ebgrid]; modelType{1} = [1 1 1 1];
% % THREE VARIABLES
% A{2} = [ cornerdistgrid hdgrid speedgrid ]; modelType{2} = [1 1 1 0];
% A{3} = [ cornerdistgrid hdgrid  ebgrid]; modelType{3} = [1 1 0 1];
% A{4} = [ cornerdistgrid  speedgrid ebgrid]; modelType{4} = [1 0 1 1];
% A{5} = [  hdgrid speedgrid ebgrid]; modelType{5} = [0 1 1 1];
% % TWO VARIABLES
% A{6} = [ cornerdistgrid hdgrid]; modelType{6} = [1 1 0 0];
% A{7} = [ cornerdistgrid  speedgrid ]; modelType{7} = [1 0 1 0];
% A{8} = [ cornerdistgrid   ebgrid]; modelType{8} = [1 0 0 1];
% A{9} = [  hdgrid speedgrid ]; modelType{9} = [0 1 1 0];
% A{10} = [  hdgrid  ebgrid]; modelType{10} = [0 1 0 1];
% A{11} = [  speedgrid ebgrid]; modelType{11} = [0 0 1 1];
% % ONE VARIABLE
% A{12} = cornerdistgrid; modelType{12} = [1 0 0 0];
% A{13} = hdgrid; modelType{13} = [0 1 0 0];
% A{14} = speedgrid; modelType{14} = [0 0 1 0];
% A{15} = ebgrid; modelType{15} = [0 0 0 1];

% % using wall bearing and distance
% A{1} = [ walldistgrid hdgrid speedgrid wallbgrid]; modelType{1} = [1 1 1 1];
% % THREE VARIABLES
% A{2} = [ walldistgrid hdgrid speedgrid ]; modelType{2} = [1 1 1 0];
% A{3} = [ walldistgrid hdgrid  wallbgrid]; modelType{3} = [1 1 0 1];
% A{4} = [ walldistgrid  speedgrid wallbgrid]; modelType{4} = [1 0 1 1];
% A{5} = [  hdgrid speedgrid wallbgrid]; modelType{5} = [0 1 1 1];
% % TWO VARIABLES
% A{6} = [ walldistgrid hdgrid]; modelType{6} = [1 1 0 0];
% A{7} = [ walldistgrid  speedgrid ]; modelType{7} = [1 0 1 0];
% A{8} = [ walldistgrid   wallbgrid]; modelType{8} = [1 0 0 1];
% A{9} = [  hdgrid speedgrid ]; modelType{9} = [0 1 1 0];
% A{10} = [  hdgrid  wallbgrid]; modelType{10} = [0 1 0 1];
% A{11} = [  speedgrid wallbgrid]; modelType{11} = [0 0 1 1];
% % ONE VARIABLE
% A{12} = walldistgrid; modelType{12} = [1 0 0 0];
% A{13} = hdgrid; modelType{13} = [0 1 0 0];
% A{14} = speedgrid; modelType{14} = [0 0 1 0];
% A{15} = wallbgrid; modelType{15} = [0 0 0 1];

% using center bearing and distance
A{1} = [ ctrdistgrid hdgrid speedgrid ctrbgrid]; modelType{1} = [1 1 1 1];
% THREE VARIABLES
A{2} = [ ctrdistgrid hdgrid speedgrid ]; modelType{2} = [1 1 1 0];
A{3} = [ ctrdistgrid hdgrid  ctrbgrid]; modelType{3} = [1 1 0 1];
A{4} = [ ctrdistgrid  speedgrid ctrbgrid]; modelType{4} = [1 0 1 1];
A{5} = [  hdgrid speedgrid ctrbgrid]; modelType{5} = [0 1 1 1];
% TWO VARIABLES
A{6} = [ ctrdistgrid hdgrid]; modelType{6} = [1 1 0 0];
A{7} = [ ctrdistgrid  speedgrid ]; modelType{7} = [1 0 1 0];
A{8} = [ ctrdistgrid   ctrbgrid]; modelType{8} = [1 0 0 1];
A{9} = [  hdgrid speedgrid ]; modelType{9} = [0 1 1 0];
A{10} = [  hdgrid  ctrbgrid]; modelType{10} = [0 1 0 1];
A{11} = [  speedgrid ctrbgrid]; modelType{11} = [0 0 1 1];
% ONE VARIABLE
A{12} = ctrdistgrid; modelType{12} = [1 0 0 0];
A{13} = hdgrid; modelType{13} = [0 1 0 0];
A{14} = speedgrid; modelType{14} = [0 0 1 0];
A{15} = ctrbgrid; modelType{15} = [0 0 0 1];
end

