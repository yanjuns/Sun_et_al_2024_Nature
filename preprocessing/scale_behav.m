function scale_behav(datapath, fscale)
%This function aims to scale the trakced behav position to reflect the true
%size of the arena. 
if ~exist('fscale','var') || isempty(fscale)
    fscale = 1.25;
end
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd;
end
cd(datapath)
%% assemble behavIndividuals and adjust the scale of position
load('sessions.mat','sessions');
behavIndividuals = cell(1,length(sessions));
for ii = 1:length(sessions)
    %get corrected tracklength
    ss = sessions{ii};
    if contains(ss,'circle')
        trackLength = 35; %cm
    elseif contains(ss, ['square','lw'])
        trackLength = 33;
    elseif contains(ss,'square')
        trackLength = 30;
    elseif contains(ss,'squareCA1')
        trackLength = 25;
    elseif contains(ss,'triangle')
        trackLength = 32;
    elseif contains(ss,'hex')
        trackLength = 37;
    elseif contains(ss,'cross')
        trackLength = 30;
    elseif contains(ss,'rect')
        trackLength = 32;
    elseif contains(ss,'irreg')
        trackLength = 30;
    elseif contains(ss,'shuttle')
        trackLength = 52;
    elseif contains(ss,['large','RT'])
        trackLength = 46;
    elseif contains(ss,'rightTri')
        trackLength = 46;
    elseif contains(ss,'trapezoid')
        trackLength = 46;
    elseif contains(ss,'large')
        trackLength = 40;
    elseif contains(ss,'icm')
        trackLength = 40;
    elseif contains(ss,'oval')
        trackLength = 36;
    else
        trackLength = 40;%cm  % cpp experiments
    end
    
    %load behav
    load([sessions{ii},'_Behav.mat'], 'behav');
    behav0 = behav;
    t = find(diff(behav0.time)<=0);
    while ~isempty(t)
        behav0.time(t+1) = behav0.time(t)+1;
        t = find(diff(behav0.time)<=0);
    end
    
    %correct position
    edgecorrect = min(behav0.position);
    behav0.position = behav0.position - edgecorrect;
    scalef = trackLength/max(max(behav0.position));
    behav0.position = behav0.position * scalef;
    %correct positionblue
    behav0.positionblue = behav0.positionblue - edgecorrect;
    behav0.positionblue = behav0.positionblue * scalef;
    %correct speed and distance
    dx = [0; diff(behav0.position(:,1))];
    dy = [0; diff(behav0.position(:,2))];
    behav0.speed = sqrt((dx).^2+(dy).^2)/behav0.dt;
    behav0.speed = smoothts(behav0.speed','b',ceil(1/behav0.dt/5));
    behav0.distancered = sum(sqrt((dx).^2+(dy).^2));
    %correct blue speed and distance
    dx = [0; diff(behav0.positionblue(:,1))];
    dy = [0; diff(behav0.positionblue(:,2))];
    behav0.speedblue = sqrt((dx).^2+(dy).^2)/behav0.dt;
    behav0.speedblue = smoothts(behav0.speedblue','b',ceil(1/behav0.dt/5));
    behav0.distanceblue = sum(sqrt((dx).^2+(dy).^2));
    
    %% multiply by a final scale for a better binning
    behav0.position = behav0.position .* fscale;
    behav0.positionblue = behav0.positionblue .* fscale;
    behav0.trackLength = trackLength * fscale;
    %assign corrected vars
    behavIndividuals{ii} = behav0;
end

save ('behavIndividuals.mat', 'behavIndividuals', '-v7.3')

end

