function get_egobearing_wall_center(datapath)
%this function calculates animals egeocentric bearings to the nearest wall
% and distance, egocentric bearings to the environmental center and
% distance.
%for detailed calculation method, see LaChance et al, 2019, A sense of
%space in postrhinal cortex, Science.
%yanjuns@stanford.edu, 1/27/2023

if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)

load ('neuronIndividualsf.mat','neuronIndividualsf');
load ('behavIndividualsf.mat','behavIndividualsf');
% load ('env_geometry.mat', 'env_coor','S','ccav_coor');
load ('env_geometry.mat', 'env_coor','S');
binsize = 2;

for ss = 1:length(neuronIndividualsf)
    neuron = neuronIndividualsf{ss};
    behav = behavIndividualsf{ss};
    %% egocentric bearing to environment center
    env_coor_ea = env_coor{ss};
    env_coor_ctr = env_coor_ea(end,:);
    xyposition = behav.position/binsize;
    behav.ctrdist = pdist2(xyposition,env_coor_ctr)*2;
    allob = atan2(env_coor_ctr(:,2) - xyposition(:,2),env_coor_ctr(:,1) - xyposition(:,1));
    behav.ctrb = angdiff(behav.hd,allob);
    neuron.ctrb = interp1(behav.time,behav.ctrb,neuron.time);
    neuron.ctrdist = interp1(behav.time,behav.ctrdist,neuron.time);
    %% egocentric bearing to nearest wall
    % egocentric bearing to the nearest wall
    if strcmp(S{ss},'largeX') || strcmp(S{ss},'largeXm1')
        env_coor_wl = [ccav_coor{ss}(1,:);env_coor{ss}(1,:);...
            ccav_coor{ss}(2:3,:);env_coor{ss}(2,:);...
            ccav_coor{ss}(4:5,:);env_coor{ss}(3,:);...
            ccav_coor{ss}(6:7,:);env_coor{ss}(4,:);...
            ccav_coor{ss}(8,:)];
    else
        env_coor_wl = env_coor_ea(1:end-1,:);
    end
    % obtaining wall coordinates for the environment
    walls = cell(1,size(env_coor_wl,1));
    for ii = 1:size(env_coor_wl,1)
        if ii < size(env_coor_wl,1)
            resolution = max([numel(env_coor_wl(ii,1):0.1:env_coor_wl(ii+1,1));...
                numel(env_coor_wl(ii,2):0.1:env_coor_wl(ii+1,2));...
                numel(env_coor_wl(ii,1):-0.1:env_coor_wl(ii+1,1));...
                numel(env_coor_wl(ii,2):-0.1:env_coor_wl(ii+1,2))]);
            walls{ii} = [linspace(env_coor_wl(ii,1),env_coor_wl(ii+1,1),resolution);...
                linspace(env_coor_wl(ii,2),env_coor_wl(ii+1,2),resolution)]';
        else
            resolution = max([numel(env_coor_wl(ii,1):0.1:env_coor_wl(1,1));...
                numel(env_coor_wl(ii,2):0.1:env_coor_wl(1,2));...
                numel(env_coor_wl(ii,1):-0.1:env_coor_wl(1,1));...
                numel(env_coor_wl(ii,2):-0.1:env_coor_wl(1,2))]);
            walls{ii} = [linspace(env_coor_wl(ii,1),env_coor_wl(1,1),resolution);...
                linspace(env_coor_wl(ii,2),env_coor_wl(1,2),resolution)]';
        end
    end
%     figure
%     hold on
%     for ii = 1:size(walls,2)
%         plot(walls{ii}(:,1), walls{ii}(:,2))
%     end
    % calculate distance and bearing to the nearest wall
    behav.walldist = NaN(length(behav.hd),1);
    behav.wallb = NaN(length(behav.hd),1);
    for ii = 1:length(behav.position)
        alldist = cellfun(@(x) pdist2(x,behav.position(ii,:)/binsize), walls, 'UniformOutput',0);
        [mindist, minid] = cellfun(@min, alldist);
        [walldist,idx] = min(mindist);
        behav.walldist(ii,1) = walldist*binsize;
        wall_vector = walls{idx}(minid(idx),:);
        ecb_vector = wall_vector - behav.position(ii,:)/binsize;
        allowall = atan2(ecb_vector(:,2),ecb_vector(:,1));
        ewb = angdiff(behav.hd(ii),allowall);
        behav.wallb(ii) = ewb;
    end
    neuron.walldist = interp1(behav.time,behav.walldist,neuron.time);
    neuron.wallb = interp1(behav.time,behav.wallb,neuron.time);

    %% output data
    behavIndividualsf{ss} = behav;
    neuronIndividualsf{ss} = neuron;
end

save('neuronIndividualsf.mat','neuronIndividualsf','-append')
save('behavIndividualsf.mat','behavIndividualsf','-append')

end

