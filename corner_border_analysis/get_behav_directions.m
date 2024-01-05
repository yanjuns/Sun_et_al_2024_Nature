function get_behav_directions(datapath)
%this function calculates animal's movement direction, allocentric bearing
%to the environmental corners, egocentric bearing to the corners, and the
%egocentric distance to the corners. 
%for detailed calculation method, see LaChance et al, 2019, A sense of
%space in postrhinal cortex, Science. 
%yanjuns@stanford.edu, 1/27/2023

if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)

load ('neuronIndividualsf.mat','neuronIndividualsf');
load ('behavIndividualsf.mat','behavIndividualsf');
load ('env_geometry.mat', 'env_coor')
figure
hold on
for ss = 1:length(neuronIndividualsf)
    %calculate head direction
    neuronIndividualsf{ss}.hd = interp1(behavIndividualsf{ss}.time,behavIndividualsf{ss}.hd,neuronIndividualsf{ss}.time);
    %calculate movement direction
    behavIndividualsf{ss}.md = atan2([0;diff(behavIndividualsf{ss}.position(:,2))],[0;diff(behavIndividualsf{ss}.position(:,1))]);
    neuronIndividualsf{ss}.md = interp1(behavIndividualsf{ss}.time,behavIndividualsf{ss}.md,neuronIndividualsf{ss}.time);
    
    %calculate allocentric bearing and egocentric bearing to corners
    env_coor_ea = env_coor{ss};
    %last row is centroid
    env_coor_ea = env_coor_ea(1:end-1,:); 
    %flip the x and y for correct calculation
%     env_coor_ea = [env_coor_ea(:,2), env_coor_ea(:,1)]; 
    %downscale xy position to match the location of the corner
    xyposition = behavIndividualsf{ss}.position/2;
    dist2corner = pdist2(xyposition,env_coor_ea);
    [D,I] = min(dist2corner,[],2);
    D = D*2;
    ab = NaN(length(dist2corner),1);
    eb = NaN(length(dist2corner),1);
    for ii = 1:length(dist2corner)
        ab(ii) = atan2(env_coor_ea(I(ii),2) - xyposition(ii,2),...
            env_coor_ea(I(ii),1) - xyposition(ii,1));
        eb(ii) = angdiff(behavIndividualsf{ss}.hd(ii),ab(ii));
    end
    behavIndividualsf{ss}.ab = ab;
    behavIndividualsf{ss}.eb = eb;
    behavIndividualsf{ss}.edist = D;
    neuronIndividualsf{ss}.ab = interp1(behavIndividualsf{ss}.time,ab,neuronIndividualsf{ss}.time);
    neuronIndividualsf{ss}.eb = interp1(behavIndividualsf{ss}.time,eb,neuronIndividualsf{ss}.time);
    neuronIndividualsf{ss}.edist = interp1(behavIndividualsf{ss}.time,D,neuronIndividualsf{ss}.time);

    %plot the allocentric bearing
    subplot(round(length(neuronIndividualsf)/2),2,ss);
    scatter(behavIndividualsf{ss}.position(:,1),behavIndividualsf{ss}.position(:,2)...
        ,20,behavIndividualsf{ss}.ab,'filled')
    colormap hsv
    axis image
end
savefig(gcf,'allocentric_bearing')
save('neuronIndividualsf.mat','neuronIndividualsf','-append')
save('behavIndividualsf.mat','behavIndividualsf','-append')



end

