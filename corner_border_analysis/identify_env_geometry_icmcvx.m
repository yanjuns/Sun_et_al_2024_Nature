function [cvex_coor,cvex_coori] = identify_env_geometry_icmcvx(datapath)
% this function is specifically for picking convex corners from the
% inserted corner manipulation (icm) experiment. 

if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd;
end
cd(datapath)

load('firingrateAll.mat','countTime')
load('behavIndividualsf.mat','behavIndividualsf')
load('env_geometry.mat','S')
numcor = 3; %number of corners
idx = strcmp(S,'icm0');
ssnum = 1:length(S);
ii = ssnum(idx);
ii = ii(1);

env = countTime{ii};
BW = env;
BW = BW > 0;
stats = regionprops(BW, {'Centroid','Extrema'});
figure;
set(gcf, 'Position', [50, 300, 1600, 800]);
subplot(1,2,1)
plot(behavIndividualsf{ii}.position(:,1),behavIndividualsf{ii}.position(:,2));
axis image
set(gca, 'YDir','reverse')
subplot(1,2,2)
imagesc(env);
hold on
plot(stats.Centroid(:,1),stats.Centroid(:,2),'r.', 'MarkerSize', 30)
plot(stats.Extrema(:,1),stats.Extrema(:,2),'r.', 'MarkerSize', 30)
axis image
set (gcf, 'WindowButtonMotionFcn', @mouseMove);
%manually select corners
[x,y] = ginput(numcor);
while any(x(:) < 0) || any(y(:) < 0) ||...
        any(x(:) > size(env,2)+0.5) || any(y(:) > size(env,1)+0.5)
    fprintf('One or more selected points exceed environment bounds,please redo selection\n')
    [x,y] = ginput(numcor);
end
prompt = 'Keep current coordinates? Y/N: ';
str = input(prompt,'s');
while strcmp(str,'N')==1 || strcmp(str,'n')==1
    close;
    figure;
    set(gcf, 'Position', [50, 300, 1600, 800]);
    subplot(1,2,1)
    plot(behavIndividualsf{ii}.position(:,1),behavIndividualsf{ii}.position(:,2));
    axis image
    set(gca, 'YDir','reverse')
    subplot(1,2,2)
    imagesc(env);
    hold on
    plot(stats.Centroid(:,1),stats.Centroid(:,2),'r.', 'MarkerSize', 30)
    plot(stats.Extrema(:,1),stats.Extrema(:,2),'r.', 'MarkerSize', 30)
    axis image
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    [x,y] = ginput(numcor);
    while any(x(:) < 0) || any(y(:) < 0) ||...
            any(x(:) > size(env,2)+0.5) || any(y(:) > size(env,1)+0.5)
        fprintf('One or more selected points exceed environment bounds,please redo selection\n')
        [x,y] = ginput(numcor);
    end
    prompt = 'Keep current coordinates? Y/N: ';
    str = input(prompt,'s');
end
close;
cvx_coor = [[x,y];[mean(x),mean(y)]];

%% find array index for each corresponding hand-picked corner
A = countTime{1,ii};
envv_ori = cvx_coor;
%get the xy coordinate of non-NaN elements in the ratemap
row = [];col = [];
for m = 1:size(A,1)
    for jj = 1:size(A,2)
        if A(m,jj) > 0
            row = [row;m];
            col = [col;jj];
        end
    end
end
coor = [col,row];
%search the nearest ratemap bin to the hand-picked corner
k = dsearchn(coor, envv_ori);
cvx_coori = coor(k,:); %new corner coordinate

%% change these to cell array for further calculation
cvex_coor = {};
cvex_coor{1} = cvx_coor;
cvex_coor = repmat(cvex_coor, 1,length(S));
cvex_coori = {};
cvex_coori{1} = cvx_coori;
cvex_coori = repmat(cvex_coori, 1,length(S));

figure
imagesc(countTime{ii})
axis image
hold on
plot(cvx_coor(:,1),cvx_coor(:,2),'.','color','r','markersize',20)

end
