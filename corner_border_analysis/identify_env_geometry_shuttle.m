function [N,env_coor,env_coori] = identify_env_geometry_shuttle(datapath)
%identify_env_geometry
%   This function aims to find the coordinates of corners and the center of
%   the environments in the experiment. It uses regionprops to
%   automatically find the centroid and corners and then manually select
%   corners, as automatic corner detections are not perfect.
%   S: session names
%   N: number of corners for each session
%   env_coor: x y coordinates of corners and the last row is the
%   coordinates of the centroid.

if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd;
end
cd(datapath)
load('NBindivLR.mat','behavIndivLR')
behavIndividualsf = behavIndivLR;
ssnum = length(behavIndivLR);

%% get session names and their number of corners
N = NaN(1,ssnum); %number of corners for all the sessions.
for ii = 1:ssnum
    N(ii) = 4; %for the condition of shuttle box
end

%% identify corners with mannual interventions
load('NBindivLR.mat','countTime')
env_coor = cell(1,ssnum);
for ii = 1:ssnum
    if N(ii) > 0
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
        [x,y] = ginput(N(ii));
        while any(x(:) < 0) || any(y(:) < 0) ||...
                any(x(:) > size(env,2)+0.5) || any(y(:) > size(env,1)+0.5)
            fprintf('One or more selected points exceed environment bounds,please redo selection\n')
            [x,y] = ginput(N(ii));
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
            [x,y] = ginput(N(ii));
            while any(x(:) < 0) || any(y(:) < 0) ||...
                    any(x(:) > size(env,2)+0.5) || any(y(:) > size(env,1)+0.5)
                fprintf('One or more selected points exceed environment bounds,please redo selection\n')
                [x,y] = ginput(N(ii));
            end
            prompt = 'Keep current coordinates? Y/N: ';
            str = input(prompt,'s');
        end
        close;
        env_coor{ii} = [[x,y];stats.Centroid];
    end
end

%% find array index for each corresponding hand-picked corner
env_coori = cell(1,length(countTime));
for n = 1:length(countTime)
    A = countTime{1,n};
    envv_ori = env_coor{1,n};
    %get the xy coordinate of non-NaN elements in the ratemap
    row = [];col = [];
    for ii = 1:size(A,1)
        for jj = 1:size(A,2)
            if A(ii,jj) > 0
                row = [row;ii];
                col = [col;jj];
            end
        end
    end
    coor = [col,row];
    %search the nearest ratemap bin to the hand-picked corner
    k = dsearchn(coor, envv_ori);
    env_coori{n} = coor(k,:); %new corner coordinate
end
%% plot to check the selected corners
figure
k = length(env_coor);
for jj = 1:k
    subplot(ceil(k/2),2,jj)
    if ~isempty(env_coor{jj})
        imagesc(countTime{jj})
        axis image
        hold on
        plot(env_coor{jj}(:,1),env_coor{jj}(:,2),'.','color','r','markersize',20)
    end
end

figure
k = length(env_coor);
for jj = 1:k
    subplot(ceil(k/2),2,jj)
    if ~isempty(env_coor{jj})
        imagesc(countTime{jj})
        axis image
        hold on
        plot(env_coori{jj}(:,1),env_coori{jj}(:,2),'.','color','r','markersize',20)
    end
end

end

