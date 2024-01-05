function behav = behavPlot(behav,speed,frames_pass)
scale=behav.trackLength/100;
% plot behavior data
figure
ha = tight_subplot(5,1,[.001 .03],[.01 .02],[.2 .2]);
axes(ha(1))
% subplot(2,2,1)
plot(behav.position(:,1),behav.position(:,2),'k-')
axis tight
axis image
axis off
title(['Original data (n=',num2str(size(behav.position,1)),')'],'FontSize',10)

% subplot(2,2,2)
axes(ha(2))
plot(behav.position(frames_pass,1),behav.position(frames_pass,2),'k-')
axis tight
axis image
axis off
title(['Filtered data (n=',num2str(length(frames_pass)),')'],'FontSize',10)

% subplot(2,2,3)
axes(ha(3))
scatter(behav.position(frames_pass,1),behav.position(frames_pass,2),2,'k','filled')
axis tight
axis image
axis off
title('Filtered data (scatter plot)','FontSize',10)

axes(ha(4))
scatter(behav.position(:,1),behav.position(:,2),2,speed/scale,'filled')
% colormap(flipud(colormap(hot)))
% colormap(hot)
newmap = hot;
idx = find(newmap(:,3)==0,1,'last');
newmap(idx:end,:) = repmat([1 1 0],size(newmap,1)-idx+1,1);        %set that position to white
colormap(newmap);                %activate it
caxis([0,max(speed(frames_pass))/scale])
% caxis([0 50])
axis tight
axis image
axis off
title('Original data (colored by speed)','FontSize',10)

axes(ha(5))
scatter(behav.position(frames_pass,1),behav.position(frames_pass,2),2,speed(frames_pass)/scale,'filled')
% colormap(hot)
colormap(newmap);  
caxis([0,max(speed(frames_pass))/scale])
% caxis([0 60])
axis tight
axis image
axis off
title('Filtered data (colored by speed)','FontSize',10)

% % compare the speed before and after filtering
% figure
% subplot(1,2,1)
% plot(behav.time/1000,speed/scale,'m-')
% xlim([0,max(behav.time/1000)]);ylim([0,max(speed/scale)+2])
% hold on
% line(get(gca,'Xlim'),[behav.filterPara.thresh_speed_trans,behav.filterPara.thresh_speed_trans],'Color','k','LineStyle','--')
% line(get(gca,'Xlim'),[behav.filterPara.thresh_speed_min,behav.filterPara.thresh_speed_min],'Color','k','LineStyle','--')
% ylabel('Speed (cm/s)','FontSize',10)
% xlabel('Time (s)','FontSize',10)
% title('Original data','FontSize',10)
% subplot(1,2,2)
% plot(behav.time(frames_pass)/1000,speed(frames_pass)/scale,'m-')
% xlim([0,max(behav.time/1000)]);ylim([0,max(speed/scale)+2])
% hold on
% line(get(gca,'Xlim'),[behav.filterPara.thresh_speed_trans,behav.filterPara.thresh_speed_trans],'Color','k','LineStyle','--')
% line(get(gca,'Xlim'),[behav.filterPara.thresh_speed_min,behav.filterPara.thresh_speed_min],'Color','k','LineStyle','--')
% ylabel('Speed (cm/s)','FontSize',10)
% xlabel('Time (s)','FontSize',10)
% title('Filtered data','FontSize',10)

% correct the data
behav.position = behav.position(frames_pass,:);
behav.speed = speed(frames_pass);
behav.time = behav.time(frames_pass,:);
if isfield(behav,'positionblue')
    behav.positionblue = behav.positionblue(frames_pass,:);
end
if isfield(behav,'speedblue')
    behav.speedblue = behav.speedblue(frames_pass);
end
if isfield(behav,'hd')
    behav.hd = behav.hd(frames_pass,:);
end
% behav.vidNum = behav.vidNum(frames_pass);
% behav.frameNum = behav.frameNum(frames_pass);


