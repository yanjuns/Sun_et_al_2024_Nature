function plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,firingrateAll,cellnumber,thresh,dataType,scale,savefigure)
%% This function is used to plot firing rate map of longitudinally tracked neurons
% Inputs:
% (1) neuronIndividualsf: a cell array of neuron data across different days, each cell is a source2D variable. Obtained by CNMF-E
% (2) behavIndividualsf: a cell array of filtered behavior data, obtained by using Tristan's code
% (3) firinrateAll: a cell array of firing rates of all the neurons across all the recording days
% (4) cellnumber: IDs of the cell that want to be plotted, e.g. 1:size(neuron.trace,1)
% (5) thresh: a vector contains cut off threshold of calcium signals for each neuron, usually defined as 2*SD of neuron.S or neuron.trace
% (6) dataType: 'S' or 'trace'
% (7) scale: max heatmap ploting scale (Hz); use 'auto' if want to plot as the max firing rate of each individual cell
if ~exist('savefigure','var')||isempty(savefigure)
    savefigure = true;
end

for jj = 1:length(cellnumber);
cellnum = cellnumber(jj);
ax = figure;
set(ax, 'Position', [0, 200, 1200, 800]);
hold on;
k = 0;

for ii = 1: length(neuronIndividualsf);
    neuron = neuronIndividualsf{ii};
    behav = behavIndividualsf{ii};
    firingRate = firingrateAll{ii};
    if ~exist('dataType','var') || isempty(dataType)
        dataType = 'S';
    end

    if strcmpi(dataType,'trace')
        dataFiring = neuron.trace;
    elseif strcmpi(dataType,'S')
        dataFiring = neuron.S;
    end

%% Calculate neuron.pos for plotting
    downsampling = length(neuron.time)/size(neuron.S,2);
    if downsampling ~= 1
        neuron.time = double(neuron.time);
        neuron.time = neuron.time(1:downsampling:end);
    end
    temp = find(diff(behav.time)<=0);
    behav.time(temp+1) = behav.time(temp)+1;
    neuron.pos = interp1(behav.time,behav.position,neuron.time); %%
    
    thresh0 = thresh(cellnum);
    idx = dataFiring(cellnum,:)>thresh0;
    
    %% plotting animal behavior trajectries
    k = k+1;
    subplot(round(length(neuronIndividualsf)/2),4,2*k-1);
    plot(neuron.pos(:,1),neuron.pos(:,2),'k');
    hold on;
    plot(neuron.pos(idx,1),neuron.pos(idx,2),'r.', 'MarkerSize', 7);
%     title(['Cell' num2str(cellnum)],'FontSize',8,'FontName','Arial')
    set(gca,'FontSize',8);
    set(gca,'xticklabel',{[]});
    axis image;
%     xlim([min(neuron.pos(:,1)) max(neuron.pos(:,1))]);
%     ylim([min(neuron.pos(:,2)) max(neuron.pos(:,2))]);
%     plot(ms.pos(idx2,1),ms.pos(idx2,2),'r.')
%     hold off
%     set(gca,'Xtick',[])
%     if mod(k,numFig) == 0 || k == max(cellnum)
%         set(gca,'Xtick',[1 ceil(neuron.num2read/2) neuron.num2read])
%     end
    
    %% plot cell firing map
%     firingRateSmoothing = filter2DMatrices(firingRate{cellnum}, 1);
    firingRateSmoothing = firingRate{cellnum};
    maxCount = scale;
    if strcmpi(scale,'auto')
        if isempty(firingRateSmoothing)
             maxCount = 1;
        else
             maxCount = max(max(firingRateSmoothing));
             if maxCount > 2
                 maxCount = 2;
             end
        end
    end
    
    subplot(round(length(neuronIndividualsf)/2),4,2*k);
    pcolor(firingRateSmoothing);
    colormap(jet);
    caxis([0,maxCount]);
    cb = colorbar('eastoutside');
    if strcmpi(scale,'auto') && maxCount > 0
        set(cb,'Ticks',[0, maxCount]);
    end
    set(gca,'xticklabel',{[]});
    shading flat;
    axis image;
%     title(['Cell', num2str(cellnum)],'FontName','Arial','FontSize',8,'FontWeight','bold')
% 
%     if mod(k,numFig) ~= 0 
%         set(gca,'Xtick',[])
%     end
%     if mod(k,numFig) == 0 
%         colorbar('eastoutside');
%     end

%     if mod(k,numFig) == 0 || k == max(cellnum)
%         saveas(gcf,fullfile(fpath,strcat(experimentName,'CellFiringBeaviorSpatial',num2str(i),'.fig')))
%     end
end
suptitle(['Cell' num2str(cellnum)]);

if savefigure
    folderName = 'RatemapFigures';
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    fpath=folderName;

    if strcmpi(scale,'auto')
        saveas(gcf,fullfile(fpath,strcat('Cell',num2str(cellnum),'_ratemap_auto','.fig')));
    else
        saveas(gcf,fullfile(fpath,strcat('Cell',num2str(cellnum),'_ratemap','.fig')));
    end
end

end
end
