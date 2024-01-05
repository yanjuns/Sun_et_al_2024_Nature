function plot_eb_map_longitudinal(neuronIndividualsf,thresh,cellnumber,dataType,savefigure)
%% This function is used to plot egocentric bearing of neurons to environment corners
% Inputs: hd, a cell array
if ~exist('savefigure','var')||isempty(savefigure)
    savefigure = true;
end
if ~exist('dataType','var')||isempty(dataType)
    dataType = 'eb';
end

for jj = 1:length(cellnumber)
    cellnum = cellnumber(jj);
    ax = figure;
    set(ax, 'Position', [0, 200, 1200, 800]);
    hold on;
    k = 0;
    for ii = 1: length(neuronIndividualsf)
        neuron = neuronIndividualsf{ii};
        dataFiring = neuron.S;
        thresh0 = thresh(cellnum);
        idx = dataFiring(cellnum,:)>thresh0;
        % color coded raster plot
        k = k+1;
        subplot(round(length(neuronIndividualsf)/2),4,2*k-1);
        plot(neuron.pos(:,1),neuron.pos(:,2),'Color',[109 110 113]/255);
        hold on;
        axis image;
        switch dataType
            case {'eb'}
                scatter(neuron.pos(idx,1),neuron.pos(idx,2),20,neuron.eb(idx),'filled')
            case {'hd'}
                scatter(neuron.pos(idx,1),neuron.pos(idx,2),20,neuron.hd(idx),'filled')
            case {'wallb'}
                scatter(neuron.pos(idx,1),neuron.pos(idx,2),20,neuron.wallb(idx),'filled')
            case {'ctrb'}
                scatter(neuron.pos(idx,1),neuron.pos(idx,2),20,neuron.ctrb(idx),'filled')
        end
        colormap hsv
        % plot egocentric bearing histogram
        subplot(round(length(neuronIndividualsf)/2),4,2*k);
        switch dataType
            case {'eb'}
                polarhistogram(neuron.eb(idx),20)
            case {'hd'}
                polarhistogram(neuron.hd(idx),20)
            case {'wallb'}
                polarhistogram(neuron.wallb(idx),20)
            case {'ctrb'}
                polarhistogram(neuron.ctrb(idx),20)
        end
    end
    suptitle(['Cell' num2str(cellnum)]);
    
    if savefigure
        folderName = 'RatemapFigures';
        if ~exist(folderName,'dir')
            mkdir(folderName);
        end
        fpath=folderName;
        saveas(gcf,fullfile(fpath,strcat('Cell',num2str(cellnum),'_ebmap','.fig')));
    end
end
end
