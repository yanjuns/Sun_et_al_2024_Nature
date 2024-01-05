function plot_ratemap(FRS, cellnumber, score, colbar, scale, nperfig, nfigcol)
% This function is used to plot firing rate map of linear track virtual
% reality data
% Inputs:
% FRS: cell array, smoothed firing rate of all the neurons
% cellID: a vector of identified good neurons
% cellnumber: range of plot, all neurons by default
% colbar: logical, plot with color bar or not
% scale: define max heatmap color scale, if 'auto', color scale align to
% max peak firing rate of each neuron
% nperfig: max number of neurons plot in each figure
% nfigcol: max number of columns plot in each figure
% @Yanjun Sun, Stanford University, 8/10/2019

if ~exist('nfigcol', 'var') || isempty(nfigcol)
    nfigcol = 8;
end
if ~exist('nperfig', 'var') || isempty(nperfig)
    nperfig = 40;
end
if ~exist('scale', 'var') || isempty(scale)
    scale = 'auto';
end
if ~exist('colbar', 'var') || isempty(colbar)
    colbar = false;
end
if ~exist('score', 'var') || isempty(score)
    score = [];
end
if ~exist('cellnumber', 'var') || isempty(cellnumber)
    cellnumber = 1:length(FRS);
end
%% create a folder to save figures
folderName = 'RatemapFigures';
if ~exist(folderName,'dir')
    mkdir(folderName);
end
fpath=folderName;
%% plot rate map figures
% get the size of monitor to define figure size
Pix_SS = get(0,'screensize');
% start to plot
n = 0;
for ii = 1:ceil(length(cellnumber)/nperfig)
k = 0;
ax = figure;
set(ax, 'Position', Pix_SS);
hold on;
    % inside each figure
    for jj = (n*nperfig+1) : (n*nperfig+nperfig)
        cellnum = cellnumber(jj);
        k = k+1;
        firingRateSmoothing = FRS{cellnum};
%         x = size(firingRateSmoothing,2);
%         y = size(firingRateSmoothing,1);
        maxCount = scale;
        if strcmpi(scale,'auto')
            if isempty(firingRateSmoothing)
                 maxCount = 1;
            else
                 maxCount = round(max(max(firingRateSmoothing)),2);
            end
        end

        subplot(nperfig/nfigcol,nfigcol,k);
%         xlim([1 x]);
%         ylim([1 y]);
        hold on
        pcolor(firingRateSmoothing);
        colormap(jet);
        caxis([0,maxCount]);
        shading flat
        axis image
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'fontname','Arial')
        h = gca;
        h.XAxis.Visible = 'off';
        h.YAxis.Visible = 'off';
        title(['cell ' num2str(cellnum)], 'FontName', 'Arial', 'Fontsize', 12);
        if colbar
            cb = colorbar('eastoutside');
            if strcmpi(scale,'auto')
                set(cb,'Ticks',[0, maxCount]);
            end
        else
            text(min(xlim), min(ylim-1), ['fr: ', num2str(maxCount), 'Hz'], 'FontName', 'Arial', 'Fontsize', 12);
        end
        if ~isempty(score)
            cscore = round(score(k),2);
           text(min(xlim+12), min(ylim-1), ['c: ', num2str(cscore)], 'FontName', 'Arial', 'Fontsize', 12);
        end
        if jj == length(cellnumber)
            break
        end
    end

% if strcmpi(scale,'auto')
%     saveas(gcf,fullfile(fpath,strcat('Cell',num2str(n*nperfig+1), 'to', num2str(n*nperfig+nperfig),'_ratemap_auto','.fig')));
% else
%     saveas(gcf,fullfile(fpath,strcat('Cell',num2str(n*nperfig+1), 'to', num2str(n*nperfig+nperfig),'_ratemap','.fig')));
% end

n = n+1;
end

end