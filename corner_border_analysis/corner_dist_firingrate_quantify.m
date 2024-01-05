function d2cornerq = corner_dist_firingrate_quantify(session_name, d2cornerb, C, datapath)
%quantify and average data for each session type
%yanjuns@stanford.edu
if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)
if ~exist('C','var')||isempty(C)
    load('corner_metricsR1.mat','C','cornercellx')
end
% %expA
% load('corner_metrics.mat','d2cornerb','C','cornercellx')
% %expC
% load('corner_metrics.mat','d2cornerb','C2','cornercellx')
% C = C2;
% %expC concave
% load('corner_metrics.mat','d2cornerb_ccav','C2','cornercellx')
% d2cornerb = d2cornerb_ccav;
% C = C2;
if ~exist('d2cornerb','var')||isempty(d2cornerb)
    load('corner_metrics.mat','d2cornerb')
end
load('env_geometry.mat','S')
d2cornerq = struct; %quantification
d2cornerq.d = cell(1,length(session_name));
d2cornerq.frcc = cell(1,length(session_name));
d2cornerq.frncc = cell(1,length(session_name));
for ii = 1:length(session_name)
    idx = strcmp(S, session_name{ii});
    % for distance
    d2c_d = d2cornerb.d(idx);
    d = [];
    for jj = 1:length(d2c_d)
        d = padconcatenation(d,d2c_d{jj},1);
    end
    d2cornerq.d{ii} = d;
    % for firing rate
    d2c_fr = d2cornerb.fr(idx);
    allcell = 1:length(d2c_fr{1,1});
    cornercell = C.cornercell(idx);
    fr = [];frncc = [];frx = [];
    for jj = 1:length(d2c_fr)
        % for corner cell
        frss = nanmean(d2c_fr{jj}(cornercell{jj},:),1);
        fr = padconcatenation(fr,frss,1);
        % for non-corner cell
        ind = ~ismember(allcell, cornercell{jj});
        ncc = allcell(ind);
        frssncc = nanmean(d2c_fr{jj}(ncc,:),1);
%         frssncc = nanmean(d2c_fr{jj}); %use all neurons
        frncc = padconcatenation(frncc,frssncc,1);
        % for tracked corner cell
        if ~exist('cornercellx','var')
            frx = [];
        else
            frssx = nanmean(d2c_fr{jj}(cornercellx,:),1);
            frx = padconcatenation(frx,frssx,1);
        end
    end
    d2cornerq.frcc{ii} = fr;
    d2cornerq.frncc{ii} = frncc;
    d2cornerq.frccx{ii} = frx;
end

end

