function [pkfrcc,pkfrccx] = find_pkfr_eachfield(datapath)
%This function aims to find the peak firing rate for each fields in each
%corner cell. 
%Yanjun Sun, yanjuns@stanford.edu, 12/9/2022

if ~exist('datapath','var')||isempty(datapath)
    datapath = pwd;
end
cd(datapath)

load('corner_metrics.mat','C','cornercellx')
fieldloc = cellfun(@(x) x.field_coor, C.cmetrics, 'uni',0);
%pkfr for each field in the corner cell
load('firingrateAll.mat','ratemap')

% for corner cells from each session
pkfrcc = cell(1,length(ratemap));
for n = 1:length(ratemap)
    cornercell = C.cornercell{n};
    if isempty(cornercell)
        pkfrcc{n} = [];
    else
        ratemapss = ratemap{n}(cornercell);
        floc = fieldloc{n}(cornercell);
        pkfr = cell(1,length(cornercell));
        for ii = 1:length(ratemapss)
            ratemap_ea = ratemapss{ii};
            x = [];
            for jj = 1:size(floc{ii},1)
                x = [x;ratemap_ea(floc{ii}(jj,2),floc{ii}(jj,1))];
            end
            pkfr{ii} = x;
        end
        pkfrcc{n} = pkfr;
    end
end

% for cornercellx
pkfrccx = cell(1,length(ratemap));
for n = 1:length(ratemap)
    cornercell = cornercellx;
    if isempty(cornercell)
        pkfrccx{n} = [];
    else
        ratemapss = ratemap{n}(cornercell);
        floc = fieldloc{n}(cornercell);
        pkfr = cell(1,length(cornercell));
        for ii = 1:length(ratemapss)
            ratemap_ea = ratemapss{ii};
            x = [];
            for jj = 1:size(floc{ii},1)
                x = [x;ratemap_ea(floc{ii}(jj,2),floc{ii}(jj,1))];
            end
            pkfr{ii} = x;
        end
        pkfrccx{n} = pkfr;
    end
end
    
end

