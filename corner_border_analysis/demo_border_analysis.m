%% identify border cells
%% calculate border score in the environment as a whole
load('D:\Miniscope_analysis_subiculum_shuttle\datapath.mat')
for n = 1:length(datapath)
    cd(datapath{n})
    load('firingrateAll.mat','firingrateAll')
    bscore = {};
    for ii = 1:length(firingrateAll)
        ratemap = firingrateAll{1,ii};
        ratemapS = cellfun(@(x) filter2DMatrices(x,1), ratemap,'uni',0);
        %remove rows that are all NaNs
        ratemapS = cellfun(@(x) x(~all(isnan(x),2),:), ratemapS, 'uni', 0);
        %remove cols that are all NaNs
        ratemapS = cellfun(@(x) x(:,~all(isnan(x),1)), ratemapS, 'uni', 0);
        %trim the edge
        ratemapS = cellfun(@(x) x(2:end-1,2:end-1), ratemapS, 'uni', 0);
        bs = cellfun(@borderScoreCalculation, ratemapS);
        bs = bs';
        bscore{ii} = bs;
    end
    save('border.mat','bscore')
end
%intersect border score with place cells to get border cell
for n = 1:length(datapath)
    cd(datapath{n})
    load('PlaceCellsLR.mat','placecellsLRm')
    load('border.mat','bscore')
    pc = cellfun(@(x,y) union(x,y), placecellsLRm(1,:),placecellsLRm(2,:),'uni',0);
    bc = cellfun(@(x) find(x > 0.45), bscore, 'uni',0);
    bcell = cellfun(@(x,y) intersect(x,y), pc, bc, 'uni',0);
    save('border.mat','bcell','-append')
end

%% calculate border score in splitted data in the shuttle box
load('D:\Miniscope_analysis_subiculum_shuttle\datapath.mat')
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat','frLR')
    tic
    bscoreLR = {};
    for ii = 1:size(frLR,1)
        for jj = 1:size(frLR,2)
            ratemap = frLR{ii,jj};
            ratemapS = cellfun(@(x) filter2DMatrices(x,1), ratemap,'uni',0);
            %remove rows that are all NaNs
            ratemapS = cellfun(@(x) x(~all(isnan(x),2),:), ratemapS, 'uni', 0);
            %remove cols that are all NaNs
            ratemapS = cellfun(@(x) x(:,~all(isnan(x),1)), ratemapS, 'uni', 0);
            %trim the edge
            ratemapS = cellfun(@(x) x(2:end-1,2:end-1), ratemapS, 'uni', 0);
            bs = cellfun(@borderScoreCalculation, ratemapS);
            bs = bs';
            bscoreLR{ii,jj} = bs;
        end
    end
    toc
    save('border.mat','bscoreLR','-append')
end
%intersect border score with place cells to get border cell
for n = 1:length(datapath)
    cd(datapath{n})
    load('PlaceCellsLR.mat','placecellsLRm')
    load('border.mat','bscoreLR')
    bc = cellfun(@(x) find(x > 0.45), bscoreLR, 'uni',0);
    bcellLR = cellfun(@(x,y) intersect(x,y), placecellsLRm, bc, 'uni', 0);
    save('border.mat','bcellLR','-append')
end

%% plot to visualize
load('neuronIndividualsf.mat');
load('behavIndividualsf.mat');
load('Thresh.mat');
load('firingrateAll.mat','firingrateAll')
load('border.mat','bcellLR')
x = bcellLR{1,1};
plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,firingrateAll,...
    x,thresh,'S',2);

%% identify corner coding cells
%% corner coding cells, maybe a specialized border cell/place cell
load('D:\Miniscope_analysis_subiculum_shuttle\datapath.mat')
for n = 1:length(datapath)
    cd(datapath{n})
    load('NBindivLR.mat', 'frLR')
    load('PlaceCellsLR.mat','placecellsLRm')
    ratemap = frLR;
    cmetrics = cell(size(ratemap,1),size(ratemap,2));
    for ii = 1:size(ratemap,1)
        for jj = 1:size(ratemap,2)
            cmetrics{ii,jj} = cellfun(@identify_corner_coding, ratemap{ii,jj});
        end
    end
    coridx = cellfun(@(x) find([x(:).acpeaks] >=7 & [x(:).acpeaks] <=9 & [x(:).atcorner] == 1),...
        cmetrics,'uni',false);
    placecell = cellfun(@(x,y) union(x,y), placecellsLRm(1,:),placecellsLRm(2,:),'uni',0);
    placecell = [placecell;placecell];
    corcell = cellfun(@(x,y) intersect(x,y), coridx, placecell, 'uni',0);
    save('corner.mat','cmetrics','corcell','-v7.3')
end
%% plot corner cells
load('corner.mat')
%plot to check selected neurons;
s2pr = 2; %session to plot, row
s2pc = 2; %session to plot, col
for  ii = 1:length(corcell{s2pr,s2pc})
    figure
    subplot(2,2,1)
    imagesc(cmetrics{s2pr, s2pc}(corcell{s2pr,s2pc}(ii)).ratemapS)
    axis image
    subplot(2,2,2)
    imagesc(cmetrics{s2pr, s2pc}(corcell{s2pr,s2pc}(ii)).ratemapS_thresh)
    axis image
    subplot(2,2,3)
    imagesc(cmetrics{s2pr, s2pc}(corcell{s2pr,s2pc}(ii)).BWS)
    axis image
    subplot(2,2,4)
    imagesc(cmetrics{s2pr, s2pc}(corcell{s2pr,s2pc}(ii)).ac)
    axis image
end
%check any neurons
cell2plot = 110;
figure
subplot(2,2,1)
imagesc(cmetrics{s2pr, s2pc}(cell2plot).ratemapS)
axis image
subplot(2,2,2)
imagesc(cmetrics{s2pr, s2pc}(cell2plot).ratemapS_thresh)
axis image
subplot(2,2,3)
imagesc(cmetrics{s2pr, s2pc}(cell2plot).BWS)
axis image
subplot(2,2,4)
imagesc(cmetrics{s2pr, s2pc}(cell2plot).ac)
axis image
% colormap jet
%% analyze corner cells
% compile data
load('D:\Miniscope_analysis_subiculum_shuttle\datapath.mat')
cor_metrics.Ctrl_fieldsc = []; cor_metrics.Ctrl_fieldsc_indiv = [];
cor_metrics.Ctrl_corboth = [];
cor_metrics.Ctrlstab = []; cor_metrics.Ctrlstab_indiv = [];
cor_metrics.CtrlcorrLR = []; cor_metrics.CtrlcorrLR_indiv = [];
for n = 1:length(Ctrl)
    cd(Ctrl{n})
    load('corner.mat')
    load('Metrics_match.mat', 'stability','corr_LRm_all')
    %number of fields in the contralateral box
    cmetricsfp = [cmetrics(2,:);cmetrics(1,:)]';
    pk = cellfun(@(x,y) [x(y).peaks], cmetricsfp, corcell, 'uni',false);
    cor_metrics.Ctrl_fieldsc_indiv = [cor_metrics.Ctrl_fieldsc_indiv;[pk{:}]'];
    cor_metrics.Ctrl_fieldsc = [cor_metrics.Ctrl_fieldsc;cellfun(@mean, pk)];
    %percent of corner cells that tune to both contexts
    cor_metrics.Ctrl_corboth = [cor_metrics.Ctrl_corboth; cellfun(@(x,y) ...
        length(intersect(x,y))/length(union(x,y)), corcell(1,:),corcell(2,:))];
    %stability
    stab = cellfun(@(x,y) x(y), stability, corcell,'uni',0);
    cor_metrics.Ctrlstab_indiv = [cor_metrics.Ctrlstab_indiv; vertcat(stab{:})];
    stab = cellfun(@mean, stab);
    cor_metrics.Ctrlstab = [cor_metrics.Ctrlstab; stab];
    %LR correlation
    corr_LRm = mat2cell(corr_LRm_all,size(corr_LRm_all,1),[1,1]);
    ccell = cellfun(@(x,y) union(x,y), corcell(1,:), corcell(2,:), 'uni', 0);
    corrLRpc = cellfun(@(x,y) x(y), corr_LRm, ccell,'uni',0);
    corrLRpc = padcat(corrLRpc{1},corrLRpc{2});
    cor_metrics.CtrlcorrLR_indiv = [cor_metrics.CtrlcorrLR_indiv;...
        corrLRpc];
    cor_metrics.CtrlcorrLR = [cor_metrics.CtrlcorrLR;nanmean(corrLRpc)];
end
cor_metrics.KO_fieldsc = [];cor_metrics.KO_fieldsc_indiv = [];
cor_metrics.KO_corboth = [];
cor_metrics.KOstab = []; cor_metrics.KOstab_indiv = [];
cor_metrics.KOcorrLR = []; cor_metrics.KOcorrLR_indiv = [];
for n = 1:length(KO)
    cd(KO{n})
    load('corner.mat')
    load('Metrics_match.mat', 'stability','corr_LRm_all')
    %number of fields in the contralateral box
    cmetricsfp = [cmetrics(2,:);cmetrics(1,:)]';
    pk = cellfun(@(x,y) [x(y).peaks], cmetricsfp, corcell, 'uni',false);
    cor_metrics.KO_fieldsc_indiv = [cor_metrics.KO_fieldsc_indiv;[pk{:}]'];
    cor_metrics.KO_fieldsc = [cor_metrics.KO_fieldsc;cellfun(@mean, pk)];
    %percent of corner cells that tune to both contexts
    cor_metrics.KO_corboth = [cor_metrics.KO_corboth; cellfun(@(x,y) ...
        length(intersect(x,y))/length(union(x,y)), corcell(1,:),corcell(2,:))];
    %stability
    stab = cellfun(@(x,y) x(y), stability, corcell,'uni',0);
    cor_metrics.KOstab_indiv = [cor_metrics.KOstab_indiv; vertcat(stab{:})];
    stab = cellfun(@mean, stab);
    cor_metrics.KOstab = [cor_metrics.KOstab; stab];
    %LR correlation
    corr_LRm = mat2cell(corr_LRm_all,size(corr_LRm_all,1),[1,1]);
    ccell = cellfun(@(x,y) union(x,y), corcell(1,:), corcell(2,:), 'uni', 0);
    corrLRpc = cellfun(@(x,y) x(y), corr_LRm, ccell,'uni',0);
    corrLRpc = padcat(corrLRpc{1},corrLRpc{2});
    cor_metrics.KOcorrLR_indiv = [cor_metrics.KOcorrLR_indiv;...
        corrLRpc];
    cor_metrics.KOcorrLR = [cor_metrics.KOcorrLR;nanmean(corrLRpc)];
end
save('D:\Miniscope_analysis_subiculum_shuttle\Results\cor_metrics.mat','cor_metrics')

%% plot corner cell results
load('D:\Miniscope_analysis_subiculum_shuttle\Results\cor_metrics.mat')
figure
%Stability
subplot(2,4,1)
histogram(cor_metrics.Ctrlstab_indiv(:),[-0.3:0.08:1],'normalization','probability')
hold on
histogram(cor_metrics.KOstab_indiv(:),[-0.3:0.08:1],'normalization','probability')
legend('Ctrl','KO')
xlabel('Stability (within session)')
ylabel('Prop of total neurons')
subplot(2,4,2)
data = [cor_metrics.Ctrlstab(:),cor_metrics.KOstab(:)];
plot(data','.','Color','k','MarkerSize',20)
xlim([0.5,2.5])
xticks([1,2])
xticklabels({'Ctrl','KO'})
ylabel('Stability')
[h,p,~,st] = ttest2(data(:,1),data(:,2));
text(0.75,0.45,['p=',num2str(p)])
text(0.75,0.4,['t=',num2str(st.tstat)])
%corrLR
subplot(2,4,3)
histogram(cor_metrics.CtrlcorrLR_indiv(:),[-0.4:0.1:1],'normalization','probability')
hold on
histogram(cor_metrics.KOcorrLR_indiv(:),[-0.4:0.1:1],'normalization','probability')
legend('Ctrl','KO')
xlabel('LR correlation')
ylabel('Proportion of total neurons')
subplot(2,4,4)
data = [cor_metrics.CtrlcorrLR(:),cor_metrics.KOcorrLR(:)];
plot(data','.','Color','k','MarkerSize',20)
xlim([0.5,2.5])
xticks([1,2])
xticklabels({'Ctrl','KO'})
ylabel('LR correlation')
ylim([0,0.6])
[h,p,~,st] = ttest2(data(:,1),data(:,2));
text(0.75,0.25,['p=',num2str(p)])
text(0.75,0.2,['t=',num2str(st.tstat)])
%number of fields in the contralateral context
subplot(2,4,5)
histogram(cor_metrics.Ctrl_fieldsc_indiv(:),'normalization','probability')
hold on
histogram(cor_metrics.KO_fieldsc_indiv(:),'normalization','probability')
legend('Ctrl','KO')
xlabel('Number of fields (contralateral)')
ylabel('Prop of total neurons')
subplot(2,4,6)
data = [cor_metrics.Ctrl_fieldsc(:),cor_metrics.KO_fieldsc(:)];
plot(data','.','Color','k','MarkerSize',20)
xlim([0.5,2.5])
xticks([1,2])
xticklabels({'Ctrl','KO'})
ylabel('Number of fields')
[h,p,~,st] = ttest2(data(:,1),data(:,2));
text(0.75,4.5,['p=',num2str(p)])
text(0.75,4.3,['t=',num2str(st.tstat)])
%proportion of place cells
subplot(2,4,7)
data = [cor_metrics.Ctrl_corboth(:),cor_metrics.KO_corboth(:)];
plot(data','.','Color','k','MarkerSize',20)
xlim([0.5,2.5])
xticks([1,2])
xticklabels({'Ctrl','KO'})
ylabel('Prop of total corner cells')
ylim([0 0.3])
[h,p,~,st] = ttest2(data(:,1),data(:,2));
text(1.2,0.25,['p=',num2str(p)])
text(1.2,0.2,['t=',num2str(st.tstat)])


