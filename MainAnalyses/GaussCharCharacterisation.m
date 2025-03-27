%% Visual speed tuning class characterisation
% Figure 1

%% get excitatory and suppressive responses

FRscale = 5;

for iunit = 1:numel(goodUnits)
    
    % stationary
    goodUnits(iunit).sigChangeP_stat = nan(4,7);
    goodUnits(iunit).sigChangeBool_stat = nan(4,7);
    goodUnits(iunit).ExcBool_stat = nan(4,7);
    goodUnits(iunit).nSig_stat = nan(4,1);
    goodUnits(iunit).nExc_stat = nan(4,1);
    goodUnits(iunit).nSupp_stat = nan(4,1);


    % find if stimulus evoked FR is higher or lower than baseline
    goodUnits(iunit).ExcBool_stat = double(cellfun(@mean, goodUnits(iunit).spikeCounts_stat)>(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_stat)*FRscale));
    goodUnits(iunit).ExcBool_stat(goodUnits(iunit).ExcBool_stat==0)=-1; % change zeros to minus 1

    % use sign-rank test to get sig difference p value for baseline vs
    % stimulus evoked


    for irow = 1:size(goodUnits(iunit).baselineSpikeCounts_stat,1)
        for icol = 1:size(goodUnits(iunit).baselineSpikeCounts_stat,2)
            try
                goodUnits(iunit).sigChangeP_stat(irow,icol) = ...
                    signrank(goodUnits(iunit).spikeCounts_stat{irow,icol},...
                    goodUnits(iunit).baselineSpikeCounts_stat{irow,icol}.*FRscale);
            catch
            end

        end
    end

    % convert to boolean significance test using bonferroni corrected p<0.05
    goodUnits(iunit).sigChangeBool_stat = goodUnits(iunit).sigChangeP_stat<0.05/7; % bonferroni correction
    goodUnits(iunit).nSig_stat = sum(goodUnits(iunit).sigChangeBool_stat,2)'; % number of significant responses
    goodUnits(iunit).nExc_stat = sum(goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==1,2)'; % number of sig excitatory responses
    goodUnits(iunit).nSupp_stat = sum(goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==-1,2)'; % number of sig suppressive responses



    % locomotion
    goodUnits(iunit).sigChangeP_run = nan(4,7);
    goodUnits(iunit).sigChangeBool_run = nan(4,7);
    goodUnits(iunit).ExcBool_run = nan(4,7);
    goodUnits(iunit).nSig_run = nan(4,1);
    goodUnits(iunit).nExc_run = nan(4,1);
    goodUnits(iunit).nSupp_run = nan(4,1);


    % find if stimulus evoked FR is higher or lower than baseline
    goodUnits(iunit).ExcBool_run = double(cellfun(@mean, goodUnits(iunit).spikeCounts_run)>(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_run)*FRscale));
    goodUnits(iunit).ExcBool_run(goodUnits(iunit).ExcBool_run==0)=-1; % change zeros to minus 1

    % use sign-rank test to get sig difference p value for baseline vs
    % stimulus evoked


    for irow = 1:size(goodUnits(iunit).baselineSpikeCounts_run,1)
        for icol = 1:size(goodUnits(iunit).baselineSpikeCounts_run,2)
            try
                goodUnits(iunit).sigChangeP_run(irow,icol) = ...
                    signrank(goodUnits(iunit).spikeCounts_run{irow,icol},...
                    goodUnits(iunit).baselineSpikeCounts_run{irow,icol}.*FRscale);
            catch
            end

        end
    end

    % convert to boolean significance test using bonferroni corrected p<0.05
    goodUnits(iunit).sigChangeBool_run = goodUnits(iunit).sigChangeP_run<0.05/7; % bonferroni correction
    goodUnits(iunit).nSig_run = sum(goodUnits(iunit).sigChangeBool_run,2)'; % number of significant responses
    goodUnits(iunit).nExc_run = sum(goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==1,2)'; % number of sig excitatory responses
    goodUnits(iunit).nSupp_run = sum(goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==-1,2)'; % number of sig suppressive responses


end


%% get basic properties

allTuning = cat(1,goodUnits.tuning_stat,goodUnits.tuning_run);
allSSI = cat(1,goodUnits.SSI_stat,goodUnits.SSI_run);

allr2 = [[goodUnits.r2_stat], [goodUnits.r2_run]]';
allr2p = [[goodUnits.r2pval_stat], [goodUnits.r2pval_run]]';
allGC = [[goodUnits.gaussChar_stat], [goodUnits.gaussChar_run]]';
allNSig = [[goodUnits.nSig_stat], [goodUnits.nSig_run]]';
allNExc = [[goodUnits.nExc_stat], [goodUnits.nExc_run]]';
allNSupp = [[goodUnits.nSupp_stat], [goodUnits.nSupp_run]]';

allPrefSpeed = [[goodUnits.prefSpeed_stat],[goodUnits.prefSpeed_run]]';

r2p_thresh = 0.05;
r2_thresh=0.1;

for igc = 1:4


    idx = find(allr2>=r2_thresh & allr2p<r2p_thresh & allGC==igc);
    gc(igc).r2vals = allr2(idx);
    gc(igc).nSigvals = allNSig(idx);
    gc(igc).nExcvals = allNExc(idx);
    gc(igc).nSuppvals = allNSupp(idx);
    gc(igc).prefSpeed = allPrefSpeed(idx);
    gc(igc).allTuning = allTuning(idx,:);
    gc(igc).allSSI = allSSI(idx,:);
    [~,gc(igc).prefSSI] = max(gc(igc).allSSI,[],2)
    

end


%% average tuning curves and SSI for each gauss char

plotOrder = [1 3 4 2];

figure
ii=0;
for igc = plotOrder
    ii=ii+1;
    subplot(2,4,ii)
    shadedErrorBar(1:7, mean(gc(igc).allTuning), sem(gc(igc).allTuning,1))
    ylim([6 18]), ax = gca; ax.XTick= 1:7; ax.YTick = 6:2:18;
    defaultAxesProperties(gca, true)
end

ii=0;
%average SSI
for igc = plotOrder
    ii=ii+1;
    subplot(2,4,ii+4)
    shadedErrorBar(1:7, mean(gc(igc).allSSI), sem(gc(igc).allSSI,1))
    ylim([0.7 0.9]), ax = gca; ax.XTick= 1:7; ax.YTick = 0.7:0.1:0.9;
    defaultAxesProperties(gca, true)
end

%% R2 vals by gc

vals = {gc(plotOrder(1)).r2vals, gc(plotOrder(2)).r2vals, gc(plotOrder(3)).r2vals, gc(plotOrder(4)).r2vals};
xvals = 1:4;
% violin plot
figure
distributionPlot(vals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', 0, 'distWidth', 1, 'Color', [0 0 0],...
    'showMM', 6)
ylim([0.1 1]); ax = gca; ax.YTick = 0.1:0.1:1;
defaultAxesProperties(gca, false);

vals =cat(1,gc.r2vals);
nVals = cellfun(@numel, {gc.r2vals});
group = repelem(1:4, nVals(1:4));
[p, anovatable, stats] = kruskalwallis(vals,group(:))
pArray = nan(3,4);
pvals = multcompare(stats,'CriticalValueType','bonferroni');
for irow = 1:6
    pArray(pvals(irow,1), pvals(irow,2)) = pvals(irow,6);
end


%% excitation and suppression by tuning class

for igc = 1:4
    SigVals_mean(igc) = mean(gc(igc).nSigvals);
    SigVals_sem(igc) = sem(gc(igc).nSigvals);

    ExcVals_mean(igc) = mean(gc(igc).nExcvals);
    ExcVals_sem(igc) = sem(gc(igc).nExcvals);

    SuppVals_mean(igc) = mean(gc(igc).nSuppvals);
    SuppVals_sem(igc) = sem(gc(igc).nSuppvals);
end

figure, hold on
% bar(1:4, SigVals_mean(plotOrder))
bar(1:4, [ExcVals_mean(plotOrder); SuppVals_mean(plotOrder)],'stacked')
ylabel('# Significant responses')

vals =cat(1,gc.nSigvals);
nVals = cellfun(@numel, {gc.nSigvals});
group = repelem(1:4, nVals(1:4));
[p, anovatable, stats] = kruskalwallis(vals,group(:))
pArray = [];
pvals = multcompare(stats,'CriticalValueType','bonferroni');
for irow = 1:6
    pArray(pvals(irow,1), pvals(irow,2)) = pvals(irow,6);
end


%% bar plot of p(gc)

gcNums = cellfun(@numel, {gc.r2vals});
allChar = gcNums./sum(gcNums);

figure
bar(1:4, allChar(plotOrder))
defaultAxesProperties(gca,false)

%% bar plot of pref speed

allPref = cat(1,gc.prefSpeed);
allChar = repelem(1:4, cellfun(@numel, {gc.prefSpeed}));

ii=0;
for ichar = plotOrder
    ii=ii+1;
    gcPrefs = allPref(allChar==ichar);
    for ispeed = 1:7
        gcSpeedVals(ispeed,ii)=sum(gcPrefs==ispeed);
    end
end


figure
bar(gcSpeedVals./sum(gcSpeedVals(:)),'stacked')
defaultAxesProperties(gca,false)
title('Preferred speeds')

%% p(SSI,given pref speed) conditional distribution

allPrefSSI = cat(1,gc.prefSSI);
allPrefSpeeds = cat(1,gc.prefSpeed);

vals = histcounts2(allPrefSSI(:), allPrefSpeeds(:));
vals_norm = vals./sum(vals,1);
figure
imagesc(vals_norm), colorbar, colormap(cmocean('tempo')), axis xy
xlabel('Pref speed'), ylabel('Pref SSI'), caxis([0, 0.5])


%% plot example units and population averages

plotOrder = [1 3 4 2];

gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));
idx=[];
% IDs = [951005818,951167247,950997203,951002872];
IDs = [951005818,951167265,950997203,951002872];

IDs = IDs(plotOrder);
for iid = 1:numel(IDs)
    idx(iid)=find(cat(1,areaUnits.ID)==IDs(iid));
end

idirs = [2 4 2 4];

figure

for iunit = 1:numel(idx)
subplot(4,4,iunit), hold on
    for ispeed = 1:7
        plot(repelem(ispeed,1,numel(areaUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed})),...
            areaUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed}, 'k.', 'MarkerSize',10)
    end
    plot(1:0.01:7, feval(gaussFun, areaUnits(idx(iunit)).gaussParams_run(idirs(iunit),:),  1:0.01:7), 'r-')
    ax = gca; ax.XTick = 1:7;
    defaultAxesProperties(gca, true)
    
    subplot(4,4,iunit+4)
    plot(1:7, areaUnits(idx(iunit)).SSI_run(idirs(iunit),:), 'k-')
    ylim([0.5 1.9]); 
    ax= gca; ax.YTick = 0.5:0.2:1.9; ax.XTick = 1:7;
    defaultAxesProperties(gca, true)
end

ii=0;
for igc = plotOrder
    ii=ii+1;
    subplot(4,4,8+ii)
    shadedErrorBar(1:7, mean(gc(igc).allTuning), sem(gc(igc).allTuning,1))
    ylim([6 18]), ax = gca; ax.XTick= 1:7; ax.YTick = 6:2:18;
    defaultAxesProperties(gca, true)
end

%average SSI
%average SSI
ii=0;
for igc = plotOrder
    ii=ii+1;
    subplot(4,4,12+ii)
    shadedErrorBar(1:7, mean(gc(igc).allSSI), sem(gc(igc).allSSI,1))
    ylim([0.7 0.9]), ax = gca; ax.XTick= 1:7; ax.YTick = 0.7:0.1:0.9;
    defaultAxesProperties(gca, true)
end



%% looking for example cells

% allChar = cat(1,areaUnits.gaussChar_run);
% allr2 = cat(1,areaUnits.r2_run);
% figure
% [idx, idirs ] = find(allChar==2 & allr2>0.65);
% 
% for iunit = 21:numel(idx)
%     figure
%     iunit
%     subplot(121), hold on
%     for ispeed = 1:7
% 
%         plot(repelem(ispeed,1,numel(areaUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed})),...
%             areaUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed}, 'k.', 'MarkerSize',10)
%     end
% 
%     plot(1:0.01:7, feval(gaussFun, areaUnits(idx(iunit)).gaussParams_run(idirs(iunit),:),  1:0.01:7), 'r-')
%     % ax = gca; ax.XTick = 1:7;
%     defaultAxesProperties(gca, true)
% 
%     subplot(122)
%     plot(1:7, areaUnits(idx(iunit)).SSI_run(idirs(iunit),:), 'k-')
%     ylim([0.5 1.9]); 
%     % ax= gca; ax.YTick = 0.5:0.2:1.9; ax.XTick = 1:7;
%     defaultAxesProperties(gca, true)
%     pause
%     close all
% end