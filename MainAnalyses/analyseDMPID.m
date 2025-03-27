%% analyse decoding

%% load the data

dataDir = 'D:\AllenMatFiles\FC';

dataFiles = dir(fullfile(dataDir, 'dmPID*'));
for ifile = 1:numel(dataFiles)
    dataFiles(ifile).sessionID = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'));
    dataFiles(ifile).filename = fullfile(dataFiles(ifile).folder, dataFiles(ifile).name);
end

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);



for isession = 1:numel(dataFiles)

    load(dataFiles(isession).filename)
    session(isession).stat = stat;
    session(isession).run=run;
end


%% get data into useable format for analysing
% each pCorrect with area, session and dir

nDirs=4;

for iarea=1:8
ta(iarea).statMeanPerf=[];
ta(iarea).statDir = [];
ta(iarea).statSesh = [];

ta(iarea).runMeanPerf=[];
ta(iarea).runDir = [];
ta(iarea).runSesh = [];
end

for isession = 1:numel(dataFiles)
    for iarea = 1:8
            for iperm = 1:numel(session(isession).stat.ta(iarea).perm)
                ta(iarea).statMeanPerf = cat(2,ta(iarea).statMeanPerf, session(isession).stat.ta(iarea).perm(iperm).dir.mean_pCorrect);
                ta(iarea).statDir = cat(2,ta(iarea).statDir, 1:nDirs);
                ta(iarea).statSesh = cat(2,ta(iarea).statSesh, repelem(isession,1,nDirs));
            end
    end
end

for isession = 1:numel(dataFiles)
    for iarea = 1:8
            for iperm = 1:numel(session(isession).run.ta(iarea).perm)
                ta(iarea).runMeanPerf = cat(2,ta(iarea).runMeanPerf, session(isession).run.ta(iarea).perm(iperm).dir.mean_pCorrect);
                ta(iarea).runDir = cat(2,ta(iarea).runDir, 1:nDirs);
                ta(iarea).runSesh = cat(2,ta(iarea).runSesh, repelem(isession,1,nDirs));
            end
    end
end

%% statsc

for iarea = 1:8
    pvals(iarea) = ranksum(ta(iarea).statMeanPerf, ta(iarea).runMeanPerf);
end


%% paired box plots

figure, hold on
nestedCellArray = {{ta.statMeanPerf},{ta.runMeanPerf}};
colorCellArray = {{areacols.*0.7}, {areacols}};
subGroupSpacing = [-0.15, +0.15];
[bp, a] = pairedBoxPlot(nestedCellArray, subGroupSpacing, colorCellArray);
ax= gca;
ax.XTick = 1:8;
ylim([0 0.72])


%% box plots for stat/run separately

figure, 
subplot(121), hold on
nestedCellArray = {{ta.statMeanPerf}};
colorCellArray = {{areacols}};
subGroupSpacing = [0];
[bp, a] = pairedBoxPlot(nestedCellArray, subGroupSpacing, colorCellArray);
ax= gca;
ax.XTick = 1:8;
ylim([0 0.72])

subplot(122), hold on
nestedCellArray = {{ta.runMeanPerf}};
colorCellArray = {{areacols}};
subGroupSpacing = [0];
[bp, a] = pairedBoxPlot(nestedCellArray, subGroupSpacing, colorCellArray);
ax= gca;
ax.XTick = 1:8;
ylim([0 0.72])


%% GLME analysis

stat_perf = [ta.statMeanPerf];
stat_dir = [ta.statDir];
stat_session = [ta.statSesh];
stat_area = repelem(areas, cellfun(@numel, {ta.statMeanPerf}));
stat_state = zeros(size(stat_perf));

run_perf = [ta.runMeanPerf];
run_dir = [ta.runDir];
run_session = [ta.runSesh];
run_area = repelem(areas, cellfun(@numel, {ta.runMeanPerf}));
run_state = ones(size(run_perf));

perfVec = [stat_perf(:); run_perf(:)];
stateVec = [stat_state(:); run_state(:)];
sessionVec = [stat_session(:); run_session(:)];
areaVec = [stat_area(:); run_area(:)];
dirVec = [stat_dir(:); run_dir(:)];

tbl = table(perfVec, stateVec, sessionVec, areaVec, dirVec,...
    'VariableNames', {'perf', 'state', 'session', 'area', 'dir'});

tbl.session = categorical(tbl.session);
tbl.state = categorical(tbl.state);
tbl.dir = categorical(tbl.dir);

f = 'perf ~ -1 + area:state + (1|dir)';
lme = fitlme(tbl, f, 'DummyVarCoding', 'full')

H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1*2-1) = -1; H_temp(iarea1*2) = 1;
    pvals(iarea1) = coefTest(lme,H_temp);
end


[psi,dispersion,stats] = covarianceParameters(lme);
stats{1}
stats{2}

vals = [lme.Coefficients.Estimate, lme.Coefficients.Lower, lme.Coefficients.Upper];

% scatjit
figure, hold on
for iarea = 1:8;
    % background bar chart
    % bar(iarea-0.2, vals(iarea*2-1,1), 'FaceColor',  areacols(iarea,:)*0.7, 'EdgeColor', 'none','BarWidth',0.35, 'FaceAlpha',0.8)
    % bar(iarea+0.2, vals(iarea*2,1), 'FaceColor',  areacols(iarea,:), 'EdgeColor', 'none','BarWidth',0.35,'FaceAlpha',0.8)

    % errorbar
    % plot([iarea iarea], [vals(iarea,2), vals(iarea,3)], '-', 'Color', 'k', 'LineWidth', 1)
    % plot(iarea, vals(iarea,1), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2)
    % errorbar
    % scatjit data
    s1 = scatJit(ta(iarea).statMeanPerf(:), 0.1, iarea-0.2, 8, areacols(iarea,:).*0.7);
    s2 = scatJit(ta(iarea).runMeanPerf(:), 0.1, iarea+0.2, 8, areacols(iarea,:));

    plot([iarea-0.2 iarea-0.2], [lme.Coefficients.Estimate(iarea*2-1)-lme.Coefficients.SE(iarea*2-1), lme.Coefficients.Estimate(iarea*2-1)+lme.Coefficients.SE(iarea*2-1)],...
        'k', 'LineWidth', 1.5)
     plot([iarea-0.2], [lme.Coefficients.Estimate(iarea*2-1)], 'wo','MarkerFaceColor', 'k')
    
    plot([iarea+0.2 iarea+0.2], [lme.Coefficients.Estimate(iarea*2)-lme.Coefficients.SE(iarea*2), lme.Coefficients.Estimate(iarea*2)+lme.Coefficients.SE(iarea*2)], ...
    'k', 'LineWidth', 1.5)
     plot([iarea+0.2], [lme.Coefficients.Estimate(iarea*2)], 'wo','MarkerFaceColor', 'k')
        
end

ax = gca; ax.XTick = 1:8; ax.XTickLabels = areas; ax.YLim = [0 1];


%% dist plot

figure, hold on
for iarea = 1:8
    plotVals = {ta(iarea).statMeanPerf};%, run.ta(iarea).pCorr(:)};
    xvals = [iarea*2-0.3];
    col2use = areacols(iarea,:).*0.7;
    
    val = distributionPlot(plotVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', false, 'distWidth', 0.4, 'Color', col2use,...
    'showMM', 6);

    plotVals = {ta(iarea).runMeanPerf};%, run.ta(iarea).pCorr(:)};
    xvals = [iarea*2+0.3];
    col2use = areacols(iarea,:);
    
    distributionPlot(plotVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', false, 'distWidth', 0.4, 'Color', col2use,...
    'showMM', 6)
end



%% stat for decoding performance between areas within a state

% stat
statP=[];
runP=[];

for iarea1 = 1:8
    for iarea2 = 1:8
    H=zeros(1,16);
    H(iarea1*2-1)=-1;
    H(iarea2*2-1)=1;
    statP(iarea1,iarea2) = coefTest(lme,H);

    H=zeros(1,16);
    H(iarea1*2)=-1;
    H(iarea2*2)=1;
    runP(iarea1,iarea2) = coefTest(lme,H);

    end
end


statP(eye(size(statP))==1) = nan;
runP(eye(size(runP))==1) = nan;

statPLogic = statP<0.05/28;
runPLogic = runP<0.05/28;

figure
subplot(121)
imagesc(statPLogic); 
ax=gca; ax.XTickLabel = areas; ax.YTickLabel = areas;
caxis([0, 1]), colormap(gray), colorbar

subplot(122)
imagesc(runPLogic); 
ax=gca; ax.XTickLabel = areas; ax.YTickLabel = areas;
caxis([0, 1]), colormap(gray), colorbar




%% imsc plot for area diffs in each state
statDiff=nan(8);
runDiff=nan(8);

statVals = cellfun(@nanmean, {ta.statMeanPerf});
runVals = cellfun(@nanmean, {ta.runMeanPerf});

% stat
for iarea1 = 1:8
    for iarea2 = 1:8
        statDiff(iarea1,iarea2) = statVals(iarea1)-statVals(iarea2);
        runDiff(iarea1,iarea2) = runVals(iarea1)-statVals(iarea2);
    end
end

figure
subplot(121)
imagesc(statDiff); colormap(crameri('vik'));colorbar
subplot(122)
imagesc(runDiff);colormap(crameri('vik'));colorbar


%% imagesc version

statVals = cellfun(@nanmean, {ta.statMeanPerf});
runVals = cellfun(@nanmean, {ta.runMeanPerf});


figure
imagesc([statVals; runVals])
colormap(cmocean('tempo'))
clim([0.2 0.5])
colorbar

figure
imagesc([runVals-statVals]);
colormap(crameri('vik'))
clim([-0.12 0.12])
colorbar
