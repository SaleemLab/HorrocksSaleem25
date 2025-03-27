%% Analyse speed tuning strength

r2_thresh = 0.1; r2p_thresh = 0.05;
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);
%% Tuning strength p(tuned)

% find fraction tuned for each area in stat/run sessions and use GLME to
% compare statistically.

% because different sessions have different directions (both direction and
% #) we treat each direction of motion independently (i.e. a neuron can
% contribute up to 4 tuning curves).


for iarea = 1:numel(areas)
    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));
    % stat
    allr2 = [areaUnits.r2_stat]';
    allr2p = [areaUnits.r2pval_stat]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),4,1);
    allDir = repmat([1:4]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,4,1);

    statUnitsr2 = cat(1,areaUnits.r2_stat);
    ta(iarea).nStatUnits = sum(any(~isnan(statUnitsr2),2));



    % remove instances where insufficient trials to calculate
    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allDir(idx)=[];
    allSigRF(idx)=[];

    ta(iarea).r2_stat = allr2;
    ta(iarea).allr2p_stat = allr2p;

    ta(iarea).areaVec_stat = allArea;
    ta(iarea).sessionVec_stat = allSession;
    ta(iarea).stateVec_stat = zeros(size(ta(iarea).sessionVec_stat));
    ta(iarea).dirVec_stat = allDir;
    ta(iarea).sigRFVec_stat = allSigRF;

    ta(iarea).nTuningCurves_stat = numel(allr2p);
    ta(iarea).nTuned_stat = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_stat = ta(iarea).nTuned_stat/ta(iarea).nTuningCurves_stat;
    ta(iarea).tunedFlag_stat = ta(iarea).r2_stat>r2_thresh & ta(iarea).allr2p_stat<r2p_thresh;


    % run
    allr2 = [areaUnits.r2_run]';
    allr2p = [areaUnits.r2pval_run]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),4,1);
    allDir = repmat([1:4]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,4,1);

    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allDir(idx)=[];
    allSigRF(idx)=[];


    ta(iarea).r2_run = allr2;
    ta(iarea).allr2p_run = allr2p;
    ta(iarea).areaVec_run = allArea;
    ta(iarea).sessionVec_run = allSession;
    ta(iarea).stateVec_run = ones(size(ta(iarea).sessionVec_run));
    ta(iarea).dirVec_run = allDir;
    ta(iarea).sigRFVec_run = allSigRF;


    runUnitsr2 = cat(1,areaUnits.r2_run);
    ta(iarea).nRunUnits = sum(any(~isnan(runUnitsr2),2));

    ta(iarea).nTuningCurves_run = numel(allr2p);
    ta(iarea).nTuned_run = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_run = ta(iarea).nTuned_run/ta(iarea).nTuningCurves_run;
    ta(iarea).tunedFlag_run = ta(iarea).r2_run>r2_thresh & ta(iarea).allr2p_run<r2p_thresh;


    ta(iarea).stateVec = categorical(cat(1,ta(iarea).stateVec_stat,ta(iarea).stateVec_run));
    ta(iarea).tunedVec = cat(1,ta(iarea).tunedFlag_stat,ta(iarea).tunedFlag_run);
    ta(iarea).sessionVec = categorical(cat(1,ta(iarea).sessionVec_stat,ta(iarea).sessionVec_run));
    ta(iarea).dirVec = categorical(cat(1,ta(iarea).dirVec_stat, ta(iarea).dirVec_run));
    ta(iarea).RFVec = categorical(cat(1,ta(iarea).sigRFVec_stat, ta(iarea).sigRFVec_run));
    ta(iarea).areaVec = categorical(cat(1,ta(iarea).areaVec_stat, ta(iarea).areaVec_run));

end


% glme on all units combined

stateVec =cat(1,ta.stateVec);
tunedVec = cat(1,ta.tunedVec);
sessionVec = cat(1,ta.sessionVec);
dirVec = cat(1,ta.dirVec);
RFVec = cat(1,ta.RFVec);
areaVec = cat(1,ta.areaVec);

tbl = table(tunedVec, stateVec, sessionVec,dirVec,RFVec,areaVec,...
    'VariableNames', {'tunedFlag', 'state', 'session','dir','RF','area'});

f = 'tunedFlag ~ -1 + area:state + (1|session) + (1|dir) + (1|RF)';

glme = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

for iarea = 1:8
    H=zeros(1,16);
    H(iarea*2-1)=-1;
    H(iarea*2)=1;
    H;
    ta(iarea).pValue2 = coefTest(glme,H);
end

[ta.pValue2]

%% plot ptuned coefficients


% stationary
xvals = [1 3 4 5 7 8 10 11];

vals = glme.Coefficients.Estimate;
vals = exp(vals);
prob_vals = vals./(1+vals);

vals = glme.Coefficients.Estimate-glme.Coefficients.SE;
vals = exp(vals);
prob_vals_lower = vals./(1+vals);

vals = glme.Coefficients.Estimate+glme.Coefficients.SE;
vals = exp(vals);
prob_vals_upper = vals./(1+vals);

% stat
figure, subplot(2,2,1), hold on
for iarea = 1:8
    icoeff = iarea*2-1;
    bar(xvals(iarea), prob_vals(icoeff), 'FaceColor',  areacols(iarea,:), 'EdgeColor', 'none','FaceAlpha',0.7)
    plot(xvals(iarea), prob_vals(icoeff), 'wo', 'MarkerFaceColor', areacols(iarea,:))
    plot([xvals(iarea), xvals(iarea)], [prob_vals_lower(icoeff),...
        prob_vals_upper(icoeff)], 'Color',areacols(iarea,:))
end
ylim([0 0.5])
ax= gca; ax.XTick = xvals; ax.XTickLabels=areas;
defaultAxesProperties(gca,false)


% compare p(tuned) between areas for stationary sessions 
statP=ones(8);

for iarea1 = 1:8
    for iarea2 = 1:8
        if iarea1~=iarea2
    H=zeros(1,16);
    H(iarea1*2-1)=-1;
    H(iarea2*2-1)=1;
    statP(iarea1,iarea2) = coefTest(glme,H);
        end
    end
end

subplot(2,2,2)
tickvals = [0 0.0001 0.001 0.01 0.05 28];
pvals_disc = discretize(statP*28,tickvals,'IncludedEdge','right');
colmap = colormap(slanCM('gray',6));
colmap=flipud(colmap);
colmap = colmap(1:4,:);
colmap = [colmap; 0 0 0];
imagesc(pvals_disc), colormap(colmap); 
c = colorbar;
c.Ticks =linspace(1,5,6); c.TickLabels = {"","0.0001","0.001", "0.01", "0.05", ">0.05"};

ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
defaultAxesProperties(gca,false)
hold on
for iarea = 1:numel(areas)
for i=1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
end


% run
subplot(2,2,3), hold on
for iarea = 1:8
    icoeff = iarea*2;
    bar(xvals(iarea), prob_vals(icoeff), 'FaceColor',  areacols(iarea,:), 'EdgeColor', 'none','FaceAlpha',0.7)
    plot(xvals(iarea), prob_vals(icoeff), 'wo', 'MarkerFaceColor', areacols(iarea,:))
    plot([xvals(iarea), xvals(iarea)], [prob_vals_lower(icoeff),...
        prob_vals_upper(icoeff)], 'Color',areacols(iarea,:))
end
ax= gca; ax.XTick = xvals; ax.XTickLabels=areas;
ylim([0 0.5])
defaultAxesProperties(gca,false)


% compare p(tuned) between areas for locomotion sessions 
runP=ones(8);

for iarea1 = 1:8
    for iarea2 = 1:8
        if iarea1~=iarea2
    H=zeros(1,16);
    H(iarea1*2)=-1;
    H(iarea2*2)=1;
    runP(iarea1,iarea2) = coefTest(glme,H);
        end
    end
end

subplot(2,2,4)
tickvals = [0 0.0001 0.001 0.01 0.05 28];
pvals_disc = discretize(runP*28,tickvals,'IncludedEdge','right');
colmap = colormap(slanCM('gray',6));
colmap=flipud(colmap);
colmap = colmap(1:4,:);
colmap = [colmap; 0 0 0];
imagesc(pvals_disc), colormap(colmap); 
c = colorbar;
c.Ticks =linspace(1,5,6); c.TickLabels = {"","0.0001","0.001", "0.01", "0.05", ">0.05"};

ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
defaultAxesProperties(gca,false)
hold on
for iarea = 1:numel(areas)
for i=1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
end



%% LME to test for difference in tuning strength

allAreaUnits = goodUnits(ismember([goodUnits.ecephys_structure_acronym], areas));

% tunedFlag, r2, stat, run, sessionID, area

stat_r2 = [allAreaUnits.r2_stat];
stat_r2p = [allAreaUnits.r2pval_stat];
stat_area = repelem([allAreaUnits.ecephys_structure_acronym],1,4);
stat_session = repelem([allAreaUnits.sessionID],1,4);
stat_dir = repmat(1:4,1,numel(allAreaUnits));
stat_state = zeros(size(stat_r2));
stat_sigRF = repelem([allAreaUnits.p_value_rf]<=0.05,1,4);

run_r2 = [allAreaUnits.r2_run];
run_r2p = [allAreaUnits.r2pval_run];
run_area = repelem([allAreaUnits.ecephys_structure_acronym],1,4);
run_session = repelem([allAreaUnits.sessionID],1,4);
run_dir = repmat(1:4,1,numel(allAreaUnits));
run_state = ones(size(run_r2));
run_sigRF = repelem([allAreaUnits.p_value_rf]<=0.05,1,4);

for iunit = 1:numel(allAreaUnits)
    allAreaUnits(iunit).fanoFactor_stat =  allAreaUnits(iunit).fanoFactor_stat(:)';
    allAreaUnits(iunit).fanoFactor_run = allAreaUnits(iunit).fanoFactor_run(:)';
    allAreaUnits(iunit).dynamicRange_stat = allAreaUnits(iunit).dynamicRange_stat(:)';
    allAreaUnits(iunit).dynamicRange_run = allAreaUnits(iunit).dynamicRange_run(:)';
end

stat_ff = [allAreaUnits.fanoFactor_stat];
stat_dr = [allAreaUnits.dynamicRange_stat];
run_ff = [allAreaUnits.fanoFactor_run];
run_dr = [allAreaUnits.dynamicRange_run];

all_r2 = [stat_r2, run_r2]; all_r2(all_r2<0)=0;
all_r2p = [stat_r2p, run_r2p];
all_area = [stat_area, run_area];
all_session = [stat_session, run_session];
all_dir = [stat_dir, run_dir];
all_state = [stat_state, run_state];
tunedFlag = all_r2>=r2_thresh & all_r2p<=r2p_thresh;
all_sigRF = [stat_sigRF, run_sigRF];
all_ff = [stat_ff, run_ff];
all_dr = [stat_dr, run_dr];

idx = find(isnan(all_r2));
all_r2(idx) = [];
all_r2p(idx) = [];
all_area(idx) = [];
all_session(idx) = [];
all_dir(idx) = [];
all_state(idx) = [];
tunedFlag(idx) = [];
all_sigRF(idx) = [];
all_ff(idx) = [];
all_dr(idx) = [];


tbl = table(tunedFlag', all_r2(:), all_area', categorical(all_session'),...
    categorical(all_dir'), categorical(all_state'), categorical(all_sigRF(:)),...
    all_ff(:), all_dr(:),...
    'VariableNames', {'tunedFlag', 'tuningStrength', 'area', 'session', 'dir', 'state', 'sigRF', 'fanoFactor', 'dynamicRange'});


%areaOrder = [7, 3, 5, 8, 1, 2, 4, 6];


%f = 'tuningStrength ~ -1 + area:state + (1|session) + (1|dir:sigRF:area)'% + (1|sigRF:area)';
f = 'tuningStrength ~ -1 + area:state + (1|session) + (1|dir) + (1|sigRF)';

lme = fitlme(tbl,f,'DummyVarCoding','full');
areaOrder = [3, 7, 5, 8, 1, 2, 4, 6];
[~, sortidx] = sort(areaOrder);
H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1) = -1; H_temp(iarea1+8) = 1;
    pvals(iarea1) = coefTest(lme,H_temp);
end

pvals = pvals(sortidx)



%% plot tuning strength cdfs, stat vs run area subplots

figure
for iarea = 1:8
    subplot(2,4,iarea), 
    hold on
    statvals = ta(iarea).r2_stat; statvals(statvals<0)=-.1;
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;

    [f, x] = ecdf(statvals);
    plot(x,1-f,'Color', areacols(iarea,:).*0.7);
    [f, x] = ecdf(runvals);
    plot(x,1-f,'Color', areacols(iarea,:));
    xlim([0 1])
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
end

%% plot tuning strength cdfs, all areas, by state

figure
    subplot(121), 
    hold on
    for iarea = 1:8
    statvals = ta(iarea).r2_stat; statvals(statvals<0)=-.1;
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;
    [f, x] = ecdf(statvals);
    plot(x,1-f,'Color', areacols(iarea,:).*0.7);
    xlim([0 1])
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
    end

        subplot(122), 
    hold on
    for iarea = 1:8
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;
    [f, x] = ecdf(runvals);
    plot(x,1-f,'Color', areacols(iarea,:));
    xlim([0 1])
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
    end


%% imagesc plot of p(tuned)
figure
imagesc([[ta.pTuned_stat]; [ta.pTuned_run]]), colormap(cmocean('tempo')), caxis([0 0.5]), colorbar
figure
imagesc([[ta.pTuned_run]-[ta.pTuned_stat]]), colormap(crameri('vik')), caxis([-0.26 0.26]), colorbar
