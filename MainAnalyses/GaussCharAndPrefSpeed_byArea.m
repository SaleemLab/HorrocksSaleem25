%% Brain area vis speed characterisation

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




%%
areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
r2_thresh = 0.1;
r2p_thresh=0.05;
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);


%original: 1 = lowpass, 2 = highpass, 3 = bandpass, 4 = inverted
%plot order: 1 = lowpass, 2 = bandpass, 3 = inverted, 4 = highpass
plotOrder = [1 3 4 2];
%% get distributions of preferred speed and shape by visual area

for iarea = 1:numel(areas)

    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));

    % tuned stat units
    idx = find([areaUnits.r2_stat]>r2_thresh & [areaUnits.r2pval_stat]<r2p_thresh);

    allPrefSpeeds = [areaUnits.prefSpeed_stat];
    allGaussChar = [areaUnits.gaussChar_stat];
    allGaussR2 = [areaUnits.gaussR2_stat];
    allr2 = [areaUnits.r2_stat];
    allr2p = [areaUnits.r2pval_stat];
    allSession = repelem([areaUnits.sessionID],1,4);
    allDir = repmat(1:4,1,numel(areaUnits));

    allnSig = [areaUnits.nSig_stat];
    allnExc(iunit).nExc_stat = [areaUnits.nExc_stat];
    allnSupp(iunit).nSupp_stat = [areaUnits.nSupp_stat];

    statta(iarea).prefSpeeds = allPrefSpeeds(idx);
    statta(iarea).gaussChar = allGaussChar(idx);
    statta(iarea).gaussR2 = allGaussR2(idx);
    statta(iarea).r2 = allr2(idx);
    statta(iarea).r2p = allr2p(idx);
    statta(iarea).nSig = allnSig(idx);
    statta(iarea).nExc = allnExc(idx);
    statta(iarea).nSupp = allnSupp(idx);
    statta(iarea).session = allSession(idx);
    statta(iarea).dir = allDir(idx);

    % tuned stat units
    idx = find([areaUnits.r2_run]>r2_thresh & [areaUnits.r2pval_run]<r2p_thresh);

    allPrefSpeeds = [areaUnits.prefSpeed_run];
    allGaussChar = [areaUnits.gaussChar_run];
    allGaussR2 = [areaUnits.gaussR2_run];
    allr2 = [areaUnits.r2_run];
    allr2p = [areaUnits.r2pval_run];

    allnSig = [areaUnits.nSig_run];
    allnExc(iunit).nExc_run = [areaUnits.nExc_run];
    allnSupp(iunit).nSupp_run = [areaUnits.nSupp_run];

    runta(iarea).prefSpeeds = allPrefSpeeds(idx);
    runta(iarea).gaussChar = allGaussChar(idx);
    runta(iarea).gaussR2 = allGaussR2(idx);
    runta(iarea).r2 = allr2(idx);
    runta(iarea).r2p = allr2p(idx);
    runta(iarea).nSig = allnSig(idx);
    runta(iarea).nExc = allnExc(idx);
    runta(iarea).nSupp = allnSupp(idx);
    runta(iarea).session = allSession(idx);
    runta(iarea).dir = allDir(idx);



end



%% mixed effects model for statistical analysis
allPrefSpeeds = cat(2,statta.prefSpeeds,runta.prefSpeeds)';
allGaussChar = categorical(cat(2,statta.gaussChar,runta.gaussChar)');
allAreaVec = categorical(cat(2,repelem(1:8,1,cellfun(@numel, {statta.prefSpeeds})),repelem(1:8,1,cellfun(@numel, {runta.prefSpeeds}))))';
allDir = categorical(cat(2,statta.dir,runta.dir)');
allSession = categorical(cat(2,statta.session,runta.session)');


tbl = table(allPrefSpeeds, allGaussChar, allAreaVec,allDir, allSession,...
    'VariableNames', {'prefSpeeds', 'gaussChar', 'areas', 'dir', 'session'});

f = 'prefSpeeds ~ -1 + areas + (1|session) + (1|dir)';
lme = fitlme(tbl,f,'dummyVarCoding','full')

for iarea1 = 1:8
    for iarea2 = 1:8
        H=zeros(1,8);
        H(iarea1)=1;H(iarea2)=-1;
        pvals(iarea1,iarea2) = coefTest(lme,H);
    end
end

pvals_bc = pvals*28; % bonferroni correction

tickvals = [0 0.001 0.01 0.05 28];
pvals_disc = discretize(pvals_bc,tickvals,'IncludedEdge','right');
colmap = colormap(slanCM('gray',4));
colmap=flipud(colmap);
colmap = colmap(1:3,:);
colmap = [colmap; 1 1 1];
imagesc(pvals_disc), colormap(colmap); 
c = colorbar;
c.Ticks =linspace(1,4,5); c.TickLabels = {"","0.001", "0.01", "0.05", ">0.05"};
ax=gca;
ax.XTickLabel = areas; ax.YTickLabel = areas;
grid on


%% gauss char glme for area

for igc = 1:4
    thisGCvec = double(allGaussChar);
    
    thisGCvec(thisGCvec==igc)=-1;
    thisGCvec(thisGCvec>0)=0;
    thisGCvec(thisGCvec==-1)=1;
    thisGCvec = thisGCvec;

    tbl = table(allPrefSpeeds, thisGCvec, allAreaVec,allDir,allSession,...
    'VariableNames', {'prefSpeeds', 'allGaussChar', 'areas','dir','session'});

f = 'allGaussChar ~ -1 + areas + (1|session) + (1|dir)';
glme = fitglme(tbl,f,'dummyVarCoding','full','Distribution','Binomial');

for iarea1 = 1:8
    for iarea2 = 1:8
        H=zeros(1,8);
        H(iarea1)=1;H(iarea2)=-1;
        pvals(iarea1,iarea2) = coefTest(glme,H);
    end
end

gc(igc).pvals_bc = pvals*28; % bonferroni correction

end

%% ks-distance for gauss char dists and pref speeds

% gauss char
for iarea = 1:8
    ta(iarea).gaussChar = [statta(iarea).gaussChar, runta(iarea).gaussChar];
end

for iarea1=1:8
    for iarea2=1:8
        [~,kspvalchar(iarea1,iarea2),kval(iarea1,iarea2)] = kstest2(ta(iarea1).gaussChar, ta(iarea2).gaussChar);
    end
end


figure
imagesc(kval)
colorbar
ax=gca; ax.XTickLabels = areas; ax.YTickLabels = areas;
colormap(slanCM('reds'))


for iarea = 1:8
    ta(iarea).prefSpeeds = [statta(iarea).prefSpeeds, runta(iarea).prefSpeeds];
end

for iarea1=1:8
    for iarea2=1:8
        [~,kspvalspd(iarea1,iarea2),kval(iarea1,iarea2)] = kstest2(ta(iarea1).prefSpeeds, ta(iarea2).prefSpeeds);
    end
end


figure
imagesc(kval)
colorbar
ax=gca; ax.XTickLabels = areas; ax.YTickLabels = areas;
colormap(slanCM('reds'))




%% stacked bar charts
gcVals= nan(4,8);

for iarea = 1:8
    for igc = 1:4
        gcVals(igc,iarea)=numel(find([[statta(iarea).gaussChar], [runta(iarea).gaussChar]]==igc));
    end
end

gcVals_norm = gcVals./sum(gcVals,1);
gcVals_norm = gcVals_norm(plotOrder,:);

figure
bar(gcVals_norm', 'stacked')
defaultAxesProperties(gca, true)

speedcols = inferno(7);

speedVals= nan(7,8);

for iarea = 1:8
    for ispeed = 1:7
        speedVals(ispeed,iarea)=numel(find([[statta(iarea).prefSpeeds], [runta(iarea).prefSpeeds]]==ispeed));
    end
end

speedVals_norm = speedVals./sum(speedVals,1);

figure
b = bar(speedVals_norm', 'stacked');
for ispeed = 1:7
    b(ispeed).FaceColor = speedcols(ispeed,:);
end
defaultAxesProperties(gca, false)




%% bar plot of pref speed histograms

figure
for iarea = 1:8

    subplot(2,4,iarea)
    ta(iarea).prefSpeeds = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    h = histcounts(ta(iarea).prefSpeeds);
    bar(1:7,h./sum(h), 'FaceColor', areacols(iarea,:))
    hold on
    %plot(mean([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]]), 0.45, 'v')
    ylim([0 0.5]); ax = gca; ax.XTick = 1:7; ax.YTick = 0:0.1:0.5;
    %title(num2str(2^(3+mean([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]])),3))
    ta(iarea).rawmean = mean(ta(iarea).prefSpeeds);
    ta(iarea).meanPref = 2^(3+mean([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]]));
    ta(iarea).varPref = var([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]]);
    meanval = mean(ta(iarea).prefSpeeds);
    plot([meanval-ta(iarea).varPref/2, meanval+ta(iarea).varPref/2],[0.45, 0.45],'Color',areacols(iarea,:))
    plot(mean(ta(iarea).prefSpeeds), 0.45,'o','MarkerFaceColor',areacols(iarea,:),'MarkerEdgeColor','w')
    defaultAxesProperties(gca, false)

end

%% colormap of pref speeds
figure,
imagesc([ta.rawmean]);
colormap(slanCM('reds'))
caxis([3.5 5])

%% error bars for mean and variance of preferred speed dists.

% mean preferred speed
figure, hold on
for iarea = 1:8
    vals = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    mean_val = mean(vals);
    sem_val = sem(vals(:)).*1.96;
    plot(iarea,mean_val, 'o', 'MarkerFaceColor', areacols(iarea,:), 'MarkerEdgeColor', 'w')
    plot([iarea, iarea], [mean_val-sem_val, mean_val+sem_val], 'Color', areacols(iarea,:));

end
ylim([3 5])
defaultAxesProperties(gca, true)


% variance of distribution
figure, hold on
for iarea = 1:8
    vals = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    var_val = var(vals);
    [h,p,ci,stats] = vartest(vals,0, 'Alpha', 0.05);
    plot(iarea,var_val, 'o', 'MarkerFaceColor', areacols(iarea,:), 'MarkerEdgeColor', 'w')
    plot([iarea, iarea], [ci(1), ci(2)], 'Color', areacols(iarea,:));

end
ylim([1 4])
defaultAxesProperties(gca, true)




%% statistical test of variance of speed distributions

for iarea = 1:8
    ta(iarea).prefSpeeds = cat(2,statta(iarea).prefSpeeds, runta(iarea).prefSpeeds);
end

group = repelem(1:8, [cellfun(@numel, {ta.prefSpeeds})]);
vals = [ta.prefSpeeds];

[p, stats] = vartestn(vals(:),group(:), 'TestType', 'LeveneQuadratic','display','off');


for iarea1 = 1:8
    for iarea2 = 1:8
        if iarea1==iarea2
            pArray(iarea1,iarea2) = nan;
        else
            group = repelem([iarea1, iarea2], [cellfun(@numel, {ta([iarea1, iarea2]).prefSpeeds})]);
            vals = [ta([iarea1, iarea2]).prefSpeeds];
            pArray(iarea1,iarea2) = vartestn(vals(:),group(:), 'TestType', 'LeveneQuadratic','display', 'off');
        end
    end
end


for iarea1 = 1:8
    for iarea2 = 1:8
        valDiff(iarea1,iarea2) = ta(iarea1).varPref-ta(iarea2).varPref;
    end
end

pArray = pArray.*nchoosek(8,2); % bonf corr
figure, hold on
imagesc(valDiff), axis ij
for iarea1 = 1:8
    for iarea2 = 1:8
        if pArray(iarea1,iarea2)<=0.001
            txt = '***';
        elseif pArray(iarea1,iarea2)<=0.01
            txt = '**';
        elseif pArray(iarea1,iarea2)<=0.05
            txt = '*';
        else
            txt = '';
        end
        text(iarea1,iarea2, txt,'HorizontalAlignment', 'Center');
    end
end

colormap(crameri('vik')), colorbar

