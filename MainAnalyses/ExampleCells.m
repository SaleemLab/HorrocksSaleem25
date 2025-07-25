%% example neuron plot  gauss char tuning curves and PSTHs/raster plots

%% get excitation/suppression

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

%% get example tuned cells

% --- 1. Initialization ---
n = 50;
% This script assumes a struct array named 'goodUnits' exists in your workspace.
% Each element of 'goodUnits' should contain the following fields:
%   - gaussChar_stat: 1x4 vector of cluster IDs for the 'stat' condition
%   - r2_stat:        1x4 vector of R-squared values for the 'stat' condition
%   - gaussChar_run:  1x4 vector of cluster IDs for the 'run' condition
%   - r2_run:         1x4 vector of R-squared values for the 'run' condition

% --- 2. Data Aggregation into a Table ---
% Pre-allocate arrays for table construction for better performance
numCurves = numel(goodUnits) * 8;
clusterIDs = zeros(numCurves, 1);
r2Values = zeros(numCurves, 1);
unitIndices = zeros(numCurves, 1);
dirIndices = zeros(numCurves, 1);
conditions = cell(numCurves, 1);

currentIndex = 1;
for iunit = 1:numel(goodUnits)
    % Process '_stat' data (4 directions)
    range = currentIndex:(currentIndex + 3);
    clusterIDs(range) = goodUnits(iunit).gaussChar_stat;
    r2Values(range) = goodUnits(iunit).r2_stat; 
    unitIndices(range) = iunit;
    dirIndices(range) = 1:4;
    conditions(range) = {'stat'};
    currentIndex = currentIndex + 4;

    % Process '_run' data (4 directions)
    range = currentIndex:(currentIndex + 3);
    clusterIDs(range) = goodUnits(iunit).gaussChar_run;
    r2Values(range) = goodUnits(iunit).r2_run;
    unitIndices(range) = iunit;
    dirIndices(range) = 1:4;
    conditions(range) = {'run'};
    currentIndex = currentIndex + 4;
end

% Create the master table
allTuningCurves = table(clusterIDs, r2Values, unitIndices, dirIndices, conditions, ...
    'VariableNames', {'Cluster', 'R2', 'Unit', 'Direction', 'Condition'});

% ADDED: Remove any rows where the R2 value is NaN before proceeding
allTuningCurves = allTuningCurves(~isnan(allTuningCurves.R2), :);

% Convert the 'Condition' column to a categorical array for efficient filtering
allTuningCurves.Condition = categorical(allTuningCurves.Condition);

% --- 3. Find and Store Top 'n' Curves Using Table Operations ---
topResults = struct();
conditionsToAnalyze = {'stat', 'run'};
numClusters = 4;

for iCond = 1:numel(conditionsToAnalyze)
    currentCond = conditionsToAnalyze{iCond};
    
    % Pre-allocate a cell array to store the results tables for this condition
    topResults.(currentCond) = cell(numClusters, 1);
    
    % Filter the table for the current condition
    condTable = allTuningCurves(allTuningCurves.Condition == currentCond, :);
    
    for clusterNum = 1:numClusters
        % Filter for the current cluster
        clusterTable = condTable(condTable.Cluster == clusterNum, :);
        
        % Sort the cluster's table by R2 value in descending order
        sortedClusterTable = sortrows(clusterTable, 'R2', 'descend');
        
        % Determine how many top results to retrieve
        numTop = min(n, height(sortedClusterTable));
        
        % Select the top 'n' rows and store the result
        topResults.(currentCond){clusterNum} = sortedClusterTable(1:numTop, :);
    end
end

% --- 4. Display the Separated Results from the Tables ---
for iCond = 1:numel(conditionsToAnalyze)
    currentCond = conditionsToAnalyze{iCond};
    fprintf('==============================================\n');
    fprintf('=== ANALYSIS FOR CONDITION: %s ===\n', upper(currentCond));
    fprintf('==============================================\n\n');
    
    for clusterNum = 1:numClusters
        fprintf('--- Top %d Tuning Curves for Cluster %d ---\n', n, clusterNum);
        
        resultsTable = topResults.(currentCond){clusterNum};
        
        if isempty(resultsTable)
            fprintf('No tuning curves found for this cluster.\n');
        else
            % Display the resulting table directly. MATLAB's table display is clean.
            disp(resultsTable);
        end
        fprintf('\n');
    end
end


%% plot examples
FRscale = 5;
state='run';
gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));
speedcols = inferno(7);
icond = 4;

for ires = 3%:n;
    % try
    iunit = topResults.(state){icond}.Unit(ires);
    idir = topResults.(state){icond}.Direction(ires);
    
    goodUnits(iunit).ID

    switch state
        case 'stat'
            scName = 'spikeCounts_stat';
            gpName = 'gaussParams_stat';
            ssiName = 'SSI_stat';
            psthName = 'statdirpsth';
            prefSpeedName = 'prefSpeed_stat';
            % ExcVec = goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==1;
            % ExcVec = ExcVec(idir,:);
            % SuppVec = goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==-1;
            % SuppVec = SuppVec(idir,:);
            blFR_mean = mean(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_stat(idir,:))*FRscale);
            blFR_sem = sem(cat(1,goodUnits(iunit).baselineSpikeCounts_stat{idir,:})*FRscale);
        case 'run'
            scName = 'spikeCounts_run';
            gpName = 'gaussParams_run';
            ssiName = 'SSI_run';
            psthName = 'rundirpsth';
            prefSpeedName = 'prefSpeed_run';
            % ExcVec = goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==1;
            % ExcVec = ExcVec(idir,:);
            % SuppVec = goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==-1;
            % SuppVec = SuppVec(idir,:);
            blFR_mean = mean(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_run(idir,:))*FRscale);
            blFR_sem = sem(cat(1,goodUnits(iunit).baselineSpikeCounts_run{idir,:})*FRscale);
    end


    figure
    subplot(2,2,3:4), % plot PSTHs
    sgtitle(sprintf('id: %.0f, area: %s idir: %.0f, istate: %s', goodUnits(iunit).ID, goodUnits(iunit).ecephys_structure_acronym{:}, idir, state));
    hold on
    for ispeed = 1:7
    plot(goodUnits(iunit).(psthName)(idir).PSTH(ispeed).psth, 'Color', speedcols(ispeed,:))
    end
    ax = gca; ax.XTick = 0:50:200;
    defaultAxesProperties(gca, true)
    ylabel('Firing Rate (Hz)')

    subplot(2,2,1), hold on % plot tuning curve and baseline FR
    for ispeed = 1:7

        plot(repelem(ispeed,1,numel(goodUnits(iunit).(scName){idir,ispeed})),...
            goodUnits(iunit).(scName){idir,ispeed}, 'k.', 'MarkerSize',10)
    end
    % errorbar(1:7, cellfun(@mean, goodUnits(iunit).(scName)(idir,:)),...
    %     cellfun(@sem, goodUnits(iunit).(scName)(idir,:)), 'Marker', 'o', 'Color','k')
    
    % plot(1:7, cellfun(@mean, goodUnits(iunit).(scName)(idir,:)),'k') %mean FR
     plot(1:0.01:7, feval(gaussFun, goodUnits(iunit).(gpName)(idir,:),  1:0.01:7), 'k-') % gaussian fit

    shadedErrorBar([1 7], [blFR_mean blFR_mean], [blFR_sem blFR_sem], 'lineProps', 'k')
    
    
    ax = gca; ax.XTick = 1:7; xlim([0.5 7.5]); 
    ylabel('Firing Rate (Hz)')
    defaultAxesProperties(gca, false)
    plot(goodUnits(iunit).(prefSpeedName)(idir), ax.YLim(1), 'kv')


    subplot(2,2,2), hold on
    plot(goodUnits(iunit).(ssiName)(idir,:), 'k.', 'LineStyle', '-')
    ax = gca; ax.XTick = 1:7; 
    xlim([0.5 7.5])
    ylabel('SSI (bits)')
    defaultAxesProperties(gca, false)
    [val, maxssi] = max(goodUnits(iunit).(ssiName)(idir,:));
    plot(maxssi, ax.YLim(1), 'kv')

    % ax.YLim = [min(goodUnits(iunit).(ssiName)(idir,:)), max(goodUnits(iunit).(ssiName)(idir,:))];




    % for ispeed = 1:7
    % 
    %     plot(repelem(ispeed,1,numel(goodUnits(iunit).(scName){idir,ispeed})),...
    %         goodUnits(iunit). (scName){idir,ispeed}, 'k.', 'MarkerSize',10)
    % end
    % 
    

    % subplot(122)
    % plot(1:7, goodUnits(iunit).(ssiName)(idir,:), 'k-')
    % % ylim([0.5 1.9]); 
    % % ax= gca; ax.XTick = 1:7;ax.YTick = 0.5:0.2:1.9; 
    % defaultAxesProperties(gca, true)

   
% --- The key part: Wait for user input ---
disp('Figure ready. Press ''s'' to save, or any other key to continue.');
k = 0;
while k ~= 1
    figure(gcf); % <-- ADD THIS LINE to ensure the figure has focus
    k = waitforbuttonpress;
end

% --- A key was pressed. Now, check if it was 's' ---
key_pressed = get(gcf, 'CurrentCharacter');
if strcmpi(key_pressed, 's')
    % Construct the filename
    filename = ['00_', num2str(icond), '_', state, '_', num2str(goodUnits(iunit).ID), '_', num2str(idir)];

    % Save the figure
    savefig(gcf, filename); % Automatically adds the .fig extension

    disp(['Figure saved as: ' filename '.png']);
else
    disp(['Continuing to next figure... (Key pressed: ''' key_pressed ''')']);
end

close 
    % % catch
    % end

end



%%

saveflag=true;
if saveflag
    % Get handles to all open figures
all_figs = findall(0, 'type', 'figure');

if isempty(all_figs)
    disp('No open figures found.');
    return;
end

% Loop through each figure
for i = 1:length(all_figs)
    fig = all_figs(length(all_figs) - i + 1);
    
    % --- Smart Filename Generation Logic ---
    
    % Priority 1: Use the figure's original FileName property if it exists
    if ~isempty(fig.FileName)
        % Extract the base name (e.g., 'lp1') from the full path
        [~, filename_base, ~] = fileparts(fig.FileName);
        
    % Priority 2: If no FileName, use the figure's Name property
    elseif ~isempty(fig.Name)
        % Sanitize the name for file system compatibility
        clean_name = strrep(fig.Name, ' ', '_');
        filename_base = regexprep(clean_name, '[^\w-]', '');
        
    % Priority 3: If all else fails, use the figure number
    else
        filename_base = ['figure_', num2str(fig.Number)];
    end
    
    % Add the .svg extension to create the final filename
    filename = [filename_base, '.svg'];
    
    % Print the figure to an SVG file
    print(fig, filename, '-dsvg', '-painters');
    
    fprintf('Saved Figure %d as %s\n', fig.Number, filename);
end

disp('--- All figures saved! ✅ ---');

end
