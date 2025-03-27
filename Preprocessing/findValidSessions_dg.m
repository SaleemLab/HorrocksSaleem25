%% Find Valid Sessions

%dataDir = 'X:\ibn-vision\USERS\Edd\AllenDataAnalysis_2022\data\FunctionalConnectivity22';
dataDir = 'D:\AllenMatFiles\BO';

dataFiles = dir(fullfile(dataDir, 'pre_session*'));

%% paramters to choose sessions to analyse
tic
statMeanThresh = 1;
statUpperThresh = 3;
runMeanThresh = 3;
runLowerThresh = 0.5;

for ifile = 1:numel(dataFiles)

    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name))

    dmtrials = trials(strcmp([trials.stimulus_name], 'drifting_gratings'));
    dmtrials_valid = dmtrials(~[dmtrials.containsInvalidTime]);

    % if numel(unique([dmtrials.Speed]==7)) % 2 sessions don't have the 7 speeds

    run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials_valid.runTrace}));
    stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials_valid.runTrace}));

    statTrials = dmtrials_valid(stat_idx);
    runTrials = dmtrials_valid(run_idx);

    % statTrials = dmtrials_valid([dmtrials_valid.meanRunSpeed]<=statMeanThresh &...
    %     [dmtrials_valid.maxRunSpeed]<=statUpperThresh);
    %
    % runTrials = dmtrials_valid([dmtrials_valid.meanRunSpeed]>=runMeanThresh &...
    %     [dmtrials_valid.minRunSpeed]>=runLowerThresh);

    uniqueDirs = unique([dmtrials.orientation]);
    uniqueDirs(isnan(uniqueDirs))=[];
    nDirs=numel(uniqueDirs);
    uniqueTF = unique([dmtrials.temporal_frequency]);
    uniqueTF(isnan(uniqueTF))=[];
    nTF = numel(uniqueTF);

    nStat = nan([nDirs,nTF]); nRun = nan([nDirs,nTF]);

    for idir = 1:numel(uniqueDirs)
        for ispeed = 1:numel(uniqueTF)

            nStat(idir,ispeed) = sum([statTrials.orientation]==uniqueDirs(idir) &...
                [statTrials.temporal_frequency]==uniqueTF(ispeed));

            nRun(idir,ispeed) = sum([runTrials.orientation]==uniqueDirs(idir) &...
                [runTrials.temporal_frequency]==uniqueTF(ispeed));
        end
    end


    session(ifile).nStat = nStat;
    session(ifile).nRun = nRun;
    %         session(ifile).statDirs = find(all(nStat>=reqTrials,2));
    %         session(ifile).runDirs = find(all(nRun>=reqTrials,2));



end

toc
%%

reqTrials = 8;

for ifile = 1:numel(dataFiles)

    session(ifile).sessionID = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'));
    session(ifile).statTF = find(all(session(ifile).nStat>=reqTrials,2));
    session(ifile).runTF = find(all(session(ifile).nRun>=reqTrials,2));
    session(ifile).statDirs = find(all(session(ifile).nStat>=reqTrials,1));
    session(ifile).runDirs = find(all(session(ifile).nRun>=reqTrials,1));
end


nStatDirs = cellfun(@numel, {session.statDirs});
nRunDirs = cellfun(@numel, {session.runDirs});

idx = find(nStatDirs~=0);
numStatSessions = numel(idx);
meanStatDirs = mean(nStatDirs(idx));

idx = find(nRunDirs~=0);
numRunSessions = numel(idx);
meanRunDirs = mean(nRunDirs(idx));

DirTuningSessions = [numStatSessions, numRunSessions; meanStatDirs, meanRunDirs]


nStatTF = cellfun(@numel, {session.statTF});
nRunTF = cellfun(@numel, {session.runTF});

idx = find(nStatTF~=0);
numStatSessions = numel(idx);
meanStatTF = mean(nStatTF(idx));

idx = find(nRunTF~=0);
numRunSessions = numel(idx);
meanRunTF = mean(nRunTF(idx));

TFTuningSessions = [numStatSessions, numRunSessions; meanStatTF, meanRunTF]

