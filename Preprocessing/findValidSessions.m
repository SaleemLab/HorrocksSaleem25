%% Find Valid Sessions

dataDir = 'D:\AllenMatFiles\FC';
dataFiles = dir(fullfile(dataDir, 'pre_session*'));

%% paramters to choose sessions to analyse


tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name))
    
    dmtrials = trials(strcmp([trials.stimulus_name], 'dot_motion'));
    dmtrials_valid = dmtrials(~[dmtrials.containsInvalidTime]);
    
    if numel(unique([dmtrials.Speed]==7)) % 2 sessions don't have the 7 speeds

        run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials_valid.runTrace}));
        stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials_valid.runTrace}));

        stat_idx = [dmtrials_valid.meanRunSpeed]<=1;
        run_idx = [dmtrials_valid.meanRunSpeed]>=3;



        statTrials = dmtrials_valid(stat_idx);
        runTrials = dmtrials_valid(run_idx);
        
        % statTrials = dmtrials_valid([dmtrials_valid.meanRunSpeed]<=statMeanThresh &...
        %     [dmtrials_valid.maxRunSpeed]<=statUpperThresh);
        % 
        % runTrials = dmtrials_valid([dmtrials_valid.meanRunSpeed]>=runMeanThresh &...
        %     [dmtrials_valid.minRunSpeed]>=runLowerThresh);
        
        uniqueDirs = unique([dmtrials.Dir]);
        uniqueSpeeds = unique([dmtrials.Speed]);
        
        nStat = nan([4,7]); nRun = nan([4,7]);
        
        for idir = 1:numel(uniqueDirs)
            for ispeed = 1:numel(uniqueSpeeds)
                
                nStat(idir,ispeed) = sum([statTrials.Dir]==uniqueDirs(idir) &...
                    [statTrials.Speed]==uniqueSpeeds(ispeed));
                
                nRun(idir,ispeed) = sum([runTrials.Dir]==uniqueDirs(idir) &...
                    [runTrials.Speed]==uniqueSpeeds(ispeed));
            end
        end
        
        
        session(ifile).nStat = nStat;
        session(ifile).nRun = nRun;
%         session(ifile).statDirs = find(all(nStat>=reqTrials,2));
%         session(ifile).runDirs = find(all(nRun>=reqTrials,2));
        
        
    else % not valid session
        
%         session(ifile).statDirs = [];
%         session(ifile).runDirs = [];
        session(ifile).nStat = 0;
        session(ifile).nRun = 0;
        
    end
    
    
    
end
toc

%%

reqTrials = 4;

for ifile = 1:numel(dataFiles)

session(ifile).sessionID = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'));
session(ifile).statDirs = find(all(session(ifile).nStat>=reqTrials,2));
session(ifile).runDirs = find(all(session(ifile).nRun>=reqTrials,2));
session(ifile).mixedDirs = find(all(session(ifile).nStat>=reqTrials,2) & all(session(ifile).nRun>=reqTrials,2));

end


nStatDirs = cellfun(@numel, {session.statDirs});
idx = find(nStatDirs~=0);
numStatSessions = numel(idx);
meanStatDirs = mean(nStatDirs(idx));

nRunDirs = cellfun(@numel, {session.runDirs});
idx = find(nRunDirs~=0);
numRunSessions = numel(idx);
meanRunDirs = mean(nRunDirs(idx));

[numStatSessions, numRunSessions; meanStatDirs, meanRunDirs]
