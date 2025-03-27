%% load sessions

dataDir = 'D:\AllenMatFiles\FC';
dataDir = 'D:\AllenMatFiles\BO';
%dataDir = 'Z:\ibn-vision\DATA\ExternalDatasets\ABI_VisualCoding_Neuropixels\BrainObservatory'
dataFiles = dir(fullfile(dataDir, 'session*.mat'));
% dataFiles_done = dir(fullfile(dataDir, 'pre*'));



% for ifile = 1:numel(dataFiles)
%     sesid2(ifile) = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'))
% end
% 
% for ifile = 1:numel(dataFiles_done)
%     sesid(ifile) = str2double(regexp(dataFiles_done(ifile).name, '\d*', 'match'))
% end
% 
% todo = find(~ismember(sesid2,sesid))

todo = 1:numel(dataFiles);

%%
savePrefix = 'pre_session_';

for ifile =1:numel(dataFiles)
    tic
    ifile
    sessionID = str2double(regexp(dataFiles(todo(ifile)).name, '\d*', 'match'));
    inputFileName = fullfile(dataDir, dataFiles(todo(ifile)).name);
    outputFileName = fullfile(dataDir, [savePrefix, num2str(sessionID), '.mat']);

    
    preprocessAllenSession(sessionID, inputFileName, outputFileName)
    toc
end