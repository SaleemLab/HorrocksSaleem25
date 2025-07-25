%% load data


 dataDir = 'D:\AllenMatFiles\FC';
  % dataDir = 'D:\AllenMatFiles\BO';

% dataDir = 'X:\ibn-vision\USERS\Edd\AllenDataAnalysis_2022\data\FunctionalConnectivity22'
 % dataFiles = dir(fullfile(dataDir, 'dgTuning_statrun*')); % ''dmTuning_statRun_8t*'dgTuning_statrun
 dataFiles = dir(fullfile(dataDir, 'dmTuning_statRun_8t_wpsth*')); % ''dmTuning_statRun_8t*'dgTuning_statrun

tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'units', 'trials');
    session(ifile).units = units;
end
toc
%% generate allUnits struct

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};


allUnits = [session.units];

nUnits = numel(allUnits)

goodUnits = allUnits([allUnits.isi_violations]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.waveform_amplitude]>=50 &...
    ismember([allUnits.ecephys_structure_acronym],areas));
nGoodUnits = numel(goodUnits)

areaUnits = goodUnits(ismember([goodUnits.ecephys_structure_acronym],areas));
nAreaUnits = numel(areaUnits)

%% get cell-type

getCellType=false;
if getCellType
tic
for iunit = 1:numel(goodUnits)
    iunit
    [ccg, t] = CCG(goodUnits(iunit).spiketimes, ones(size(goodUnits(iunit).spiketimes)),...
    'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
    goodUnits(iunit).acg = ccg;
    fit_params_out = fit_ACG(ccg,false);

    goodUnits(iunit).tau_rise = fit_params_out.acg_tau_rise;
end
toc
narrow_idx = find([goodUnits.waveform_duration]<=0.45);
wide_idx = find([goodUnits.waveform_duration]>0.45 & [goodUnits.tau_rise]>6);
pyr_idx = find(~ismember(1:numel(goodUnits), [narrow_idx,wide_idx]));

[goodUnits.cellType] = deal(nan);
[goodUnits(pyr_idx).cellType] = deal(1);
[goodUnits(narrow_idx).cellType] = deal(2);
[goodUnits(wide_idx).cellType] = deal(3);
end