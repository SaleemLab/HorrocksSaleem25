# required imports
import os
import scipy.io as sio
import ast

import numpy as np
import xarray as xr
import pandas as pd

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession,
    removed_unused_stimulus_presentation_columns
)

# tell pandas to show all columns when we display a DataFrame
pd.set_option("display.max_columns", None)

# set path to maifest and data
manifest_path = os.path.join("D:\AllenSDK", "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

# load a specific session - best is to download before running
session_id = 778240327 # func connec
session = cache.get_session_data(session_id)

# get stimulus table and process required parameters
stimtable = session.get_stimulus_table();

stimtable["pos"].replace(to_replace="null", value="None", inplace=True)
stimtable["phase"].replace(to_replace="null", value="None", inplace=True)
stimtable["size"].replace(to_replace="null", value="None", inplace=True)
stimtable["spatial_frequency"].replace(to_replace="null", value="None", inplace=True)

pos = stimtable["pos"].tolist()
phase = stimtable["phase"].tolist()
size = stimtable["size"].tolist()
spatial_frequency = stimtable["spatial_frequency"].tolist()
length = len(pos)

for i in range(length):
    pos[i] = ast.literal_eval(pos[i])
    phase[i] = ast.literal_eval(phase[i])
    size[i] = ast.literal_eval(size[i])
    spatial_frequency[i] = ast.literal_eval(spatial_frequency[i])


stimtable["pos"] = pos
stimtable["phase"] = phase
stimtable["size"] = size
stimtable["spatial_frequency"] = spatial_frequency
stimtable.fillna(value=pd.np.nan, inplace=True)
stimtable.replace(to_replace="null", value=pd.np.nan, inplace=True)

# get generic unit info, e.g. unit location
unitInfo = session.units
unitIndex = np.array(unitInfo.index)

# dict of spike times for each unit, make matlab compatible
spikeTimes = session.spike_times
newSpikeTimes = {}
for key in spikeTimes.keys() :
    newSpikeTimes["c_" + str(key)] = [spikeTimes[key]]

# running speed info + pupil info
run_times = session.running_speed["start_time"] + \
    (session.running_speed["end_time"] - session.running_speed["start_time"]) / 2
run_speed = session.running_speed["velocity"]
run_times = np.array(run_times)
run_speed = np.array(run_speed)
runInfo = np.column_stack((run_times, run_speed))

pupilTable = session.get_pupil_data()
pupilTime = np.array(pupilTable.index)
pupilWidth = np.array(pupilTable["pupil_width"])
pupilHeight = np.array(pupilTable["pupil_height"])
pupil_x = np.array(pupilTable["pupil_center_x"])
pupil_y = np.array(pupilTable["pupil_center_y"])

# running speed info
pupilInfo = {'pupilTime' : pupilTime,
              'pupilWidth': pupilWidth,
              'pupilHeight' : pupilHeight,
              'pupil_x': pupil_x,
              'pupil_y' : pupil_y
             }

outputVars = {'stim_table': {name: col.values for name, col in stimtable.items()},
              'unitInformation' : {name: col.values for name, col in unitInfo.items()},
              'unitIndex': unitIndex,
              'spikeTimes' : newSpikeTimes,
              'runInfo' : runInfo,
              'pupilInfo' : pupilInfo
            }


sio.savemat('session_778240327.mat', outputVars,  long_field_names=True)
