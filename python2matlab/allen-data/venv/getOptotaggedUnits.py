import os

import numpy as np
import pandas as pd
import xarray as xr
import scipy.io as sio


import matplotlib.pyplot as plt

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache


########## PARSE ARGS AND ACCESS SESSION ##########
##ap = argparse.ArgumentParser()
##ap.add_argument("-i", "--input", required=True,
##type=int, help="session ID number")
##ap.add_argument("-o", "--output", required=True,
##help="path to output file")

##args = ap.parse_args()

##sessionNumber = args.input

# set path to maifest and data
manifest_path = os.path.join("D:\AllenSDK", "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

sessions = cache.get_session_table()
pvalb_sessions = sessions[(sessions.full_genotype.str.match('Vip')) &
                         (sessions.session_type == 'functional_connectivity')]

sessionValues = pvalb_sessions.index.values
taggedUnits = np.empty(shape=1,dtype=np.int64)

for sessionid in sessionValues:
    session = cache.get_session_data(sessionid)
    print(sessionid)
    trials = session.optogenetic_stimulation_epochs[(session.optogenetic_stimulation_epochs.duration > 0.004) & \
                                                (session.optogenetic_stimulation_epochs.duration < 0.02)]

    units = session.units  # [session.units.ecephys_structure_acronym.str.match('VIS')]
    time_resolution = 0.0005  # 0.5 ms bins
    bin_edges = np.arange(-0.01, 0.025, time_resolution)


    def optotagging_spike_counts(bin_edges, trials, units):
        time_resolution = np.mean(np.diff(bin_edges))

        spike_matrix = np.zeros((len(trials), len(bin_edges), len(units)))

        for unit_idx, unit_id in enumerate(units.index.values):

            spike_times = session.spike_times[unit_id]

            for trial_idx, trial_start in enumerate(trials.start_time.values):
                in_range = (spike_times > (trial_start + bin_edges[0])) * \
                           (spike_times < (trial_start + bin_edges[-1]))

                binned_times = ((spike_times[in_range] - (trial_start + bin_edges[0])) / time_resolution).astype('int')
                spike_matrix[trial_idx, binned_times, unit_idx] = 1

        return xr.DataArray(
            name='spike_counts',
            data=spike_matrix,
            coords={
                'trial_id': trials.index.values,
                'time_relative_to_stimulus_onset': bin_edges,
                'unit_id': units.index.values
            },
            dims=['trial_id', 'time_relative_to_stimulus_onset', 'unit_id']
        )

    da = optotagging_spike_counts(bin_edges, trials, units)

    baseline = da.sel(time_relative_to_stimulus_onset=slice(-0.01, -0.002))

    baseline_rate = baseline.sum(dim='time_relative_to_stimulus_onset').mean(dim='trial_id') / 0.008

    evoked = da.sel(time_relative_to_stimulus_onset=slice(0.001, 0.009))

    evoked_rate = evoked.sum(dim='time_relative_to_stimulus_onset').mean(dim='trial_id') / 0.008

    cre_pos_units = da.unit_id[(evoked_rate / (baseline_rate + 1)) > 2].values
    cre_pos_units = np.array(cre_pos_units)
    print(cre_pos_units.dtype)
    print(taggedUnits.dtype)

    taggedUnits = np.concatenate([taggedUnits, cre_pos_units])


print(taggedUnits)



########## GENERATE OUTPUT VARS DICT AND SAVE ##########

outputVars = {'VIPtags' : taggedUnits}

sio.savemat('VIPtags.mat', outputVars,  long_field_names=True)