########## IMPORTS ##########
import ast
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import os
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession,
    removed_unused_stimulus_presentation_columns
)
from allensdk.brain_observatory.ecephys.visualization import plot_mean_waveforms, plot_spike_counts, raster_plot
from allensdk.brain_observatory.visualization import plot_running_speed

# tell pandas to show all columns when we display a DataFrame
pd.set_option("display.max_columns", None)

manifest_path = os.path.join("D:\AllenSDK", "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
type=int, help="session ID number")
#ap.add_argument("-o", "--output", required=True,
#help="path to output file")

args = ap.parse_args()

sessionNumber = args.input
session_id = sessionNumber # for example
session = cache.get_session_data(session_id)

#EcephysSession?
#session.probes["has_lfp_data"]
df = session.probes
df = df.loc[df['has_lfp_data'] == True]
probe_id = np.array(df.index)
print(probe_id)

from scipy.ndimage.filters import gaussian_filter


## save CSD figure
for id in range(len(probe_id)):
    csd = session.get_current_source_density(probe_id[id])

    _ = plt.figure(figsize=(10, 10))
    filtered_csd = gaussian_filter(csd.data, sigma=(5, 1))
    maxVal = np.max(filtered_csd)
    minVal = np.abs(np.min(filtered_csd))
    absMax = np.max([minVal, maxVal])
    print(absMax)

    fig, ax = plt.subplots(figsize=(10, 20))
    c = ax.pcolor(csd["time"], csd["vertical_position"], filtered_csd, vmin=-absMax, vmax=absMax, cmap='jet')
    _ = ax.set_xlabel("time relative to stimulus onset (s)")
    _ = ax.set_ylabel("vertical position (um)")
    ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.set_xlim((0, 0.3))
    fig.colorbar(c, ax=ax)
    fig.suptitle('session:' + str(session_id) + ', ' + 'probe id:' + str(probe_id[id]))
    plt.savefig('session_' + str(session_id) + '_' + 'probe id_' + str(probe_id[id]) + '.png')

########## GENERATE OUTPUT VARS DICT AND SAVE ##########
#outputVars = {'csd' : csd,
#            }


#sio.savemat('csdtest.mat', outputVars,  long_field_names=True)

