import numpy as np
from pynwb import NWBHDF5IO

io = NWBHDF5IO('session_821695405.nwb', 'r')
nwbfile_in = io.read()

