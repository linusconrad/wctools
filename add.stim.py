# This script creates a stimulus waveform of an .abf file with the pyabf API

import pyabf
import pandas

abf = pyabf.ABF("2020_12_15 sgn P7 2DIV 1B cell5_0004.abf")

print(abf.headerText)

# Make a container for the data
stimulus = []

# iterate over the sweeps to generate all the stimulus waveform
for i in range(abf.sweepCount):
    abf.setSweep(i)
    
