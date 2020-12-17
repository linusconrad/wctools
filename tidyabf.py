# This script creates a stimulus waveform of an .abf file with the pyabf API

import pyabf
import pandas as pd

def tidyabf(file):
    abf = pyabf.ABF(file)
# Make a container for the data
    abfcontents = []
# iterate over the sweeps to generate all the stimulus waveform
    for i in range(abf.sweepCount):
        abf.setSweep(i)
        loopdf = pd.DataFrame({'t': pd.Series(abf.sweepX),
                              'value': pd.Series(abf.sweepY),
                              'command' : pd.Series(abf.sweepC)})
        abfcontents.append(loopdf)
    # collect the data to one frame
    abfdata = pd.concat(abfcontents)
    abf.sweepPointCount
    # add sweep information
    sweep = pd.Series(abf.sweepList).repeat(abf.sweepPointCount)

    # need to use "values" here else pandas get problems with indices
    abfdata['sweep'] = sweep.values

    return(abfdata)
