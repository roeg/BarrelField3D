import reader
from celltable import *

def load_synapse_table(fname):
    if fname.endswith('.csv') or fname.endswith('.CSV'):
        return reader.load_csv_synapse_table(fname)
    elif fname.endswith('.am') or fname.endswith('.AM'):
        return reader.load_amira_synapse_table(fname)
    else:
        errstr = 'Unknown file type: %s' % fname
        raise RuntimeError(errstr)

def load_synapses_per_celltype_per_column_table(fname):
    if fname.endswith('.csv') or fname.endswith('.CSV'):
        return reader.load_csv_avg_synapses_table(fname)
    else:
        errstr = 'Unknown file type: %s' % fname
        raise RuntimeError(errstr)