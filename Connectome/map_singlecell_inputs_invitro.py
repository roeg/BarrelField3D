#===============================================================================
# SingleCellInputMapper
# Tool for estimating connectivity (inputs) of individual neuron morphologies
# registered into standard barrel cortex model.
# Based on methods and data presented in:
# Egger, Dercksen et al., Frontiers Neuroanatomy 2014
# 
# Inputs:
# - single neuron morphology
# - CellType 'L5tt'
# - e.g. map_singlecell_inputs_invitro(
#   '//nas1/nas1/Data_daniel/Network/L5/L5_InVitro/DendriteHoc/cell_0_20_shiftX_-25.hoc',
#   'L5tt')
# Outputs:
# - .csv file with number of synapses onto Apical and Basal Dendrite for each bouton density
#
# Modified by: Daniel Udvary
# Author: Robert Egger
#         Computational Neuroanatomy
#         Max Planck Institute for Biological Cybernetics
#         Tuebingen, Germany
#         Email: robert.egger@tuebingen.mpg.de
#===============================================================================

# Python standard library imports
import sys
import os.path
import glob
import time
import numpy as np
import singlecell_input_mapper as sim

exTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct')
inhTypes = ('SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
            'L1','L23Trans','L45Sym','L45Peak','L56Trans')

def map_singlecell_inputs_invitro(cellName, cellTypeName):

    if not (cellTypeName in exTypes) and not (cellTypeName in inhTypes):
        errstr = 'Unknown cell type %s!' % cellTypeName
        print 'ERROR! %s' % errstr
        raise TypeError(errstr)

    if cellName[-4:] != '.hoc':
        errstr = 'Can only read in hoc morphologies: %s!' % cellName
        print 'ERROR! %s' % errstr
        raise TypeError(errstr)
    
    startTime = time.time()

    # Main Path
    if sys.platform == 'win32':
        prefix = '//nas1/nas1/Data_daniel/Network/L5/L5_InVitro/'
    else:
        prefix = '/nas1/Data_daniel/Network/L5/L5_InVitro/'

    # Load Postsynaptic hoc Morphology
    print 'Loading cell morphology...'
    parser = sim.CellParser(cellName)
    parser.spatialgraph_to_cell()
    singleCell = parser.get_cell()

    # Position of Postsynaptic Morphology and get path to correct bouton densities
    strtmp = cellName.replace('_invivo','')
    s = strtmp.split('_')
    postsynapse_y = s[-3]
    postsynapse_x = s[-4]
    if postsynapse_y == '0':
        postsynapse_y = s[-1][:-4].strip('.am')
        boutonDensityFolderName = os.path.join(prefix,'BoutonDensity_y0')
    elif postsynapse_x == '0' :
        postsynapse_x = s[-1][:-4].strip('.am')
        boutonDensityFolderName = os.path.join(prefix,'BoutonDensity_x0')
    else:
        errstr = 'Postsynaptic morphology lacks information about sliced dimension ( %s _ %s )!' % (postsynapse_x,postsynapse_y)
        print 'ERROR! %s' % errstr
        raise TypeError(errstr)

    # Required data for calculation synapse density
    connectionsSpreadsheetName = os.path.join(prefix,'data/CellTypeConnectionSpreadsheet.csv')
    ExPSTDensityName = os.path.join(prefix,'data/EXNormalizationPSTs.am')
    connectionsSpreadsheet = sim.read_connections_spreadsheet(connectionsSpreadsheetName)
    ExPSTDensity = sim.read_scalar_field(ExPSTDensityName)
    ExPSTDensity.resize_mesh()
    InhPSTDensity = []
    outputfname = cellName.replace('.hoc','.csv')

    synapseDensityComputation = sim.SynapseDensity(singleCell, cellTypeName, connectionsSpreadsheet,\
                                               exTypes, inhTypes, ExPSTDensity, InhPSTDensity)

    boutonDensityList = os.path.join(boutonDensityFolderName,'list.txt')
    outputstr = 'Presynapse,Postsynapse,Presynapse_x,Presynapse_y,Postsynapse_x,Postsynapse_y,synapsesBasal,synapsesApical\n'

     # Load each bouton density and compute Innervation
    with open(boutonDensityList, 'r') as spreadsheet:
        for line in spreadsheet:

            # Load Bouton Density
            tmpLine = line.strip('\n\r')
            densityName = os.path.join(boutonDensityFolderName,tmpLine)
            boutonDensity = sim.read_scalar_field(densityName)
            boutonDensity.resize_mesh()

            # Compute Innervation
            syn = synapseDensityComputation.compute_synapse_density(boutonDensity, cellTypeName)

            # Extract Innervation, if tmp == None, no innervation
            synApical = 0
            synBasal = 0

            if syn is not None:
                if 'Dendrite' in syn.keys():
                    synBasal = np.sum(syn['Dendrite'].mesh)
                if 'ApicalDendrite' in syn.keys():
                    synApical = np.sum(syn['ApicalDendrite'].mesh)

            # Position of Presynapse/BoutonDensity
            tmpLine = tmpLine.replace('_invivo','')
            s = tmpLine.split('_')
            presynapse_y = s[-2]
            presynapse_x = s[-3]

            tmpStr = '%s,%s,%s,%s,%s,%s,%f,%f\n' % (tmpLine,cellName.split('/')[-1],
                                                    presynapse_x,presynapse_y,postsynapse_x,postsynapse_y,
                                                    synBasal,synApical)
            outputstr += tmpStr

    endTime = time.time()
    duration = (endTime - startTime)/60.0
    print 'Runtime: %.1f minutes' % duration

    with open(outputfname, 'w') as csvOutput:
        csvOutput.write(outputstr)

        print 'Saved data as csv in %s' % outputfname

if __name__ == '__main__':
    if len(sys.argv) == 3:
        fname = sys.argv[1]
        cellTypeName = sys.argv[2]
        map_singlecell_inputs_invitro(fname, cellTypeName)
    else:
        print 'Usage: python map_singlecell_inputs.py [morphology filename] [postsynaptic cell type name]'