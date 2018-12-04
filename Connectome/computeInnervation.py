#===============================================================================
# computeInnervation
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

def computeInnervation(dendriteListPath, dendriteFolderName, cellTypeName, boutonListPath, boutonFolderName, outputCSVName):

    if not (cellTypeName in exTypes) and not (cellTypeName in inhTypes):
        errstr = 'Unknown cell type %s!' % cellTypeName
        print 'ERROR! %s' % errstr
        raise TypeError(errstr)

    startTime = time.time()

    # Main Path
    if sys.platform == 'win32':
        prefix = '//nas1/nas1/Data_daniel/Network/L5/L5_InVitro/'
    else:
        prefix = '/nas1/Data_daniel/Network/L5/L5_InVitro/'

    # Required data for calculation synapse density
    connectionsSpreadsheetName = os.path.join(prefix,'data/CellTypeConnectionSpreadsheet.csv')
    ExPSTDensityName = os.path.join(prefix,'data/EXNormalizationPSTs.am')
    connectionsSpreadsheet = sim.read_connections_spreadsheet(connectionsSpreadsheetName)
    ExPSTDensity = sim.read_scalar_field(ExPSTDensityName)
    ExPSTDensity.resize_mesh()
    InhPSTDensity = []
    boutonDensityList = []
    boutonDensityListName = []

    # Load each bouton density and store it
    with open(boutonListPath, 'r') as spreadsheet:
        for line in spreadsheet:

            # Load Bouton Density
            tmpLine = line.strip('\n\r')
            densityName = os.path.join(boutonFolderName,tmpLine)
            boutonDensity = sim.read_scalar_field(densityName)
            boutonDensity.resize_mesh()
            boutonDensityList.append(boutonDensity)
            boutonDensityListName.append(tmpLine)

    if len(boutonDensityList) == 0:
        errstr = 'No Bouton Densities were found'
        print 'ERROR! %s' % errstr
        raise TypeError(errstr)

    outputstr = 'Presynapse,Postsynapse,synapsesBasal,synapsesApical\n'

    # Go through all dendrites
    with open(dendriteListPath, 'r') as spreadsheet:
        for line in spreadsheet:

            # Read in Dendrite
            tmpLine = line.strip('\n\r')
            cellNameDend = os.path.join(dendriteFolderName,tmpLine)
            if cellNameDend[-4:] != '.hoc':
                errstr = 'Can only read in hoc morphologies: %s!' % cellNameDend
                print 'ERROR! %s' % errstr
                raise TypeError(errstr)

            # Load Postsynaptic hoc Morphology
            print 'Loading cell morphology... %s' % tmpLine
            parser = sim.CellParser(cellNameDend)
            parser.spatialgraph_to_cell()
            singleCell = parser.get_cell()

            synapseDensityComputation = sim.SynapseDensity(singleCell, cellTypeName, connectionsSpreadsheet,\
                                               exTypes, inhTypes, ExPSTDensity, InhPSTDensity)

            # Compute Innervation
            for b in range(len(boutonDensityList)):

                syn = synapseDensityComputation.compute_synapse_density(boutonDensityList[b], cellTypeName)

                # Extract Innervation, if tmp == None, no innervation
                synApical = 0
                synBasal = 0

                if syn is not None:
                    if 'Dendrite' in syn.keys():
                        synBasal = np.sum(syn['Dendrite'].mesh)
                    if 'ApicalDendrite' in syn.keys():
                        synApical = np.sum(syn['ApicalDendrite'].mesh)

                # Name of Pre- and Postsynapse
                dendName = cellNameDend.split('/')[-1]
                dendName = dendName.replace(".am","")
                dendName = dendName.replace(".hoc","")
                boutonName = boutonDensityListName[b].replace(".am","")
                tmpStr = '%s,%s,%f,%f\n' % (boutonName,dendName,synBasal,synApical)
                outputstr += tmpStr

    # Display runtime
    endTime = time.time()
    duration = (endTime - startTime)/60.0
    print 'Runtime: %.1f minutes' % duration

    # Store Innervation as .csv file
    with open(outputCSVName, 'w') as csvOutput:
        csvOutput.write(outputstr)

        print 'Saved data as csv in %s' % outputCSVName

if __name__ == '__main__':
    if len(sys.argv) == 7:
        dendriteListPath = sys.argv[1]
        dendriteFolderName = sys.argv[2]
        cellTypeName = sys.argv[3]
        boutonListPath = sys.argv[4]
        boutonFolderName = sys.argv[5]
        outputFolderName = sys.argv[6]

        computeInnervation(dendriteListPath, dendriteFolderName, cellTypeName, boutonListPath, boutonFolderName, outputFolderName)
    else:
        print 'Usage: python computeInnervation.py dendriteListPath dendriteFolderName cellTypeName boutonListPath boutonFolderName outputFolderName'