'''
Created on Apr 21, 2014

I/O for NeuroNet output files generated with Amira

@author: regger
'''

import numpy as np
from celltable import *

def load_csv_synapse_table(fname):
    header, raw_data, totalFlag = _load_raw_file(fname, transpose=1)
    
    positions = []
    for i in range(len(raw_data[3])):
        pos = (raw_data[3][i], raw_data[4][i], raw_data[5][i])
        positions.append(pos)
    positions = np.array(positions)
#    synPerCell = []
#    for i in range(len(raw_data[7])):
#        syns = [raw_data[j][i] for j in range(8, len(raw_data))]
#        synPerCell.append(syns)
#    synPerCell = np.array(synPerCell)
    synPerCellOffset = 7
    if totalFlag:
        synPerCellOffset = 8
    synPerCell = []
    for i in range(len(raw_data[0])):
        syns = [raw_data[j][i] for j in range(synPerCellOffset, len(raw_data))]
        synPerCell.append(syns)
    synPerCell = np.array(synPerCell)
    
    data = {}
    data['CELLID'] = raw_data[0]
    data['CELLTYPE'] = raw_data[1]
    data['COLUMN'] = raw_data[2]
    data['INSIDE_COLUMN'] = raw_data[6]
    if totalFlag:
        data['TOTAL_SYNAPSES'] = np.array(raw_data[7])
    data['SOMA_POS'] = positions
    data['SYN_PER_CONNECTED_CELL'] = synPerCell
    return CellTable(header, data)

def _load_raw_file(fname, skiprows=0, transpose=0):
    print 'Reading csv synapse table %s ...' % fname
    data = []
    header = []
    totalSynFlag = False
    foundHeader = False
#    if skiprows:
#        headerLine = skiprows-1
#    else:
#        headerLine = 0
    with open(fname, 'r') as dataFile:
#        lineCnt = 0
        for line in dataFile:
            if not foundHeader:
#            if lineCnt == headerLine:
                splitHeader = []
                if line.find('\t') > -1:
                    splitHeader = line.strip().split('\t')
                elif line.find(',') > -1:
                    splitHeader = line.strip().split(',')
                else:
                    splitHeader = [line.strip()]
                if splitHeader[0] == 'CELLID':
                    foundHeader = True
                    if 'TOTAL' in splitHeader[7]:
                        totalSynFlag = True
                    for i in range(7+totalSynFlag, len(splitHeader)):
                        splitInfo = splitHeader[i].split('_')
                        infoDict = {}
                        infoDict['CONNECTED_COLUMN'] = splitInfo[0]
                        if 'axon' in splitInfo[1]:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1][:-4]
                        else:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1]
                        header.append(infoDict)
#                    column = splitInfo[0]
#                    preCelltype = splitInfo[1][:-4]
#                    header.append((column, preCelltype))
#            lineCnt += 1
                continue
#            if lineCnt <= skiprows:
#                continue
            splitLine = []
            if line.find('\t') > -1:
                splitLine = line.strip().split('\t')
            elif line.find(',') > -1:
                splitLine = line.strip().split(',')
            else:
                splitLine = [line.strip()]
            cellID = int(splitLine[0])
            cellType = splitLine[1]
            column = splitLine[2]
            somax = float(splitLine[3])
            somay = float(splitLine[4])
            somaz = float(splitLine[5])
            insideColumn = int(splitLine[6])
            if totalSynFlag:
                totalSyn = float(splitLine[7])
            lineData = []
            lineData.append(cellID)
            lineData.append(cellType)
            lineData.append(column)
            lineData.append(somax)
            lineData.append(somay)
            lineData.append(somaz)
            lineData.append(insideColumn)
            if totalSynFlag:
                lineData.append(totalSyn)
            for i in range(7+totalSynFlag, len(splitLine)):
                lineData.append(float(splitLine[i]))
            data.append(lineData)
    if not transpose:
        return header, data, totalSynFlag
    else:
        dataTrans = [[] for n in range(len(data[0]))]
        for i in range(len(data)):
            for j in range(len(data[i])):
                dataTrans[j].append(data[i][j])
        return header, dataTrans, totalSynFlag

def load_amira_synapse_table(fname):
    print 'Reading AmiraMesh synapse table %s ...' % fname
    header = []
    data = {}
    dataTrans = []
    data['CELLID'] = []
    data['CELLTYPE'] = []
    data['COLUMN'] = []
    data['INSIDE_COLUMN'] = []
    somaX = []
    somaY = []
    somaZ = []
    
    foundHeader = False
    foundReferences = False
    dataReferences = {}
    tmpStr = ''
    totalSynFlag = False
    currentColReadFlag = -1
    currentRefFlag = -1
    with open(fname, 'r') as dataFile:
        for line in dataFile:
            line = line.strip()
            if not line:
                continue
            if 'define' in line:
                continue
            if not foundHeader:
#                implement header and data references
                if '}' in line:
#                    only checking TOTAL_SYNAPSES table,
#                    not AVERAGE_SYNAPSES table
                    foundHeader = True
                if 'Column' in line:
                    splitLine = line.strip().split(' ')
                    colNr = int(splitLine[0][-4:])
                    if colNr == 7:
                        if 'TOTAL' in splitLine[1]:
                            totalSynFlag = True
                            data['TOTAL_SYNAPSES'] = []
                        if not totalSynFlag:
                            splitInfo = splitLine[1].strip('\",').split('_')
                            infoDict = {}
                            infoDict['CONNECTED_COLUMN'] = splitInfo[0]
                            if 'axon' in splitInfo[1]:
                                infoDict['CONNECTED_CELLTYPE'] = splitInfo[1][:-4]
                            else:
                                infoDict['CONNECTED_CELLTYPE'] = splitInfo[1]
                            header.append(infoDict)
                            dataTrans.append([])
                    if colNr > 7:
                        splitInfo = splitLine[1].strip('\",').split('_')
                        infoDict = {}
                        infoDict['CONNECTED_COLUMN'] = splitInfo[0]
                        if 'axon' in splitInfo[1]:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1][:-4]
                        else:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1]
                        header.append(infoDict)
                        dataTrans.append([])
                continue
            if not foundReferences:
                if 'Table0000Column' in line:
                    splitLine = line.strip().split(' ')
                    colNr = int(splitLine[0][-4:])
                    refNr = int(splitLine[-1][1:])
                    dataReferences[refNr] = colNr
                if len(dataReferences) == 7 + totalSynFlag + len(header):
                    foundReferences = True
                continue
#            implement actual data reading
            if line[0] == '@':
                currentRefFlag = int(line[1:])
                if currentRefFlag in dataReferences.keys():
                    currentColReadFlag = dataReferences[currentRefFlag]
                continue
            if currentRefFlag not in dataReferences.keys():
                continue
            if currentColReadFlag == 0:
#                CELLID (int)
                data['CELLID'].append(int(line))
            if currentColReadFlag == 1:
#                CELLTYPE (ASCII 0-terminated)
                if int(line) == 0:
                    data['CELLTYPE'].append(tmpStr)
                    tmpStr = ''
                else:
                    tmpStr += chr(int(line))
            if currentColReadFlag == 2:
#                COLUMN (ASCII 0-terminated)
                if int(line) == 0:
                    data['COLUMN'].append(tmpStr)
                    tmpStr = ''
                else:
                    tmpStr += chr(int(line))
            if currentColReadFlag == 3:
#                SOMA_X (float)
                somaX.append(float(line))
            if currentColReadFlag == 4:
#                SOMA_Y (float)
                somaY.append(float(line))
            if currentColReadFlag == 5:
#                SOMA_Z (float)
                somaZ.append(float(line))
            if currentColReadFlag == 6:
#                INSIDE_COLUMN (int)
                data['INSIDE_COLUMN'].append(int(line))
            if currentColReadFlag == 7 and totalSynFlag:
#                TOTAL_SYNAPSES (float)
                data['TOTAL_SYNAPSES'].append(float(line))
            if currentColReadFlag == 7 and not totalSynFlag:
#                synapses (float)
                colIndex = currentColReadFlag - totalSynFlag - 7
                dataTrans[colIndex].append(float(line))
            if currentColReadFlag > 7:
                colIndex = currentColReadFlag - totalSynFlag - 7
                dataTrans[colIndex].append(float(line))
    
    positions = [(somaX[i], somaY[i], somaZ[i]) for i in range(len(somaX))]
    data['SOMA_POS'] = np.array(positions)
#    synPerCell = [[dataTrans[j][i] for j in range(len(dataTrans))] for i in range(len(data['CELLID']))]
#    data['SYN_PER_CONNECTED_CELL'] = np.array(synPerCell)
    data['SYN_PER_CONNECTED_CELL'] = np.array(dataTrans).transpose()
    
    return CellTable(header, data)

def load_csv_avg_synapses_table(fname):
    header, raw_data = _load_raw_avg_file(fname, transpose=1)
    
#    skipped E/I type during reading already...
    synPerCellOffset = 2
    synPerCell = []
    for i in range(len(raw_data[0])):
        syns = [raw_data[j][i] for j in range(synPerCellOffset, len(raw_data))]
        synPerCell.append(syns)
    synPerCell = np.array(synPerCell)
    
    data = {}
    data['CELLID'] = None
    data['COLUMN'] = raw_data[0]
    data['CELLTYPE'] = raw_data[1]
    data['INSIDE_COLUMN'] = None
    data['TOTAL_SYNAPSES'] = None
    data['SOMA_POS'] = None
    data['SYN_PER_CONNECTED_CELL'] = synPerCell
    return CellTable(header, data)

def _load_raw_avg_file(fname, skiprows=0, transpose=0):
    print 'Reading csv average synapse table %s ...' % fname
    data = []
    header = []
    foundHeader = False
#    if skiprows:
#        headerLine = skiprows-1
#    else:
#        headerLine = 0
    with open(fname, 'r') as dataFile:
#        lineCnt = 0
        for line in dataFile:
            if not foundHeader:
#            if lineCnt == headerLine:
                splitHeader = line.strip().split('\t')
                if splitHeader[0] == 'Pre column':
                    foundHeader = True
                    for i in range(3, len(splitHeader)):
                        splitInfo = splitHeader[i].split('_')
                        infoDict = {}
                        infoDict['CONNECTED_COLUMN'] = splitInfo[0]
                        if 'axon' in splitInfo[1]:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1][:-4]
                        else:
                            infoDict['CONNECTED_CELLTYPE'] = splitInfo[1]
                        header.append(infoDict)
#                    column = splitInfo[0]
#                    preCelltype = splitInfo[1][:-4]
#                    header.append((column, preCelltype))
#            lineCnt += 1
                continue
#            if lineCnt <= skiprows:
#                continue
            splitLine = line.strip().split('\t')
            column = splitLine[0]
            cellType = splitLine[1]
            lineData = []
            lineData.append(column)
            lineData.append(cellType)
            for i in range(3, len(splitLine)):
                lineData.append(float(splitLine[i]))
            data.append(lineData)
            
    if not transpose:
        return header, data
    else:
        dataTrans = [[] for n in range(len(data[0]))]
        for i in range(len(data)):
            for j in range(len(data[i])):
                dataTrans[j].append(data[i][j])
        return header, dataTrans

