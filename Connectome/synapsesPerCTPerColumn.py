import sys
import numpy as np
import connectome_analysis as ca

columns = ('Alpha','A1','A2','A3','A4','Beta','B1','B2','B3','B4',\
           'Gamma','C1','C2','C3','C4','Delta','D1','D2','D3','D4',\
           'E1','E2','E3','E4')
columns3x3 = ('B1','B2','B3','C1','C2','C3','D1','D2','D3')
allTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
            'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
            'L1','L23Trans','L45Sym','L45Peak','L56Trans')
#            'SymLocal','L1','L23Trans','L45Sym','L45Peak','L56Trans')
allJoinedTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
            'Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
#VPMType = ('VPM',)
#allL4Types = ('VPM','L34','L4py','L4sp','L4ss',\
#            'SymLocal','L23Trans','L45Sym','L45Peak')
#exTypes = ('VPM','L34','L4py','L4sp','L4ss')
#inhTypes = ('SymLocal','L23Trans','L45Sym','L45Peak')
exTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct')
#inhTypes = ('SymLocal','L1','L23Trans','L45Sym','L45Peak','L56Trans')
inhTypes = ('SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
            'L1','L23Trans','L45Sym','L45Peak','L56Trans')
inhJoinedTypes = ('Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')

joinCellTypes = ('SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
joinName = 'Local'

D2 = ('D2',)

def compute_synapsesPerCellTypePerColumn(fname):
    synapseTable = ca.load_synapse_table(fname)
    
    preIDs = {}
    summaryData = {}
    summaryPerPreTypeColumnData = {}
    for preColumn in columns:
        preIDs[preColumn] = {}
        summaryPerPreTypeColumnData[preColumn] = {}
        for preType in allTypes:
            tmpIDs = synapseTable.get_connected_type_ids([preType], [preColumn])
            if not tmpIDs:
                continue
            tmpPreType = preType
            if tmpPreType in joinCellTypes:
                tmpPreType = joinName
                if not preIDs[preColumn].has_key(tmpPreType):
                    preIDs[preColumn][tmpPreType] = [tmpIDs[0]]
                else:
                    preIDs[preColumn][tmpPreType].append(tmpIDs[0])
            else:
                preIDs[preColumn][tmpPreType] = tmpIDs[0]
            if not summaryPerPreTypeColumnData[preColumn].has_key(tmpPreType):
                summaryPerPreTypeColumnData[preColumn][tmpPreType] = {}
                for postColumn in columns3x3:
                    summaryPerPreTypeColumnData[preColumn][tmpPreType][postColumn] = {}
                    for postCellType in allTypes:
                        tmpPostType = postCellType
                        if tmpPostType in joinCellTypes:
                            tmpPostType = joinName
                        if not summaryPerPreTypeColumnData[preColumn].has_key(tmpPostType):
                            summaryPerPreTypeColumnData[preColumn][tmpPreType][postColumn][tmpPostType] = []
    for postColumn in columns3x3:
        summaryData[postColumn] = {}
        for postCellType in allTypes:
            tmpPostType = postCellType
            if tmpPostType in joinCellTypes:
                tmpPostType = joinName
            if not summaryData[postColumn].has_key(tmpPostType):
                summaryData[postColumn][tmpPostType] = {}
    
    for postColumn in columns3x3:
        analyzedJoinCellTypes = False
        for postCellType in allTypes:
            postCellTypeName = postCellType
            print 'Analyzing postsynaptic column/cell type %s/%s ...' % (postColumn, postCellType)
            postCellIDs = synapseTable.get_cell_ids([postCellType], [postColumn])
            if postCellType in joinCellTypes and not analyzedJoinCellTypes:
                postCellTypeName = joinName
                analyzedJoinCellTypes = True
                postCellIDs = []
                for tmpPostType in joinCellTypes:
                    print 'Analyzing postsynaptic column/cell type %s/%s ...' % (postColumn, tmpPostType)
                    postCellIDs += synapseTable.get_cell_ids([tmpPostType], [postColumn])
            elif postCellType in joinCellTypes and analyzedJoinCellTypes:
                print '\tAlready analyzed; skipping to next column/cell type combination!'
                continue
            
            nCells = 0
            totalSyns = []
            preTypeColumnSyns = {}
            for preColumn in columns:
                preTypeColumnSyns[preColumn] = {}
#                for preCellType in allTypes:
                for preCellType in allJoinedTypes:
                    preTypeColumnSyns[preColumn][preCellType] = []
            for postID in postCellIDs:
                nCells += 1
                totalSynsTmp = synapseTable.data['TOTAL_SYNAPSES'][postID]
                totalSyns.append(totalSynsTmp)
                for preColumn in columns:
#                    for preCellType in allTypes:
                    for preCellType in allJoinedTypes:
                        if preCellType not in preIDs[preColumn].keys():
                            continue
                        preID = preIDs[preColumn][preCellType]
#                        if joined types: iterate over pre columns and add;
#                        otherwise just read value from column
                        try:
                            preColumnTypeSynsTmp = 0.0
                            for col in preID:
                                preColumnTypeSynsTmp += synapseTable.data['SYN_PER_CONNECTED_CELL'][postID][col]
                            preTypeColumnSyns[preColumn][preCellType].append(preColumnTypeSynsTmp)
                        except TypeError:
                            preColumnTypeSynsTmp = synapseTable.data['SYN_PER_CONNECTED_CELL'][postID][preID]
                            preTypeColumnSyns[preColumn][preCellType].append(preColumnTypeSynsTmp)
            
            if nCells:
                summaryData[postColumn][postCellTypeName]['Number of cells'] = nCells
                summaryData[postColumn][postCellTypeName]['Synapses per cell'] = [np.mean(totalSyns), np.std(totalSyns), np.min(totalSyns), np.max(totalSyns)]
                for preColumn in preTypeColumnSyns.keys():
                    for preCellType in preTypeColumnSyns[preColumn].keys():
                        if not preTypeColumnSyns[preColumn][preCellType]:
                            continue
                        avgSyn = np.mean(preTypeColumnSyns[preColumn][preCellType])
                        STDSyn = np.std(preTypeColumnSyns[preColumn][preCellType])
                        minSyn = np.min(preTypeColumnSyns[preColumn][preCellType])
                        maxSyn = np.max(preTypeColumnSyns[preColumn][preCellType])
                        summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellTypeName] = [avgSyn, STDSyn, minSyn, maxSyn]
    
    prefix = ''
    if fname.endswith('.csv'):
        prefix = fname[:-4]
    elif fname.endswith('.am'):
        prefix = fname[:-3]
    outName = prefix + '_total_synapses_per_cell.csv'
    with open(outName, 'w') as outFile:
        header = 'Column\tcell type\tE/I\tquantity\tmean\tSTD\tmin\tmax\n'
        outFile.write(header)
        for column in summaryData.keys():
            for cellType in summaryData[column].keys():
                for quantity in summaryData[column][cellType].keys():
                    funcType = 'E'
#                    if cellType in inhTypes:
                    if cellType in inhJoinedTypes:
                        funcType = 'I'
                    if quantity == 'Number of cells':
                        line = column + '\t' + cellType + '\t' + funcType + '\t' + quantity + '\t' + str(summaryData[column][cellType][quantity]) + '\n'
                        outFile.write(line)
                    else:
                        line = column + '\t' + cellType + '\t' + funcType + '\t' + quantity
                        for i in range(len(summaryData[column][cellType][quantity])):
                            line += '\t' + str(summaryData[column][cellType][quantity][i])
                        line += '\n'
                        outFile.write(line)
    
    outName = prefix + '_synapses_per_pre_celltype_column.csv'
    with open(outName, 'w') as outFile:
        header = 'Pre column\tpre cell type\t pre E/I\tpost column\tpost cell type\tpost E/I\tavg\tSTD\tmin\tmax\n'
        outFile.write(header)
        for preColumn in summaryPerPreTypeColumnData.keys():
            for preCellType in summaryPerPreTypeColumnData[preColumn].keys():
                for postColumn in summaryPerPreTypeColumnData[preColumn][preCellType].keys():
                    for postCellType in summaryPerPreTypeColumnData[preColumn][preCellType][postColumn].keys():
                        preFuncType = 'E'
#                        if preCellType in inhTypes:
                        if preCellType in inhJoinedTypes:
                            preFuncType = 'I'
                        postFuncType = 'E'
#                        if postCellType in inhTypes:
                        if postCellType in inhJoinedTypes:
                            postFuncType = 'I'
                        line = preColumn + '\t' + preCellType + '\t' + preFuncType + '\t' + postColumn + '\t' + postCellType + '\t' + postFuncType
                        tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
                        if not len(tmpData):
                            continue
                        for i in range(len(tmpData)):
                            line += '\t' + str(tmpData[i])
                        line += '\n'
                        outFile.write(line)
    
    outName = prefix + '_synapses_per_pre_celltype_column_avg_matrix.csv'
    with open(outName, 'w') as outFile:
        header = 'Pre column\tpre cell type\t pre E/I'
        for postColumn in columns3x3:
            for postCellType in exTypes:
                header += '\t' + postColumn + '_' + postCellType
#            for postCellType in inhTypes:
            for postCellType in inhJoinedTypes:
                header += '\t' + postColumn + '_' + postCellType
        header += '\n'
        outFile.write(header)
        for preColumn in summaryPerPreTypeColumnData.keys():
            for preCellType in summaryPerPreTypeColumnData[preColumn].keys():
                preFuncType = 'E'
#                if preCellType in inhTypes:
                if preCellType in inhJoinedTypes:
                    preFuncType = 'I'
                line = preColumn + '\t' + preCellType + '\t' + preFuncType
                for postColumn in columns3x3:
                    for postCellType in exTypes:
                        line += '\t'
                        try:
                            tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
                            if not len(tmpData):
                                line += '0'
                                continue
                            line += str(tmpData[0])
                        except KeyError:
                            line += '0'
                            continue
#                    for postCellType in inhTypes:
                    for postCellType in inhJoinedTypes:
                        line += '\t'
                        try:
                            tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
                            if not len(tmpData):
                                line += '0'
                                continue
                            line += str(tmpData[0])
                        except KeyError:
                            line += '0'
                            continue
                line += '\n'
                outFile.write(line)
    
##    BUG: nCells not available for presynaptic cell type VPM 
#    outName = prefix + '_innervation_per_pre_cell_column_3x3.csv'
#    with open(outName, 'w') as outFile:
#        header = 'Pre column\tpre cell type\t pre E/I\tpost column\tpost cell type\tpost E/I\tavg\tSTD\tmin\tmax\n'
#        outFile.write(header)
#        for preColumn in summaryPerPreTypeColumnData.keys():
#            if preColumn not in columns3x3:
#                continue
#            for preCellType in summaryPerPreTypeColumnData[preColumn].keys():
#                for postColumn in summaryPerPreTypeColumnData[preColumn][preCellType].keys():
#                    for postCellType in summaryPerPreTypeColumnData[preColumn][preCellType][postColumn].keys():
#                        preFuncType = 'E'
#                        if preCellType in inhTypes:
#                            preFuncType = 'I'
#                        postFuncType = 'E'
#                        if postCellType in inhTypes:
#                            postFuncType = 'I'
#                        line = preColumn + '\t' + preCellType + '\t' + preFuncType + '\t' + postColumn + '\t' + postCellType + '\t' + postFuncType
#                        tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
#                        if not len(tmpData):
#                            continue
#                        for i in range(len(tmpData)):
#                            nCells = summaryData[preColumn][preCellType]['Number of cells']
#                            avgInnervation = tmpData[i]/nCells
#                            line += '\t' + str(avgInnervation)
#                        line += '\n'
#                        outFile.write(line)
#    
#    outName = prefix + '_innervation_per_pre_cell_column_3x3_avg_matrix.csv'
#    with open(outName, 'w') as outFile:
#        header = 'Pre column\tpre cell type\t pre E/I'
#        for postColumn in columns3x3:
#            for postCellType in exTypes:
#                header += '\t' + postColumn + '_' + postCellType
#            for postCellType in inhTypes:
#                header += '\t' + postColumn + '_' + postCellType
#        header += '\n'
#        outFile.write(header)
#        for preColumn in summaryPerPreTypeColumnData.keys():
#            for preCellType in summaryPerPreTypeColumnData[preColumn].keys():
#                if preColumn not in columns3x3:
#                    continue
#                preFuncType = 'E'
#                if preCellType in inhTypes:
#                    preFuncType = 'I'
#                line = preColumn + '\t' + preCellType + '\t' + preFuncType
#                for postColumn in columns3x3:
#                    for postCellType in exTypes:
#                        line += '\t'
#                        try:
#                            tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
#                            if not len(tmpData):
#                                line += '0'
#                                continue
#                            nCells = summaryData[preColumn][preCellType]['Number of cells']
#                            avgInnervation = tmpData[0]/nCells
#                            line += str(avgInnervation)
#                        except KeyError:
#                            line += '0'
#                            continue
#                    for postCellType in inhTypes:
#                        line += '\t'
#                        try:
#                            tmpData = summaryPerPreTypeColumnData[preColumn][preCellType][postColumn][postCellType]
#                            if not len(tmpData):
#                                line += '0'
#                                continue
#                            nCells = summaryData[preColumn][preCellType]['Number of cells']
#                            avgInnervation = tmpData[0]/nCells
#                            line += str(avgInnervation)
#                        except KeyError:
#                            line += '0'
#                            continue
#                line += '\n'
#                outFile.write(line)

if __name__=='__main__':
    if len(sys.argv) == 2:
        fname = sys.argv[1]
        compute_synapsesPerCellTypePerColumn(fname)
        