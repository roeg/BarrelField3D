import sys
import numpy as np
import connectome_analysis as ca
import matplotlib.pyplot as plt

columns = ('A1','A2','A3','A4','Alpha','B1','B2','B3','B4','Beta',\
           'C1','C2','C3','C4','Gamma','D1','D2','D3','D4','Delta',\
           'E1','E2','E3','E4')
columns3x3 = ('B1','B2','B3','C1','C2','C3','D1','D2','D3')
PC = ('C2',)
row = ('C1','C3')
arc = ('B2','D2')
diag = ('B1','B3','D1','D3')
otherColumns = ('A1','A2','A3','A4','Alpha','B4','Beta',\
         'C4','Gamma','D4','Delta','E1','E2','E3','E4')
allTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
            'Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
#allTypes = ('VPM','L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
#            'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
#            'L1','L23Trans','L45Sym','L45Peak','L56Trans')
allPostTypes = ('L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
            'Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
allPreTypes = ('L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct','VPM',\
            'Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
#allPostTypes = ('L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct',\
#            'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
#            'L1','L23Trans','L45Sym','L45Peak','L56Trans')
#VPMType = ('VPM',)
#allL4Types = ('VPM','L34','L4py','L4sp','L4ss',\
#            'SymLocal','L23Trans','L45Sym','L45Peak')
#exTypes = ('VPM','L34','L4py','L4sp','L4ss')
#inhTypes = ('SymLocal','L23Trans','L45Sym','L45Peak')
exTypes = ('L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct','VPM')
exPostTypes = ('L2','L34','L4py','L4sp','L4ss','L5st','L5tt','L6cc','L6ccinv','L6ct')
inhTypes = ('Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
#inhTypes = ('SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
#            'L1','L23Trans','L45Sym','L45Peak','L56Trans')
inhPostTypes = ('Local','L1','L23Trans','L45Sym','L45Peak','L56Trans')
#inhPostTypes = ('SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6',\
#            'L1','L23Trans','L45Sym','L45Peak','L56Trans')

columnColors = ('red','gray','darkgray','#d3d3d3','white')
exInhColors = ('black','red')
exCellTypeColors = ('cyan','blue','palegreen','darkgreen','lime','yellow',\
                    'orange','indigo','lavender','deeppink','black')
#inhCellTypeColors = ('blue','blue','blue','blue','blue','blue',\
#                     'black','cyan','palegreen','darkgreen','orange')
inhCellTypeColors = ('blue','black','cyan','palegreen','darkgreen','orange')

def create_synapsesPerCellTypePerColumn_pie_charts(fname):
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    figCount = 0
    for postCellType in allPostTypes:
        print 'Creating pie charts for postsynaptic cell type %s in PC %s ...' % (postCellType, PC[0])
        postCellIDs = synapseTable.get_connected_type_ids([postCellType], PC)
        postID = postCellIDs[0]
        if len(postCellIDs) > 1:
            errstr = 'Oops! there should not be more than 1 postsynaptic cell type selected...'
            raise RuntimeError(errstr)
        
#        PC/row/arc etc
        PCSyns = 0
        rowSyns = 0
        arcSyns = 0
        diagSyns = 0
        otherSyns = 0
        PCCellIDs = synapseTable.get_cell_ids(allTypes, PC)
        rowCellIDs = synapseTable.get_cell_ids(allTypes, row)
        arcCellIDs = synapseTable.get_cell_ids(allTypes, arc)
        diagCellIDs = synapseTable.get_cell_ids(allTypes, diag)
        otherColumnCellIDs = synapseTable.get_cell_ids(allTypes, otherColumns)
        for preID in PCCellIDs:
                PCSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in rowCellIDs:
                rowSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in arcCellIDs:
                arcSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in diagCellIDs:
                diagSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in otherColumnCellIDs:
                otherSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapses = np.array([PCSyns, rowSyns, arcSyns, diagSyns, otherSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapses, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_column_inputs.pdf'
        plt.savefig(pltName)
        
#        E/I ratio
        exSyns = 0
        inhSyns = 0
        exCellIDs = synapseTable.get_cell_ids(exTypes, columns)
        inhCellIDs = synapseTable.get_cell_ids(inhTypes, columns)
        for preID in exCellIDs:
                exSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in inhCellIDs:
                inhSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        exInhSynapses = np.array([exSyns, inhSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(exInhSynapses, colors=exInhColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exInhSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_ex_inh_inputs.pdf'
        plt.savefig(pltName)
        
#        EX PC/row/arc etc
        PCSynsEx = 0
        rowSynsEx = 0
        arcSynsEx = 0
        diagSynsEx = 0
        otherSynsEx = 0
        PCCellIDsEx = synapseTable.get_cell_ids(exTypes, PC)
        rowCellIDsEx = synapseTable.get_cell_ids(exTypes, row)
        arcCellIDsEx = synapseTable.get_cell_ids(exTypes, arc)
        diagCellIDsEx = synapseTable.get_cell_ids(exTypes, diag)
        otherColumnCellIDsEx = synapseTable.get_cell_ids(exTypes, otherColumns)
        for preID in PCCellIDsEx:
                PCSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in rowCellIDsEx:
                rowSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in arcCellIDsEx:
                arcSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in diagCellIDsEx:
                diagSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in otherColumnCellIDsEx:
                otherSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapsesEx = np.array([PCSynsEx, rowSynsEx, arcSynsEx, diagSynsEx, otherSynsEx])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapsesEx, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapsesEx.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_ex_column_inputs.pdf'
        plt.savefig(pltName)
        
#        INH PC/row/arc etc
        PCSynsInh = 0
        rowSynsInh = 0
        arcSynsInh = 0
        diagSynsInh = 0
        otherSynsInh = 0
        PCCellIDsInh = synapseTable.get_cell_ids(inhTypes, PC)
        rowCellIDsInh = synapseTable.get_cell_ids(inhTypes, row)
        arcCellIDsInh = synapseTable.get_cell_ids(inhTypes, arc)
        diagCellIDsInh = synapseTable.get_cell_ids(inhTypes, diag)
        otherColumnCellIDsInh = synapseTable.get_cell_ids(inhTypes, otherColumns)
        for preID in PCCellIDsInh:
                PCSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in rowCellIDsInh:
                rowSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in arcCellIDsInh:
                arcSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in diagCellIDsInh:
                diagSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for preID in otherColumnCellIDsInh:
                otherSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapsesInh = np.array([PCSynsInh, rowSynsInh, arcSynsInh, diagSynsInh, otherSynsInh])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapsesInh, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapsesInh.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_inh_column_inputs.pdf'
        plt.savefig(pltName)
        
#        all EX types
        exTypeSyns = []
        for tmpType in exTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], columns)
            for preID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            exTypeSyns.append(tmpSyns)
        
        exTypeSynapses = np.array(exTypeSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(exTypeSynapses, colors=exCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exTypeSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_ex_type_all_inputs.pdf'
        plt.savefig(pltName)
        
#        PC EX types
        exTypePCSyns = []
        for tmpType in exTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], PC)
            for preID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            exTypePCSyns.append(tmpSyns)
        
        exTypePCSynapses = np.array(exTypePCSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(exTypePCSynapses, colors=exCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exTypePCSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_ex_type_PC_inputs.pdf'
        plt.savefig(pltName)
        
#        all INH types
        inhTypeSyns = []
        for tmpType in inhTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], columns)
            for preID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            inhTypeSyns.append(tmpSyns)
        
        inhTypeSynapses = np.array(inhTypeSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(inhTypeSynapses, colors=inhCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % inhTypeSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_inh_type_all_inputs.pdf'
        plt.savefig(pltName)
        
#        PC INH types
        inhTypePCSyns = []
        for tmpType in inhTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], PC)
            for preID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            inhTypePCSyns.append(tmpSyns)
        
        inhTypePCSynapses = np.array(inhTypePCSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(inhTypePCSynapses, colors=inhCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % inhTypePCSynapses.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_inh_type_PC_inputs.pdf'
        plt.savefig(pltName)

def create_synapsesPerCellTypePerColumn_output_pie_charts(fname):
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    figCount = 0
    for preCellType in allPreTypes:
        print 'Creating pie charts for presynaptic cell type %s in PC %s ...' % (preCellType, PC[0])
        preCellIDs = synapseTable.get_cell_ids([preCellType], PC)
        preID = preCellIDs[0]
        if len(preCellIDs) > 1:
            errstr = 'Oops! there should not be more than 1 presynaptic cell type selected...'
            raise RuntimeError(errstr)
        
#        PC/row/arc etc
        PCSyns = 0
        rowSyns = 0
        arcSyns = 0
        diagSyns = 0
        PCCellIDs = synapseTable.get_connected_type_ids(allTypes, PC)
        rowCellIDs = synapseTable.get_connected_type_ids(allTypes, row)
        arcCellIDs = synapseTable.get_connected_type_ids(allTypes, arc)
        diagCellIDs = synapseTable.get_connected_type_ids(allTypes, diag)
        for postID in PCCellIDs:
                PCSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in rowCellIDs:
                rowSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in arcCellIDs:
                arcSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in diagCellIDs:
                diagSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapses = np.array([PCSyns, rowSyns, arcSyns, diagSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapses, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_column_outputs.pdf'
        plt.savefig(pltName)
        
#        E/I ratio
        exSyns = 0
        inhSyns = 0
        exCellIDs = synapseTable.get_connected_type_ids(exTypes, columns3x3)
        inhCellIDs = synapseTable.get_connected_type_ids(inhTypes, columns3x3)
        for postID in exCellIDs:
                exSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in inhCellIDs:
                inhSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        exInhSynapses = np.array([exSyns, inhSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(exInhSynapses, colors=exInhColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exInhSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_ex_inh_outputs.pdf'
        plt.savefig(pltName)
        
#        EX PC/row/arc etc
        PCSynsEx = 0
        rowSynsEx = 0
        arcSynsEx = 0
        diagSynsEx = 0
        PCCellIDsEx = synapseTable.get_connected_type_ids(exTypes, PC)
        rowCellIDsEx = synapseTable.get_connected_type_ids(exTypes, row)
        arcCellIDsEx = synapseTable.get_connected_type_ids(exTypes, arc)
        diagCellIDsEx = synapseTable.get_connected_type_ids(exTypes, diag)
        for postID in PCCellIDsEx:
                PCSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in rowCellIDsEx:
                rowSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in arcCellIDsEx:
                arcSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in diagCellIDsEx:
                diagSynsEx += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapsesEx = np.array([PCSynsEx, rowSynsEx, arcSynsEx, diagSynsEx])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapsesEx, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapsesEx.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_ex_column_outputs.pdf'
        plt.savefig(pltName)
        
#        INH PC/row/arc etc
        PCSynsInh = 0
        rowSynsInh = 0
        arcSynsInh = 0
        diagSynsInh = 0
        PCCellIDsInh = synapseTable.get_connected_type_ids(inhTypes, PC)
        rowCellIDsInh = synapseTable.get_connected_type_ids(inhTypes, row)
        arcCellIDsInh = synapseTable.get_connected_type_ids(inhTypes, arc)
        diagCellIDsInh = synapseTable.get_connected_type_ids(inhTypes, diag)
        for postID in PCCellIDsInh:
                PCSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in rowCellIDsInh:
                rowSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in arcCellIDsInh:
                arcSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for postID in diagCellIDsInh:
                diagSynsInh += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        columnSynapsesInh = np.array([PCSynsInh, rowSynsInh, arcSynsInh, diagSynsInh])
        figCount += 1
        plt.figure(figCount)
        plt.pie(columnSynapsesInh, colors=columnColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % columnSynapsesInh.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_inh_column_outputs.pdf'
        plt.savefig(pltName)
        
#        all EX types
        exTypeSyns = []
        for tmpType in exTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_connected_type_ids([tmpType], columns3x3)
            for postID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            exTypeSyns.append(tmpSyns)
        
        exTypeSynapses = np.array(exTypeSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(exTypeSynapses, colors=exCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exTypeSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_ex_type_all_outputs.pdf'
        plt.savefig(pltName)
        
#        PC EX types
        exTypePCSyns = []
        for tmpType in exTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_connected_type_ids([tmpType], PC)
            for postID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            exTypePCSyns.append(tmpSyns)
        
        exTypePCSynapses = np.array(exTypePCSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(exTypePCSynapses, colors=exCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % exTypePCSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_ex_type_PC_outputs.pdf'
        plt.savefig(pltName)
        
#        all INH types
        inhTypeSyns = []
        for tmpType in inhTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_connected_type_ids([tmpType], columns3x3)
            for postID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            inhTypeSyns.append(tmpSyns)
        
        inhTypeSynapses = np.array(inhTypeSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(inhTypeSynapses, colors=inhCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % inhTypeSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_inh_type_all_outputs.pdf'
        plt.savefig(pltName)
        
#        PC INH types
        inhTypePCSyns = []
        for tmpType in inhTypes:
            tmpSyns = 0
            tmpCellIDs = synapseTable.get_connected_type_ids([tmpType], PC)
            for postID in tmpCellIDs:
                tmpSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
            inhTypePCSyns.append(tmpSyns)
        
        inhTypePCSynapses = np.array(inhTypePCSyns)
        figCount += 1
        plt.figure(figCount)
        plt.pie(inhTypePCSynapses, colors=inhCellTypeColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % inhTypePCSynapses.sum()
        pltTitle = preCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + preCellType + '_inh_type_PC_outputs.pdf'
        plt.savefig(pltName)

def create_supraGran_infra_pie_charts(fname):
#    exSupraGranTypes = ('L2','L34','L4py','L4sp','L4ss')
    exSupraTypes = ('L2','L34')
    exGranTypes = ('L4py','L4sp','L4ss')
    exInfraTypes = ('L5st','L5tt','L6cc','L6ccinv','L6ct')
    supraGranInfraColors = ('Blue','Green','Red')
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    figCount = 0
    for postCellType in exPostTypes:
        print 'Creating pie charts for postsynaptic cell type %s in PC %s ...' % (postCellType, PC[0])
        postCellIDs = synapseTable.get_connected_type_ids([postCellType], PC)
        postID = postCellIDs[0]
        if len(postCellIDs) > 1:
            errstr = 'Oops! there should not be more than 1 postsynaptic cell type selected...'
            raise RuntimeError(errstr)
        
#        PC S/G/I
        supraSyns = 0
        granSyns = 0
        infraSyns = 0
        for tmpType in exSupraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], PC)
            for preID in tmpCellIDs:
                supraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exGranTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], PC)
            for preID in tmpCellIDs:
                granSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exInfraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], PC)
            for preID in tmpCellIDs:
                infraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        supraGranInfraTypeSyns = np.array([supraSyns, granSyns, infraSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(supraGranInfraTypeSyns, colors=supraGranInfraColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % supraGranInfraTypeSyns.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_SGI_PC_inputs.pdf'
        plt.savefig(pltName)
        
#        row S/G/I
        supraSyns = 0
        granSyns = 0
        infraSyns = 0
        for tmpType in exSupraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], row)
            for preID in tmpCellIDs:
                supraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exGranTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], row)
            for preID in tmpCellIDs:
                granSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exInfraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], row)
            for preID in tmpCellIDs:
                infraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        supraGranInfraTypeSyns = np.array([supraSyns, granSyns, infraSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(supraGranInfraTypeSyns, colors=supraGranInfraColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % supraGranInfraTypeSyns.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_SGI_row_inputs.pdf'
        plt.savefig(pltName)
        
#        arc S/G/I
        supraSyns = 0
        granSyns = 0
        infraSyns = 0
        for tmpType in exSupraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], arc)
            for preID in tmpCellIDs:
                supraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exGranTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], arc)
            for preID in tmpCellIDs:
                granSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exInfraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], arc)
            for preID in tmpCellIDs:
                infraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        supraGranInfraTypeSyns = np.array([supraSyns, granSyns, infraSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(supraGranInfraTypeSyns, colors=supraGranInfraColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % supraGranInfraTypeSyns.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_SGI_arc_inputs.pdf'
        plt.savefig(pltName)
        
#        diagonal S/G/I
        supraSyns = 0
        granSyns = 0
        infraSyns = 0
        for tmpType in exSupraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], diag)
            for preID in tmpCellIDs:
                supraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exGranTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], diag)
            for preID in tmpCellIDs:
                granSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        for tmpType in exInfraTypes:
            tmpCellIDs = synapseTable.get_cell_ids([tmpType], diag)
            for preID in tmpCellIDs:
                infraSyns += synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
        
        supraGranInfraTypeSyns = np.array([supraSyns, granSyns, infraSyns])
        figCount += 1
        plt.figure(figCount)
        plt.pie(supraGranInfraTypeSyns, colors=supraGranInfraColors, autopct='%d%%')
        plt.axis('equal')
        totalSynStr = '%d' % supraGranInfraTypeSyns.sum()
        pltTitle = postCellType + ' - total: ' + totalSynStr
        plt.title(pltTitle)
        pltName = fname[:-4] + '_' + postCellType + '_SGI_diag_inputs.pdf'
        plt.savefig(pltName)

def compute_celltype_integrator_distributor_index(fname):
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    cellTypeIntegratorEntropies = []
    for postCellType in exPostTypes:
        print 'Computing integrator index for postsynaptic cell type %s in PC %s ...' % (postCellType, PC[0])
        postCellIDs = synapseTable.get_connected_type_ids([postCellType], PC)
        postID = postCellIDs[0]
        if len(postCellIDs) > 1:
            errstr = 'Oops! there should not be more than 1 postsynaptic cell type selected...'
            raise RuntimeError(errstr)
        
        synapses = []
        for preCellType in exTypes:
            tmpCellIDs = synapseTable.get_cell_ids([preCellType], PC)
            preID = tmpCellIDs[0]
            synapses.append(synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID])
        synapses = np.array(synapses)
        totalSyn = synapses.sum()
        cellTypeEntropy = 0.0
        for i in range(len(synapses)):
            syn = synapses[i]
            preCellType = exTypes[i]
            p = syn/totalSyn
            cellTypeEntropy += -p*np.log(p+1e-6)
            print '    p(%s) = %.2f' % (preCellType, p)
        cellTypeIntegratorEntropies.append(cellTypeEntropy)
        print '  integrator entropy = %.2f' % cellTypeEntropy
    
    print '----------------------------'
    
    cellTypeDistributorEntropies = []
    for preCellType in exTypes:
        print 'Computing distributor index for presynaptic cell type %s in PC %s ...' % (preCellType, PC[0])
        preCellIDs = synapseTable.get_cell_ids([preCellType], PC)
        preID = preCellIDs[0]
        if len(preCellIDs) > 1:
            errstr = 'Oops! there should not be more than 1 presynaptic cell type selected...'
            raise RuntimeError(errstr)
        
        synapses = []
        for postCellType in exPostTypes:
            tmpCellIDs = synapseTable.get_connected_type_ids([postCellType], PC)
            postID = tmpCellIDs[0]
            synapses.append(synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID])
        synapses = np.array(synapses)
        totalSyn = synapses.sum()
        cellTypeEntropy = 0.0
        for i in range(len(synapses)):
            syn = synapses[i]
            postCellType = exPostTypes[i]
            p = syn/totalSyn
            cellTypeEntropy += -p*np.log(p+1e-6)
            print '    p(%s) = %.2f' % (postCellType, p)
        cellTypeDistributorEntropies.append(cellTypeEntropy)
        print '  distributor entropy = %.2f' % cellTypeEntropy
    
    print '----------------------------'
    
    x1 = np.arange(1, len(cellTypeIntegratorEntropies)+1, 1)
    plt.figure(1)
    plt.plot(x1, cellTypeIntegratorEntropies, 'bo', label='Inputs')
#    plt.title('Integrator index')
    
    x2 = np.arange(1, len(cellTypeDistributorEntropies)+1, 1)
#    plt.figure(2)
    plt.plot(x2, cellTypeDistributorEntropies, 'ro', label='Outputs')
    plt.title('low = specific; high = unspecific')
    plt.xticks(x2, exTypes, rotation='vertical')
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.legend()
    plt.show()
    
def create_synapsesPerCellTypePerColumn_connection_matrix(fname):
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    outName = fname[:-4] + '_connection_matrix_image.txt'
    with open(outName, 'w') as outFile:
        preTypesList = []
        preTypesList.extend(exTypes)
        preTypesList.extend(inhTypes)
        postTypesList = []
        postTypesList.extend(exPostTypes)
        postTypesList.extend(inhPostTypes)
        for preColumn in columns:
            for preType in preTypesList:
                tmpLine = ''
                preIDs = synapseTable.get_cell_ids([preType], [preColumn])
                if len(preIDs) > 1:
                    errstr = 'Oops! there should not be more than 1 presynaptic cell type selected...'
                    raise RuntimeError(errstr)
                preID = preIDs[0]
                for postColumn in columns3x3:
                    for postType in postTypesList:
                        postIDs = synapseTable.get_connected_type_ids([postType], [postColumn])
                        if len(postIDs) > 1:
                            errstr = 'Oops! there should not be more than 1 postsynaptic cell type selected...'
                            raise RuntimeError(errstr)
                        postID = postIDs[0]
                        nrSyn = synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
                        tmpLine += str(nrSyn) + '\t'
                line = tmpLine[:-1] + '\n'
                outFile.write(line)
    
def create_synapsesPerCellTypePerColumn_connection_matrix_3x3(fname):
    synapseTable = ca.load_synapses_per_celltype_per_column_table(fname)
    
    outName = fname[:-4] + '_connection_matrix_3x3_image.txt'
    with open(outName, 'w') as outFile:
        preTypesList = []
        preTypesList.extend(exTypes)
        preTypesList.extend(inhTypes)
        postTypesList = []
        postTypesList.extend(exPostTypes)
        postTypesList.extend(inhPostTypes)
        for preColumn in columns3x3:
            for preType in preTypesList:
                tmpLine = ''
                preIDs = synapseTable.get_cell_ids([preType], [preColumn])
                if len(preIDs) > 1:
                    errstr = 'Oops! there should not be more than 1 presynaptic cell type selected...'
                    raise RuntimeError(errstr)
                preID = preIDs[0]
                for postColumn in columns3x3:
                    for postType in postTypesList:
                        postIDs = synapseTable.get_connected_type_ids([postType], [postColumn])
                        if len(postIDs) > 1:
                            errstr = 'Oops! there should not be more than 1 postsynaptic cell type selected...'
                            raise RuntimeError(errstr)
                        postID = postIDs[0]
                        nrSyn = synapseTable.data['SYN_PER_CONNECTED_CELL'][preID][postID]
                        tmpLine += str(nrSyn) + '\t'
                line = tmpLine[:-1] + '\n'
                outFile.write(line)
        

if __name__=='__main__':
    if len(sys.argv) == 2:
        fname = sys.argv[1]
#        create_synapsesPerCellTypePerColumn_pie_charts(fname)
        create_synapsesPerCellTypePerColumn_output_pie_charts(fname)
#        create_synapsesPerCellTypePerColumn_connection_matrix(fname)
#        create_synapsesPerCellTypePerColumn_connection_matrix_3x3(fname)
#        create_supraGran_infra_pie_charts(fname)
#        compute_celltype_integrator_distributor_index(fname)
        