'''
Created on Apr 21, 2014

@author: regger
'''

class CellTable(object):
    '''
    table of connectivity computation results per cell
    
    attributes:
        header
            - list of (CONNECTED_COLUMN,CONNECTED_CELLTYPE) dictionaries
              containing presynaptic columns and cell types in
              the same order as in the data table
        data
            - table containing one row per postsynaptic cell
            - for each postsynaptic cell, it contains a dictionary
              with the following information:
              CELLID    (NeuroNet cell ID corresponding to Morphologies.am)
              CELLTYPE
              COLUMN
              INSIDE_COLUMN
              TOTAL_SYNAPSES    (if present in input data)
              SOMA_POS
              SYN_PER_CONNECTED_CELL
                  list with number of synapses per presynaptic
                  cell type and column (ordered as in the header list)
    
    '''
    
    def __init__(self, header=None, data=None):
        '''
        Constructor
        '''
        self.header = header
        self.data = data
    
    def get_cell_ids(self, cellTypes=None, columns=None):
        cellTypeIDs = set()
        columnIDs = set()
        if cellTypes is not None:
            for i in range(len(self.data['CELLTYPE'])):
                postCellType = self.data['CELLTYPE'][i]
                if postCellType in cellTypes:
                    cellTypeIDs.add(i)
        if columns is not None:
            for i in range(len(self.data['COLUMN'])):
                postColumn = self.data['COLUMN'][i]
                if postColumn in columns:
                    columnIDs.add(i)
        IDs = []
        if len(cellTypeIDs) and len(columnIDs):
            cellTypeIDs.intersection_update(columnIDs)
            IDs = list(cellTypeIDs)
        elif len(cellTypeIDs) and columns is None:
            IDs = list(cellTypeIDs)
        elif len(columnIDs) and cellTypes is None:
            IDs = list(columnIDs)
        IDs.sort()
        return IDs
    
    def get_connected_type_ids(self, cellTypes=None, columns=None):
        cellTypeIDs = set()
        columnIDs = set()
        if cellTypes is not None:
            for i in range(len(self.header)):
                postCellType = self.header[i]['CONNECTED_CELLTYPE']
                if postCellType in cellTypes:
                    cellTypeIDs.add(i)
        if columns is not None:
            for i in range(len(self.header)):
                postColumn = self.header[i]['CONNECTED_COLUMN']
                if postColumn in columns:
                    columnIDs.add(i)
        IDs = []
        if len(cellTypeIDs) and len(columnIDs):
            cellTypeIDs.intersection_update(columnIDs)
            IDs = list(cellTypeIDs)
        elif len(cellTypeIDs) and columns is None:
            IDs = list(cellTypeIDs)
        elif len(columnIDs) and cellTypes is None:
            IDs = list(columnIDs)
        IDs.sort()
        return IDs
    
    