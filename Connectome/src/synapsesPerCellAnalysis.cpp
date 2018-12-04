/****************************************************************************/
/*                                                                          */
/* Program:   ConnectionMatrixVisualization                                 */
/*                                                                          */
/* File:      uniqueConnectionMatrixRows.cpp                                */
/*                                                                          */
/* Purpose:   load NeuroNet connection matrix and save only unique rows     */
/*            (i.e. individual presynaptic axons)                           */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max Planck Institute for Biological Cybernetics               */
/*            Spemannstr. 38-44                                             */
/*            72076 Tuebingen                                               */
/*            Germany                                                       */
/*                                                                          */
/* EMail:     robert.egger@tuebingen.mpg.de                                 */
/*                                                                          */
/* History:   16.06.2014                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "matrixanalyzer.h"

int main( int argc , char * argv[])
{
	if(argc == 5)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		int argOffset = 3;
		
		CellTable * table = new CellTable;
		Reader * cellTableReader = new Reader(inputFilename, inputFilename);
		cellTableReader->readSynapsesPerCellTable(table);
		
		MatrixAnalyzer * matrixAnalysis = new MatrixAnalyzer();
		std::vector< unsigned int > params = matrixAnalysis->parseInputParameters(argc, argv, argOffset);
		if(params.size() == 2)
		{
			unsigned int postColumn, postType;
			postColumn = params[0];
			postType = params[1];
			std::list< unsigned int > postColumns, postTypes;
			postColumns.push_back(postColumn);
			postTypes.push_back(postType);
			matrixAnalysis->writeSynapsesPerCellRows(outputFilename, postColumns, postTypes, table);
		}
		
		delete cellTableReader, delete table;
	}
	else
	{
		std::cout << "Parameters: [inputFilename] [outputFilename] [post column] [post cell type] " << std::endl;
	}
	
	return 0;
}


