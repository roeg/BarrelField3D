/****************************************************************************/
/*                                                                          */
/* Program:   ConnectionMatrixVisualization                                 */
/*                                                                          */
/* File:      connectionMatrixVisualization.cpp                             */
/*                                                                          */
/* Purpose:   load NeuroNet connection matrix and save specified            */
/*            cell types and columns as image for visualization             */
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
#include "../../common/amiraReader.h"
#include "matrixanalyzer.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		ConnectionMatrix * connectome = new ConnectionMatrix;
		Reader * matrixReader = new Reader(inputFilename);
		matrixReader->readConnectionMatrix(connectome);
		delete matrixReader, delete connectome;
	}
	if(argc == 7)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		int argOffset = 3;
		ConnectionMatrix * connectome = new ConnectionMatrix;
		Reader * matrixReader = new Reader(inputFilename);
		matrixReader->readConnectionMatrix(connectome);
		
		MatrixAnalyzer * matrixAnalysis = new MatrixAnalyzer(connectome);
		std::vector< unsigned int > params = matrixAnalysis->parseInputParameters(argc, argv, argOffset);
		if(params.size() == 4)
		{
			unsigned int preColumn, preType, postColumn, postType;
			preColumn = params[0];
			preType = params[1];
			postColumn = params[2];
			postType = params[3];
			matrixAnalysis->writeConnectionMatrix(outputFilename, preColumn, preType, postColumn, postType);
		}
		
		delete matrixReader, delete connectome;
	}
	else
	{
		std::cout << "Parameters: [inputFilename] [outputFilename] [preColumn] [preType] [postColumn] [postType]" << std::endl;
	}
	
	return 0;
}
