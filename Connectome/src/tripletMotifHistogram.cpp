/****************************************************************************/
/*                                                                          */
/* Program:   TripletMotifHistogram                                         */
/*                                                                          */
/* File:      tripletMotifHistogram.cpp                                     */
/*                                                                          */
/* Purpose:   load NeuroNet connection matrix and compute triplet motif     */
/*            histogram of specified cell types and columns                 */
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
	if(argc == 8)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		unsigned int nrOfTriplets = atoi(argv[3]);
		int argOffset = 4;
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
// 			double binSize = atof(argv[3]);
			matrixAnalysis->analyzeTripletMotifs(outputFilename, preColumn, preType, postColumn, postType, nrOfTriplets);
// 			Profile * tripletMotifHist = matrixAnalysis->computeTripletMotifDistribution(preColumn, preType, postColumn, postType, nrOfTriplets);
// 			std::string motifOutName(outputFilename);
// 			motifOutName += "_motifHist.csv";
// 			tripletMotifHist->writeProfile(motifOutName.c_str());
// 			delete tripletMotifHist;
		}
		delete matrixReader, delete matrixAnalysis, delete connectome;
	}
	else
	{
		std::cout << "Parameters: [inputFilename] [outputFilename] [preColumn] [preType] [postColumn] [postType]" << std::endl;
	}
	
	return 0;
}
