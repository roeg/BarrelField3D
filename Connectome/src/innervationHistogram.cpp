/****************************************************************************/
/*                                                                          */
/* Program:   InnervationHistogram                                          */
/*                                                                          */
/* File:      innervationHistogram.cpp                                      */
/*                                                                          */
/* Purpose:   load NeuroNet connection matrix and save innervation          */
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
// 			double binSize = atof(argv[3]);
			double innervationBinSize = 0.1;
			double probabilityBinSize = 0.025;
			Profile * innervationHist = matrixAnalysis->computeInnervationHistogram(preColumn, preType, postColumn, postType, innervationBinSize);
			std::string innervationOutName(outputFilename);
			innervationOutName += "_innervation.csv";
			innervationHist->writeProfile(innervationOutName.c_str(), 0.5*innervationBinSize);
			Profile * probabilityHist = matrixAnalysis->computeProbabilityHistogram(preColumn, preType, postColumn, postType, probabilityBinSize);
			std::string probabilityOutName(outputFilename);
			probabilityOutName += "_prob.csv";
			probabilityHist->writeProfile(probabilityOutName.c_str(), 0.5*probabilityBinSize);
			Profile * synapseNumberHist = matrixAnalysis->computeSynapseNumberHistogram(preColumn, preType, postColumn, postType);
			std::string synapseOutName(outputFilename);
			synapseOutName += "_nSyn.csv";
			synapseNumberHist->writeProfile(synapseOutName.c_str());
			delete innervationHist;
		}
		delete matrixReader, delete matrixAnalysis, delete connectome;
	}
	else
	{
		std::cout << "Parameters: [inputFilename] [outputFilename] [preColumn] [preType] [postColumn] [postType]" << std::endl;
	}
	
	return 0;
}
