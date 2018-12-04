/****************************************************************************/
/*                                                                          */
/* Program:   MorphAnalyzer                                                 */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   Program for analysis of registered neuron morphologies with   */
/*            respect to columns, septa and layers in standardized barrel   */
/*            cortex. E.g., computes distribution of axon in different      */
/*            columns/layers, and computes z-profiles of axons taking local */
/*            orientation into account                                      */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "morph_analyzer.h"
#include "../../common/inputcheckpoint.h"

int main( int argc , char * argv[])
{
	if(argc == 2 || argc == 3)
	{
		const char * inputFilename = argv[1];
		double binSize = 50;
		if(argc == 3)
		{
			binSize = atof(argv[2]);
		}
		
		std::string ifName(inputFilename);
		size_t suffix;
		
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".hoc") != std::string::npos)
		{
			hocReader->readHocFile();
			suffix = ifName.find(".hoc");
		}
		else if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(0);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			return 0;
		}
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		if(checkPoint->getParameters().axonFlag)
		{
			Analyzer * morphAnalyzer = new Analyzer(hocReader->getSpatialGraph(), checkPoint->getParameters());
			morphAnalyzer->computeAxonClusterParameters(ofName.c_str(), binSize);
			
// 			std::vector< double > zProfile = morphAnalyzer->compute1DProfileLocally("Axon", 50);
// 			ofName += "_axon_1D_profile.csv";
// 			if(zProfile.size())
// 			{
// 				std::ofstream ProfileWriter;
// 				ProfileWriter.open(ofName.c_str());
// 				ProfileWriter << "# 1D Profile of axon" << std::endl;
// 				ProfileWriter << "Depth [um]\tlength per bin [um]" << std::endl;
// 				for(int ii = 0; ii < zProfile.size(); ++ii)
// 					ProfileWriter << ii*50+25 << "\t" << zProfile[ii] << std::endl;
// 				ProfileWriter.close();
// 			}
// 			else
// 				std::cout << "Error! Z Profile is empty!" << std::endl;
			
			delete morphAnalyzer;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	if(argc == 6)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		const char * label = argv[3];
		double binSize = atof(argv[4]);
		int mode = atoi(argv[5]);
		
		std::string ofName(outputFilename);
		std::string ifName(inputFilename);
		
		Reader * hocReader = new Reader(inputFilename, outputFilename);
		if(ifName.find(".hoc") != std::string::npos)
			hocReader->readHocFile();
		else if(ifName.find(".am") != std::string::npos)
			hocReader->readSpatialGraphFile(0);
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			return 0;
		}
		
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		Analyzer * morphAnalyzer = new Analyzer(hocReader->getSpatialGraph(), checkPoint->getParameters());
		std::vector< double >* zProfile;
		
		if(mode == 1)
		{
			zProfile = morphAnalyzer->compute1DProfileGlobally(label, binSize);
			ofName += "_1DProfile_global.csv";
		}
		if(mode == 2)
		{
			zProfile = morphAnalyzer->compute1DProfileLocally(label, binSize);
			ofName += "_1DProfile_local.csv";

		}
		
		if(zProfile->size())
		{
			std::ofstream ProfileWriter;
			ProfileWriter.open(ofName.c_str());
			ProfileWriter << "# 1D Profile of " << label << std::endl;
			ProfileWriter << "Depth [um]\tlength per bin [um]" << std::endl;
			for(int ii = 0; ii < zProfile->size(); ++ii)
				ProfileWriter << ii*binSize << "\t" << (*zProfile)[ii] << std::endl;
			ProfileWriter.close();
			
 			//debug
 			std::cout << "1D Profile of " << label << std::endl;
 			std::cout << "Depth [um]\tlength per bin [um]" << std::endl;
 			for(int ii = 0; ii < zProfile->size(); ++ii)
				std::cout << ii*binSize << "\t" << (*zProfile)[ii] << std::endl;
		}
		else
			std::cout << "Error! Z Profile is empty!" << std::endl;
		
		delete morphAnalyzer, delete hocReader;
	}
	
// 	else
// 	{
// 		std::cout << "Usage: MorphAnalysis [Input filename] [Output filename] [Label] [Binsize in micron] [Mode: Global - 1 ; Local - 2]" << std::endl;
// 	}
	return 0;
}
