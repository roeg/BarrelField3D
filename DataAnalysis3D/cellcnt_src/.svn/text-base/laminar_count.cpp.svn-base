/****************************************************************************/
/*                                                                          */
/* Program:   LaminarCountAnalysis                                          */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            one layer (barrel columns reconstructed in 3D are passed      */
/*            as an argument)                                               */
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
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/inputparameters.h"
#include "cortexreconstruction.h"
#include "landmarkanalyzer.h"

#define USECELLS
// #define USEREGISTRATION

int main( int argc , char * argv[])
{
	if(argc == 5 || argc == 6 || argc == 9)
	{
		const char * cellFilename = argv[1];
		const char * barrelColumnsFilename = argv[2];
		const char * outputFilename = argv[3];
		const char * layerColumnsFilename = NULL;
		const char * supraFilename = NULL;
		const char * granFilename = NULL;
		const char * infraFilename = NULL;
		// mode: classic = 0, MO layers = 1
		int mode = atoi(argv[4]);
		bool IN = 0;
		if(argc == 6)
		{
			IN = 1;
			layerColumnsFilename = argv[5];
		}
		if(argc == 9)
		{
			supraFilename = argv[5];
			granFilename = argv[6];
			infraFilename = argv[7];
			IN = atoi(argv[8]);
		}
		
		PointsPointerType cellSomata;
		AmiraSpatialGraph * barrelColumns, * layerColumns;
		
		#ifdef USECELLS
		std::string cellStr(cellFilename);
		Reader * cellFileReader = new Reader(cellFilename, outputFilename);
		if(cellStr.find(".landmarkAscii") != std::string::npos)
			cellSomata = cellFileReader->readLandmarkFile();
		else
		{
			std::cout << "Error! Cell landmark file has to be Amira '.landmarkAscii' file!" << std::endl;
			return 0;
		}
		#endif
		
		std::string landmarkStr(barrelColumnsFilename);
		Reader * landmarkFileReader = new Reader(barrelColumnsFilename, outputFilename);
		if(landmarkStr.find(".am") != std::string::npos)
		{
			landmarkFileReader->readSpatialGraphFile(0);
			barrelColumns = landmarkFileReader->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			return 0;
		}
		
		// check if input is ok and set flags
		InputCheckpoint * checkpoint = new InputCheckpoint(barrelColumns);
		checkpoint->checkBarrelField();
		if(!checkpoint->getInputOK())
		{
			std::cout << "Error! Landmark file corrupt. Aborting..." << std::endl;
			delete checkpoint;
			#ifdef USECELLS
			delete cellFileReader;
			#endif
			delete landmarkFileReader;
			delete barrelColumns;
			return 0;
		}
		
		if(layerColumnsFilename)
		{
			std::string landmarkStr2(layerColumnsFilename);
			Reader * landmarkFileReader2 = new Reader(layerColumnsFilename, outputFilename);
			if(landmarkStr.find(".am") != std::string::npos)
			{
				landmarkFileReader2->readSpatialGraphFile(0);
				layerColumns = landmarkFileReader2->getSpatialGraph();
			}
			else
			{
				std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
				return 0;
			}
			
			// check if input is ok and set flags
			InputCheckpoint * checkpoint2 = new InputCheckpoint(layerColumns);
			checkpoint->checkBarrelField();
			if(!checkpoint->getInputOK())
			{
				std::cout << "Error! Landmark file corrupt. Aborting..." << std::endl;
				delete checkpoint;
				#ifdef USECELLS
				delete cellFileReader;
				#endif
				delete landmarkFileReader2;
				delete layerColumns;
				return 0;
			}
			delete checkpoint2;
		}
		
		// create barrel field
		CortexReconstruction * bfRecon = new CortexReconstruction(barrelColumns);
		bfRecon->convertSGToBarrelField();
		std::map< int, Column * > barrelField = bfRecon->getBarrelField();
		std::map< int, Column * > barrelFieldLayer;
		
		// and finally, count cells...
		#ifdef USECELLS
		LandmarkAnalyzer * cellAnalyzer = new LandmarkAnalyzer;
		cellAnalyzer->setLandmarkSet(cellSomata);
		cellAnalyzer->setBarrelField(barrelField);
		if(layerColumnsFilename)
		{
			// create barrel field
			CortexReconstruction * bfReconLayer = new CortexReconstruction(layerColumns);
			bfReconLayer->convertSGToBarrelField();
			barrelFieldLayer = bfReconLayer->getBarrelField();
			cellAnalyzer->setBarrelFieldLayer(barrelFieldLayer);
		}
		if(mode == 0)
		{
			if(!IN)
			{
				cellAnalyzer->computeColumnProfiles(outputFilename);
			}
			if(IN)
			{
				cellAnalyzer->computeCorrectedColumnProfilesInLayer(outputFilename);
			}
		}
		if(mode == 1 && argc == 9)
		{
			AmiraSpatialGraph * supraColumns;
			std::string supraStr(supraFilename);
			Reader * supraReader = new Reader(supraFilename, outputFilename);
			if(supraStr.find(".am") != std::string::npos)
			{
				supraReader->readSpatialGraphFile(0);
				supraColumns = supraReader->getSpatialGraph();
			}
			else
			{
				std::cout << "Error! Supra landmark file has to be Amira '.am' file!" << std::endl;
				return 0;
			}
			
			// check if input is ok and set flags
			InputCheckpoint * supraCheckpoint = new InputCheckpoint(supraColumns);
			supraCheckpoint->checkBarrelField();
			if(!supraCheckpoint->getInputOK())
			{
				std::cout << "Error! Supra landmark file corrupt. Aborting..." << std::endl;
				delete supraCheckpoint;
				delete supraReader;
				delete supraColumns;
				return 0;
			}
			delete supraCheckpoint;
			CortexReconstruction * supraRecon = new CortexReconstruction(supraColumns);
			supraRecon->convertSGToBarrelField();
			cellAnalyzer->setSupragranularLayer(supraRecon->getBarrelField());
			
			AmiraSpatialGraph * granColumns;
			std::string granStr(granFilename);
			Reader * granReader = new Reader(granFilename, outputFilename);
			if(granStr.find(".am") != std::string::npos)
			{
				granReader->readSpatialGraphFile(0);
				granColumns = granReader->getSpatialGraph();
			}
			else
			{
				std::cout << "Error! Granular landmark file has to be Amira '.am' file!" << std::endl;
				return 0;
			}
			
			// check if input is ok and set flags
			InputCheckpoint * granCheckpoint = new InputCheckpoint(granColumns);
			granCheckpoint->checkBarrelField();
			if(!granCheckpoint->getInputOK())
			{
				std::cout << "Error! Granular landmark file corrupt. Aborting..." << std::endl;
				delete granCheckpoint;
				delete granReader;
				delete granColumns;
				return 0;
			}
			delete granCheckpoint;
			CortexReconstruction * granRecon = new CortexReconstruction(granColumns);
			granRecon->convertSGToBarrelField();
			cellAnalyzer->setGranularLayer(granRecon->getBarrelField());
			
			AmiraSpatialGraph * infraColumns;
			std::string infraStr(infraFilename);
			Reader * infraReader = new Reader(infraFilename, outputFilename);
			if(infraStr.find(".am") != std::string::npos)
			{
				infraReader->readSpatialGraphFile(0);
				infraColumns = infraReader->getSpatialGraph();
			}
			else
			{
				std::cout << "Error! Granular landmark file has to be Amira '.am' file!" << std::endl;
				return 0;
			}
			
			// check if input is ok and set flags
			InputCheckpoint * infraCheckpoint = new InputCheckpoint(infraColumns);
			infraCheckpoint->checkBarrelField();
			if(!infraCheckpoint->getInputOK())
			{
				std::cout << "Error! Infra landmark file corrupt. Aborting..." << std::endl;
				delete infraCheckpoint;
				delete infraReader;
				delete infraColumns;
				return 0;
			}
			delete infraCheckpoint;
			CortexReconstruction * infraRecon = new CortexReconstruction(infraColumns);
			infraRecon->convertSGToBarrelField();
			cellAnalyzer->setInfragranularLayer(infraRecon->getBarrelField());
			
			cellAnalyzer->countLaminarNeuronNumbers(outputFilename, IN);
			
			delete supraReader;
			delete supraColumns;
			delete supraRecon;
			delete granReader;
			delete granColumns;
			delete granRecon;
			delete infraReader;
			delete infraColumns;
			delete infraRecon;
		}
		#endif
		
		delete checkpoint;
		#ifdef USECELLS
		delete cellAnalyzer;
		delete cellFileReader;
		#endif
		delete bfRecon;
		delete landmarkFileReader, delete barrelColumns;
	}
	
	return 0;
}
