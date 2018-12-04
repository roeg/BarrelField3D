/****************************************************************************/
/*                                                                          */
/* Program:   NeuroCountAnalysis                                            */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            barrel columns reconstructed in 3D                            */
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
#include "landmarkregistration.h"
#include "landmarkanalyzer.h"

#define USECELLS
// #define USEREGISTRATION

int main( int argc , char * argv[])
{
	if(argc == 4 || argc == 5)
	{
		const char * cellFilename = argv[1];
		const char * landmarkFilename = argv[2];
		const char * outputFilename = argv[3];
		bool IN = 0;
		if(argc == 5) IN = atoi(argv[4]);
		
		PointsPointerType cellSomata;
		AmiraSpatialGraph * landmarkSG;
		
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
		
		std::string landmarkStr(landmarkFilename);
		Reader * landmarkFileReader = new Reader(landmarkFilename, outputFilename);
		if(landmarkStr.find(".am") != std::string::npos)
		{
			landmarkFileReader->readSpatialGraphFile(0);
			landmarkSG = landmarkFileReader->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			return 0;
		}
		
		// check if input is ok and set flags
		InputCheckpoint * checkpoint = new InputCheckpoint(landmarkSG);
		checkpoint->run();
		if(!checkpoint->getInputOK())
		{
			std::cout << "Error! Landmark file corrupt. Aborting..." << std::endl;
			delete checkpoint;
			#ifdef USECELLS
			delete cellFileReader;
			#endif
			delete landmarkFileReader;
			delete landmarkSG;
			return 0;
		}
		// otherwise, everything fine; start 3D reconstruction
		CortexReconstruction * recon3D = new CortexReconstruction(landmarkSG, checkpoint->getParameters());
		recon3D->startReconstruction(outputFilename);
		Surface * pia = recon3D->getPiaSurface();
		Surface * WM = recon3D->getWMSurface();
		std::map< int, Column * > barrelField = recon3D->getBarrelField();
		
		// registration to Standard BF
		#ifdef USEREGISTRATION
		LandmarkRegistration * bfReg = new LandmarkRegistration(barrelField);
		bfReg->startRegistration();
		bfReg->writeTransform(outputFilename);
		TransformPointerType regTransform = bfReg->getTransform();
		#endif
		
		// and finally, count cells...
		#ifdef USECELLS
		LandmarkAnalyzer * cellAnalyzer = new LandmarkAnalyzer;
		cellAnalyzer->setBarrelField(barrelField);
		cellAnalyzer->setPiaSurface(pia);
		cellAnalyzer->setWMSurface(WM);
		cellAnalyzer->setLandmarkSet(cellSomata);
		if(IN)
		{
			cellAnalyzer->computeCorrectedColumnProfiles(outputFilename);
		}
		else
		{
			cellAnalyzer->computeColumnProfiles(outputFilename);
			cellAnalyzer->computeSeparateColumnProfiles(outputFilename);
		}
		#endif
		
		delete checkpoint;
		#ifdef USECELLS
		delete cellAnalyzer;
		delete cellFileReader;
		#endif
		#ifdef USEREGISTRATION
		delete bfReg;
		#endif
		delete recon3D,delete pia, delete WM;
		delete landmarkFileReader, delete landmarkSG;
	}
	
	if(argc == 3)
	{
		const char * cellFilename = argv[1];
		const char * outputFilename = argv[2];
		
		PointsPointerType cellSomata;
		
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
		
		// and finally, count cells...
		#ifdef USECELLS
		LandmarkAnalyzer * cellAnalyzer = new LandmarkAnalyzer;
		cellAnalyzer->standardBFAnalysis();
		cellAnalyzer->setLandmarkSet(cellSomata);
		cellAnalyzer->computeColumnProfiles(outputFilename);
		cellAnalyzer->computeSeparateColumnProfiles(outputFilename);
		#endif
		
		#ifdef USECELLS
		delete cellAnalyzer;
		delete cellFileReader;
		#endif
	}
	
	return 0;
}
