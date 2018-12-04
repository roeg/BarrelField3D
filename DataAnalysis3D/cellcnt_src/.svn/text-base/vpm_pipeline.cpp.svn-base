/****************************************************************************/
/*                                                                          */
/* Program:   3DCellCountAnalysis                                           */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            barreloids in VPM reconstructed in 3D                         */
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
#include "vpmreconstruction.h"
#include "barreloidsurfacereconstruct.h"
#include "landmarkanalyzer.h"

#define USECELLS

int main( int argc , char * argv[])
{
	if(argc == 4)
	{
		const char * cellFilename = argv[1];
		const char * landmarkFilename = argv[2];
		const char * outputFilename = argv[3];
		
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
		
// 		// test barreloid surface reconstruction:
// 		PolyDataPointerType C1barreloid = PolyDataPointerType::New();
// 		if(landmarkSG->extractLandmark(E3, C1barreloid))
// 		{
// 			BarreloidSurfaceReconstruct * barreloidReconstruct = new BarreloidSurfaceReconstruct(checkpoint->getParameters(), C1barreloid);
// 			PolyDataPointerType C1barreloidSurface = barreloidReconstruct->surfaceReconstruction(E3);
// 			Reader * surfWriter = new Reader(outputFilename, outputFilename);
// 			surfWriter->writeAmiraSurfaceFile(C1barreloidSurface);
// 			delete surfWriter;
// 			delete barreloidReconstruct;
// 		}
// 		// test ellipsoid approximation:
// 		PolyDataPointerType C1barreloid = PolyDataPointerType::New();
// 		if(landmarkSG->extractLandmark(C2, C1barreloid))
// 		{
// 			BarreloidSurfaceReconstruct * barreloidReconstruct = new BarreloidSurfaceReconstruct(checkpoint->getParameters(), C1barreloid);
// 			PolyDataPointerType barreloidSurface = barreloidReconstruct->surfaceReconstruction(C2);
// 			barreloidReconstruct->getPrincipalAxes();
// 			delete barreloidReconstruct;
// 		}
		
		VPMReconstruction * recon3D = new VPMReconstruction(landmarkSG, checkpoint->getParameters());
		recon3D->startReconstruction(outputFilename);
		std::map< int, ClosedSurface * > barreloidField = recon3D->getBarreloids();
		
		// and finally, count cells...
		#ifdef USECELLS
		LandmarkAnalyzer * cellAnalyzer = new LandmarkAnalyzer;
		cellAnalyzer->setBarreloidField(barreloidField);
		cellAnalyzer->setLandmarkSet(cellSomata);
		cellAnalyzer->countCellsInBarreloids(outputFilename);
		#endif
		
		delete checkpoint;
		#ifdef USECELLS
		delete cellAnalyzer;
		delete cellFileReader;
		#endif
// 		delete recon3D;
		delete landmarkFileReader, delete landmarkSG;
	}
	
	return 0;
}
