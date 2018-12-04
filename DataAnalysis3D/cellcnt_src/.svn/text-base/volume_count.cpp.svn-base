/****************************************************************************/
/*                                                                          */
/* Program:   VolumeCountAnalysis                                           */
/*                                                                          */
/* File:      volume_count.cpp                                              */
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
#include "landmarkanalyzer.h"

#define USECELLS
#define USEREGISTRATION

int main( int argc , char * argv[])
{
	if(argc == 7 || argc == 8)
	{
		const char * cellFilename = argv[1];
		const char * barrelColumnsFilename = argv[2];
		const char * S1Filename = argv[3];
		const char * piaFilename = argv[4];
		const char * wmFilename = argv[5];
		const char * outputFilename = argv[6];
		const char * layerColumnsFilename = NULL;
		if(argc == 8)
			layerColumnsFilename = argv[7];
		
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
		
		// create barrel field
		// and read surfaces
		CortexReconstruction * bfRecon = new CortexReconstruction(barrelColumns);
		bfRecon->convertSGToBarrelField();
		std::map< int, Column * > barrelField = bfRecon->getBarrelField();
		Reader * S1Reader = new Reader(S1Filename, S1Filename);
		Reader * piaReader = new Reader(piaFilename, piaFilename);
		Reader * wmReader = new Reader(wmFilename, wmFilename);
		Surface * S1ConvexHull = new Surface(S1Reader->readAmiraSurfaceFile());
		Surface * pia = new Surface(piaReader->readAmiraSurfaceFile());
		Surface * WM = new Surface(wmReader->readAmiraSurfaceFile());
		
		if(argc == 8)
		{
			std::string landmarkStr2(layerColumnsFilename);
			Reader * landmarkFileReader2 = new Reader(layerColumnsFilename, outputFilename);
			if(landmarkStr2.find(".am") != std::string::npos)
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
			checkpoint2->checkBarrelField();
			if(!checkpoint2->getInputOK())
			{
				std::cout << "Error! Landmark file corrupt. Aborting..." << std::endl;
				delete checkpoint2;
				delete landmarkFileReader2;
				delete layerColumns;
				return 0;
			}
		}
		
		// and finally, count cells...
		#ifdef USECELLS
		LandmarkAnalyzer * cellAnalyzer = new LandmarkAnalyzer;
		cellAnalyzer->setBarrelField(barrelField);
		cellAnalyzer->setPiaSurface(pia);
		cellAnalyzer->setWMSurface(WM);
		cellAnalyzer->setLandmarkSet(cellSomata);
		if(argc == 7)
		{
			cellAnalyzer->countCellsInSeptum(S1ConvexHull, outputFilename);
		}
		if(argc == 8)
		{
			CortexReconstruction * bfRecon2 = new CortexReconstruction(layerColumns);
			bfRecon2->convertSGToBarrelField();
			std::map< int, Column * > barrelFieldLayer = bfRecon2->getBarrelField();
			cellAnalyzer->setBarrelFieldLayer(barrelFieldLayer);
			cellAnalyzer->countInhCellsInSeptumInLayer(S1ConvexHull, outputFilename);
		}
		#endif
		
		delete checkpoint;
		#ifdef USECELLS
		delete cellAnalyzer;
		delete cellFileReader;
		#endif
		delete bfRecon, delete pia, delete WM, delete S1ConvexHull;
		delete S1Reader, delete piaReader, delete wmReader;
		delete landmarkFileReader, delete barrelColumns;
	}
	
	return 0;
}
