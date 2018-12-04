/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   program for processing of contour data obtained from the      */
/*            SurfaceExtraction image processing pipeline                   */
/*            -Surfaces are calculated for Pia and White Matter from the    */
/*            raw contour data                                              */
/*            -blood vessels are connected in 3D by a greedy algorithm      */
/*            -barrel contours are smoothed in the stack-z direction        */
/*            -for each barrel, a new z-axis is calculated based on the     */
/*            distance to Pia, orientation of that axis w.r.t. Pia at the   */
/*            intersection point and orientation of blood vessels in the    */
/*            neighborhood of the barrel                                    */
/*            -based on these axes, barrel anatomy is calculated            */
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
#include "geometry.h"
/*old Values*/
float XYSAMPLING = 0;
float ZSAMPLING = 0;

/*float XYSAMPLING = 0.311;
float ZSAMPLING = 0.6;

float XYSAMPLING = 0.395879;
float ZSAMPLING = 1;*/

float averageSomaRadius = 6.0;
float zScale = 1.5;
unsigned long BINSIZE = 0;
unsigned long BOXSIZE = 150;

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		std::string ofName(inputFilename);
		ofName += "_S1_volume";
		Reader * amiraReader = new Reader(inputFilename, ofName.c_str());
		amiraReader->readSpatialGraphFile(0);
		Geometry * spatialGraphProcessor = new Geometry(amiraReader->getSpatialGraph());
		spatialGraphProcessor->computeTotalVolumes(inputFilename);
		amiraReader->setSpatialGraph(spatialGraphProcessor->getSpatialGraph());
		amiraReader->writeSpatialGraphFile();
		delete spatialGraphProcessor;
		delete amiraReader;
	}
// 	if(argc == 3)
// 	{
// 		const char * inputFilename = argv[1];
// 		const char * outputFilename = argv[2];
// 		
// 		Reader * amiraReader = new Reader(inputFilename, outputFilename);
// 		bool applyTransform = 1;
// 		amiraReader->readSpatialGraphFile(applyTransform);
// 		amiraReader->writeSpatialGraphFile();
// 		delete amiraReader;
// 	}
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		Reader * amiraReader = new Reader(inputFilename, outputFilename);
		amiraReader->readHocFile();
		
		Geometry * spatialGraphProcessor = new Geometry(amiraReader->getSpatialGraph());
		spatialGraphProcessor->reconstructionError(outputFilename);
		
		amiraReader->setSpatialGraph(spatialGraphProcessor->getSpatialGraph2());
		amiraReader->writeSpatialGraphFile();
		delete amiraReader;
	}
	if(argc == 4)
	{
		const char * spatialGraphFilename = argv[1];
		const char * surfaceFilename = argv[2];
		const char * outputFilename = argv[3];
		
		float binsize = 5/*atof(argv[4])*/;
		
		Geometry * spatialGraphProcessor = new Geometry();
		spatialGraphProcessor->computeVesselToNormalAngleHistogram(spatialGraphFilename, surfaceFilename, outputFilename, binsize);
		
		Reader * amiraReader = new Reader(spatialGraphFilename, outputFilename);
		amiraReader->setSpatialGraph(spatialGraphProcessor->getSpatialGraph());
		amiraReader->writeSpatialGraphFileFromEdges();
		
		delete spatialGraphProcessor;
	}
	if(argc == 5)
	{
		const char * spatialGraphFilename1 = argv[1];
		const char * spatialGraphFilename2 = argv[2];
		const char * outputFilename = argv[3];
		double alpha = atof(argv[4]);
		alpha = 0.5;
		
		bool applyTransform = 0;
		
		Reader * piaReader = new Reader(spatialGraphFilename1, outputFilename);
		piaReader->readSpatialGraphFile(applyTransform);
		Reader * barrelReader = new Reader(spatialGraphFilename2, outputFilename);
		barrelReader->readSpatialGraphFile(applyTransform);
		
		AmiraSpatialGraph * sg = piaReader->getSpatialGraph();
		sg->mergeSpatialGraph(barrelReader->getSpatialGraph());
		delete barrelReader;
		
		Geometry * spatialGraphProcessor = new Geometry(sg);
		spatialGraphProcessor->computeBarrelSurfaces(outputFilename, alpha);
		
		piaReader->setSpatialGraph(spatialGraphProcessor->getSpatialGraph2());
		piaReader->writeSpatialGraphFileFromEdges();
		
		delete spatialGraphProcessor;
		delete piaReader;
	}
	
	return 0;
}