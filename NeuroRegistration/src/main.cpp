/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
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
#include "morph_reg.h"

// #define REG_ACCURACY
#ifdef REG_ACCURACY
double var_alpha;
double var_gamma;
#endif

int main( int argc , char * argv[])
{
// 	if(argc == 3)
// 	{
// 		const char * inputFilename = argv[1];
// 		const char * outputFilename = argv[2];
// 		Reader * hocReader = new Reader(inputFilename, outputFilename);
// 		hocReader->readHocFile();
// 		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph());
// 		morphRegistration->neuronMorphologyZAxis();
// 		hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
// 		hocReader->writeSpatialGraphFile();
// 		delete hocReader;
// 	}
	
	#ifdef REG_ACCURACY
	if(argc == 9)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		int neuronAxis;
		int cellType;
		int goodApicalDend;
		bool noBarrels;
		//hard-core mode for e.g. VPM axons
		bool onlyAxon;
		
		neuronAxis = atoi(argv[3]);
		cellType = atoi(argv[4]);
		goodApicalDend = atoi(argv[5]);
		noBarrels = atoi(argv[6]);
		onlyAxon = 0;
		
		var_alpha = atof(argv[7]);
		var_gamma = atof(argv[8]);
		
		Reader * hocReader = new Reader(inputFilename, outputFilename);
// 		hocReader->readSpatialGraphFile(0);
		hocReader->readHocFile();
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 1);
		int success = morphRegistration->startNeuronRegistration(outputFilename, neuronAxis, cellType, goodApicalDend, noBarrels, onlyAxon);
		if(success)
		{
			hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
			hocReader->writeSpatialGraphFile();
			hocReader->writeHocFile();
		}
		else
			std::cout << "Could not register neuron morphology!" << std::endl;
		delete morphRegistration, delete hocReader;
	}
	#endif
	
	#ifndef REG_ACCURACY
	if(argc >= 7 && argc <= 10)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		int neuronAxis;
		int cellType;
		int goodApicalDend;
		int landmarkMode;
		bool piaCorrection;
		//hard-core mode for e.g. VPM axons
		bool onlyAxon;
		std::string manualHomeBarrel("0");
		if(argc == 7)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			landmarkMode = atoi(argv[6]);
			piaCorrection = 0;
			onlyAxon = 0;
		}
		else if(argc == 8)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			landmarkMode = atoi(argv[6]);
			piaCorrection = atoi(argv[7]);
			onlyAxon = 0;
		}
		else if(argc == 9)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			landmarkMode = atoi(argv[6]);
			piaCorrection = atoi(argv[7]);
			onlyAxon = atoi(argv[8]);
		}
		else if(argc == 10)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			landmarkMode = atoi(argv[6]);
			piaCorrection = atoi(argv[7]);
			onlyAxon = atoi(argv[8]);
			manualHomeBarrel = std::string(argv[9]);
		}
		
		Reader * hocReader = new Reader(inputFilename, outputFilename);
// 		hocReader->readSpatialGraphFile(0);
		hocReader->readHocFile();
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 1, landmarkMode);
		int success = morphRegistration->startNeuronRegistration(outputFilename, neuronAxis, cellType, goodApicalDend,
																 piaCorrection, onlyAxon, manualHomeBarrel);
		if(success)
		{
			hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
			hocReader->writeSpatialGraphFile();
// 			hocReader->writeHocFile();
			hocReader->writeSeparateHocFiles();
		}
		else
			std::cout << "Could not register neuron morphology!" << std::endl;
		delete morphRegistration, delete hocReader;
	}
	if(argc == 4 || argc == 5)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		const char * newColumnLabel = argv[3];
		//hard-core mode for e.g. VPM axons
		bool onlyAxon = 0;
		if(argc == 5)
			onlyAxon = atoi(argv[4]);
		
		Reader * hocReader = new Reader(inputFilename, outputFilename);
		hocReader->readHocFile();
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 0, 0);
		int success = morphRegistration->registerToDifferentColumn(outputFilename, newColumnLabel, onlyAxon);
		if(success)
		{
			hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
			hocReader->writeHocFile();
		}
		else
			std::cout << "Could not register neuron morphology to new column " << newColumnLabel << std::endl;
		delete morphRegistration, delete hocReader;
	}
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		Reader * hocReader = new Reader(inputFilename, outputFilename);
		std::string inputNameStr(inputFilename);
		if(inputNameStr.find(".hoc") == inputNameStr.size()-4)
		{
			hocReader->readHocFile();
		}
		else if(inputNameStr.find(".am") == inputNameStr.size()-3)
		{
			hocReader->readSpatialGraphFile(0);
		}
		else
		{
			std::cout << "Error: file has to be AmiraMesh or hoc format!" << std::endl;
			return 0;
		}
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 1, 0);
		morphRegistration->reconstructionError(outputFilename);
		hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
		hocReader->writeSpatialGraphFile();
		delete morphRegistration, delete hocReader;
	}
	#endif
	
	return 0;
}
