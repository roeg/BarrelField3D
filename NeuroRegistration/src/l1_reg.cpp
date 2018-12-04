/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      l1_reg.cpp                                                    */
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

int main( int argc , char * argv[])
{
	if(argc >= 6 && argc <= 8)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		int neuronAxis;
		int cellType;
		bool goodApicalDend;
		bool piaCorrection;
		//hard-core mode for e.g. VPM axons
		bool onlyAxon;
		std::string manualHomeBarrel("0");
		if(argc == 6)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			piaCorrection = 0;
			onlyAxon = 0;
		}
		else if(argc == 7)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			piaCorrection = atoi(argv[6]);
			onlyAxon = 0;
		}
		else if(argc == 8)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			piaCorrection = atoi(argv[6]);
			onlyAxon = atoi(argv[7]);
		}
		else if(argc == 9)
		{
			neuronAxis = atoi(argv[3]);
			cellType = atoi(argv[4]);
			goodApicalDend = atoi(argv[5]);
			piaCorrection = atoi(argv[6]);
			onlyAxon = atoi(argv[7]);
			manualHomeBarrel = std::string(argv[8]);
		}
		
		Reader * hocReader = new Reader(inputFilename, outputFilename);
// 		hocReader->readSpatialGraphFile(0);
		hocReader->readHocFile();
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 1);
		morphRegistration->setL1flag(1);
		int success = morphRegistration->startNeuronRegistration(outputFilename, neuronAxis, cellType, goodApicalDend,
																 piaCorrection, onlyAxon, manualHomeBarrel);
		if(success)
		{
			hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
			hocReader->writeSpatialGraphFile();
			hocReader->writeSeparateHocFiles();
// 			hocReader->writeHocFile();
		}
		else
			std::cout << "Could not register neuron morphology!" << std::endl;
		delete morphRegistration, delete hocReader;
	}
	
	return 0;
}