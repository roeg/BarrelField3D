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
// #define USEREGISTRATION
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/inputparameters.h"
#include "cortexreconstruction.h"
#ifdef USEREGISTRATION
#include "landmarkregistration.h"
#endif

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * vesselFilename = argv[1];
		const char * outputFilename = argv[2];
		
		AmiraSpatialGraph * landmarkSG;
		
		std::string landmarkStr(vesselFilename);
		Reader * landmarkFileReader = new Reader(vesselFilename, outputFilename);
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
		recon3D->vesselReconstruction(outputFilename);
		
		// registration to Standard BF
		#ifdef USEREGISTRATION
		LandmarkRegistration * bfReg = new LandmarkRegistration(barrelField);
		bfReg->startRegistration();
		bfReg->writeTransform(outputFilename);
		TransformPointerType regTransform = bfReg->getTransform();
		#endif
		
		delete checkpoint;
		#ifdef USEREGISTRATION
		delete bfReg;
		#endif
		delete recon3D;
		delete landmarkFileReader, delete landmarkSG;
	}
	
	return 0;
}
