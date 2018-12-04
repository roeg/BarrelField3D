/****************************************************************************/
/*                                                                          */
/* Program:   3DCellCountAnalysis                                           */
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

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * landmarkFilename = argv[1];
		const char * outputFilename = argv[2];
		
		AmiraSpatialGraph * landmarkSG;
		
		std::string landmarkStr(landmarkFilename);
		Reader * landmarkFileReader = new Reader(landmarkFilename, outputFilename);
		if(landmarkStr.find(".am") != std::string::npos)
		{
			landmarkFileReader->readSpatialGraphFile(0);
			landmarkSG = landmarkFileReader->getSpatialGraph();
		}
		else if(landmarkStr.find(".hoc") != std::string::npos)
		{
			landmarkFileReader->readHocFile();
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
			delete landmarkFileReader;
			delete landmarkSG;
			return 0;
		}
		// otherwise, everything fine; start 3D reconstruction
		CortexReconstruction * recon3D = new CortexReconstruction(landmarkSG, checkpoint->getParameters());
		recon3D->startReconstruction(outputFilename);
		
		delete checkpoint;
		delete recon3D;
		delete landmarkFileReader, delete landmarkSG;
	}
	
	return 0;
}
