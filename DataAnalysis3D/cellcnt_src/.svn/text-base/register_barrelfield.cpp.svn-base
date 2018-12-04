/****************************************************************************/
/*                                                                          */
/* Program:   RegisterBarrelfield                                           */
/*                                                                          */
/* File:      register_barrelfield.cpp                                      */
/*                                                                          */
/* Purpose:   tool for registration of reconstructed barrel field to        */
/*            standard barrel field. input has to be in the form of a       */
/*            barrel field consisting only of barrel top/bottom contours    */
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
#include "landmarkregistration.h"

int main( int argc , char * argv[])
{
	if(argc == 3 || argc == 4)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		bool scaling = 0;
		if(argc == 4) scaling = atoi(argv[3]);
		
		AmiraSpatialGraph * landmarkSG;
		
		std::string landmarkStr(inputFilename);
		Reader * landmarkFileReader = new Reader(inputFilename, outputFilename);
		if(landmarkStr.find(".am") != std::string::npos)
		{
			landmarkFileReader->readSpatialGraphFile(0);
			landmarkSG = landmarkFileReader->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! Barrel field file has to be Amira '.am' file!" << std::endl;
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
		
		// registration to Standard BF
		LandmarkRegistration * bfReg = new LandmarkRegistration(landmarkSG);
		if(!scaling)
		{
			bfReg->startRegistration();
// 			bfReg->startColumnRegistration();
		}
		else
		{
			bfReg->startRegAnisotropicScale();
		}
		bfReg->writeTransform(outputFilename);
		
		delete checkpoint;
		delete bfReg;
		delete landmarkFileReader, delete landmarkSG;
	}
	
	return 0;
}
