/****************************************************************************/
/*                                                                          */
/* Program:   DensityCorrection                                             */
/*                                                                          */
/* File:      density_correction.cpp                                        */
/*                                                                          */
/* Purpose:   correct 3D IN density                                         */
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
#include "landmarkanalyzer.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * densityFilename = argv[1];
		
		ImageDataPointerType totalDensity;
		
		std::string densityStr(densityFilename);
		Reader * densityFileReader = new Reader(densityFilename, densityFilename);
		if(densityStr.find(".am") != std::string::npos)
		{
			totalDensity = densityFileReader->readScalarField();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			delete densityFileReader;
			return 0;
		}
		
		std::string outStr = densityStr.substr(0, densityStr.size()-3);
		LandmarkAnalyzer * densityAnalyzer = new LandmarkAnalyzer;
// 		densityAnalyzer->cropDensityLayers(totalDensity, outStr.c_str());
		densityAnalyzer->cropINDensity(totalDensity, outStr.c_str());
// 		densityAnalyzer->cropDensityCRowArc2(totalDensity, outStr.c_str());
		
		delete densityAnalyzer;
		delete densityFileReader;
	}
	
	return 0;
}
