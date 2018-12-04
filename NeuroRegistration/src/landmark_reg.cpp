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
	if(argc == 4 || argc == 5)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		int mode = atoi(argv[3]);
// 		const char * refBarrel = (mode == 2) ? argv[4] : NULL;
		const char * refBarrel = argv[4];
		
		std::string inStr(inputFilename);
		Reader * hocReader = new Reader(inputFilename, outputFilename);
		if(inStr.find(".hoc") != std::string::npos)
			hocReader->readHocFile();
		else if(inStr.find(".am") != std::string::npos)
			hocReader->readSpatialGraphFile(0);
		else
		{
			std::cout << "Error! Can only read .hoc or .am files!" << std::endl;
			return 0;
		}
		Registration * morphRegistration = new Registration(hocReader->getSpatialGraph(), 1);
		morphRegistration->barrelFieldRegistration(outputFilename, mode, refBarrel);
		
		hocReader->setSpatialGraph(morphRegistration->getSpatialGraph());
		hocReader->writeSpatialGraphFile();
		delete morphRegistration, delete hocReader;
	}
	
	return 0;
}
