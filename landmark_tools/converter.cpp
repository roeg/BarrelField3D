/****************************************************************************/
/*                                                                          */
/* Program:   AverageBarrelField                                            */
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
#include "../common/typedefs.h"
#include "../common/amiraReader.h"
#include "../common/basics.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		Reader * converter = new Reader(inputFilename, inputFilename);
		converter->readLandmarkFile();
		converter->writeSpatialGraphFile();
		delete converter;
	}
	else
		std::cout << "Usage: ./Converter [input filename]" << std::endl;
	return 0;
}