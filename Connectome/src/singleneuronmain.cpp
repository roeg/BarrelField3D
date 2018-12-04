/****************************************************************************/
/*                                                                          */
/* Program:   SingleNeuronSynapseMapper                                     */
/*                                                                          */
/* File:      singleneuronmain.cpp                                          */
/*                                                                          */
/* Purpose:   single neuron synapse mapping                                 */
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
#include "../../common/spatialgraphset.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		SpatialGraphSet * sgs = new SpatialGraphSet(inputFilename);
		delete sgs;
	}
	
	return 0;
}
