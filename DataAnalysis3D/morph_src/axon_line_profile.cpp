/****************************************************************************/
/*                                                                          */
/* Program:   MorphAnalyzer                                                 */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   Program for analysis of registered neuron morphologies with   */
/*            respect to columns, septa and layers in standardized barrel   */
/*            cortex. E.g., computes distribution of axon in different      */
/*            columns/layers, and computes z-profiles of axons taking local */
/*            orientation into account                                      */
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
// #include "morph_analyzer.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/amiraReader.h"
#include "../../common/basics.h"
#include "../../common/profile.h"

#define HORIZONTAL 0
#define CORONAL 1
#define SAGITTAL 2

// #define DEBUG

void computeAxonLineProfile(AmiraSpatialGraph * neuronMorphology, int orientation, double binSize, const char * outputFilename);

int main( int argc , char * argv[])
{
	if(argc == 4)
	{
		const char * inputFilename = argv[1];
		int orientation = atoi(argv[2]);
		int binSize = atof(argv[3]);
		
		std::string ifName(inputFilename);
		size_t suffix;
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(0);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .am files!" << std::endl;
			delete hocReader;
			return 0;
		}
		
#ifdef DEBUG
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		if(checkPoint->getParameters().axonFlag)
		{
			computeAxonLineProfile(hocReader->getSpatialGraph(), orientation, binSize, inputFilename);
		}
#endif
#ifndef DEBUG
		computeAxonLineProfile(hocReader->getSpatialGraph(), orientation, binSize, inputFilename);
#endif
	}
	
	else
	{
		std::cout << "Usage: MorphAnalysis [Input filename] [Output filename] [Label] [Binsize in micron] [Mode: Global - 1 ; Local - 2]" << std::endl;
	}
	return 0;
}

void computeAxonLineProfile(AmiraSpatialGraph* neuronMorphology, int orientation, double binSize, const char* outputFilename)
{
	std::vector< double > lengthBins;
	std::vector< std::pair< double, double> > bins;
	
	double bounds[6];
	neuronMorphology->getBoundingBox(bounds);
	double profileBounds[2];
	if(orientation == HORIZONTAL) //  M-L
	{
		profileBounds[0] = bounds[0];
		profileBounds[1] = bounds[1];
	}
	if(orientation == CORONAL) // D-V
	{
		profileBounds[0] = bounds[4];
		profileBounds[1] = bounds[5];
	}
	if(orientation == SAGITTAL) // A-P
	{
		profileBounds[0] = bounds[2];
		profileBounds[1] = bounds[3];
	}
	
	for(double boundary = profileBounds[0]; boundary <= profileBounds[1]; boundary += binSize)
	{
		double clippingBox[6];
		if(orientation == HORIZONTAL) //  M-L
		{
			clippingBox[0] = boundary;
			clippingBox[1] = boundary + binSize;
			clippingBox[2] = -1e5;
			clippingBox[3] = 1e5;
			clippingBox[4] = -1e5;
			clippingBox[5] = 1e5;
		}
		if(orientation == CORONAL) //  D-V
		{
			clippingBox[0] = -1e5;
			clippingBox[1] = 1e5;
			clippingBox[2] = -1e5;
			clippingBox[3] = 1e5;
			clippingBox[4] = boundary;
			clippingBox[5] = boundary + binSize;
		}
		if(orientation == SAGITTAL) //  A-P
		{
			clippingBox[0] = -1e5;
			clippingBox[1] = 1e5;
			clippingBox[2] = boundary;
			clippingBox[3] = boundary + binSize;
			clippingBox[4] = -1e5;
			clippingBox[5] = 1e5;
		}
		
		AmiraSpatialGraph * clippedSG = neuronMorphology->clipSpatialGraph(clippingBox);
		double tmpLength = 0;
		if(clippedSG)
		{
			for(int i = 0; i < clippedSG->edgesPointer()->size(); ++i)
			{
				Edge * tmpEdge = clippedSG->edgesPointer()->at(i);
				if(tmpEdge->label == Axon)
				{
					tmpLength += tmpEdge->segmentLength();
				}
			}
		}
		lengthBins.push_back(tmpLength);
		bins.push_back(std::pair< double, double >(boundary, boundary+binSize));
	}
#ifdef DEBUG
		std::cout << "\nDone processing axon length!" << std::endl;
#endif
	
	#ifndef DEBUG
	// write output files
	
	char * outNameTemplate = new char[256];
	if(orientation == HORIZONTAL)
	{
		sprintf(outNameTemplate, "%s_line_profile_M-L_binSize_%.0fmu.csv", outputFilename, binSize);
	}
	if(orientation == CORONAL)
	{
		sprintf(outNameTemplate, "%s_line_profile_D-V_binSize_%.0fmu.csv", outputFilename, binSize);
	}
	if(orientation == SAGITTAL)
	{
		sprintf(outNameTemplate, "%s_line_profile_A-P_binSize_%.0fmu.csv", outputFilename, binSize);
	}
	std::string axonOutName(outNameTemplate);
	std::ofstream AxonFile;
	AxonFile.open(axonOutName.c_str());
	AxonFile << "Bin begin (mu m)\tbin end(mu m)\tLength (mu m)" << std::endl;
	for(int i = 0; i < lengthBins.size(); ++i)
	{
		AxonFile << bins[i].first << "\t" << bins[i].second << "\t" << lengthBins[i] << std::endl;
	}
	AxonFile.close();
	#endif
}







