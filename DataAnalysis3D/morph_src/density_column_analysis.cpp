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
#include "../../common/barrel_field.h"

void computeDensityColumnParameters(const char * outputFilename, ImageDataPointerType density);

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * densityFilename = argv[1];
		
		ImageDataPointerType density;
		
		std::string densityStr(densityFilename);
		Reader * densityFileReader = new Reader(densityFilename, densityFilename);
		if(densityStr.find(".am") != std::string::npos)
		{
			density = densityFileReader->readScalarField();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			delete densityFileReader;
			return 0;
		}
		
		std::string outStr = densityStr.substr(0, densityStr.size()-3);
		computeDensityColumnParameters(outStr.c_str(), density);
		delete densityFileReader;
	}
	
	return 0;
}


// computes integrated density per column per layer
void computeDensityColumnParameters(const char * outputFilename, ImageDataPointerType density)
{
	// all density integrals supra/gran/infra
	// same for branch points
	// for all columns and septum
	BarrelField * SBF = new BarrelField(false);
	std::map< int, std::map< int, double > > columnIntegrals;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		// ignore Arcs 5 & 6 b/c we don't have standard barrels there
		if(*labelIt == C5 || *labelIt == C6
			|| *labelIt == D5 || *labelIt == D6
			|| *labelIt == E5 || *labelIt == E6)
		{
			continue;
		}
		std::map< int, double > layerIntegrals;
		for(int layer = 0; layer <= INFRA; ++layer)
		{
			layerIntegrals.insert(std::pair< int, double >(layer, 0));
		}
		columnIntegrals.insert(std::pair< int, std::map< int, double > >(*labelIt, layerIntegrals));
	}
	std::map< int, double > septumLayerIntegrals;
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		septumLayerIntegrals.insert(std::pair< int, double >(layer, 0));
	}
	columnIntegrals.insert(std::pair< int, std::map< int, double > >(Septum, septumLayerIntegrals));
	double outS1Integral = 0;
	
	double origin[3], voxelSpacing[3];
	int dimensions[3];
	density->GetOrigin(origin);
	density->GetSpacing(voxelSpacing);
	density->GetDimensions(dimensions);
	
	for(int kk = 0; kk < dimensions[0]; ++kk)
		for(int ll = 0; ll < dimensions[1]; ++ll)
			for(int mm = 0; mm < dimensions[2]; ++mm)
			{
				int klm[3];
				double xyz[3];
				klm[0] = kk, klm[1] = ll, klm[2] = mm;
				xyz[0] = origin[0] + kk*voxelSpacing[0];
				xyz[1] = origin[1] + ll*voxelSpacing[1];
				xyz[2] = origin[2] + mm*voxelSpacing[2];
				
				double * densVal = static_cast< double * >(density->GetScalarPointer(klm));
				
				if(SBF->isInsideS1(xyz))
				{
					int layer = SBF->laminarPosition(xyz);
					int column = SBF->insideColumn(xyz);
					if(columnIntegrals.find(column) != columnIntegrals.end())
					{
						columnIntegrals[column][layer] += *densVal;
					}
					else
					{
						std::cout << "Warning! Invalid column ID " << column << std::endl;
					}
				}
				else
				{
					outS1Integral += *densVal;
				}
			}
	
	std::string axonOutName(outputFilename);
	axonOutName += "_density_column_layer_integrals.csv";
	std::ofstream AxonFile;
	AxonFile.open(axonOutName.c_str());
	AxonFile << "Layer\t";
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		// ignore Arcs 5 & 6 for now
		// b/c we don't have standard barrels there
		if(*labelIt == C5 || *labelIt == C6
			|| *labelIt == D5 || *labelIt == D6
			|| *labelIt == E5 || *labelIt == E6)
		{
			continue;
		}
		AxonFile << SBF->int2Labels[*labelIt] << "\t";
	}
	AxonFile << "Septum\t";
	AxonFile << "Outside S1";
	AxonFile << std::endl;
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		if(!layer)
			AxonFile << "Other\t";
		if(layer == SUPRA)
			AxonFile << "Supra\t";
		if(layer == GRAN)
			AxonFile << "Gran\t";
		if(layer == INFRA)
			AxonFile << "Infra\t";
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			// ignore Arcs 5 & 6 for now
			// b/c we don't have standard barrels there
			if(*labelIt == C5 || *labelIt == C6
				|| *labelIt == D5 || *labelIt == D6
				|| *labelIt == E5 || *labelIt == E6)
			{
				continue;
			}
			
			AxonFile << columnIntegrals[*labelIt][layer] << "\t";
		}
		
		AxonFile << columnIntegrals[Septum][layer] << "\t";
		
		if(!layer)
		{
			AxonFile << outS1Integral;
		}
		AxonFile << std::endl;
	}
	AxonFile.close();
}
