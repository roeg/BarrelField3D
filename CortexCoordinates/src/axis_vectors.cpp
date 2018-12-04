/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
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
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "../../common/basics.h"
#include "../../common/barrel_field.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		std::string ofName(inputFilename);
		ofName += "_column_orientations.csv";
		Reader * amiraReader = new Reader(inputFilename, ofName.c_str());
		amiraReader->readSpatialGraphFile(0);
		AmiraSpatialGraph * inputSG = amiraReader->getSpatialGraph();
		BarrelField *SBF = new BarrelField(0);
		
		std::map< int, std::vector< double > > colOrientations;
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			PolyDataPointerType column = PolyDataPointerType::New();
			if(inputSG->extractLandmark(ID, column))
			{
				double pt1[3], pt2[3], vec[3];
				double pCoords[3];
				int subID;
				double * weights0 = new double[column->GetCell(0)->GetNumberOfPoints()];
				double * weights1 = new double[column->GetCell(1)->GetNumberOfPoints()];
				column->GetCell(0)->GetParametricCenter(pCoords);
				column->GetCell(0)->EvaluateLocation(subID, pCoords, pt1, weights0);
				column->GetCell(1)->GetParametricCenter(pCoords);
				column->GetCell(1)->EvaluateLocation(subID, pCoords, pt2, weights1);
				delete [] weights0, delete [] weights1;
				
				vtkMath::Subtract(pt1, pt2, vec);
				vtkMath::Normalize(vec);
				std::vector< double > orientation;
				orientation.push_back(vec[0]);
				orientation.push_back(vec[1]);
				orientation.push_back(vec[2]);
				colOrientations.insert(std::pair< int, std::vector< double > >(ID, orientation));
			}
		}
		
		std::ofstream orientationFile;
		orientationFile.open(ofName.c_str());
		orientationFile << "Barrel\torientation x\torientation y\torientation z" << std::endl;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			orientationFile << SBF->int2Labels[ID];
			if(colOrientations.find(ID) != colOrientations.end())
			{
				for(int ii = 0; ii < 3; ++ii)
					orientationFile << "\t" << colOrientations[ID][ii];
			}
			orientationFile << std::endl;
		}
		orientationFile.close();
		
		delete SBF, delete inputSG, delete amiraReader;
	}
	
	return 0;
}