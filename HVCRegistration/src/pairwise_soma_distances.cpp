/****************************************************************************/
/*                                                                          */
/* Program:   3DCellCountAnalysis                                           */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            barreloids in VPM reconstructed in 3D                         */
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

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * sgListFilename = argv[1];
		const char * outputFilename = argv[2];
		
		std::vector< char *  > individualSpatialGraphNames;
		std::vector< double * > somaLocations;
		std::ifstream inputStream(sgListFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading SpatialGraph files from " << sgListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					individualSpatialGraphNames.push_back(tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * sgReader = new Reader(tmpChar, outputFilename);
					if(sgStr.find(".am") != std::string::npos)
					{
						std::cout << "Loading SpatialGraph file " << tmpChar << std::endl;
						sgReader->readSpatialGraphFile(1);
						AmiraSpatialGraph * thisSG = sgReader->getSpatialGraph();
						PolyDataPointerType thisSoma = PolyDataPointerType::New();
						if(thisSG->extractLandmark(Soma, thisSoma))
						{
							int subId;
							double * somaLoc = new double[3];
							double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
							thisSoma->GetCell(0)->GetParametricCenter(pCenter);
							thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
							somaLocations.push_back(somaLoc);
							delete [] weights;
						}
						else
						{
							std::cout << "Error! SpatialGraph file " << tmpChar << " does not contain a soma!" << std::endl;
							delete sgReader;
							return 0;
						}
					}
					else
					{
						std::cout << "Error! SpatialGraph file " << tmpChar << " has to be Amira '.am' ascii file!" << std::endl;
						delete sgReader;
						return 0;
					}
					
					delete sgReader;
				}
		}
		
		std::vector< std::map< int, double > > pairWiseDistances;
		for(int i = 0; i < somaLocations.size(); ++i)
		{
			std::map< int, double > tmpDistances;
			for(int j = i+1; j < somaLocations.size(); ++j)
			{
				double * soma1 = somaLocations[i];
				double * soma2 = somaLocations[j];
				double dist = sqrt(vtkMath::Distance2BetweenPoints(soma1, soma2));
				tmpDistances[j] = dist;
			}
			pairWiseDistances.push_back(tmpDistances);
		}
		
		std::string outName(outputFilename);
		if(outName.find(".csv") == std::string::npos)
		{
			outName += ".csv";
		}
		
		std::ofstream outFile(outName.c_str());
		for(int i = 0; i < individualSpatialGraphNames.size(); ++i)
		{
			outFile << "\t" << individualSpatialGraphNames[i];
		}
		outFile << std::endl;
		for(int i = 0; i < individualSpatialGraphNames.size(); ++i)
		{
			outFile << individualSpatialGraphNames[i];
			for(int j = 0; j < individualSpatialGraphNames.size(); ++j)
			{
				outFile << "\t";
				if(pairWiseDistances[i].find(j) != pairWiseDistances[i].end())
				{
					outFile << pairWiseDistances[i][j];
				}
				else
				{
					outFile << "---";
				}
			}
			outFile << std::endl;
		}
		outFile.close();
		
		std::string outName2(outputFilename);
		outName2 += "_soma_locations.csv";
	
		std::ofstream outFile2(outName2.c_str());
		outFile2 << "Cell ID\tM-L\tA-P\tD-V" << std::endl;
		for(int i = 0; i < somaLocations.size(); ++i)
		{
			outFile2 << individualSpatialGraphNames[i] << "\t" << somaLocations[i][0] << "\t" << somaLocations[i][1] << "\t" << somaLocations[i][2] << std::endl;
		}
		outFile2.close();
	}
	
	if(argc == 4)
	{
		const char * sgListFilename = argv[1];
		const char * surfaceListFilename = argv[2];
		const char * outputFilename = argv[3];

		std::vector< char *  > individualSpatialGraphNames;
		std::vector< char *  > individualSurfaceNames;
		std::vector< double * > somaLocations;
		std::vector< PolyDataPointerType > HVCSurfaces;
		
		std::ifstream inputStream1(sgListFilename);
		if(!inputStream1.fail())
		{
			std::cout << "Loading SpatialGraph files from " << sgListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream1, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					individualSpatialGraphNames.push_back(tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * sgReader = new Reader(tmpChar, outputFilename);
					if(sgStr.find(".am") != std::string::npos)
					{
						std::cout << "Loading SpatialGraph file " << tmpChar << std::endl;
						sgReader->readSpatialGraphFile(1);
						AmiraSpatialGraph * thisSG = sgReader->getSpatialGraph();
						PolyDataPointerType thisSoma = PolyDataPointerType::New();
						if(thisSG->extractLandmark(Soma, thisSoma))
						{
							int subId;
							double * somaLoc = new double[3];
							double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
							thisSoma->GetCell(0)->GetParametricCenter(pCenter);
							thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
							somaLocations.push_back(somaLoc);
							delete [] weights;
						}
						else
						{
							std::cout << "Error! SpatialGraph file " << tmpChar << " does not contain a soma!" << std::endl;
							delete sgReader;
							return 0;
						}
					}
					else
					{
						std::cout << "Error! SpatialGraph file " << tmpChar << " has to be Amira '.am' ascii file!" << std::endl;
						delete sgReader;
						return 0;
					}
					
					delete sgReader;
				}
		}
		
		std::ifstream inputStream2(surfaceListFilename);
		if(!inputStream2.fail())
		{
			std::cout << "Loading surface files from " << surfaceListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream2, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * surfaceReader = new Reader(tmpChar, outputFilename);
					if(sgStr.find(".surf") != std::string::npos)
					{
						std::cout << "Loading surface file " << tmpChar << std::endl;
						PolyDataPointerType thisSurface = PolyDataPointerType::New();
						thisSurface = surfaceReader->readAmiraSurfaceFile();
						HVCSurfaces.push_back(thisSurface);
					}
					else
					{
						std::cout << "Error! Surface file " << tmpChar << " has to be Amira '.surf' ascii file!" << std::endl;
						delete surfaceReader;
						return 0;
					}
					
					delete surfaceReader;
				}
		}
		
		std::vector< double > boundaryDistances;
		for(int i = 0; i < HVCSurfaces.size(); ++i)
		{
			double minDist = 1e6;
			for(int j = 0; j < HVCSurfaces[i]->GetNumberOfCells(); ++j)
			{
				CellPointerType currentCell = HVCSurfaces[i]->GetCell(j);
				if(currentCell->GetNumberOfPoints() != 3)
				{
					std::cout << "Warning! Surface data is not triangles!" << std::endl;
				}
				PointsPointerType tmpPoints = currentCell->GetPoints();
				double centerPoint[] = {0, 0, 0};
				int nrOfPoints = tmpPoints->GetNumberOfPoints();
				for(int k = 0; k < nrOfPoints; ++k)
				{
					double tmpCoords[3];
					tmpPoints->GetPoint(k, tmpCoords);
					for(int n = 0; n < 3; ++n)
						centerPoint[n] += tmpCoords[n];
				}
				for(int k = 0; k < 3; ++k)
					centerPoint[k] = centerPoint[k]/(double)nrOfPoints;
				
				double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(centerPoint, somaLocations[i]));
				if(tmpDist < minDist)
				{
					minDist = tmpDist;
				}
			}
			boundaryDistances.push_back(minDist);
		}
		
		std::string outName(outputFilename);
		if(outName.find(".csv") == std::string::npos)
		{
			outName += ".csv";
		}
		
		std::ofstream outFile(outName.c_str());
		outFile << "SpatialGraphName\tboundary distance (microns)" << std::endl;
		for(int i = 0; i < individualSpatialGraphNames.size(); ++i)
		{
			outFile << individualSpatialGraphNames[i] << "\t" << boundaryDistances[i] << std::endl;
		}
		outFile.close();
	}
	
	return 0;
}



