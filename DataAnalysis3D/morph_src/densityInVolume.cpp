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
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
		std::string ifName1(inputFilename);
		
		Reader * sgReader1 = new Reader(inputFilename, inputFilename);
		AmiraSpatialGraph * sg1;
		if(ifName1.find(".am") != std::string::npos)
		{
			sgReader1->readSpatialGraphFile(0);
			sg1 = sgReader1->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete sgReader1;
			return 0;
		}
		
		double axonCOM[] = {0,0,0};
		double nrOfPoints = 0;
		for(int i = 0; i < sg1->edgesPointer()->size(); ++i)
		{
			if(sg1->edgesPointer()->at(i)->label != Axon)
			{
				continue;
			}
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = sg1->edgesPointer()->at(i)->edgePointCoordinates.begin(); edgePtIt != sg1->edgesPointer()->at(i)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * tmp = *edgePtIt;
				axonCOM[0] += tmp[0];
				axonCOM[1] += tmp[1];
				axonCOM[2] += tmp[2];
				nrOfPoints += 1.0;
			}
		}
		
		axonCOM[0] /= nrOfPoints;
		axonCOM[1] /= nrOfPoints;
		axonCOM[2] /= nrOfPoints;
		
		std::cout << inputFilename << std::endl;
		std::cout << "1 0 0 0 0 1 0 0 0 0 1 0 " << -axonCOM[0] << " " << -axonCOM[1] << " " << -axonCOM[2] << " 1" << std::endl;
		
		delete sgReader1;
	}
	else if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * surfaceFilename = argv[2];
		std::string ifName(inputFilename);
		std::string surfName(surfaceFilename);
		
		Reader * densityFileReader = new Reader(inputFilename, inputFilename);
		ImageDataPointerType densityField;
		if(ifName.find(".am") != std::string::npos)
		{
			densityField = densityFileReader->readScalarField();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			delete densityFileReader;
			return 0;
		}
		
		Reader * surfReader = new Reader(surfaceFilename, surfaceFilename);
		PolyDataPointerType surface;
		if(surfName.find(".surf") != std::string::npos)
		{
			surface = surfReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
			delete surfReader;
			return 0;
		}
		ClosedSurface * convexHull = new ClosedSurface(surface);
		
		
		std::vector< double > insideVoxelValues;
		double origin[3], voxelSpacing[3];
		int dimensions[3];
		densityField->GetOrigin(origin);
		densityField->GetSpacing(voxelSpacing);
		densityField->GetDimensions(dimensions);
		
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
					
					double * densVal = static_cast< double * >(densityField->GetScalarPointer(klm));
					
					if(convexHull->isPointInsideSurface(xyz))
					{
						insideVoxelValues.push_back(*densVal);
					}
				}
		
		std::string outName(inputFilename);
		outName += "_inside_convex_hull.csv";
		std::ofstream outFile;
		outFile.open(outName.c_str());
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			outFile << insideVoxelValues[i] << std::endl;
		}
		outFile.close();
		
		double touchedVoxels = 0;
		double averageDensity = 0;
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			if(insideVoxelValues[i] > 0)
			{
				++touchedVoxels;
			}
			averageDensity += insideVoxelValues[i];
		}
		touchedVoxels /= insideVoxelValues.size();
		averageDensity /= insideVoxelValues.size();
		double densitySTD = 0;
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			densitySTD += (insideVoxelValues[i] - averageDensity)*(insideVoxelValues[i] - averageDensity);
		}
		densitySTD = sqrt(densitySTD/insideVoxelValues.size());
		
		std::string outName2(inputFilename);
		outName2 += "_inside_convex_hull_summary.csv";
		std::ofstream outFile2;
		outFile2.open(outName2.c_str());
		outFile2 << "Touched voxel fraction\t" << touchedVoxels << std::endl;
		outFile2 << "Density avg\t" << averageDensity << std::endl;
		outFile2 << "Density STD\t" << densitySTD << std::endl;
		outFile2.close();
		
		delete densityFileReader, delete surfReader, delete convexHull;
	}
	
	else if(argc == 4)
	{
		const char * inputFilename = argv[1];
		const char * surfaceFilename = argv[2];
		double voxelSize = atof(argv[3]);
		std::string ifName(inputFilename);
		std::string surfName(surfaceFilename);
		
		Reader * SGREader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".am") != std::string::npos)
		{
			SGREader->readSpatialGraphFile(0);
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete SGREader;
			return 0;
		}
		AmiraSpatialGraph * neuronMorphology = SGREader->getSpatialGraph();
		
		Reader * surfReader = new Reader(surfaceFilename, surfaceFilename);
		PolyDataPointerType surface;
		if(surfName.find(".surf") != std::string::npos)
		{
			surface = surfReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
			delete surfReader;
			return 0;
		}
		ClosedSurface * HVCSurface = new ClosedSurface(surface);
		
		std::vector< double > insideVoxelValues;
		double boundingBox[6];
		HVCSurface->ptr()->GetBounds(boundingBox);
		for(int x = boundingBox[0]; x <= boundingBox[1]; x += voxelSize)
			for(int y = boundingBox[2]; y <= boundingBox[3]; y += voxelSize)
				for(int z = boundingBox[4]; z <= boundingBox[5]; z += voxelSize)
				{
					double xyz[3];
					xyz[0] = x + 0.5*voxelSize;
					xyz[1] = y + 0.5*voxelSize;
					xyz[2] = z + 0.5*voxelSize;
					if(HVCSurface->isPointInsideSurface(xyz))
					{
						double box[6];
						box[0] = x;
						box[1] = x + voxelSize;
						box[2] = y;
						box[3] = y + voxelSize;
						box[4] = z;
						box[5] = z + voxelSize;
						
						AmiraSpatialGraph * tmpSG = neuronMorphology->clipSpatialGraph(box);
						double tmpLength = 0;
						if(tmpSG)
						{
							for(int i = 0; i < tmpSG->edgesPointer()->size(); ++i)
							{
								Edge * tmpEdge = tmpSG->edgesPointer()->at(i);
								if(tmpEdge->label == Axon || tmpEdge->label == ProjectionAxon)
								{
									tmpLength += tmpEdge->segmentLength();
								}
							}
						}
						insideVoxelValues.push_back(tmpLength);
					}
				}
		
		char * outNameTemplate = new char[128];
		sprintf(outNameTemplate, "%s_voxels_%.0fmu", inputFilename, voxelSize);
		std::string outName(outNameTemplate);
		outName += "_inside_HVC.csv";
		std::ofstream outFile;
		outFile.open(outName.c_str());
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			outFile << insideVoxelValues[i] << std::endl;
		}
		outFile.close();
		
		double touchedVoxels = 0;
		double averageDensity = 0;
		double voxelVol = voxelSize*voxelSize*voxelSize;
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			if(insideVoxelValues[i] > 0)
			{
				++touchedVoxels;
			}
			averageDensity += insideVoxelValues[i]/voxelVol;
		}
		touchedVoxels /= insideVoxelValues.size();
		averageDensity /= insideVoxelValues.size();
		double densitySTD = 0;
		for(int i = 0; i < insideVoxelValues.size(); ++i)
		{
			densitySTD += (insideVoxelValues[i]/voxelVol - averageDensity)*(insideVoxelValues[i]/voxelVol - averageDensity);
		}
		densitySTD = sqrt(densitySTD/insideVoxelValues.size());
		
		std::string outName2(outNameTemplate);
		outName2 += "_inside_HVC_summary.csv";
		std::ofstream outFile2;
		outFile2.open(outName2.c_str());
		outFile2 << "Touched voxel fraction\t" << touchedVoxels << std::endl;
		outFile2 << "Density avg\t" << averageDensity << std::endl;
		outFile2 << "Density STD\t" << densitySTD << std::endl;
		outFile2.close();
		
		delete SGREader, delete surfReader, delete HVCSurface;
	}
	
// 	else if(argc == 4)
// 	{
// 		const char * inputFilename1 = argv[1];
// 		const char * inputFilename2 = argv[2];
// 		const char * surfaceFilename = argv[3];
// 		std::string ifName1(inputFilename1);
// 		std::string ifName2(inputFilename2);
// 		std::string surfName(surfaceFilename);
// 		
// 		Reader * densityFileReader1 = new Reader(inputFilename1, inputFilename1);
// 		ImageDataPointerType densityField1;
// 		if(ifName1.find(".am") != std::string::npos)
// 		{
// 			densityField1 = densityFileReader1->readScalarField();
// 		}
// 		else
// 		{
// 			std::cout << "Error! Density file has to be Amira '.am' file!" << std::endl;
// 			delete densityFileReader1;
// 			return 0;
// 		}
// 		
// 		Reader * densityFileReader2 = new Reader(inputFilename2, inputFilename2);
// 		ImageDataPointerType densityField2;
// 		if(ifName2.find(".am") != std::string::npos)
// 		{
// 			densityField2 = densityFileReader2->readScalarField();
// 		}
// 		else
// 		{
// 			std::cout << "Error! Density file has to be Amira '.am' file!" << std::endl;
// 			delete densityFileReader2;
// 			return 0;
// 		}
// 		
// 		Reader * surfReader = new Reader(surfaceFilename, surfaceFilename);
// 		PolyDataPointerType surface;
// 		if(surfName.find(".surf") != std::string::npos)
// 		{
// 			surface = surfReader->readAmiraSurfaceFile();
// 		}
// 		else
// 		{
// 			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
// 			delete surfReader;
// 			return 0;
// 		}
// 		ClosedSurface * convexHull = new ClosedSurface(surface);
// 		
// 		int overlap = 0, dens1Only = 0, dens2Only = 0, insideConvexHull = 0;
// 		std::vector< double > insideVoxelValues1;
// 		std::vector< double > insideVoxelValues2;
// 		
// 		double origin[3], voxelSpacing[3];
// 		int dimensions[3];
// 		densityField1->GetOrigin(origin);
// 		densityField1->GetSpacing(voxelSpacing);
// 		densityField1->GetDimensions(dimensions);
// 		
// 		for(int kk = 0; kk < dimensions[0]; ++kk)
// 			for(int ll = 0; ll < dimensions[1]; ++ll)
// 				for(int mm = 0; mm < dimensions[2]; ++mm)
// 				{
// 					int klm[3];
// 					klm[0] = kk, klm[1] = ll, klm[2] = mm;
// 					
// 					double * densVal1 = static_cast< double * >(densityField1->GetScalarPointer(klm));
// 					double * densVal2 = static_cast< double * >(densityField2->GetScalarPointer(klm));
// 					
// 					bool cornerInsideSurface = false;
// 					double xyz[3];
// 					xyz[0] = origin[0] + kk*voxelSpacing[0];
// 					xyz[1] = origin[1] + ll*voxelSpacing[1];
// 					xyz[2] = origin[2] + mm*voxelSpacing[2];
// 					if(convexHull->isPointInsideSurface(xyz))
// 					{
// 						cornerInsideSurface = true;
// 					}
// 					
// 					for(int i = -1; i < 3; i += 2)
// 						for(int j = -1; j < 3; j += 2)
// 							for(int k = -1; k < 3; k += 2)
// 							{
// 								double pt[3];
// 								pt[0] = xyz[0] + i*0.5*voxelSpacing[0];
// 								pt[1] = xyz[1] + j*0.5*voxelSpacing[1];
// 								pt[2] = xyz[2] + k*0.5*voxelSpacing[2];
// 								if(convexHull->isPointInsideSurface(pt))
// 								{
// 									cornerInsideSurface = true;
// 								}
// 							}
// 					
// 					if(!cornerInsideSurface)
// 					{
// 						*densVal1 = 0;
// 						*densVal2 = 0;
// 					}
// 					else
// 					{
// 						++insideConvexHull;
// 						insideVoxelValues1.push_back(*densVal1);
// 						insideVoxelValues2.push_back(*densVal2);
// 						if(*densVal1 && *densVal2)
// 						{
// 							*densVal1 = 3;
// 							*densVal2 = 3;
// 							++overlap;
// 						}
// 						else if(*densVal1 && *densVal2 == 0)
// 						{
// 							*densVal1 = 2;
// 							*densVal2 = 1;
// 							++dens1Only;
// 						}
// 						else if(*densVal1 == 0 && *densVal2)
// 						{
// 							*densVal1 = 1;
// 							*densVal2 = 2;
// 							++dens2Only;
// 						}
// 					}
// 				}
// 		
// 		std::string outName1(inputFilename1);
// 		outName1 += "_";
// 		outName1 += inputFilename2;
// 		outName1 += "_overlap_density";
// // 		std::string outName2(inputFilename2);
// // 		outName1 += "_";
// // 		outName1 += inputFilename1;
// // 		outName2 += "_overlap_density";
// 		
// 		Reader * densityWriter1 = new Reader(outName1.c_str(), outName1.c_str());
// 		densityWriter1->writeScalarField(densityField1);
// // 		Reader * densityWriter2 = new Reader(outName2.c_str(), outName2.c_str());
// // 		densityWriter2->writeScalarField(densityField2);
// 		
// 		std::string summaryOutName(inputFilename1);
// 		summaryOutName += "_";
// 		summaryOutName += inputFilename2;
// 		summaryOutName += "_overlap_summary.csv";
// 		std::ofstream SummaryFile;
// 		SummaryFile.open(summaryOutName.c_str());
// 		SummaryFile << "File\t" << inputFilename1 << "\t" << inputFilename2 << std::endl;
// 		SummaryFile << "Overlap voxels\t" << overlap << "\t" << overlap << std::endl;
// 		SummaryFile << "Total voxels\t" << overlap+dens1Only << "\t" << overlap+dens2Only << std::endl;
// 		SummaryFile << "Fraction common\t" << float(overlap)/(overlap+dens1Only+dens2Only) << std::endl;
// 		SummaryFile << "Convex hull voxels\t" << insideConvexHull << std::endl;
// 		SummaryFile << "Fraction common of hull\t" << float(overlap)/insideConvexHull << std::endl;
// 		SummaryFile.close();
// 		
// 		std::string summaryOutName2(inputFilename1);
// 		summaryOutName2 += "_";
// 		summaryOutName2 += inputFilename2;
// 		summaryOutName2 += "_overlap_correlation.csv";
// 		std::ofstream SummaryFile2;
// 		SummaryFile2.open(summaryOutName2.c_str());
// 		for(int i = 0; i < insideVoxelValues1.size(); ++i)
// 		{
// 			SummaryFile2 << insideVoxelValues1[i] << "\t" << insideVoxelValues2[i] << std::endl;
// 		}
// 		SummaryFile2.close();
// 		
// 		delete densityFileReader1, delete densityFileReader2, delete densityWriter1, /*delete densityWriter2,*/ delete surfReader, delete convexHull;
// 	}
	
	// pair-wise convex hull comparison
	else if(argc == 5)
	{
		const char * inputFilename1 = argv[1];
		const char * inputFilename2 = argv[2];
		const char * spatialGraphName1 = argv[3];
		const char * spatialGraphName2 = argv[4];
		std::string ifName1(inputFilename1);
		std::string ifName2(inputFilename2);
		std::string sgName1(spatialGraphName1);
		std::string sgName2(spatialGraphName2);
		
		Reader * densityFileReader1 = new Reader(inputFilename1, inputFilename1);
		ImageDataPointerType densityField1;
		if(ifName1.find(".am") != std::string::npos)
		{
			densityField1 = densityFileReader1->readScalarField();
		}
		else
		{
			std::cout << "Error! Density file has to be Amira '.am' file!" << std::endl;
			delete densityFileReader1;
			return 0;
		}
		
		Reader * densityFileReader2 = new Reader(inputFilename2, inputFilename2);
		ImageDataPointerType densityField2;
		if(ifName2.find(".am") != std::string::npos)
		{
			densityField2 = densityFileReader2->readScalarField();
		}
		else
		{
			std::cout << "Error! Density file has to be Amira '.am' file!" << std::endl;
			delete densityFileReader2;
			return 0;
		}
		
		Reader * sgReader1 = new Reader(spatialGraphName1, spatialGraphName1);
		AmiraSpatialGraph * sg1;
		if(sgName1.find(".am") != std::string::npos)
		{
			sgReader1->readSpatialGraphFile(0);
			sg1 = sgReader1->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete sgReader1;
			return 0;
		}
		
		Reader * sgReader2 = new Reader(spatialGraphName2, spatialGraphName2);
		AmiraSpatialGraph * sg2;
		if(sgName2.find(".am") != std::string::npos)
		{
			sgReader2->readSpatialGraphFile(0);
			sg2 = sgReader2->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete sgReader2;
			return 0;
		}
		
		PolyDataPointerType axon1 = PolyDataPointerType::New();
		sg1->extractLandmark(Axon, axon1);
		PolyDataPointerType axon2 = PolyDataPointerType::New();
		sg2->extractLandmark(Axon, axon2);
		
		AppendPolyDataPointerType combineAxonsFilter = AppendPolyDataPointerType::New();
		combineAxonsFilter->AddInput(axon1);
		combineAxonsFilter->AddInput(axon2);
		combineAxonsFilter->Update();
		PolyDataPointerType allAxons = combineAxonsFilter->GetOutput();
		
		Delaunay3DFilterPointerType delaunay3DFilter = Delaunay3DFilterPointerType::New();
// 		delaunay3DFilter->SetInput(allSGPoints);
		delaunay3DFilter->SetInput(allAxons);
		delaunay3DFilter->SetTolerance(1e-6);
		delaunay3DFilter->Update();
		DataSetSurfaceFilterPointerType surfaceFilter = DataSetSurfaceFilterPointerType::New();
		surfaceFilter->SetInput(delaunay3DFilter->GetOutput());
		surfaceFilter->Update();
		PolyDataPointerType hullData = surfaceFilter->GetOutput();
// 		std::string hullOutName(sgName1);
// 		hullOutName += "_";
// 		hullOutName += sgName2;
// 		Reader * hullWriter = new Reader(hullOutName.c_str(), hullOutName.c_str());
// 		hullWriter->writeAmiraSurfaceFile(convexHull);
		
		ClosedSurface * convexHull = new ClosedSurface(hullData);
		
		int overlap = 0, dens1Only = 0, dens2Only = 0, insideConvexHull = 0;
		double overlapLength = 0, totalLength = 0;
		std::vector< double > insideVoxelValues1;
		std::vector< double > insideVoxelValues2;
		
		double origin[3], voxelSpacing[3];
		int dimensions[3];
		densityField1->GetOrigin(origin);
		densityField1->GetSpacing(voxelSpacing);
		densityField1->GetDimensions(dimensions);
		
		for(int kk = 0; kk < dimensions[0]; ++kk)
			for(int ll = 0; ll < dimensions[1]; ++ll)
				for(int mm = 0; mm < dimensions[2]; ++mm)
				{
					int klm[3];
					klm[0] = kk, klm[1] = ll, klm[2] = mm;
					
					double * densVal1 = static_cast< double * >(densityField1->GetScalarPointer(klm));
					double * densVal2 = static_cast< double * >(densityField2->GetScalarPointer(klm));
					
					bool cornerInsideSurface = false;
					double xyz[3];
					xyz[0] = origin[0] + kk*voxelSpacing[0];
					xyz[1] = origin[1] + ll*voxelSpacing[1];
					xyz[2] = origin[2] + mm*voxelSpacing[2];
					if(convexHull->isPointInsideSurface(xyz))
					{
						cornerInsideSurface = true;
					}
					
					for(int i = -1; i < 3; i += 2)
						for(int j = -1; j < 3; j += 2)
							for(int k = -1; k < 3; k += 2)
							{
								double pt[3];
								pt[0] = xyz[0] + i*0.5*voxelSpacing[0];
								pt[1] = xyz[1] + j*0.5*voxelSpacing[1];
								pt[2] = xyz[2] + k*0.5*voxelSpacing[2];
								if(convexHull->isPointInsideSurface(pt))
								{
									cornerInsideSurface = true;
								}
							}
					
					if(!cornerInsideSurface)
					{
						*densVal1 = 0;
						*densVal2 = 0;
					}
					else
					{
						++insideConvexHull;
						insideVoxelValues1.push_back(*densVal1);
						insideVoxelValues2.push_back(*densVal2);
						if(*densVal1 && *densVal2)
						{
							overlapLength += *densVal1;
							overlapLength += *densVal2;
							totalLength += *densVal1;
							totalLength += *densVal2;
							*densVal1 = 3;
							*densVal2 = 3;
							++overlap;
						}
						else if(*densVal1 && *densVal2 == 0)
						{
							totalLength += *densVal1;
							*densVal1 = 2;
							*densVal2 = 1;
							++dens1Only;
						}
						else if(*densVal1 == 0 && *densVal2)
						{
							totalLength += *densVal2;
							*densVal1 = 1;
							*densVal2 = 2;
							++dens2Only;
						}
					}
				}
		
		std::string summaryOutName(inputFilename1);
		summaryOutName += "_";
		summaryOutName += inputFilename2;
		summaryOutName += "_overlap_voxel_length_summary.csv";
		std::ofstream SummaryFile;
		SummaryFile.open(summaryOutName.c_str());
		SummaryFile << "File\t" << inputFilename1 << "\t" << inputFilename2 << std::endl;
		SummaryFile << "Overlap voxels\t" << overlap << "\t" << overlap << std::endl;
		SummaryFile << "Total voxels\t" << overlap+dens1Only << "\t" << overlap+dens2Only << std::endl;
		SummaryFile << "Fraction common\t" << float(overlap)/(overlap+dens1Only+dens2Only) << std::endl;
		SummaryFile << "Convex hull voxels\t" << insideConvexHull << std::endl;
		SummaryFile << "Fraction common of hull\t" << float(overlap)/insideConvexHull << std::endl;
		SummaryFile << std::endl;
		SummaryFile << "Overlap length\t" << overlapLength << std::endl;
		SummaryFile << "Total length\t" << totalLength << std::endl;
		SummaryFile << "Fraction common length\t" << overlapLength/totalLength << std::endl;
		SummaryFile.close();
		
		std::string summaryOutName2(inputFilename1);
		summaryOutName2 += "_";
		summaryOutName2 += inputFilename2;
		summaryOutName2 += "_overlap_correlation.csv";
		std::ofstream SummaryFile2;
		SummaryFile2.open(summaryOutName2.c_str());
		for(int i = 0; i < insideVoxelValues1.size(); ++i)
		{
			SummaryFile2 << insideVoxelValues1[i] << "\t" << insideVoxelValues2[i] << std::endl;
		}
		SummaryFile2.close();
		
		delete densityFileReader1, delete densityFileReader2, delete sgReader1, delete sgReader2/*, delete hullWriter*/;
	}
	
	// pair-wise comparison using axon COM
	// computes density on global grid on the fly
// 	else if(argc == 4)
// 	{
// 		const char * spatialGraphName1 = argv[1];
// 		const char * spatialGraphName2 = argv[2];
// 		double spacing = atof(argv[3]);
// 		std::string sgName1(spatialGraphName1);
// 		std::string sgName2(spatialGraphName2);
// 		
// 		Reader * sgReader1 = new Reader(spatialGraphName1, spatialGraphName1);
// 		AmiraSpatialGraph * sg1;
// 		if(sgName1.find(".am") != std::string::npos)
// 		{
// 			sgReader1->readSpatialGraphFile(0);
// 			sg1 = sgReader1->getSpatialGraph();
// 		}
// 		else
// 		{
// 			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
// 			delete sgReader1;
// 			return 0;
// 		}
// 		
// 		Reader * sgReader2 = new Reader(spatialGraphName2, spatialGraphName2);
// 		AmiraSpatialGraph * sg2;
// 		if(sgName2.find(".am") != std::string::npos)
// 		{
// 			sgReader2->readSpatialGraphFile(0);
// 			sg2 = sgReader2->getSpatialGraph();
// 		}
// 		else
// 		{
// 			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
// 			delete sgReader2;
// 			return 0;
// 		}
// 		
// 		double axon1COM[] = {0,0,0};
// 		double nrOfPoints = 0;
// 		for(int i = 0; i < sg1->edgesPointer()->size(); ++i)
// 		{
// 			if(sg1->edgesPointer()->at(i)->label != Axon)
// 			{
// 				continue;
// 			}
// 			
// 			std::list< double * >::const_iterator edgePtIt;
// 			for(edgePtIt = sg1->edgesPointer()->at(i)->edgePointCoordinates.begin(); edgePtIt != sg1->edgesPointer()->at(i)->edgePointCoordinates.end(); ++edgePtIt)
// 			{
// 				double * tmp = *edgePtIt;
// 				axon1COM[0] += tmp[0];
// 				axon1COM[1] += tmp[1];
// 				axon1COM[2] += tmp[2];
// 				nrOfPoints += 1.0;
// 			}
// 		}
// 		axon1COM[0] /= -nrOfPoints;
// 		axon1COM[1] /= -nrOfPoints;
// 		axon1COM[2] /= -nrOfPoints;
// 		
// 		double axon2COM[] = {0,0,0};
// 		nrOfPoints = 0;
// 		for(int i = 0; i < sg2->edgesPointer()->size(); ++i)
// 		{
// 			if(sg2->edgesPointer()->at(i)->label != Axon)
// 			{
// 				continue;
// 			}
// 			
// 			std::list< double * >::const_iterator edgePtIt;
// 			for(edgePtIt = sg2->edgesPointer()->at(i)->edgePointCoordinates.begin(); edgePtIt != sg2->edgesPointer()->at(i)->edgePointCoordinates.end(); ++edgePtIt)
// 			{
// 				double * tmp = *edgePtIt;
// 				axon2COM[0] += tmp[0];
// 				axon2COM[1] += tmp[1];
// 				axon2COM[2] += tmp[2];
// 				nrOfPoints += 1.0;
// 			}
// 		}
// 		axon2COM[0] /= -nrOfPoints;
// 		axon2COM[1] /= -nrOfPoints;
// 		axon2COM[2] /= -nrOfPoints;
// 		
// // 		TransformPointerType axon1Translate = TransformPointerType::New();
// // 		axon1Translate->Translate(axon1COM);
// // 		sg1->setTransformation(axon1Translate);
// // 		sg1->applyTransformation();
// // 		TransformPointerType axon2Translate = TransformPointerType::New();
// // 		axon2Translate->Translate(axon2COM);
// // 		sg2->setTransformation(axon2Translate);
// // 		sg2->applyTransformation();
// 		
// 		PolyDataPointerType axon1 = PolyDataPointerType::New();
// 		sg1->extractLandmark(Axon, axon1);
// 		PolyDataPointerType axon2 = PolyDataPointerType::New();
// 		sg2->extractLandmark(Axon, axon2);
// 		
// 		AppendPolyDataPointerType combineAxonsFilter = AppendPolyDataPointerType::New();
// 		combineAxonsFilter->AddInput(axon1);
// 		combineAxonsFilter->AddInput(axon2);
// 		combineAxonsFilter->Update();
// 		PolyDataPointerType allAxons = combineAxonsFilter->GetOutput();
// 		
// 		Delaunay3DFilterPointerType delaunay3DFilter = Delaunay3DFilterPointerType::New();
// // 		delaunay3DFilter->SetInput(allSGPoints);
// 		delaunay3DFilter->SetInput(allAxons);
// 		delaunay3DFilter->SetTolerance(1e-6);
// 		delaunay3DFilter->Update();
// 		DataSetSurfaceFilterPointerType surfaceFilter = DataSetSurfaceFilterPointerType::New();
// 		surfaceFilter->SetInput(delaunay3DFilter->GetOutput());
// 		surfaceFilter->Update();
// 		PolyDataPointerType hullData = surfaceFilter->GetOutput();
// // 		std::string hullOutName(sgName1);
// // 		hullOutName += "_";
// // 		hullOutName += sgName2;
// // 		Reader * hullWriter = new Reader(hullOutName.c_str(), hullOutName.c_str());
// // 		hullWriter->writeAmiraSurfaceFile(convexHull);
// 		
// 		ClosedSurface * convexHull = new ClosedSurface(hullData);
// 		
// 		int overlap = 0, dens1Only = 0, dens2Only = 0, insideConvexHull = 0;
// 		double overlapLength = 0, totalLength = 0;
// 		std::vector< double > insideVoxelValues1;
// 		std::vector< double > insideVoxelValues2;
// 		
// 		double hullBounds[6];
// 		convexHull->ptr()->GetBounds(hullBounds);
// 		
// 		double origin[3], voxelSpacing[3];
// 		int dimensions[3];
// 		origin[0] = (int(hullBounds[0]/spacing) - 0.5)*spacing;
// 		origin[1] = (int(hullBounds[2]/spacing) - 0.5)*spacing;
// 		origin[2] = (int(hullBounds[4]/spacing) - 0.5)*spacing;
// 		voxelSpacing[0] = spacing;
// 		voxelSpacing[1] = spacing;
// 		voxelSpacing[2] = spacing;
// 		dimensions[0] = int((hullBounds[1]-hullBounds[0])/spacing) + 1;
// 		dimensions[1] = int((hullBounds[3]-hullBounds[2])/spacing) + 1;
// 		dimensions[2] = int((hullBounds[5]-hullBounds[4])/spacing) + 1;
// 		
// 		ImageDataPointerType overlapDensity = ImageDataPointerType::New();
// 		overlapDensity->SetScalarTypeToDouble();
// 		overlapDensity->SetOrigin(origin);
// 		overlapDensity->SetSpacing(voxelSpacing);
// 		overlapDensity->SetDimensions(dimensions);
// 		overlapDensity->AllocateScalars();
// 		
// 		for(int kk = 0; kk < dimensions[0]; ++kk)
// 			for(int ll = 0; ll < dimensions[1]; ++ll)
// 				for(int mm = 0; mm < dimensions[2]; ++mm)
// 				{
// 					int klm[3];
// 					klm[0] = kk, klm[1] = ll, klm[2] = mm;
// 					
// 					double * densVal = static_cast< double * >(overlapDensity->GetScalarPointer(klm));
// 					
// 					bool cornerInsideSurface = false;
// 					double xyz[3];
// 					xyz[0] = origin[0] + kk*voxelSpacing[0];
// 					xyz[1] = origin[1] + ll*voxelSpacing[1];
// 					xyz[2] = origin[2] + mm*voxelSpacing[2];
// 					if(convexHull->isPointInsideSurface(xyz))
// 					{
// 						cornerInsideSurface = true;
// 					}
// 					
// 					for(int i = -1; i < 3; i += 2)
// 						for(int j = -1; j < 3; j += 2)
// 							for(int k = -1; k < 3; k += 2)
// 							{
// 								double pt[3];
// 								pt[0] = xyz[0] + i*0.5*voxelSpacing[0];
// 								pt[1] = xyz[1] + j*0.5*voxelSpacing[1];
// 								pt[2] = xyz[2] + k*0.5*voxelSpacing[2];
// 								if(convexHull->isPointInsideSurface(pt))
// 								{
// 									cornerInsideSurface = true;
// 								}
// 							}
// 					
// 					if(!cornerInsideSurface)
// 					{
// 						*densVal = 0;
// 					}
// 					else
// 					{
// 						++insideConvexHull;
// 					
// 						double tmpBox[6];
// 						tmpBox[0] = xyz[0] - 0.5*spacing;
// 						tmpBox[1] = xyz[0] + 0.5*spacing;
// 						tmpBox[2] = xyz[1] - 0.5*spacing;
// 						tmpBox[3] = xyz[1] + 0.5*spacing;
// 						tmpBox[4] = xyz[2] - 0.5*spacing;
// 						tmpBox[5] = xyz[2] + 0.5*spacing;
// 						AmiraSpatialGraph * tmpSG1 = sg1->clipSpatialGraph(tmpBox);
// 						AmiraSpatialGraph * tmpSG2 = sg2->clipSpatialGraph(tmpBox);
// 						
// 						double lengthSG1 = 0;
// 						double lengthSG2 = 0;
// 						if(tmpSG1)
// 						{
// 							for(int i = 0; i < tmpSG1->edgesPointer()->size(); ++i)
// 							{
// 								if(tmpSG1->edgesPointer()->at(i)->label == Axon)
// 								{
// 									lengthSG1 += tmpSG1->edgesPointer()->at(i)->segmentLength();
// 								}
// 							}
// 						}
// 						if(tmpSG2)
// 						{
// 							for(int i = 0; i < tmpSG2->edgesPointer()->size(); ++i)
// 							{
// 								if(tmpSG2->edgesPointer()->at(i)->label == Axon)
// 								{
// 									lengthSG2 += tmpSG2->edgesPointer()->at(i)->segmentLength();
// 								}
// 							}
// 						}
// 						
// 						insideVoxelValues1.push_back(lengthSG1);
// 						insideVoxelValues2.push_back(lengthSG2);
// 						if(lengthSG1 && lengthSG2)
// 						{
// 							overlapLength += lengthSG1;
// 							overlapLength += lengthSG2;
// 							totalLength += lengthSG1;
// 							totalLength += lengthSG2;
// 							*densVal = 3;
// 							++overlap;
// 						}
// 						else if(lengthSG1 && lengthSG2 == 0)
// 						{
// 							totalLength += lengthSG1;
// 							*densVal = 2;
// 							++dens1Only;
// 						}
// 						else if(lengthSG1 == 0 && lengthSG2)
// 						{
// 							totalLength += lengthSG2;
// 							*densVal = 1;
// 							++dens2Only;
// 						}
// 					}
// 				}
// 		
// 		char * voxelSizeChar = new char[64];
// 		sprintf(voxelSizeChar, "_voxel_size_%.0f_", spacing);
// 		
// 		std::string overlapName(spatialGraphName1);
// 		overlapName += "_";
// 		overlapName += spatialGraphName2;
// 		overlapName += voxelSizeChar;
// // 		overlapName += "axon_COM_overlap_density";
// 		overlapName += "soma_COM_overlap_density";
// 		Reader * overlapDensityWriter = new Reader(overlapName.c_str(), overlapName.c_str());
// 		overlapDensityWriter->writeScalarField(overlapDensity);
// 		
// 		std::string summaryOutName(spatialGraphName1);
// 		summaryOutName += "_";
// 		summaryOutName += spatialGraphName2;
// 		summaryOutName += voxelSizeChar;
// // 		summaryOutName += "axon_COM_overlap_voxel_length_summary.csv";
// 		summaryOutName += "soma_COM_overlap_voxel_length_summary.csv";
// 		std::ofstream SummaryFile;
// 		SummaryFile.open(summaryOutName.c_str());
// 		SummaryFile << "File\t" << spatialGraphName1 << "\t" << spatialGraphName2 << std::endl;
// 		SummaryFile << "Overlap voxels\t" << overlap << "\t" << overlap << std::endl;
// 		SummaryFile << "Total voxels\t" << overlap+dens1Only << "\t" << overlap+dens2Only << std::endl;
// 		SummaryFile << "Fraction common\t" << float(overlap)/(overlap+dens1Only+dens2Only) << std::endl;
// 		SummaryFile << "Convex hull voxels\t" << insideConvexHull << std::endl;
// 		SummaryFile << "Fraction common of hull\t" << float(overlap)/insideConvexHull << std::endl;
// 		SummaryFile << std::endl;
// 		SummaryFile << "Overlap length\t" << overlapLength << std::endl;
// 		SummaryFile << "Total length\t" << totalLength << std::endl;
// 		SummaryFile << "Fraction common length\t" << overlapLength/totalLength << std::endl;
// 		SummaryFile.close();
// 		
// 		std::string summaryOutName2(spatialGraphName1);
// 		summaryOutName2 += "_";
// 		summaryOutName2 += spatialGraphName2;
// 		summaryOutName2 += voxelSizeChar;
// // 		summaryOutName2 += "axon_COM_overlap_correlation.csv";
// 		summaryOutName2 += "soma_COM_overlap_correlation.csv";
// 		std::ofstream SummaryFile2;
// 		SummaryFile2.open(summaryOutName2.c_str());
// 		for(int i = 0; i < insideVoxelValues1.size(); ++i)
// 		{
// 			SummaryFile2 << insideVoxelValues1[i] << "\t" << insideVoxelValues2[i] << std::endl;
// 		}
// 		SummaryFile2.close();
// 		
// 		delete sgReader1, delete sgReader2/*, delete hullWriter*/, delete voxelSizeChar;
// 	}
	
	// 2D version
	else if(argc == 4)
	{
#define AXONCOM
// #define ZAXIS
// #define XAXIS
#define YAXIS
		
		const char * spatialGraphName1 = argv[1];
		const char * spatialGraphName2 = argv[2];
		double spacing = atof(argv[3]);
		std::string sgName1(spatialGraphName1);
		std::string sgName2(spatialGraphName2);
		
		Reader * sgReader1 = new Reader(spatialGraphName1, spatialGraphName1);
		AmiraSpatialGraph * sg1;
		if(sgName1.find(".am") != std::string::npos)
		{
			sgReader1->readSpatialGraphFile(0);
			sg1 = sgReader1->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete sgReader1;
			return 0;
		}
		
		Reader * sgReader2 = new Reader(spatialGraphName2, spatialGraphName2);
		AmiraSpatialGraph * sg2;
		if(sgName2.find(".am") != std::string::npos)
		{
			sgReader2->readSpatialGraphFile(0);
			sg2 = sgReader2->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! SpatialGraph file has to be Amira '.am' file!" << std::endl;
			delete sgReader2;
			return 0;
		}
		
		ConvexHull2DFilterPointerType convexHull2DFilter = ConvexHull2DFilterPointerType::New();
		convexHull2DFilter->SetDataTypeToFloat();
		
#ifdef AXONCOM
		double axon1COM[] = {0,0,0};
		double nrOfPoints = 0;
		for(int i = 0; i < sg1->edgesPointer()->size(); ++i)
		{
			if(sg1->edgesPointer()->at(i)->label != Axon)
			{
				continue;
			}
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = sg1->edgesPointer()->at(i)->edgePointCoordinates.begin(); edgePtIt != sg1->edgesPointer()->at(i)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * tmp = *edgePtIt;
				axon1COM[0] += tmp[0];
				axon1COM[1] += tmp[1];
				axon1COM[2] += tmp[2];
				nrOfPoints += 1.0;
				convexHull2DFilter->InsertNextPoint(tmp);
			}
		}
		axon1COM[0] /= -nrOfPoints;
		axon1COM[1] /= -nrOfPoints;
		axon1COM[2] /= -nrOfPoints;
		
		double axon2COM[] = {0,0,0};
		nrOfPoints = 0;
		for(int i = 0; i < sg2->edgesPointer()->size(); ++i)
		{
			if(sg2->edgesPointer()->at(i)->label != Axon)
			{
				continue;
			}
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = sg2->edgesPointer()->at(i)->edgePointCoordinates.begin(); edgePtIt != sg2->edgesPointer()->at(i)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * tmp = *edgePtIt;
				axon2COM[0] += tmp[0];
				axon2COM[1] += tmp[1];
				axon2COM[2] += tmp[2];
				nrOfPoints += 1.0;
				convexHull2DFilter->InsertNextPoint(tmp);
			}
		}
		axon2COM[0] /= -nrOfPoints;
		axon2COM[1] /= -nrOfPoints;
		axon2COM[2] /= -nrOfPoints;
		
		TransformPointerType axon1Translate = TransformPointerType::New();
		axon1Translate->Translate(axon1COM);
		sg1->setTransformation(axon1Translate);
		sg1->applyTransformation();
		TransformPointerType axon2Translate = TransformPointerType::New();
		axon2Translate->Translate(axon2COM);
		sg2->setTransformation(axon2Translate);
		sg2->applyTransformation();
#endif
		
		PolyDataPointerType axon1 = PolyDataPointerType::New();
		sg1->extractLandmark(Axon, axon1);
		PolyDataPointerType axon2 = PolyDataPointerType::New();
		sg2->extractLandmark(Axon, axon2);
		
		AppendPolyDataPointerType combineAxonsFilter = AppendPolyDataPointerType::New();
		combineAxonsFilter->AddInput(axon1);
		combineAxonsFilter->AddInput(axon2);
		combineAxonsFilter->Update();
		PolyDataPointerType allAxons = combineAxonsFilter->GetOutput();
		
		TransformPointerType xToZAxis = TransformPointerType::New();
		xToZAxis->RotateY(-90.0);
		TransformPointerType zToXAxis = TransformPointerType::New();
		zToXAxis->RotateY(90.0);
		TransformPointerType yToZAxis = TransformPointerType::New();
		yToZAxis->RotateX(-90.0);
		TransformPointerType zToYAxis = TransformPointerType::New();
		zToYAxis->RotateX(90.0);
		
		ConvexHull2DFilterPointerType2 convexHull2DFilterZ = ConvexHull2DFilterPointerType2::New();
		convexHull2DFilterZ->SetInput(allAxons);
		convexHull2DFilterZ->Update();
		PolyDataPointerType convexHullZ = convexHull2DFilterZ->GetOutput();
// 		convexHullZ->Print(std::cout);
		
		TransformFilterType allAxonsYToZ = TransformFilterType::New();
		allAxonsYToZ->SetTransform(yToZAxis);
		allAxonsYToZ->SetInput(allAxons);
		allAxonsYToZ->Update();
		ConvexHull2DFilterPointerType2 convexHull2DFilterY = ConvexHull2DFilterPointerType2::New();
		convexHull2DFilterY->SetInput(allAxonsYToZ->GetOutput());
		convexHull2DFilterY->Update();
		TransformFilterType allAxonsHullZToY = TransformFilterType::New();
		allAxonsHullZToY->SetTransform(zToYAxis);
		allAxonsHullZToY->SetInput(convexHull2DFilterY->GetOutput());
		allAxonsHullZToY->Update();
		PolyDataPointerType convexHullY = allAxonsHullZToY->GetOutput();
// 		convexHullY->Print(std::cout);
		
		TransformFilterType allAxonsXToZ = TransformFilterType::New();
		allAxonsXToZ->SetTransform(xToZAxis);
		allAxonsXToZ->SetInput(allAxons);
		allAxonsXToZ->Update();
		ConvexHull2DFilterPointerType2 convexHull2DFilterX = ConvexHull2DFilterPointerType2::New();
		convexHull2DFilterX->SetInput(allAxonsXToZ->GetOutput());
		convexHull2DFilterX->Update();
		TransformFilterType allAxonsHullZToX = TransformFilterType::New();
		allAxonsHullZToX->SetTransform(zToXAxis);
		allAxonsHullZToX->SetInput(convexHull2DFilterX->GetOutput());
		allAxonsHullZToX->Update();
		PolyDataPointerType convexHullX = allAxonsHullZToX->GetOutput();
// 		convexHullX->Print(std::cout);
		
// 		std::cout << "X-projection convex hull nr. of points: " << convexHull2DFilter->GetSizeCCWHullX() << std::endl;
// 		std::cout << "Y-projection convex hull nr. of points: " << convexHull2DFilter->GetSizeCCWHullY() << std::endl;
// 		std::cout << "Z-projection convex hull nr. of points: " << convexHull2DFilter->GetSizeCCWHullZ() << std::endl;
		
// 		Delaunay2DFilterPointerType delaunay2DFilter = Delaunay2DFilterPointerType::New();
// // 		delaunay3DFilter->SetInput(allSGPoints);
// 		delaunay2DFilter->SetInput(allAxons);
// 		delaunay2DFilter->SetTolerance(1e-6);
// 		delaunay2DFilter->Update();
// 		PolyDataPointerType hullData = delaunay2DFilter->GetOutput();
// 		hullData->Print(std::cout);
		
// 		DataSetSurfaceFilterPointerType surfaceFilter = DataSetSurfaceFilterPointerType::New();
// 		surfaceFilter->SetInput(delaunay3DFilter->GetOutput());
// 		surfaceFilter->Update();
// 		PolyDataPointerType hullData = surfaceFilter->GetOutput();
// 		std::string hullOutName(sgName1);
// 		hullOutName += "_";
// 		hullOutName += sgName2;
// 		Reader * hullWriter = new Reader(hullOutName.c_str(), hullOutName.c_str());
// 		hullWriter->writeAmiraSurfaceFile(convexHull);
		
// 		ClosedSurface * convexHull = new ClosedSurface(hullData);
		
		int overlap = 0, dens1Only = 0, dens2Only = 0, insideConvexHull = 0;
		double overlapLength = 0, totalLength = 0;
		std::vector< double > insideVoxelValues1;
		std::vector< double > insideVoxelValues2;
		
		double hullBoundsX[6];
		double hullBoundsY[6];
		double hullBoundsZ[6];
		convexHullX->GetBounds(hullBoundsX);
		convexHullY->GetBounds(hullBoundsY);
		convexHullZ->GetBounds(hullBoundsZ);
		
		double originX[3], voxelSpacingX[3];
		int dimensionsX[3];
#ifdef ZAXIS
		originX[0] = (int(hullBoundsZ[0]/spacing) - 0.5)*spacing;
		originX[1] = (int(hullBoundsZ[2]/spacing) - 0.5)*spacing;
		originX[2] = -0.5e4;
#endif
#ifdef XAXIS
		originX[0] = -0.5e4;
		originX[1] = (int(hullBoundsX[2]/spacing) - 0.5)*spacing;
		originX[2] = (int(hullBoundsX[4]/spacing) - 0.5)*spacing;
#endif
#ifdef YAXIS
		originX[0] = (int(hullBoundsY[0]/spacing) - 0.5)*spacing;
		originX[1] = -0.5e4;
		originX[2] = (int(hullBoundsY[4]/spacing) - 0.5)*spacing;
#endif
#ifdef ZAXIS
		voxelSpacingX[0] = spacing;
		voxelSpacingX[1] = spacing;
		voxelSpacingX[2] = 1e4;
#endif
#ifdef XAXIS
		voxelSpacingX[0] = 1e4;
		voxelSpacingX[1] = spacing;
		voxelSpacingX[2] = spacing;
#endif
#ifdef YAXIS
		voxelSpacingX[0] = spacing;
		voxelSpacingX[1] = 1e4;
		voxelSpacingX[2] = spacing;
#endif
#ifdef ZAXIS
		dimensionsX[0] = int((hullBoundsZ[1]-hullBoundsZ[0])/voxelSpacingX[0]) + 1;
		dimensionsX[1] = int((hullBoundsZ[3]-hullBoundsZ[2])/voxelSpacingX[1]) + 1;
		dimensionsX[2] = int((hullBoundsZ[5]-hullBoundsZ[4])/voxelSpacingX[2]) + 1;
#endif
#ifdef XAXIS
		dimensionsX[0] = int((hullBoundsX[1]-hullBoundsX[0])/voxelSpacingX[0]) + 1;
		dimensionsX[1] = int((hullBoundsX[3]-hullBoundsX[2])/voxelSpacingX[1]) + 1;
		dimensionsX[2] = int((hullBoundsX[5]-hullBoundsX[4])/voxelSpacingX[2]) + 1;
#endif
#ifdef YAXIS
		dimensionsX[0] = int((hullBoundsY[1]-hullBoundsY[0])/voxelSpacingX[0]) + 1;
		dimensionsX[1] = int((hullBoundsY[3]-hullBoundsY[2])/voxelSpacingX[1]) + 1;
		dimensionsX[2] = int((hullBoundsY[5]-hullBoundsY[4])/voxelSpacingX[2]) + 1;
#endif
		
		ImageDataPointerType overlapDensityX = ImageDataPointerType::New();
		overlapDensityX->SetScalarTypeToDouble();
		overlapDensityX->SetOrigin(originX);
		overlapDensityX->SetSpacing(voxelSpacingX);
		overlapDensityX->SetDimensions(dimensionsX);
		overlapDensityX->AllocateScalars();
		
		for(int kk = 0; kk < dimensionsX[0]; ++kk)
			for(int ll = 0; ll < dimensionsX[1]; ++ll)
				for(int mm = 0; mm < dimensionsX[2]; ++mm)
				{
					int klm[3];
					klm[0] = kk, klm[1] = ll, klm[2] = mm;
					double xyz[3];
					xyz[0] = originX[0] + kk*voxelSpacingX[0];
					xyz[1] = originX[1] + ll*voxelSpacingX[1];
					xyz[2] = originX[2] + mm*voxelSpacingX[2];
					
					double * densVal = static_cast< double * >(overlapDensityX->GetScalarPointer(klm));
					
					double * closestPoint;
					int subId;
					double pcoords[3];
					double dist2;
#ifdef ZAXIS
					double * weights = new double[convexHullZ->GetCell(0)->GetNumberOfPoints()];
#endif
#ifdef XAXIS
					double * weights = new double[convexHullX->GetCell(0)->GetNumberOfPoints()];
#endif
#ifdef YAXIS
					double * weights = new double[convexHullY->GetCell(0)->GetNumberOfPoints()];
#endif
					
					bool cornerInsideSurface = false;
#ifdef ZAXIS
					if(convexHullZ->GetCell(0)->EvaluatePosition(xyz, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
#ifdef XAXIS
					if(convexHullX->GetCell(0)->EvaluatePosition(xyz, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
#ifdef YAXIS
					if(convexHullY->GetCell(0)->EvaluatePosition(xyz, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
					{
						cornerInsideSurface = true;
					}
					for(int i = -1; i < 3; i += 2)
						for(int j = -1; j < 3; j += 2)
							for(int k = -1; k < 3; k += 2)
							{
								double pt[3];
								pt[0] = xyz[0] + i*0.5*voxelSpacingX[0];
								pt[1] = xyz[1] + j*0.5*voxelSpacingX[1];
								pt[2] = xyz[2] + k*0.5*voxelSpacingX[2];
// 								pt[0] = xyz[0];
// 								pt[1] = xyz[1] + j*0.5*voxelSpacingX[1];
// 								pt[2] = xyz[2] + k*0.5*voxelSpacingX[2];
// 								pt[0] = xyz[0] + i*0.5*voxelSpacingX[0];
// 								pt[1] = xyz[1];
// 								pt[2] = xyz[2] + k*0.5*voxelSpacingX[2];
#ifdef ZAXIS
								if(convexHullZ->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
#ifdef XAXIS
								if(convexHullX->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
#ifdef YAXIS
								if(convexHullY->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pcoords, dist2, weights) > 0)
#endif
								{
									cornerInsideSurface = true;
								}
							}
					
					if(!cornerInsideSurface)
					{
						*densVal = 0;
					}
					else
					{
						++insideConvexHull;
						
						double tmpBox[6];
#ifdef ZAXIS
						tmpBox[0] = xyz[0] - 0.5*voxelSpacingX[0];
						tmpBox[1] = xyz[0] + 0.5*voxelSpacingX[0];
						tmpBox[2] = xyz[1] - 0.5*voxelSpacingX[1];
						tmpBox[3] = xyz[1] + 0.5*voxelSpacingX[1];
						tmpBox[4] = -5e3;
						tmpBox[5] = 5e3;
#endif
#ifdef XAXIS
						tmpBox[0] = -5e3;
						tmpBox[1] = 5e3;
						tmpBox[2] = xyz[1] - 0.5*voxelSpacingX[1];
						tmpBox[3] = xyz[1] + 0.5*voxelSpacingX[1];
						tmpBox[4] = xyz[2] - 0.5*voxelSpacingX[2];
						tmpBox[5] = xyz[2] + 0.5*voxelSpacingX[2];
#endif
#ifdef YAXIS
						tmpBox[0] = xyz[0] - 0.5*voxelSpacingX[0];
						tmpBox[1] = xyz[0] + 0.5*voxelSpacingX[0];
						tmpBox[2] = -5e3;
						tmpBox[3] = 5e3;
						tmpBox[4] = xyz[2] - 0.5*voxelSpacingX[2];
						tmpBox[5] = xyz[2] + 0.5*voxelSpacingX[2];
#endif
						AmiraSpatialGraph * tmpSG1 = sg1->clipSpatialGraph(tmpBox);
						AmiraSpatialGraph * tmpSG2 = sg2->clipSpatialGraph(tmpBox);
						
						double lengthSG1 = 0;
						double lengthSG2 = 0;
						if(tmpSG1)
						{
							for(int i = 0; i < tmpSG1->edgesPointer()->size(); ++i)
							{
								if(tmpSG1->edgesPointer()->at(i)->label == Axon)
								{
									lengthSG1 += tmpSG1->edgesPointer()->at(i)->segmentLength();
								}
							}
						}
						if(tmpSG2)
						{
							for(int i = 0; i < tmpSG2->edgesPointer()->size(); ++i)
							{
								if(tmpSG2->edgesPointer()->at(i)->label == Axon)
								{
									lengthSG2 += tmpSG2->edgesPointer()->at(i)->segmentLength();
								}
							}
						}
						
						insideVoxelValues1.push_back(lengthSG1);
						insideVoxelValues2.push_back(lengthSG2);
						if(lengthSG1 && lengthSG2)
						{
							overlapLength += lengthSG1;
							overlapLength += lengthSG2;
							totalLength += lengthSG1;
							totalLength += lengthSG2;
							*densVal = 3;
							++overlap;
						}
						else if(lengthSG1 && lengthSG2 == 0)
						{
							totalLength += lengthSG1;
							*densVal = 2;
							++dens1Only;
						}
						else if(lengthSG1 == 0 && lengthSG2)
						{
							totalLength += lengthSG2;
							*densVal = 1;
							++dens2Only;
						}
						else if(lengthSG1 == 0 && lengthSG2 == 0)
						{
							*densVal = 0;
						}
					}
					
					delete [] weights;
				}
		
		char * voxelSizeChar = new char[64];
		sprintf(voxelSizeChar, "_voxel_size_%.0f_", spacing);
		
		std::string overlapName(spatialGraphName1);
		overlapName += "_";
		overlapName += spatialGraphName2;
		overlapName += voxelSizeChar;
// 		overlapName += "axon_COM_overlap_density";
#ifdef ZAXIS
#ifndef AXONCOM
		overlapName += "soma_COM_overlap_density_X-Y_plane";
#endif
#ifdef AXONCOM
		overlapName += "axon_COM_overlap_density_X-Y_plane";
#endif
#endif
#ifdef XAXIS
#ifndef AXONCOM
		overlapName += "soma_COM_overlap_density_Y-Z_plane";
#endif
#ifdef AXONCOM
		overlapName += "axon_COM_overlap_density_Y-Z_plane";
#endif
#endif
#ifdef YAXIS
#ifndef AXONCOM
		overlapName += "soma_COM_overlap_density_X-Z_plane";
#endif
#ifdef AXONCOM
		overlapName += "axon_COM_overlap_density_X-Z_plane";
#endif
#endif
		Reader * overlapDensityWriter = new Reader(overlapName.c_str(), overlapName.c_str());
		overlapDensityWriter->writeScalarField(overlapDensityX);
		
		std::string summaryOutName(spatialGraphName1);
		summaryOutName += "_";
		summaryOutName += spatialGraphName2;
		summaryOutName += voxelSizeChar;
// 		summaryOutName += "axon_COM_overlap_voxel_length_summary.csv";
#ifdef ZAXIS
#ifndef AXONCOM
		summaryOutName += "soma_COM_overlap_voxel_length_X-Y_plane_summary.csv";
#endif
#ifdef AXONCOM
		summaryOutName += "axon_COM_overlap_voxel_length_X-Y_plane_summary.csv";
#endif
#endif
#ifdef XAXIS
#ifndef AXONCOM
		summaryOutName += "soma_COM_overlap_voxel_length_Y-Z_plane_summary.csv";
#endif
#ifdef AXONCOM
		summaryOutName += "axon_COM_overlap_voxel_length_Y-Z_plane_summary.csv";
#endif
#endif
#ifdef YAXIS
#ifndef AXONCOM
		summaryOutName += "soma_COM_overlap_voxel_length_X-Z_plane_summary.csv";
#endif
#ifdef AXONCOM
		summaryOutName += "axon_COM_overlap_voxel_length_X-Z_plane_summary.csv";
#endif
#endif
		std::ofstream SummaryFile;
		SummaryFile.open(summaryOutName.c_str());
		SummaryFile << "File\t" << spatialGraphName1 << "\t" << spatialGraphName2 << std::endl;
		SummaryFile << "Overlap voxels\t" << overlap << "\t" << overlap << std::endl;
		SummaryFile << "Total voxels\t" << overlap+dens1Only << "\t" << overlap+dens2Only << std::endl;
		SummaryFile << "Fraction common\t" << float(overlap)/(overlap+dens1Only+dens2Only) << std::endl;
		SummaryFile << "Convex hull voxels\t" << insideConvexHull << std::endl;
		SummaryFile << "Fraction common of hull\t" << float(overlap)/insideConvexHull << std::endl;
		SummaryFile << std::endl;
		SummaryFile << "Overlap length\t" << overlapLength << std::endl;
		SummaryFile << "Total length\t" << totalLength << std::endl;
		SummaryFile << "Fraction common length\t" << overlapLength/totalLength << std::endl;
		SummaryFile.close();
		
		std::string summaryOutName2(spatialGraphName1);
		summaryOutName2 += "_";
		summaryOutName2 += spatialGraphName2;
		summaryOutName2 += voxelSizeChar;
// 		summaryOutName2 += "axon_COM_overlap_correlation.csv";
#ifdef ZAXIS
#ifndef AXONCOM
		summaryOutName2 += "soma_COM_overlap_X-Y_plane_correlation.csv";
#endif
#ifdef AXONCOM
		summaryOutName2 += "axon_COM_overlap_X-Y_plane_correlation.csv";
#endif
#endif
#ifdef XAXIS
#ifndef AXONCOM
		summaryOutName2 += "soma_COM_overlap_Y-Z_plane_correlation.csv";
#endif
#ifdef AXONCOM
		summaryOutName2 += "axon_COM_overlap_Y-Z_plane_correlation.csv";
#endif
#endif
#ifdef YAXIS
#ifndef AXONCOM
		summaryOutName2 += "soma_COM_overlap_X-Z_plane_correlation.csv";
#endif
#ifdef AXONCOM
		summaryOutName2 += "axon_COM_overlap_X-Z_plane_correlation.csv";
#endif
#endif
		std::ofstream SummaryFile2;
		SummaryFile2.open(summaryOutName2.c_str());
		for(int i = 0; i < insideVoxelValues1.size(); ++i)
		{
			SummaryFile2 << insideVoxelValues1[i] << "\t" << insideVoxelValues2[i] << std::endl;
		}
		SummaryFile2.close();
		
		delete sgReader1, delete sgReader2/*, delete hullWriter*/, delete voxelSizeChar;
	}
	
 	else
 	{
 		std::cout << "Usage: DensityInVolume [Density filename] [Closed surface filename]" << std::endl;
 	}
	return 0;
}
