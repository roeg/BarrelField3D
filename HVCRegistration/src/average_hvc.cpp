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
#include "hvcreconstruction.h"

#define DEBUG

// measure length along "cardinal" anatomical axes (in aligned coordinates)
std::vector< double > getCardinalAxesLengths(PolyDataPointerType alignedSurface, double measuringPoint[3]);
std::vector< int > getOrderedContourIndices(PolyDataPointerType contours);
UnstructuredGridPointerType triangulateAdjacentContours ( PolyDataPointerType contours, int index1, int index2 );
UnstructuredGridPointerType triangulateEndContours ( PolyDataPointerType contours, int index, bool lastContour );

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * surfaceListFilename = argv[1];
		const char * outputFilename = argv[2];
		
		std::vector< PolyDataPointerType > individualSurfaces;
		std::ifstream inputStream(surfaceListFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading surface files from " << surfaceListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					
					std::string surfaceStr(tmpChar);
					Reader * hvcSurfaceReader = new Reader(tmpChar, outputFilename);
					if(surfaceStr.find(".surf") != std::string::npos)
					{
						std::cout << "Loading surface file " << tmpChar << std::endl;
						PolyDataPointerType hvcSurface = hvcSurfaceReader->readAmiraSurfaceFile();
						individualSurfaces.push_back(hvcSurface);
					}
					else
					{
						std::cout << "Error! Surface file" << tmpChar << " has to be Amira '.surf' ascii file!" << std::endl;
						return 0;
					}
					
					delete hvcSurfaceReader;
					delete [] tmpChar;
				}
		}
		
		// average M-L extent defines sampling area ?
		double MLExtent[] = {-800,800};
		double samplingStep = 50;
		double samplingNormal[] = {1,0,0};
		int nrAngles = 36;
		
		std::cout << "Computing triangulation of average surface..." << std::endl;
		Delaunay2DFilterPointerType delTriangulation1 = Delaunay2DFilterPointerType::New();
		PointsPointerType avgHVCPoints = PointsPointerType::New();
		PolyDataPointerType avgHVCMesh = PolyDataPointerType::New(); // Delaunay filter needs DataObject as input...
		avgHVCPoints->SetDataTypeToFloat();
		avgHVCMesh->Allocate(1);
		int avgHVCPointIDCount = 0;
		double DVMeshExtent;
		
		for(double samplingLocation = MLExtent[0]; samplingLocation <= MLExtent[1]; samplingLocation += samplingStep)
// 		for(double samplingLocation = 0; ; )
		{
			double samplingOrigin[] = {samplingLocation,0,0};
			PlanePointerType sagittalPlane = PlanePointerType::New();
			sagittalPlane->SetOrigin(samplingOrigin);
			sagittalPlane->SetNormal(samplingNormal);
			std::vector< std::map< int, double * > > thisPlaneSamplingPoints;	// indices: surface/angle
			
			std::flush(std::cout << "Sampling " << individualSurfaces.size() << " surfaces at location M-L = " << samplingLocation << std::endl);
			
			for(int i = 0; i < individualSurfaces.size(); ++i)
// 			for(int i = 18; ; )
			{
				Surface * currentSurface = new Surface(individualSurfaces[i]);
				CutterPointerType cutSurfaceCoronalFilter = CutterPointerType::New();
				cutSurfaceCoronalFilter->SetCutFunction(sagittalPlane);
				cutSurfaceCoronalFilter->SetInput(individualSurfaces[i]);
				cutSurfaceCoronalFilter->Update();
				PolyDataPointerType coronalSurfaceCut = cutSurfaceCoronalFilter->GetOutput();
				
// 				coronalSurfaceCut->Print(std::cout);
// 				
// 				CellArrayPointerType lineArray = coronalSurfaceCut->GetLines();
// 				coronalSurfaceCut->GetLines()->InitTraversal();
// 				coronalSurfaceCut->GetLines()->Print(std::cout);
// 				IdListPointerType cellIDs = IdListPointerType::New();
// 				int cellCnt = 0;
// 				while(lineArray->GetNextCell(cellIDs))
// 				{
// 					std::flush(std::cout << "line " << cellCnt << ":" << std::endl);
// 					std::flush(std::cout << "nrPts = " << cellIDs->GetNumberOfIds() << std::endl);
// 					for(int ptID = 0; ptID < cellIDs->GetNumberOfIds(); ++ptID)
// 					{
// 						double pt[3];
// 						coronalSurfaceCut->GetPoint(cellIDs->GetId(ptID), pt);
// 						std::flush(std::cout << "pt " << ptID << " = [" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl);
// 					}
// 					++cellCnt;
// 				}
// 				break;
				
				std::map< int, double * > angleSamplingVec;
				thisPlaneSamplingPoints.push_back(angleSamplingVec);
				
				if(!coronalSurfaceCut->GetNumberOfPoints())
				{
// 					for(int ii = 0; ii < nrAngles; ++ii)
// 					{
// 						double * avgPt = new double[3];
// 						avgPt[0] = samplingLocation;
// 						avgPt[1] = 0;
// 						avgPt[2] = 0;
// 						thisPlaneSamplingPoints[i].push_back(avgPt);
// 					}
					continue;
				}
				
				// use different angle sampling origin: has to be inside of cut surface
				double cutOrigin[] = {0,0,0};
				for(int j = 0; j < coronalSurfaceCut->GetNumberOfPoints(); ++j)
				{
					double tmpPt[3];
					coronalSurfaceCut->GetPoint(j, tmpPt);
					cutOrigin[0] += tmpPt[0];
					cutOrigin[1] += tmpPt[1];
					cutOrigin[2] += tmpPt[2];
				}
				cutOrigin[0] /= coronalSurfaceCut->GetNumberOfPoints();
				cutOrigin[1] /= coronalSurfaceCut->GetNumberOfPoints();
				cutOrigin[2] /= coronalSurfaceCut->GetNumberOfPoints();
				std::flush(std::cout << "cutOrigin = [" << cutOrigin[0] << "," << cutOrigin[1] << "," << cutOrigin[2] << "]" << std::endl);
				
// 				// identify polygon normals so each contour can be traversed in same direction
// 				double normal[3], xUnitVec[3] = {1,0,0};
// 				double rotAngle = 0;
// // 				vtkPolygon::ComputeNormal(coronalSurfaceCut->GetCell(0)->GetPoints(), normal);
// 				vtkPolygon::ComputeNormal(coronalSurfaceCut->GetPoints(), normal);
// // 				std::flush(std::cout << "normal = [" << normal[0] << "," << normal[1] << "," << normal[2] << "]" << std::endl);
// 				bool reverseDirection;
// 				if(normal[0] < 0)
// 				{
// 					reverseDirection = 1;
// 					xUnitVec[0] *= -1;
// 				}
// 				else
// 				{
// 					reverseDirection = 0;
// 				}
				
// 				PlanePointerType thisPlane = PlanePointerType::New();
// 				PointsPointerType thisCellPoints = coronalSurfaceCut->GetCell(0)->GetPoints();
// 				PointsPointerType thisCellPoints = coronalSurfaceCut->GetPoints();
// 				thisCellPoints->Print(std::cout);
				for(int jj = 0; jj < nrAngles; ++jj)	// has to be unique b/c order of the points is not clear
				{
					std::flush(std::cout << "Sampling surface " << i << " at angle " << jj << std::endl);
					double angle = jj*2*PI/nrAngles;
// 					double planeNormal[3];
// 					if(reverseDirection)
// 					{
// 						planeNormal[1] = -sin(angle);
// 						planeNormal[2] = -cos(angle);
// 					}
// 					else
// 					{
// 						planeNormal[1] = sin(angle);
// 						planeNormal[2] = cos(angle);
// 					}
// 					planeNormal[0] = 0;
// 					
// 					thisPlane->SetOrigin(cutOrigin);
// 					thisPlane->SetNormal(planeNormal);
					
					double directionVec[3];
					directionVec[0] = 0;
					directionVec[1] = sin(angle);
					directionVec[2] = cos(angle);
					
					currentSurface->intersectLineInDirection(directionVec, cutOrigin);
					if(currentSurface->isIntersectionFound())
					{
						double * intersectPt = new double[3];
						currentSurface->getLastIntersectPoint(intersectPt);
						thisPlaneSamplingPoints[i][jj] = intersectPt;
					}
					
// // 					double firstPt[3];
// // 					double lastVal;
// // 					thisCellPoints->GetPoint(0, firstPt);
// // 					double lastVal = thisPlane->EvaluateFunction(firstPt);
// // 					for(int kk = 1; kk <= thisCellPoints->GetNumberOfPoints(); ++kk)
// 					CellArrayPointerType lineArray = coronalSurfaceCut->GetLines();
// 					coronalSurfaceCut->GetLines()->InitTraversal();
// // 					coronalSurfaceCut->GetLines()->Print(std::cout);
// 					IdListPointerType cellIDs = IdListPointerType::New();
// 					int cellCnt = 0;
// 					while(lineArray->GetNextCell(cellIDs))
// 					{
// // 						std::flush(std::cout << "cell " << 0 << " pt " << kk << std::endl);
// 						
// 						double thisPt[3], lastPt[3];
// // 						thisCellPoints->GetPoint(kk%thisCellPoints->GetNumberOfPoints(), thisPt);
// 						
// 	// 					std::flush(std::cout << "line " << cellCnt << ":" << std::endl);
// // 						std::flush(std::cout << "nrPts = " << cellIDs->GetNumberOfIds() << std::endl);
// // 						for(int ptID = 0; ptID < cellIDs->GetNumberOfIds(); ++ptID)
// // 						{
// // 							double pt[3];
// // 							coronalSurfaceCut->GetPoint(cellIDs->GetId(ptID), pt);
// // 	// 						std::flush(std::cout << "pt " << ptID << " = [" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl);
// // 						}
// 						coronalSurfaceCut->GetPoint(cellIDs->GetId(0), lastPt);
// 						coronalSurfaceCut->GetPoint(cellIDs->GetId(1), thisPt);
// 						++cellCnt;
// 						
// 						double lastVal = thisPlane->EvaluateFunction(lastPt);
// 						double ptVal = thisPlane->EvaluateFunction(thisPt);
// // 						std::flush(std::cout << "lastVal = " << lastVal << std::endl);
// // 						std::flush(std::cout << "ptVal = " << ptVal << std::endl);
// 						if((ptVal > 0 && lastVal < 0))	// has to be unique b/c order of the points is not clear
// 						{
// 							// regular case
// 							double * samplePt =  new double[3];
// // 							double * lastPt = new double[3];
// // 							thisCellPoints->GetPoint(kk-1, lastPt);
// 							double direction[3];
// 							double normD = 0;
// 							for(int ll = 0; ll < 3; ++ll)
// 							{
// 								direction[ll] = lastPt[ll] - thisPt[ll];
// 								normD += direction[ll]*direction[ll];
// 							}
// 							normD = sqrt(normD);
// 							if(normD)
// 							{
// 								double distAlongLine = std::abs(ptVal);
// 								double correctionAngle = std::abs(direction[0]*planeNormal[0]/normD + direction[1]*planeNormal[1]/normD);
// 								if(correctionAngle)
// 									distAlongLine /= correctionAngle;
// 								for(int ll = 0; ll < 3; ++ll)
// 									samplePt[ll] = thisPt[ll] + direction[ll]*distAlongLine/normD;
// 							}
// 							else
// 							{
// 								for(int ll = 0; ll < 3; ++ll)
// 								{
// 									samplePt[ll] = thisPt[ll];
// 								}
// 							}
// 							thisPlaneSamplingPoints[i].push_back(samplePt);
// // 							delete [] lastPt;
// 							std::flush(std::cout << "regular order" << std::endl);
// 							std::flush(std::cout << "pt @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl);
// 						}
// 						else if(lastVal == 0 && ptVal > 0)
// 						{
// 							// case: this point is directly on plane
// 							double * samplePt =  new double[3];
// 							for(int ll = 0; ll < 3; ++ll)
// 							{
// 								samplePt[ll] = thisPt[ll];
// 							}
// 							thisPlaneSamplingPoints[i].push_back(samplePt);
// 							std::flush(std::cout << "pt on plane" << std::endl);
// 							std::flush(std::cout << "pt @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl);
// 						}
// 		// 				else if(ptVal != 0 && lastVal == 0)
// 		// 				{
// 		// 					;// do nothing
// 		// 				}
// 		// 				else if(ptVal == 0 && lastVal == 0)
// 		// 				{
// 		// 					;// either an error or both points lie directly on the plane (extremely unlikely)
// 		// 				}
// // 						lastVal = ptVal;
// // 						delete [] thisPt;
// 					} // all cell points
				} // all angles
			}
			// avg. points here
			std::flush(std::cout << "Averaging points at sample location " << samplingLocation << std::endl);
			std::flush(std::cout << "thisPlaneSamplingPoints.size() = " << thisPlaneSamplingPoints.size() << std::endl);
			for(int i = 0; i < individualSurfaces.size(); ++i)
			{
				std::flush(std::cout << "thisPlaneSamplingPoints[" << i << "].size() = " << thisPlaneSamplingPoints[i].size() << std::endl);
			}
			
			IdListPointerType currentCellPointIDs = IdListPointerType::New();
			double coronalDVBounds[] = {1e6, -1e6};
			for(int j = 0; j < nrAngles; ++j)
			{
				double avgPt[3];
				avgPt[0] = 0, avgPt[1] = 0, avgPt[2] = 0;
				int avgCount = 0;
				for(int i = 0; i < individualSurfaces.size(); ++i)
				{
					if(thisPlaneSamplingPoints[i].find(j) == thisPlaneSamplingPoints[i].end())
					{
						continue;
					}
					double * tmpPt = thisPlaneSamplingPoints[i][j];
					for(int kk = 0; kk < 3; ++kk)
					{
						avgPt[kk] += tmpPt[kk];
					}
					++avgCount;
				}
				for(int jj = 0; jj < 3; ++jj)
				{
					avgPt[jj] = avgPt[jj]/avgCount;
				}
				avgHVCPoints->InsertNextPoint(avgPt);
				currentCellPointIDs->InsertNextId(avgHVCPointIDCount);
				++avgHVCPointIDCount;
				
				if(fabs(samplingLocation) < 1e-6)
				{
					if(avgPt[2] < coronalDVBounds[0])
					{
						coronalDVBounds[0] = avgPt[2];
					}
					if(avgPt[2] > coronalDVBounds[1])
					{
						coronalDVBounds[1] = avgPt[2];
					}
				}
			}
			avgHVCMesh->InsertNextCell(VTK_POLYGON, currentCellPointIDs);
			
			if(fabs(samplingLocation) < 1e-6)
			{
				DVMeshExtent = coronalDVBounds[1] - coronalDVBounds[0];
			}
		}
		
		avgHVCMesh->SetPoints(avgHVCPoints);
		avgHVCMesh->Update();
// 		delTriangulation1->SetInput(avgHVCMesh);
// 		delTriangulation1->Update();
// 		PolyDataPointerType avgHVCSurface = delTriangulation1->GetOutput();
		
		// for some reason, this surface conversion doesn't work out of the box for these contours...
		// simply export contours and create surface in Amira instead...
		AppendFilterPointerType mergeTetrasFilter = AppendFilterPointerType::New();
		
		std::vector< int > contourIndices = getOrderedContourIndices(avgHVCMesh);
		mergeTetrasFilter->AddInput(triangulateEndContours(avgHVCMesh, contourIndices[0], 0));
		for(int ii = 0; ii < contourIndices.size()-1; ++ii)
		{
			mergeTetrasFilter->AddInput(triangulateAdjacentContours(avgHVCMesh, contourIndices[ii], contourIndices[ii+1]));
		}
		mergeTetrasFilter->AddInput(triangulateEndContours(avgHVCMesh, contourIndices.back(), 1));
		mergeTetrasFilter->Update();
		
	// 	return smoothSurface(mergeSurfacesFilter->GetOutput());
// 		DataSetSurfaceFilterPointerType surfaceExtractor = DataSetSurfaceFilterPointerType::New();
// 		surfaceExtractor->SetInput(mergeTetrasFilter->GetOutput());
// 		surfaceExtractor->Update();
	// 	#ifdef DEBUG
	// 	surfaceExtractor->GetOutput()->Print(std::cout);
	// 	#endif
		GeometryFilterPointerType surfaceExtractor = GeometryFilterPointerType::New();
		surfaceExtractor->SetInput(mergeTetrasFilter->GetOutput());
		surfaceExtractor->MergingOn();
		surfaceExtractor->Update();
		
		PolyDataPointerType avgHVCSurface = surfaceExtractor->GetOutput();
		
		// measure extent along standardized axes and scale mesh if necessary
		double registeredOrigin[] = {0,0,0};
		double horizontalNormal[] = {0,0,1};
// 		double coronalNormal[] = {0,1,0};
// 		double sagittalNormal[] = {1,0,0};
		PlanePointerType horizontalPlane = PlanePointerType::New();
		horizontalPlane->SetOrigin(registeredOrigin);
		horizontalPlane->SetNormal(horizontalNormal);
// 		PlanePointerType coronalPlane = PlanePointerType::New();
// 		coronalPlane->SetOrigin(registeredOrigin);
// 		coronalPlane->SetNormal(coronalNormal);
// 		PlanePointerType sagittalPlane = PlanePointerType::New();
// 		sagittalPlane->SetOrigin(registeredOrigin);
// 		sagittalPlane->SetNormal(sagittalNormal);
// 		CutterPointerType cutSurfaceCoronalFilter2 = CutterPointerType::New();
// 		cutSurfaceCoronalFilter2->SetCutFunction(coronalPlane);
// 		cutSurfaceCoronalFilter2->SetInput(avgHVCMesh);
// 		cutSurfaceCoronalFilter2->Update();
// 		PolyDataPointerType coronalSurfaceCut2 = cutSurfaceCoronalFilter2->GetOutput();
		CutterPointerType cutSurfaceHorizontalFilter = CutterPointerType::New();
		cutSurfaceHorizontalFilter->SetCutFunction(horizontalPlane);
		cutSurfaceHorizontalFilter->SetInput(avgHVCMesh);
		cutSurfaceHorizontalFilter->Update();
		PolyDataPointerType horizontalSurfaceCut = cutSurfaceHorizontalFilter->GetOutput();
// 		CutterPointerType cutSurfaceSagittalFilter = CutterPointerType::New();
// 		cutSurfaceSagittalFilter->SetCutFunction(sagittalPlane);
// 		cutSurfaceSagittalFilter->SetInput(avgHVCMesh);
// 		cutSurfaceSagittalFilter->Update();
// 		PolyDataPointerType sagittalSurfaceCut = cutSurfaceSagittalFilter->GetOutput();
		double coronalBB[6], horizontalBB[6], sagittalBB[6];
// 		coronalSurfaceCut2->GetBounds(coronalBB);
		horizontalSurfaceCut->GetBounds(horizontalBB);
// 		sagittalSurfaceCut->GetBounds(sagittalBB);
		
		double MLMeshExtent = horizontalBB[1] - horizontalBB[0];
		double APMeshExtent = horizontalBB[3] - horizontalBB[2];
// 		double DVMeshExtent = sagittalBB[5] - sagittalBB[4];
		const double MLAverage = 1611, APAverage = 821, DVAverage = 455;
		double MLScale = MLAverage/MLMeshExtent;
		double APScale = APAverage/APMeshExtent;
		double DVScale = DVAverage/DVMeshExtent;
		std::cout << "ML mesh = " << MLMeshExtent << "; scale factor = " << MLScale << std::endl;
		std::cout << "AP mesh = " << APMeshExtent << "; scale factor = " << APScale << std::endl;
		std::cout << "DV mesh = " << DVMeshExtent << "; scale factor = " << DVScale << std::endl;
		
		TransformPointerType meshScaleTransform = TransformPointerType::New();
		meshScaleTransform->Scale(MLScale, APScale, DVScale);
		TransformFilterType meshTransformFilter = TransformFilterType::New();
		meshTransformFilter->SetTransform(meshScaleTransform);
		meshTransformFilter->SetInput(avgHVCMesh);
		meshTransformFilter->Update();
		PolyDataPointerType avgHVCMeshScaled = PolyDataPointerType::New();
		avgHVCMeshScaled = meshTransformFilter->GetOutput();
// 		
// 		Reader * avgSurfaceWriter = new Reader(outputFilename, outputFilename);
// 		avgSurfaceWriter->writeAmiraSurfaceFile(avgHVCSurface);
// 		delete avgSurfaceWriter;
		
		AmiraSpatialGraph * averageHVCSG = new AmiraSpatialGraph;
		for(int i = 0; i < avgHVCMeshScaled->GetNumberOfCells(); ++i)
		{
			int edgeConnectivity[2];
			int numEdgePoints;
			int label = D2;
			std::list< double * > edgePointCoordinates;
			IdListPointerType edgePtIDs = avgHVCMeshScaled->GetCell(i)->GetPointIds();
			numEdgePoints = edgePtIDs->GetNumberOfIds();
			if(i == 0 || i == avgHVCMesh->GetNumberOfCells()-1)
			{	
				int subId;
				double center[3], pCenter[3], * weights = new double[numEdgePoints];
				avgHVCMeshScaled->GetCell(i)->GetParametricCenter(pCenter);
				avgHVCMeshScaled->GetCell(i)->EvaluateLocation(subId, pCenter, center, weights);
				delete [] weights;
				
// 				// put center pt at end of section
// 				// to close surface in 3D
// 				if(i == 0)
// 					center[0] -= 50;
// 				else
// 					center[0] += 50;
				
				numEdgePoints = 4;
				double * pt1 = new double[3];
				double * pt2 = new double[3];
				double * pt3 = new double[3];
				double * pt4 = new double[3];
// 				pt1 = avgHVCMeshScaled->GetCell(i)->GetPoints()->GetPoint(0);
// 				pt2 = avgHVCMeshScaled->GetCell(i)->GetPoints()->GetPoint(0);
				for(int j = 0; j < 3; ++j)
				{
					pt1[j] = center[j];
					pt2[j] = center[j];
					pt3[j] = center[j];
					pt4[j] = center[j];
				}
				pt2[1] += 10;
				pt3[2] += 10;
				edgePointCoordinates.push_back(pt1);
				edgePointCoordinates.push_back(pt2);
				edgePointCoordinates.push_back(pt3);
				edgePointCoordinates.push_back(pt4);
				edgeConnectivity[0] = averageHVCSG->getNumberOfVertices();
				edgeConnectivity[1] = averageHVCSG->getNumberOfVertices() + 1;
				Edge * newEdge = new Edge(edgeConnectivity, numEdgePoints, label, edgePointCoordinates);
				Vertex * newVert1 = new Vertex(pt1, label);
				Vertex * newVert2 = new Vertex(pt4, label);
				averageHVCSG->addVertex(newVert1);
				averageHVCSG->addVertex(newVert2);
				averageHVCSG->addEdge(newEdge);
			}
			else
			{
				for(int j = 0; j <= numEdgePoints; ++j)
				{
					double * pt = new double[3];
					avgHVCMeshScaled->GetPoints()->GetPoint(edgePtIDs->GetId(j%numEdgePoints), pt);
					edgePointCoordinates.push_back(pt);
				}
				edgeConnectivity[0] = averageHVCSG->getNumberOfVertices();
				edgeConnectivity[1] = averageHVCSG->getNumberOfVertices() + 1;
				Edge * newEdge = new Edge(edgeConnectivity, numEdgePoints+1, label, edgePointCoordinates);
				Vertex * newVert1 = new Vertex(edgePointCoordinates.front(), label);
				Vertex * newVert2 = new Vertex(edgePointCoordinates.back(), label);
				averageHVCSG->addVertex(newVert1);
				averageHVCSG->addVertex(newVert2);
				averageHVCSG->addEdge(newEdge);
			}
		}
		Reader * avgContourWriter = new Reader(outputFilename, outputFilename);
		avgContourWriter->setSpatialGraph(averageHVCSG);
		avgContourWriter->writeSpatialGraphFile();
		delete avgContourWriter/*, delete averageHVCSG*/;
	}
	
	return 0;
}
// order contour indices by their z coordinate
std::vector< int > getOrderedContourIndices(PolyDataPointerType contours)
{
	std::map< double, int > contourIndexMap;
	for(int ii = 0; ii < contours->GetNumberOfCells(); ++ii)
	{
		double pt[3];
		contours->GetCell(ii)->GetPoints()->GetPoint(0, pt);
		contourIndexMap.insert(std::pair< double, int >(pt[0], ii));
	}
	
	#ifdef DEBUG
	std::cout << " Barreloid contour order:" << std::endl;
	#endif
	std::vector< int > orderedIndices;
	std::map< double, int >::const_iterator contourIndexMapIt;
	for(contourIndexMapIt = contourIndexMap.begin(); contourIndexMapIt != contourIndexMap.end(); ++contourIndexMapIt)
	{
		orderedIndices.push_back(contourIndexMapIt->second);
		#ifdef DEBUG
		std::cout << " x = " << contourIndexMapIt->first << "\tindex = " << contourIndexMapIt->second << std::endl;
		#endif
	}
	
	return orderedIndices;
}

UnstructuredGridPointerType triangulateAdjacentContours ( PolyDataPointerType contours, int index1, int index2 )
{
	Delaunay3DFilterPointerType del3DTriangulation = Delaunay3DFilterPointerType::New();
	PolyDataPointerType contourSubset = PolyDataPointerType::New();
	IdListPointerType subsetIDs = IdListPointerType::New();
	
	subsetIDs->InsertId(0, index1);
	subsetIDs->InsertId(1, index2);
	contourSubset->Allocate(1);
	contourSubset->CopyCells(contours, subsetIDs);
	contourSubset->Update();
	
	del3DTriangulation->SetInput(contourSubset);
	del3DTriangulation->Update();
	
	return del3DTriangulation->GetOutput();
}

// triangulate first/last polygon including center point
// to ensure a closed surface
UnstructuredGridPointerType triangulateEndContours ( PolyDataPointerType contours, int index, bool lastContour )
{
	PolyDataPointerType endPolyData = PolyDataPointerType::New();
	PointsPointerType endPolyPts = PointsPointerType::New();
	PolygonPointerType endPoly = PolygonPointerType::New();
	endPolyData->Allocate(1);
	endPolyPts->SetDataTypeToFloat();
	endPoly->GetPointIds()->SetNumberOfIds(contours->GetCell(index)->GetNumberOfPoints()+1);
	
	PointsPointerType cellPts = contours->GetCell(index)->GetPoints();
	unsigned int cellPtNr = cellPts->GetNumberOfPoints();
	for(int ii = 0; ii < cellPtNr; ++ii)
	{
		double pt[3];
		cellPts->GetPoint(ii, pt);
		endPolyPts->InsertPoint(ii, pt);
		endPoly->GetPointIds()->InsertId(ii, ii);
	}
	
	int subId;
	double center[3], pCenter[3], * weights = new double[cellPtNr];
	contours->GetCell(index)->GetParametricCenter(pCenter);
	contours->GetCell(index)->EvaluateLocation(subId, pCenter, center, weights);
	delete [] weights;
	
	// put center pt at end of section
	// to close surface in 3D
	if(!lastContour)
		center[0] -= 25;
	else
		center[0] += 25;
	#ifdef DEBUG
	if(lastContour)
	{
		std::cout << " Last contour! Center @ [" << center[0] << "," << center[1] << "," << center[2] << "]" << std::endl;
	}
	else if(!lastContour)
	{
		std::cout << " First contour! Center @ [" << center[0] << "," << center[1] << "," << center[2] << "]" << std::endl;
	}
	#endif
	
	endPolyPts->InsertPoint(cellPtNr, center);
	endPoly->GetPointIds()->InsertId(cellPtNr, cellPtNr);
	endPolyData->InsertNextCell(endPoly->GetCellType(), endPoly->GetPointIds());
	endPolyData->SetPoints(endPolyPts);
	
	Delaunay3DFilterPointerType del3DTriangulation = Delaunay3DFilterPointerType::New();
	del3DTriangulation->SetInput(endPolyData);
	del3DTriangulation->Update();
	
	return del3DTriangulation->GetOutput();
}


