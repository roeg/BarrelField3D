#include "cortexreconstruction.h"
// #define DEBUG
#define MAXCONTOUR
// #define AVGCONTOUR
#define USEVESSELS

CortexReconstruction::CortexReconstruction ( AmiraSpatialGraph* inputSG, InputParameters parameters )
{
	spatialGraph = new AmiraSpatialGraph;
	spatialGraph->mergeSpatialGraph(inputSG);
	this->parameters = parameters;
	SBF = new BarrelField(false);
}

CortexReconstruction::CortexReconstruction ( AmiraSpatialGraph* inputSG )
{
	spatialGraph = new AmiraSpatialGraph;
	spatialGraph->mergeSpatialGraph(inputSG);
	SBF = new BarrelField(false);
}

CortexReconstruction::~CortexReconstruction()
{
	if(spatialGraph) delete spatialGraph;
	if(SBF) delete SBF;
}

/****************************************************************************/
/* public interface                                                         */
/****************************************************************************/

std::map< int, Column* > CortexReconstruction::getBarrelField()
{
	return finalBarrels;
}

Surface* CortexReconstruction::getPiaSurface()
{
	Surface * pia = new Surface(piaSurface);
	return pia;
}

Surface* CortexReconstruction::getWMSurface()
{
	Surface * wm = new Surface(wmSurface);
	return wm;
}

void CortexReconstruction::startReconstruction(const char* outputFilename)
{
// 	#ifdef DEBUG
	std::cout << "Starting 3D reconstruction!" << std::endl;
// 	#endif
	if(parameters.piaFlag)
	{
		// new surface reconstruction:
		PolyDataPointerType piaStructure = PolyDataPointerType::New();
		if(!spatialGraph->extractLandmark(Pia, piaStructure))
		{
			std::cout << "Error! Pia not in SpatialGraph!" << std::endl;
			return;
		}
		CortexSurfaceReconstruct * pia3DReconstruct = new CortexSurfaceReconstruct(parameters, piaStructure);
		piaSurface = pia3DReconstruct->surfaceReconstruction(Pia);
		
		// old surface reconstruction:
// 		piaSurface = surfaceReconstruction(Pia);
		barrelReconstruction();
		
		#ifdef DEBUG
		std::cout << "Writing Pia '.surf' file" << std::endl;
		#endif
		std::string piaName(outputFilename);
		piaName += "_pia";
		Reader * piaWriter = new Reader(piaName.c_str(), piaName.c_str());
		piaWriter->writeAmiraSurfaceFile(piaSurface);
		delete piaWriter;
	}
	else	//should never need to branch here IF used correctly
	{
		std::cout << "Error! Cannot reconstruct 3D barrels without Pia! This should not happen..." << std::endl;
	}
	if(parameters.wmFlag)
	{
		// new surface reconstruction:
		PolyDataPointerType wmStructure = PolyDataPointerType::New();
		if(!spatialGraph->extractLandmark(WhiteMatter, wmStructure))
		{
			std::cout << "Error! Pia not in SpatialGraph!" << std::endl;
			return;
		}
		CortexSurfaceReconstruct * wm3DReconstruct = new CortexSurfaceReconstruct(parameters, wmStructure);
		wmSurface = wm3DReconstruct->surfaceReconstruction(WhiteMatter);
		
		// old surface reconstruction:
// 		wmSurface = surfaceReconstruction(WhiteMatter);
		
		#ifdef DEBUG
		std::cout << "Writing WM '.surf' file" << std::endl;
		#endif
		std::string wmName(outputFilename);
		wmName += "_WM";
		Reader * wmWriter = new Reader(wmName.c_str(), wmName.c_str());
		wmWriter->writeAmiraSurfaceFile(wmSurface);
		delete wmWriter;
	}
	if(/*parameters.wmFlag &&*/ parameters.piaFlag)
	{
		#ifdef DEBUG
		std::cout << "Writing barrel column parameters" << std::endl;
		#endif
		writeBarrelParameters(outputFilename);
		
		#ifdef DEBUG
		std::cout << "Writing landmark '.am' file" << std::endl;
		#endif
		std::string barrelName(outputFilename);
		barrelName += "_barrels";
		Reader * barrelWriter = new Reader(barrelName.c_str(), barrelName.c_str());
		barrelWriter->setSpatialGraph(spatialGraph);
		barrelWriter->writeSpatialGraphFile();
		delete barrelWriter;
	}
}

void CortexReconstruction::convertSGToBarrelField()
{
	#ifdef DEBUG
	std::cout << "Converting input SpatialGraph to BarrelField" << std::endl;
	#endif
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		PolyDataPointerType barrel = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(ID, barrel))
		{
			double barrelTop[3], barrelBottom[3], pCenter1[3], pCenter2[3];
			int subID1, subID2;
			double * weights1 = new double[barrel->GetCell(0)->GetNumberOfPoints()];
			double * weights2 = new double[barrel->GetCell(1)->GetNumberOfPoints()];
			barrel->GetCell(0)->GetParametricCenter(pCenter1);
			barrel->GetCell(1)->GetParametricCenter(pCenter2);
			barrel->GetCell(0)->EvaluateLocation(subID1, pCenter1, barrelTop, weights1);
			barrel->GetCell(1)->EvaluateLocation(subID2, pCenter2, barrelBottom, weights2);
			delete [] weights1, delete [] weights2;
			
#ifdef DEBUG
			double topAvg[3] = {0,0,0}, bottomAvg[3] = {0,0,0};
			for(int ii = 0; ii < barrel->GetCell(0)->GetNumberOfPoints(); ++ii)
			{
				double topPt[3], bottomPt[3];
				barrel->GetCell(0)->GetPoints()->GetPoint(ii, topPt);
				barrel->GetCell(1)->GetPoints()->GetPoint(ii, bottomPt);
				for(int jj = 0; jj < 3; ++jj)
				{
					topAvg[jj] += topPt[jj];
					bottomAvg[jj] += bottomPt[jj];
				}
			}
			for(int ii = 0; ii < 3; ++ii)
			{
				topAvg[ii] = topAvg[ii]/barrel->GetCell(0)->GetNumberOfPoints();
				bottomAvg[ii] = bottomAvg[ii]/barrel->GetCell(0)->GetNumberOfPoints();
			}
			double pAxis[3];
			vtkMath::Subtract(barrelTop, barrelBottom, pAxis);
			vtkMath::Normalize(pAxis);
			double avgAxis[3];
			vtkMath::Subtract(topAvg, bottomAvg, avgAxis);
			vtkMath::Normalize(avgAxis);
			
			std::cout << "------------------------------" << std::endl;
			std::cout << "Barrel " << SBF->int2Labels[ID] << std::endl;
			std::cout << "pTop = [" << barrelTop[0] << "," << barrelTop[1] << "," << barrelTop[2] << "]" << std::endl;
			std::cout << "pBottom = [" << barrelBottom[0] << "," << barrelBottom[1] << "," << barrelBottom[2] << "]" << std::endl;
			std::cout << "avgTop = [" << topAvg[0] << "," << topAvg[1] << "," << topAvg[2] << "]" << std::endl;
			std::cout << "avgBottom = [" << bottomAvg[0] << "," << bottomAvg[1] << "," << bottomAvg[2] << "]" << std::endl;
			std::cout << "pAxis = [" << pAxis[0] << "," << pAxis[1] << "," << pAxis[2] << "]" << std::endl;
			std::cout << "avgAxis = [" << avgAxis[0] << "," << avgAxis[1] << "," << avgAxis[2] << "]" << std::endl;
			std::cout << "------------------------------" << std::endl;
#endif
			
			Column * newBarrel = new Column(barrel, barrelTop, barrelBottom);
			finalBarrels.insert(std::pair< int, Column * >(ID, newBarrel));
		}
	}
}

void CortexReconstruction::vesselReconstruction(const char* outputFilename)
{
// 	#ifdef DEBUG
	std::cout << "Starting 3D reconstruction!" << std::endl;
// 	#endif
		computeBloodVessels(1, 5);
		#ifdef DEBUG
		std::cout << "Writing landmark '.am' file" << std::endl;
		#endif
		std::string barrelName(outputFilename);
		barrelName += "_3Dvessels";
		Reader * barrelWriter = new Reader(barrelName.c_str(), barrelName.c_str());
		barrelWriter->setSpatialGraph(spatialGraph);
		barrelWriter->writeSpatialGraphFile();
		delete barrelWriter;
}

/****************************************************************************/
/* private methods                                                          */
/****************************************************************************/

void CortexReconstruction::writeBarrelParameters ( const char* ofName )
{
	std::string outFile(ofName);
	outFile += "_parameters.csv";
	std::ofstream ParameterFile(outFile.c_str());
	ParameterFile << "Barrel\tbottom\ttop\tbarrel height\tcolumn height\tarea" << std::endl;
	
	Surface * tmpPia = new Surface(piaSurface);
	Surface * tmpWM;
	if(parameters.wmFlag)
		tmpWM = new Surface(wmSurface);
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		ParameterFile << SBF->int2Labels[ID];
		if(finalBarrels.find(ID) != finalBarrels.end())
		{
			double height, area, topPiaDist, bottomPiaDist, piaWMDist;
			double topPt[3], bottomPt[3], axis[3], center[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				topPt[ii] = finalBarrels[ID]->top[ii];
				bottomPt[ii] = finalBarrels[ID]->bottom[ii];
				axis[ii] = topPt[ii] - bottomPt[ii];
				center[ii] = 0.5*(topPt[ii] + bottomPt[ii]);
			}
			vtkMath::Normalize(axis);
			
			height = finalBarrels[ID]->getHeight();
			
			PolyDataPointerType avgContours = finalBarrels[ID]->contours;
			double normal[3];
			// don't use GetCell(0)->GetPoints b/c points in PointIDList are not sorted in strictly ascending order for overlap-corrected contours!
			area = vtkPolygon::ComputeArea(avgContours->GetPoints(), avgContours->GetCell(0)->GetNumberOfPoints(), avgContours->GetCell(0)->GetPointIds()->GetPointer(0), normal);
			
			double piaIntersection[3];
			tmpPia->intersectLine(axis, center);
			if(tmpPia->isValid())
			{
				tmpPia->getLastIntersectPoint(piaIntersection);
				bottomPiaDist = sqrt(vtkMath::Distance2BetweenPoints(piaIntersection, bottomPt));
				topPiaDist = sqrt(vtkMath::Distance2BetweenPoints(piaIntersection, topPt));
			}
			else
			{
				topPiaDist = -1;
				bottomPiaDist = -1;
				piaWMDist = -1;
			}
			
			if(parameters.wmFlag)
			{
				double wmIntersection[3];
				tmpWM->intersectLine(axis, center);
				if(tmpWM->isValid() && tmpPia->isValid())
				{
					tmpWM->getLastIntersectPoint(wmIntersection);
					piaWMDist = sqrt(vtkMath::Distance2BetweenPoints(piaIntersection, wmIntersection));
					std::vector< double * > topBottomPts;
					topBottomPts.push_back(piaIntersection);
					topBottomPts.push_back(wmIntersection);
					PolyDataPointerType oneContour = PolyDataPointerType::New();
					IdListPointerType IDList = IdListPointerType::New();
					oneContour->Allocate();
					IDList->InsertId(0, 0);
					oneContour->CopyCells(avgContours, IDList);
					spatialGraph->addPolyDataObject(createTopBottomContours(oneContour, topBottomPts), ID);
				}
			}
			else
				piaWMDist = -1;
			
			ParameterFile << "\t" << bottomPiaDist << "\t" << topPiaDist << "\t";
			ParameterFile << height << "\t" << piaWMDist << "\t" << area;
		}
		ParameterFile << std::endl;
	}
	
	ParameterFile.close();
	delete tmpPia, delete tmpWM;
}

void CortexReconstruction::barrelReconstruction()
{
	double alpha = 0.5;
	#ifdef REG_ACCURACY
	alpha = var_alpha;
	#endif
// 	std::cout << "alpha = ";
// 	std::cin >> alpha;
	std::map< int, PolyDataPointerType > barrels;
	std::map< int, double * > barrelAxes;
	std::map< int, double * > barrelCenters;
	std::list< unsigned int > goodVessels;
	#ifdef USEVESSELS
	goodVessels = computeConstrainingVessels(piaSurface);
	#endif
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		#ifdef DEBUG
		std::flush(std::cout << "Reconstructing 3D barrel " << SBF->int2Labels[ID] << std::endl);
		#endif
		PolyDataPointerType barrel = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(ID, barrel))
		{
			// reconstruct barrel centroid and barrel axis
			double * barrelCentroid = new double[3], * barrelAxis = new double[3];
			calculateBarrelCentroid(barrel, barrelCentroid);
			#ifdef DEBUG
			std::cout << "Barrel centroid @ [" << barrelCentroid[0] << "," << barrelCentroid[1] << "," << barrelCentroid[2] << "]" << std::endl;
			#endif
			if(parameters.piaFlag)
			{
				newBarrelAxis(barrel, piaSurface, goodVessels, alpha, barrelAxis);
				vtkMath::Normalize(barrelAxis);
			}
			else
				barrelAxis[0] = 0, barrelAxis[1] = 0, barrelAxis[2] = 1;
			#ifdef DEBUG
			std::cout << "Barrel axis = [" << barrelAxis[0] << "," << barrelAxis[1] << "," << barrelAxis[2] << "]" << std::endl;
			#endif
			
			#ifdef MAXCONTOUR
			#ifdef DEBUG
			std::cout << "Sampling barrel contours at equal angular intervals" << std::endl;
			#endif
			PolyDataPointerType sampledBarrel = sampleBarrelContour(barrel);
			barrel->DeepCopy(sampledBarrel);
			#endif
			
			barrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, barrel));
			barrelCenters.insert(std::pair< int, double * >(*labelIt, barrelCentroid));
			barrelAxes.insert(std::pair< int, double * >(*labelIt, barrelAxis));
		}
	}
	// apply divergence constraint on barrel axis vector field
	// to avoid systematic errors that may result from
	// pia reconstruction
	enforceAxisDivergence(barrelAxes, barrelCenters);
	
	#ifdef MAXCONTOUR
	// calculate avg barrel contours and ensure mutual non-overlap
	// of those in barrel field (they can overlap in deeper layers...)
	std::map< int, std::vector< double * > > endPointMap;
	std::map< int, PolyDataPointerType > maxBarrelContours;
	computeMaxBarrelContours(barrels, barrelAxes, barrelCenters, maxBarrelContours, endPointMap);
	#endif
	
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels.find(ID) != barrels.end())
		{
			#ifdef AVGCONTOUR
			//reconstruct avg barrel contour
			PolyDataPointerType avgContour = computeAverageBarrelContour(barrels[ID], barrelCenters[ID], barrelAxes[ID]);
			// reconstruct barrel top/bottom points
			std::vector< double * > endPts;
			closeBarrelAlongNewAxis(barrelAxes[ID], barrelCenters[ID], barrels[ID], endPts);
			Column * newBarrel = new Column(createTopBottomContours(avgContour, endPts), endPts[0], endPts[1]);
			#endif
			#ifdef MAXCONTOUR
			Column * newBarrel = new Column(maxBarrelContours[ID], endPointMap[ID][0], endPointMap[ID][1]);
			#endif
			
			finalBarrels.insert(std::pair< int, Column * >(ID, newBarrel));
			
			#ifdef DEBUG
			double boundingBox[6], indirectAxis[3];
			for(int ii = 0; ii < 3; ++ii)
				indirectAxis[ii] = finalBarrels[ID]->top[ii] - finalBarrels[ID]->bottom[ii];
			vtkMath::Normalize(indirectAxis);
			finalBarrels[ID]->contours->GetBounds(boundingBox);
			std::cout << "*************************************" << std::endl;
			std::cout << "Barrel " << SBF->int2Labels[ID] << ":" << std::endl;
			std::cout << "Top @ [" << finalBarrels[ID]->top[0] << "," << finalBarrels[ID]->top[1] << "," << finalBarrels[ID]->top[2] << "]" << std::endl;
			std::cout << "Bottom @ [" << finalBarrels[ID]->bottom[0] << "," << finalBarrels[ID]->bottom[1] << "," << finalBarrels[ID]->bottom[2] << "]" << std::endl;
			std::cout << "Center @ [" << barrelCenters[ID][0] << "," << barrelCenters[ID][1] << "," << barrelCenters[ID][2] << "]" << std::endl;
			std::cout << "Axis = [" << barrelAxes[ID][0] << "," << barrelAxes[ID][1] << "," << barrelAxes[ID][2] << "]" << std::endl;
			std::cout << "Indirect axis = [" << indirectAxis[0] << "," << indirectAxis[1] << "," << indirectAxis[2] << "]" << std::endl;
			std::cout << "Bounding box:" << std::endl;
			std::cout << "[" << boundingBox[0] << "," << boundingBox[1] << "]" << std::endl;
			std::cout << "[" << boundingBox[2] << "," << boundingBox[3] << "]" << std::endl;
			std::cout << "[" << boundingBox[4] << "," << boundingBox[5] << "]" << std::endl;
			#endif
			
			delete [] barrelAxes[ID], delete [] barrelCenters[ID];
		}
	}
}

PolyDataPointerType CortexReconstruction::computeAverageBarrelContour ( PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3] )
{
	int nrAngles = 36;
	std::vector< std::vector< double * > > thisPlaneSamplingPoints;	// one vector for each sampling angle
	for(int ii = 0; ii < nrAngles; ++ii)
	{
		std::vector< double * > angleSamplingVec;
		thisPlaneSamplingPoints.push_back(angleSamplingVec);
	}
// 	std::cout << std::endl;
// 	std::cout << "nr. of home barrel contours: " << completeBarrel->GetNumberOfCells() << std::endl;
	for(int ii = 0; ii < completeBarrel->GetNumberOfCells(); ++ii)
	{
		double centerPoint[3];
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[completeBarrel->GetCell(ii)->GetNumberOfPoints()];
		completeBarrel->GetCell(ii)->GetParametricCenter(pCoords);
		completeBarrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		delete [] weights;
		
		// identify polygon normals so each contour can be traversed in same direction
		double normal[3], rotAxis[3], zUnitVec[3] = {0,0,parameters.zReversed ? -1 : 1};
		double rotAngle = 0;
		vtkPolygon::ComputeNormal(completeBarrel->GetCell(ii)->GetPoints(), normal);
		bool reverseDirection;
		if(normal[2] < 0)
		{
			reverseDirection = 1;
			zUnitVec[2] *= -1;
		}
		else
			reverseDirection = 0;
		
		// rotate sampling plane such that the 
		// orientation is always the same
		// independent of polygon normal direction
		vtkMath::Cross(normal, zUnitVec, rotAxis);
// 		vtkMath::Cross(barrelAxis, zUnitVec, rotAxis);
		vtkMath::Normalize(rotAxis);
		rotAngle = std::acos(vtkMath::Dot(normal, zUnitVec))*180/PI;
// 		rotAngle = std::acos(vtkMath::Dot(barrelAxis, zUnitVec))*180/PI;
		TransformPointerType rot = TransformPointerType::New();
		rot->RotateWXYZ(rotAngle, rotAxis);
		
// 		std::cout << std::endl;
// 		std::cout << "home barrel cell " << ii << std::endl;
// 		std::cout << "normal = [" << normal[0] << "," << normal[1] << "," << normal[2] << "]" << std::endl;
// 		std::cout << "rotAxis = [" << rotAxis[0] << "," << rotAxis[1] << "," << rotAxis[2] << "]" << std::endl;
// 		std::cout << "zUnitVec = [" << zUnitVec[0] << "," << zUnitVec[1] << "," << zUnitVec[2] << "]" << std::endl;
// 		std::cout << "rotAngle = " << rotAngle << std::endl;
		
		PlanePointerType thisPlane = PlanePointerType::New();
		PointsPointerType thisCellPoints = completeBarrel->GetCell(ii)->GetPoints();
		for(int jj = 0; jj < nrAngles; ++jj)	// has to be unique b/c order of the points is not clear
		{
// 			std::cout << std::endl;
// 			std::cout << "sampling angle " << jj << std::endl;
			
			double angle = jj*2*PI/nrAngles;
			double planeNormal[3], hPlaneNormal[4];
			if(reverseDirection)
			{
				hPlaneNormal[0] = -sin(angle);
				hPlaneNormal[1] = -cos(angle);
			}
			else
			{
				hPlaneNormal[0] = sin(angle);
				hPlaneNormal[1] = cos(angle);
			}
			hPlaneNormal[2] = 0, hPlaneNormal[3] = 1;
			rot->MultiplyPoint(hPlaneNormal, hPlaneNormal);
			planeNormal[0] = hPlaneNormal[0];
			planeNormal[1] = hPlaneNormal[1];
			planeNormal[2] = hPlaneNormal[2];
			
			thisPlane->SetOrigin(centerPoint);
// 			thisPlane->SetOrigin(barrelCentroid);
			thisPlane->SetNormal(planeNormal);
			
			double * firstPt = new double[3];
			thisCellPoints->GetPoint(0, firstPt);
			double lastVal = thisPlane->EvaluateFunction(firstPt);
			delete [] firstPt;
			for(int kk = 1; kk <= thisCellPoints->GetNumberOfPoints(); ++kk)
			{
// 				std::cout << "cell " << ii << " pt " << kk << std::endl;
				
				double * thisPt = new double[3];
				thisCellPoints->GetPoint(kk%thisCellPoints->GetNumberOfPoints(), thisPt);
				double ptVal = thisPlane->EvaluateFunction(thisPt);
				if((ptVal < 0 && lastVal > 0))	// has to be unique b/c order of the points is not clear
				{
					// regular case
					double * samplePt =  new double[3];
					double * lastPt = new double[3];
					thisCellPoints->GetPoint(kk-1, lastPt);
					double direction[3];
					double normD = 0;
					for(int ll = 0; ll < 3; ++ll)
					{
						direction[ll] = lastPt[ll] - thisPt[ll];
						normD += direction[ll]*direction[ll];
					}
					normD = sqrt(normD);
					if(normD)
					{
						double distAlongLine = std::abs(ptVal);
						double correctionAngle = std::abs(direction[0]*planeNormal[0]/normD + direction[1]*planeNormal[1]/normD);
						if(correctionAngle)
							distAlongLine /= correctionAngle;
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll] + direction[ll]*distAlongLine/normD;
					}
					else
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll];
						thisPlaneSamplingPoints[jj].push_back(samplePt);
					delete [] lastPt;
// 					std::cout << "regular order" << std::endl;
// 					std::cout << "pt " << kk << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
				}
				else if(ptVal == 0 && lastVal > 0)
				{
					// case: this point is directly on plane
					double * samplePt =  new double[3];
					for(int ll = 0; ll < 3; ++ll)
						samplePt[ll] = thisPt[ll];
					thisPlaneSamplingPoints[jj].push_back(samplePt);
// 					std::cout << "pt on plane" << std::endl;
// 					std::cout << "pt " << kk << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
				}
// 				else if(ptVal != 0 && lastVal == 0)
// 				{
// 					;// do nothing
// 				}
// 				else if(ptVal == 0 && lastVal == 0)
// 				{
// 					;// either an error or both points lie directly on the plane (extremely unlikely)
// 				}
				lastVal = ptVal;
				delete [] thisPt;
			} // all cell points
		} // all angles
	} // all barrel cells
	
	PolyDataPointerType avgHBContour = PolyDataPointerType::New();
	PointsPointerType contourPts = PointsPointerType::New();
	PolygonPointerType contourPoly = PolygonPointerType::New();
	avgHBContour->Allocate(1);
	contourPts->SetDataTypeToFloat();
	contourPts->SetNumberOfPoints(nrAngles);
	contourPoly->GetPointIds()->SetNumberOfIds(nrAngles);
	
	for(int ii = 0; ii < nrAngles; ++ii)
	{
		if(thisPlaneSamplingPoints[ii].size())
		{
// 			std::cout << std::endl;
// 			std::cout << "sampling angle " << ii << std::endl;
			
			double * avgPt = new double[3];
			avgPt[0] = 0, avgPt[1] = 0, avgPt[2] = 0;
			for(int jj = 0; jj < thisPlaneSamplingPoints[ii].size(); ++jj)
			{
				double * tmpPt = thisPlaneSamplingPoints[ii][jj];
				for(int kk = 0; kk < 3; ++kk)
					avgPt[kk] += tmpPt[kk];
				
// 				std::cout << "contour " << jj << std::endl;
// 				std::cout << "pt = [" << tmpPt[0] << "," << tmpPt[1] << "," << tmpPt[2] << "]" << std::endl;
			}
			for(int jj = 0; jj < 3; ++jj)
				avgPt[jj] = avgPt[jj]/double(thisPlaneSamplingPoints[ii].size());
			
			contourPts->InsertPoint(ii, avgPt);
			contourPoly->GetPointIds()->SetId(ii, ii);
		}
	}
	avgHBContour->InsertNextCell(contourPoly->GetCellType(), contourPoly->GetPointIds());
	avgHBContour->SetPoints(contourPts);
	avgHBContour->Update();
	
	#ifdef DEBUG
	if(spatialGraph)
		spatialGraph->addPolyDataObject(avgHBContour, Barrel);
	#endif
	
	return avgHBContour;
}

PolyDataPointerType CortexReconstruction::sampleBarrelContour ( PolyDataPointerType completeBarrel )
{
	int nrAngles = 36;
	std::vector< std::vector< double * > > contourSamplingPoints;	// one vector for each contour
	for(int ii = 0; ii < completeBarrel->GetNumberOfCells(); ++ii)
	{
		std::vector< double * > contourSamplingVec;
		contourSamplingPoints.push_back(contourSamplingVec);
	}
	#ifdef DEBUG
	std::cout << "Nr. of barrel contours: " << completeBarrel->GetNumberOfCells() << std::endl;
	#endif
	for(int ii = 0; ii < completeBarrel->GetNumberOfCells(); ++ii)
	{
		#ifdef DEBUG
		std::cout << "\tSampling contour #" << ii << std::endl;
		#endif
		double centerPoint[3];
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[completeBarrel->GetCell(ii)->GetNumberOfPoints()];
		completeBarrel->GetCell(ii)->GetParametricCenter(pCoords);
		completeBarrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		delete [] weights;
		
		// identify polygon normals so each contour can be traversed in same direction
		double normal[3], rotAxis[3], zUnitVec[3] = {0,0,parameters.zReversed ? -1 : 1};
// 		double rotAngle = 0;
		vtkPolygon::ComputeNormal(completeBarrel->GetCell(ii)->GetPoints(), normal);
		bool reverseDirection;
		if(normal[2] < 0)
		{
			reverseDirection = 1;
			zUnitVec[2] *= -1;
		}
		else
			reverseDirection = 0;
		
		#ifdef DEBUG
		std::cout << "\tCenter @ [" << centerPoint[0] << "," << centerPoint[1] << "," << centerPoint[2] << "]" << std::endl;
		std::cout << "\tnormal = [" << normal[0] << "," << normal[1] << "," << normal[2] << "]" << std::endl;
		std::cout << "\treverseDirection = " << reverseDirection << std::endl;
		#endif
		
		// rotate sampling plane such that the 
		// orientation is always the same
		// independent of polygon normal direction
// 		vtkMath::Cross(normal, zUnitVec, rotAxis);
// 		vtkMath::Normalize(rotAxis);
// 		rotAngle = std::acos(vtkMath::Dot(normal, zUnitVec))*180/PI;
// 		TransformPointerType rot = TransformPointerType::New();
// 		rot->RotateWXYZ(rotAngle, rotAxis);
		
		PlanePointerType thisPlane = PlanePointerType::New();
		PointsPointerType thisCellPoints = completeBarrel->GetCell(ii)->GetPoints();
		for(int jj = 0; jj < nrAngles; ++jj)	// has to be unique b/c order of the points is not clear
		{
			double angle = jj*2*PI/nrAngles;
			double planeNormal[3], hPlaneNormal[4];
// 			if(reverseDirection)
// 			{
// 				hPlaneNormal[0] = -sin(angle);
// 				hPlaneNormal[1] = -cos(angle);
// 			}
// 			else
// 			{
// 				hPlaneNormal[0] = sin(angle);
// 				hPlaneNormal[1] = cos(angle);
// 			}
// 			hPlaneNormal[2] = 0, hPlaneNormal[3] = 1;
// 			rot->MultiplyPoint(hPlaneNormal, hPlaneNormal);
// 			planeNormal[0] = hPlaneNormal[0];
// 			planeNormal[1] = hPlaneNormal[1];
// 			planeNormal[2] = hPlaneNormal[2];
			if(reverseDirection)
			{
				planeNormal[0] = -sin(angle);
				planeNormal[1] = -cos(angle);
				planeNormal[2] = 0;
			}
			else
			{
				planeNormal[0] = sin(angle);
				planeNormal[1] = cos(angle);
				planeNormal[2] = 0;
			}
			
			thisPlane->SetOrigin(centerPoint);
// 			thisPlane->SetOrigin(barrelCentroid);
			thisPlane->SetNormal(planeNormal);
			
			double firstPt[3];
			thisCellPoints->GetPoint(0, firstPt);
			double lastVal = thisPlane->EvaluateFunction(firstPt);
			for(int kk = 1; kk <= thisCellPoints->GetNumberOfPoints(); ++kk)
			{
// 				#ifdef DEBUG
// 				std::flush(std::cout << "\tcell " << ii << " pt " << kk << " @ [");
// 				#endif
				
				double thisPt[3];
				thisCellPoints->GetPoint(kk%thisCellPoints->GetNumberOfPoints(), thisPt);
				double ptVal = thisPlane->EvaluateFunction(thisPt);
				
// 				#ifdef DEBUG
// 				std::cout << thisPt[0] << "," << thisPt[1] << "," << thisPt[2] << "]" << std::endl;
// 				std::cout << "\tptVal = " << ptVal << " ---  lastVal = " << lastVal << std::endl;
// 				#endif
				
				if((ptVal < 0 && lastVal > 0))	// has to be unique b/c order of the points is not clear
				{
					// regular case
					double * samplePt =  new double[3];
					double lastPt[3];
					thisCellPoints->GetPoint(kk-1, lastPt);
					double direction[3];
					double normD = 0;
					for(int ll = 0; ll < 3; ++ll)
					{
						direction[ll] = lastPt[ll] - thisPt[ll];
						normD += direction[ll]*direction[ll];
					}
					normD = sqrt(normD);
					if(normD)
					{
						double distAlongLine = std::abs(ptVal);
						double correctionAngle = std::abs(direction[0]*planeNormal[0]/normD + direction[1]*planeNormal[1]/normD);
						if(correctionAngle)
							distAlongLine /= correctionAngle;
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll] + direction[ll]*distAlongLine/normD;
					}
					else
					{
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll];
					}
					contourSamplingPoints[ii].push_back(samplePt);
					#ifdef DEBUG
					std::cout << "\tRegular order!" << std::endl;
					std::cout << "\tSampled pt " << contourSamplingPoints[ii].size() << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
					#endif
				}
				else if(ptVal == 0 && lastVal > 0)
				{
					// case: this point is directly on plane
					double * samplePt =  new double[3];
					for(int ll = 0; ll < 3; ++ll)
						samplePt[ll] = thisPt[ll];
					contourSamplingPoints[ii].push_back(samplePt);
					#ifdef DEBUG
					std::cout << "\tPoint on plane!" << std::endl;
					std::cout << "\tSampled pt " << contourSamplingPoints[ii].size() << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
					#endif
				}
// 				else if(ptVal != 0 && lastVal == 0)
// 				{
// 					;// do nothing
// 				}
// 				else if(ptVal == 0 && lastVal == 0)
// 				{
// 					;// either an error or both points lie directly on the plane (extremely unlikely)
// 				}
				lastVal = ptVal;
			} // all cell points
		} // all angles
	} // all barrel cells
	
	#ifdef DEBUG
	std::cout << "\tCreating sampled barrel contours!" << std::endl;
	#endif
	
	PolyDataPointerType avgHBContour = PolyDataPointerType::New();
	PointsPointerType contourPts = PointsPointerType::New();
	avgHBContour->Allocate(1);
	contourPts->SetDataTypeToFloat();
	contourPts->SetNumberOfPoints(nrAngles*completeBarrel->GetNumberOfCells());
	for(int ii = 0; ii < completeBarrel->GetNumberOfCells(); ++ii)
	{
		#ifdef DEBUG
		std::flush(std::cout << "\tCell " << ii << " # of sampled points = " << contourSamplingPoints[ii].size() << std::endl);
		#endif
		PolygonPointerType contourPoly = PolygonPointerType::New();
		contourPoly->GetPointIds()->SetNumberOfIds(nrAngles);
		for(int jj = 0; jj < nrAngles && jj < contourSamplingPoints[ii].size(); ++jj)
		{
			contourPts->InsertPoint(ii*nrAngles + jj, contourSamplingPoints[ii][jj]);
			contourPoly->GetPointIds()->SetId(jj, ii*nrAngles + jj);
		}
		avgHBContour->InsertNextCell(contourPoly->GetCellType(), contourPoly->GetPointIds());
	}
	avgHBContour->SetPoints(contourPts);
	avgHBContour->Update();
	
	return avgHBContour;
}

void CortexReconstruction::computeMaxBarrelContours ( std::map< int, PolyDataPointerType > barrels, std::map< int, double* > barrelAxes, std::map< int, double* > barrelCenters, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double* > >& endPointMap )
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		if(barrels.find(*labelIt) != barrels.end())
		{
			std::vector< double * > endPoints;
			closeBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPoints);
			endPointMap.insert(std::pair< int, std::vector< double * > >(*labelIt, endPoints));
			avgBarrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, smoothBarrelAlongAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPointMap[*labelIt])));
			#ifdef DEBUG
			double contourNormal[3];
			vtkPolygon::ComputeNormal(avgBarrels[*labelIt]->GetCell(0)->GetPoints(), contourNormal);
			std::cout << "Barrel " << SBF->int2Labels[*labelIt] << " max contour normal = [";
			std::cout << contourNormal[0] << "," << contourNormal[1] << "," << contourNormal[2] << "]" << std::endl;
			#endif
		}
	
	enforceOverlapConstraint(barrelAxes, avgBarrels, endPointMap);
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgBarrels.find(ID) != avgBarrels.end() && spatialGraph)
		{
			spatialGraph->addPolyDataObject(avgBarrels[ID], Barrel);
			
			// update barrel top/bottom points for consistent calculations
			// later on when loading reconstructed barrels/columns
			// from Amira files and computing top/bottom points from contours
			PolyDataPointerType barrel = avgBarrels[ID];
			double barrelTop[3], barrelBottom[3], pCenter1[3], pCenter2[3];
			int subID1, subID2;
			double * weights1 = new double[barrel->GetCell(0)->GetNumberOfPoints()];
			double * weights2 = new double[barrel->GetCell(1)->GetNumberOfPoints()];
			barrel->GetCell(0)->GetParametricCenter(pCenter1);
			barrel->GetCell(1)->GetParametricCenter(pCenter2);
			barrel->GetCell(0)->EvaluateLocation(subID1, pCenter1, barrelTop, weights1);
			barrel->GetCell(1)->EvaluateLocation(subID2, pCenter2, barrelBottom, weights2);
			for(int ii = 0; ii < 3; ++ii)
			{
				endPointMap[ID][0][ii] = barrelTop[ii];
				endPointMap[ID][1][ii] = barrelBottom[ii];
			}
			delete [] weights1, delete [] weights2;
			
			#ifdef DEBUG
			spatialGraph->addLine(endPointMap[ID][0], endPointMap[ID][1], Barrel);
			#endif
		}
		#ifdef DEBUG
		double contourNormal[3];
		vtkPolygon::ComputeNormal(avgBarrels[*labelIt]->GetCell(0)->GetPoints(), contourNormal);
		std::cout << "Barrel " << SBF->int2Labels[*labelIt] << " corrected contour normal = [";
		std::cout << contourNormal[0] << "," << contourNormal[1] << "," << contourNormal[2] << "]" << std::endl;
		#endif
	}
}

PolyDataPointerType CortexReconstruction::smoothBarrelAlongAxis ( double newAxis[3], double barrelCentroid[3], PolyDataPointerType barrel, std::vector< double * > endPoints )
{
	if(barrel->GetNumberOfPoints())
	{
// 		barrel->Print(std::cout);
// 		std::flush(std::cout << "Calculate avg barrel contour" << std::endl);
		if(barrel->GetCell(0)->GetNumberOfPoints() != 36)
		{
			std::cout << "Warning! Barrel has not been sampled correctly! max barrel contour may be corrupted..." << std::endl;
		}
		double zAxis[3], barrelCenter[3];
		zAxis[0] = newAxis[0], zAxis[1] = newAxis[1], zAxis[2] = newAxis[2];
		vtkMath::Normalize(zAxis);
		double axisPt1[3], axisPt2[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			barrelCenter[ii] = 0.5*(endPoints[0][ii] + endPoints[1][ii]);
			axisPt1[ii] = barrelCenter[ii] + 1000*zAxis[ii];
			axisPt2[ii] = barrelCenter[ii] - 1000*zAxis[ii];
		}
		std::vector< std::list< double * > > avgPoints;
		std::vector< std::list< double > > avgDistances;
		for(int ii = 0; ii < 36; ++ii)
		{
			std::list< double * > ptList;
			avgPoints.push_back(ptList);
			std::list< double > distList;
			avgDistances.push_back(distList);
		}
		
		// sort all points @ each sampling angle by their distance to barrel axis
		// insert vector pointing from axis to that point for later convenience
		for(int ii = 0; ii < barrel->GetNumberOfCells(); ++ii)
		{
			PointsPointerType cellPts = barrel->GetCell(ii)->GetPoints();
			for(int jj = 0; jj < 36 && jj < cellPts->GetNumberOfPoints(); ++jj)
			{
				double t;
				double * closestPt = new double[3];
				double pt[3];
				cellPts->GetPoint(jj, pt);
				double dist = vtkLine::DistanceToLine(pt, axisPt1, axisPt2, t, closestPt);
				dist = sqrt(dist);
				for(int kk = 0; kk < 3; ++kk)
					closestPt[kk] = pt[kk] - closestPt[kk];
				avgPoints[jj].push_back(closestPt);
				avgDistances[jj].push_back(dist);
			}
		}
		
// 		for(int ii = 0; ii < 36; ++ii)
// 			std::flush(std::cout << "avgPoints[" << ii << "].size() = " << avgPoints[ii].size() << std::endl);
		// create top & bottom max contours @ extreme points 
		PolyDataPointerType avgContours = PolyDataPointerType::New();
		PointsPointerType avgContourPoints = PointsPointerType::New();
		PolygonPointerType topPoly = PolygonPointerType::New();
		PolygonPointerType bottomPoly = PolygonPointerType::New();
		avgContours->Allocate();
		avgContourPoints->SetDataTypeToFloat();
		topPoly->GetPointIds()->SetNumberOfIds(36);
		bottomPoly->GetPointIds()->SetNumberOfIds(36);
		for(int ii = 0; ii < 36; ++ii)
		{
			if(avgPoints[ii].size())
			{
				double * topPt = new double[3];
// 				double avgVec[] = {0, 0, 0};
				double * avgVec = new double[3];
				avgVec[0] = 0, avgVec[1] = 0, avgVec[2] = 0;
				topPt[0] = endPoints[0][0], topPt[1] = endPoints[0][1], topPt[2] = endPoints[0][2];
				std::list< double * >::const_iterator avgPointListIt;
				for(avgPointListIt = avgPoints[ii].begin(); avgPointListIt != avgPoints[ii].end(); ++avgPointListIt)
				{
					for(int jj = 0; jj < 3; ++jj)
						avgVec[jj] += (*avgPointListIt)[jj];
				}
// 				for(int jj = 0; jj < 3; ++jj)
// 				{
// 					avgVec[jj] = avgVec[jj]/double(avgPoints[ii].size());
// 					topPt[jj] += avgVec[jj];
// 				}
				//new version:
				vtkMath::Normalize(avgVec);
// 				std::flush(std::cout << "avgVec = " << "[" << avgVec[0] << "," << avgVec[1] << "," << avgVec[2] << "]" << std::endl);
				double avgRadius = 0, stdRadius = 0;
				std::list< double >::const_iterator radiusIt;
				for(radiusIt = avgDistances[ii].begin(); radiusIt != avgDistances[ii].end(); ++radiusIt)
					avgRadius += *radiusIt;
				if(avgDistances[ii].size())
					avgRadius = avgRadius/(double)avgDistances[ii].size();
				for(radiusIt = avgDistances[ii].begin(); radiusIt != avgDistances[ii].end(); ++radiusIt)
					stdRadius += (*radiusIt - avgRadius)*(*radiusIt - avgRadius);
				if(avgDistances[ii].size())
					stdRadius = stdRadius/(double)avgDistances[ii].size();
				stdRadius = sqrt(stdRadius);
// 				std::flush(std::cout << "this angle avg dist = " << avgRadius << " +- " << stdRadius << std::endl);
				for(int jj = 0; jj < 3; ++jj)
				{
					topPt[jj] += (avgRadius + 1.0*stdRadius)*avgVec[jj];
				}
				
				avgContourPoints->InsertNextPoint(topPt);
				topPoly->GetPointIds()->SetId(ii, ii);
			}
		}
		for(int ii = 0; ii < 36; ++ii)
		{
			if(avgPoints[ii].size())
			{
				double * bottomPt = new double[3];
// 				double avgVec[] = {0, 0, 0};
				double * avgVec = new double[3];
				avgVec[0] = 0, avgVec[1] = 0, avgVec[2] = 0;
				bottomPt[0] = endPoints[1][0], bottomPt[1] = endPoints[1][1], bottomPt[2] = endPoints[1][2];
				std::list< double * >::const_iterator avgPointListIt;
				for(avgPointListIt = avgPoints[ii].begin(); avgPointListIt != avgPoints[ii].end(); ++avgPointListIt)
				{
					for(int jj = 0; jj < 3; ++jj)
						avgVec[jj] += (*avgPointListIt)[jj];
				}
// 				for(int jj = 0; jj < 3; ++jj)
// 				{
// 					avgVec[jj] = avgVec[jj]/double(avgPoints[ii].size());
// 					bottomPt[jj] += avgVec[jj];
// 				}
				//new version:
				vtkMath::Normalize(avgVec);
				double avgRadius = 0, stdRadius = 0;
				std::list< double >::const_iterator radiusIt;
				for(radiusIt = avgDistances[ii].begin(); radiusIt != avgDistances[ii].end(); ++radiusIt)
					avgRadius += *radiusIt;
				if(avgDistances[ii].size())
					avgRadius = avgRadius/(double)avgDistances[ii].size();
				for(radiusIt = avgDistances[ii].begin(); radiusIt != avgDistances[ii].end(); ++radiusIt)
					stdRadius += (*radiusIt - avgRadius)*(*radiusIt - avgRadius);
				if(avgDistances[ii].size())
					stdRadius = stdRadius/(double)avgDistances[ii].size();
				stdRadius = sqrt(stdRadius);
				for(int jj = 0; jj < 3; ++jj)
				{
					bottomPt[jj] += (avgRadius + 1.0*stdRadius)*avgVec[jj];
				}
				
				avgContourPoints->InsertNextPoint(bottomPt);
				bottomPoly->GetPointIds()->SetId(ii, ii + 36);
			}
		}
		avgContours->InsertNextCell(topPoly->GetCellType(), topPoly->GetPointIds());
		avgContours->InsertNextCell(bottomPoly->GetCellType(), bottomPoly->GetPointIds());
		avgContours->SetPoints(avgContourPoints);
		avgContours->Update();
		
		return avgContours;
	}
	else
	{
		PolyDataPointerType avgContours = PolyDataPointerType::New();
		avgContours->Allocate();
		return avgContours;
	}
}

void CortexReconstruction::enforceOverlapConstraint ( std::map< int, double* > barrelAxes, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double* > >& endPointMap )
{
	std::map< int, std::list< int > > barrelGrid = createBarrelGrid(barrelAxes);
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			std::list< int >::const_iterator neighborIt;
			for(neighborIt = barrelGrid[*labelIt].begin(); neighborIt != barrelGrid[*labelIt].end(); ++neighborIt)
			{
				int centerBarrel = *labelIt;
				int neighborBarrel = *neighborIt;
				double centerAxis[3], neighborAxis[3];
				double centerBottom[3], neighborBottom[3];
				double projectedNeighborBottom[3];
				double centerTop[3], neighborTop[3];
				for(int ii = 0; ii < 3; ++ii)
				{
					centerAxis[ii] = barrelAxes[centerBarrel][ii];
					neighborAxis[ii] = barrelAxes[neighborBarrel][ii];
					centerBottom[ii] = endPointMap[centerBarrel][1][ii];
					neighborBottom[ii] = endPointMap[neighborBarrel][1][ii];
					projectedNeighborBottom[ii] = endPointMap[neighborBarrel][1][ii];
					centerTop[ii] = endPointMap[centerBarrel][0][ii];
					neighborTop[ii] = endPointMap[neighborBarrel][0][ii];
				}
				
				// barrel axes point upwards towards Pia by construction
				// project neighbor pts on bottom plane of center barrel
				vtkMath::Normalize(centerAxis);
				vtkMath::Normalize(neighborAxis);
				double angle = 0;
				for(int ii = 0; ii < 3; ++ii)
					angle += neighborAxis[ii]*centerAxis[ii];
				PointsPointerType neighborPoints = avgBarrels[neighborBarrel]->GetCell(1)->GetPoints();
				PolyDataPointerType projectedPoly = PolyDataPointerType::New();
				PointsPointerType projectedPoints = PointsPointerType::New();
				PolygonPointerType poly = PolygonPointerType::New();
				projectedPoly->Allocate(1);
				projectedPoints->SetDataTypeToFloat();
				poly->GetPointIds()->SetNumberOfIds(neighborPoints->GetNumberOfPoints());
				PlanePointerType centerPlane = PlanePointerType::New();
				centerPlane->SetNormal(centerAxis);
				centerPlane->SetOrigin(centerBottom);
				for(int ii = 0; ii < neighborPoints->GetNumberOfPoints(); ++ii)
				{
					double pt[3];
					neighborPoints->GetPoint(ii, pt);
					double dist = centerPlane->EvaluateFunction(pt);	// distance has to be signed -> is pt above/below plane?
					if(angle)
						dist /= angle;
					for(int jj = 0; jj < 3; ++jj)
						pt[jj] -= dist*neighborAxis[jj];
					poly->GetPointIds()->SetId(ii, ii);
					projectedPoints->InsertNextPoint(pt);
					if(!ii)
					{
						double centerDist = centerPlane->EvaluateFunction(projectedNeighborBottom);
						if(angle)
							centerDist /= angle;
						for(int jj = 0; jj < 3; ++jj)
							projectedNeighborBottom[jj] -= centerDist*neighborAxis[jj];
					}
				}
				projectedPoly->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
				projectedPoly->SetPoints(projectedPoints);
				
				// compute overlap (if any) and adjust accordingly
				double intersectPt0[3], intersectPt1[3];
				double tol = 1E-3;
				if(intersectConvexCellsInPlane(avgBarrels[centerBarrel]->GetCell(1), projectedPoly->GetCell(0), tol, intersectPt0, intersectPt1) == 2)
				{
					// first determine which points need to be adjusted
					// they need to intersect the line formed by the intersection pts
					// more robust: determine their sign wrt intersection plane
					// plane normal perpendicular to center axis and intersection line
					double pNormal[3], intersectVec[3];
					for(int ii = 0; ii < 3; ++ii)
						intersectVec[ii] = intersectPt0[ii] - intersectPt1[ii];
					vtkMath::Normalize(intersectVec);
					pNormal[0] = centerAxis[1]*intersectVec[2] - centerAxis[2]*intersectVec[1];
					pNormal[1] = centerAxis[2]*intersectVec[0] - centerAxis[0]*intersectVec[2];
					pNormal[2] = centerAxis[0]*intersectVec[1] - centerAxis[1]*intersectVec[0];
					vtkMath::Normalize(pNormal);
					PlanePointerType intersectPlane = PlanePointerType::New();
					intersectPlane->SetNormal(pNormal);
					intersectPlane->SetOrigin(intersectPt0);
					
					// now, find the points and push them away from intersection line
					// towards their home barrel center
					PointsPointerType centerPoints = avgBarrels[centerBarrel]->GetCell(1)->GetPoints();
					PointsPointerType newCenterPoints = PointsPointerType::New();
					newCenterPoints->SetDataTypeToFloat();
					bool updateCenterContours = 0;
					int centerInside;
					if(intersectPlane->EvaluateFunction(centerBottom) < 0)
						centerInside = -1;
					else
						centerInside = 1;
					for(int ii = 0; ii < centerPoints->GetNumberOfPoints(); ++ii)
					{
						double pt[3];
						centerPoints->GetPoint(ii, pt);
						if(centerInside != (intersectPlane->EvaluateFunction(pt) < 0 ? -1 : 1))
						{
							updateCenterContours = 1;
							double newDist = 0;
							for(int jj = 0; jj < 3; ++jj)
								newDist += (pt[jj] - centerBottom[jj])*(pt[jj] - centerBottom[jj]);
							newDist = sqrt(newDist);
							newDist = std::max(newDist - intersectPlane->DistanceToPlane(pt) - 8.0, 1.0);
							for(int jj = 0; jj < 3; ++jj)
								pt[jj] = pt[jj] - centerBottom[jj];
							vtkMath::Normalize(pt);
							double newTopPt[3];
							for(int jj = 0; jj < 3; ++jj)
							{
								pt[jj] = newDist*pt[jj];
								pt[jj] += centerBottom[jj];
							}
							newCenterPoints->InsertNextPoint(pt);
						}
						else
							newCenterPoints->InsertNextPoint(pt);
					}
					// check projected neighbor points
					// don't re-insert them right away, they first need to be projected back
					PointsPointerType newNeighborPoints = PointsPointerType::New();
					newNeighborPoints->SetDataTypeToFloat();
					bool updateNeighborContours = 0;
					int neighborInside = -1*centerInside;	// on opposite side of intersect plane
					for(int ii = 0; ii < projectedPoints->GetNumberOfPoints(); ++ii)
					{
						double pt[3];
						projectedPoints->GetPoint(ii, pt);
						if(neighborInside != (intersectPlane->EvaluateFunction(pt) < 0 ? -1 : 1))
						{
							updateNeighborContours = 1;
							double newDist = 0;
							for(int jj = 0; jj < 3; ++jj)
								newDist += (pt[jj] - projectedNeighborBottom[jj])*(pt[jj] - projectedNeighborBottom[jj]);
							newDist = sqrt(newDist);
							newDist = std::max(newDist - intersectPlane->DistanceToPlane(pt) - 8.0, 1.0);
							for(int jj = 0; jj < 3; ++jj)
								pt[jj] = pt[jj] - projectedNeighborBottom[jj];
							vtkMath::Normalize(pt);
							for(int jj = 0; jj < 3; ++jj)
								pt[jj] = newDist*pt[jj] + projectedNeighborBottom[jj];
							newNeighborPoints->InsertNextPoint(pt);
						}
						else
							newNeighborPoints->InsertNextPoint(pt);
					}
					// only update if actually changes occured
					if(updateCenterContours)
					{
						// project point on axis to compute radial version of pt
						// then shift it to bottom contour again
						double longZ0[3], longZ1[3];
						for(int ii = 0; ii < 3; ++ii)
						{
							longZ0[ii] = centerBottom[ii] + 1000*centerAxis[ii];
							longZ1[ii] = centerBottom[ii] - 1000*centerAxis[ii];
						}
						PolyDataPointerType newContours = PolyDataPointerType::New();
						PointsPointerType newContourPts = PointsPointerType::New();
						PolygonPointerType bottomPoly = PolygonPointerType::New();
						PolygonPointerType topPoly = PolygonPointerType::New();
						newContours->Allocate(1);
						newContourPts->SetDataTypeToFloat();
						bottomPoly->GetPointIds()->SetNumberOfIds(newCenterPoints->GetNumberOfPoints());
						topPoly->GetPointIds()->SetNumberOfIds(newCenterPoints->GetNumberOfPoints());
						for(int ii = 0; ii < newCenterPoints->GetNumberOfPoints(); ++ii)
						{
							double * pt = new double[3], axisPt[3], * newTopPt = new double[3], t;
							newCenterPoints->GetPoint(ii, pt);
							vtkLine::DistanceToLine(pt, longZ0, longZ1, t, axisPt);
							for(int jj = 0; jj < 3; ++jj)
							{
								pt[jj] = pt[jj] - axisPt[jj];
								newTopPt[jj] = pt[jj] + centerTop[jj];
								pt[jj] += centerBottom[jj];
							}
							newContourPts->InsertNextPoint(newTopPt);
							newContourPts->InsertNextPoint(pt);
							topPoly->GetPointIds()->InsertId(ii, 2*ii);
							bottomPoly->GetPointIds()->InsertId(ii, 2*ii+1);
						}
						newContours->InsertNextCell(topPoly->GetCellType(), topPoly->GetPointIds());
						newContours->InsertNextCell(bottomPoly->GetCellType(), bottomPoly->GetPointIds());
						newContours->SetPoints(newContourPts);
						newContours->Update();
						avgBarrels[centerBarrel]->DeepCopy(newContours);
					}
					if(updateNeighborContours)
					{
						// project point on axis to compute radial version of pt
						// then shift it to bottom contour again
						double longZ0[3], longZ1[3];
						for(int ii = 0; ii < 3; ++ii)
						{
							longZ0[ii] = neighborBottom[ii] + 1000*neighborAxis[ii];
							longZ1[ii] = neighborBottom[ii] - 1000*neighborAxis[ii];
						}
						PolyDataPointerType newContours = PolyDataPointerType::New();
						PointsPointerType newContourPts = PointsPointerType::New();
						PolygonPointerType bottomPoly = PolygonPointerType::New();
						PolygonPointerType topPoly = PolygonPointerType::New();
						newContours->Allocate(1);
						newContourPts->SetDataTypeToFloat();
						bottomPoly->GetPointIds()->SetNumberOfIds(newNeighborPoints->GetNumberOfPoints());
						topPoly->GetPointIds()->SetNumberOfIds(newNeighborPoints->GetNumberOfPoints());
						for(int ii = 0; ii < newNeighborPoints->GetNumberOfPoints(); ++ii)
						{
							double * pt = new double[3], axisPt[3], * newTopPt = new double[3], t;
							newNeighborPoints->GetPoint(ii, pt);
							vtkLine::DistanceToLine(pt, longZ0, longZ1, t, axisPt);
							for(int jj = 0; jj < 3; ++jj)
							{
								pt[jj] = pt[jj] - axisPt[jj];
								newTopPt[jj] = pt[jj] + neighborTop[jj];
								pt[jj] += neighborBottom[jj];
							}
							newContourPts->InsertNextPoint(newTopPt);
							newContourPts->InsertNextPoint(pt);
							topPoly->GetPointIds()->InsertId(ii, 2*ii);
							bottomPoly->GetPointIds()->InsertId(ii, 2*ii+1);
						}
						newContours->InsertNextCell(topPoly->GetCellType(), topPoly->GetPointIds());
						newContours->InsertNextCell(bottomPoly->GetCellType(), bottomPoly->GetPointIds());
						newContours->SetPoints(newContourPts);
						newContours->Update();
						avgBarrels[neighborBarrel]->DeepCopy(newContours);
					}
				}
			}
		}
	
// 	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
// 			spatialGraph->addPolyDataObject(avgBarrels[*labelIt], *labelIt);
}

PolyDataPointerType CortexReconstruction::createTopBottomContours ( PolyDataPointerType avgContour, std::vector< double* > endPoints )
{
	std::vector< std::vector< double > > radialContour;
	PointsPointerType contourPts = avgContour->GetPoints();
	unsigned int nrOfPoints = contourPts->GetNumberOfPoints();
	for(int ii = 0; ii < nrOfPoints; ++ii)
	{
		double contourPt[3], closestPt[3], t;
		std::vector< double > radialPt;
		contourPts->GetPoint(ii, contourPt);
		vtkLine::DistanceToLine(contourPt, endPoints[0], endPoints[1], t, closestPt);
		for(int jj = 0; jj < 3; ++jj)
			radialPt.push_back(contourPt[jj] - closestPt[jj]);
		radialContour.push_back(radialPt);
	}
	
	PolyDataPointerType topBottomContours = PolyDataPointerType::New();
	PointsPointerType topBottomPts = PointsPointerType::New();
	topBottomContours->Allocate(1);
	topBottomPts->SetDataTypeToFloat();
	topBottomPts->SetNumberOfPoints(2*nrOfPoints);
	for(int contourID = 0; contourID < 2; ++contourID)
	{
		PolygonPointerType contourPoly = PolygonPointerType::New();
		contourPoly->GetPointIds()->SetNumberOfIds(nrOfPoints);
		for(int ii = 0; ii < nrOfPoints; ++ii)
		{
			double newPt[3];
			for(int jj = 0; jj < 3; ++jj)
				newPt[jj] = endPoints[contourID][jj] + radialContour[ii][jj];
			
			int IDoffset = contourID*nrOfPoints;
			topBottomPts->InsertPoint(ii + IDoffset, newPt);
			contourPoly->GetPointIds()->SetId(ii, ii + IDoffset);
		}
		topBottomContours->InsertNextCell(contourPoly->GetCellType(), contourPoly->GetPointIds());
	}
	topBottomContours->SetPoints(topBottomPts);
	topBottomContours->Update();
	
	#ifdef DEBUG
	if(spatialGraph)
		spatialGraph->addPolyDataObject(topBottomContours, Barrel);
	#endif
	
	return topBottomContours;
}

void CortexReconstruction::newBarrelAxis ( PolyDataPointerType barrel, PolyDataPointerType piaSurface, std::list< unsigned int > vessels, double alpha, double axis[3] )
{
	double barrelCOM[3];
	calculateBarrelCentroid(barrel, barrelCOM);
	std::multimap< double, double * > scores = barrelAxisScores(piaSurface, barrelCOM, 2000, alpha);
	
	#ifdef USEVESSELS
	unsigned int directionVessel = vesselsAroundBarrel(barrelCOM, vessels);
	
	std::multimap< double, double * >::reverse_iterator scoreIt;
	std::multimap< double, double * >::iterator vesselIt;
	for(scoreIt = scores.rbegin(); scoreIt != scores.rend(); ++scoreIt)
	{
		bool keepAxis = 1;
		double * axis = scoreIt->second;
		if(directionVessel)
		{
			//front = top, back = bottom
			double * vesselPt1 = (*(spatialGraph->edgesPointer()))[directionVessel]->edgePointCoordinates.front();
			double * vesselPt2 = (*(spatialGraph->edgesPointer()))[directionVessel]->edgePointCoordinates.back();
			double vessel[3];
			for(int ii = 0; ii < 3; ++ii)
				vessel[ii] = vesselPt1[ii] - vesselPt2[ii];
			double aNorm = 0, vNorm = 0, angle = 0;
			for(int ii = 0; ii < 3; ++ii)
			{
				aNorm += axis[ii]*axis[ii];
				vNorm += vessel[ii]*vessel[ii];
				angle += axis[ii]*vessel[ii];
			}
			aNorm = sqrt(aNorm);
			vNorm = sqrt(vNorm);
			angle = acos(std::abs(angle)/(aNorm*vNorm))*180.0/PI;
			if(angle > 5)
				keepAxis = 0;
		}
		
		if(keepAxis)
			break;
	}
	if(scoreIt == scores.rend())
	{
		std::cout << "Warning! Could not find barrel axis satisfying the vessel constraints. Selecting axis with highest score." << std::endl;
		scoreIt = scores.rbegin();
	}
	axis[0] = scoreIt->second[0], axis[1] = scoreIt->second[1], axis[2] = scoreIt->second[2];
	#endif
	#ifndef USEVESSELS
	axis[0] = scores.rbegin()->second[0], axis[1] = scores.rbegin()->second[1], axis[2] = scores.rbegin()->second[2];
	#endif
	/*** BEGIN MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
// 		std::cout << "Writing 5 highest scoring axes out of " << scores.size() << std::endl;
// 		std::multimap< double, double * >::const_reverse_iterator axesIt;
// 		int count = 0;
// 		int nrOfAxisPoints = 2;
// 		for(axesIt = scores.rbegin(); axesIt != scores.rend() && count < 5; ++axesIt, ++count)
// 		{
// 			double * endPoint = new double[3];
// 			double * bottomPoint = new double[3];
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				endPoint[ii] = barrelCOM[ii] + axesIt->second[ii];
// 				bottomPoint[ii] = barrelCOM[ii] - 3*axesIt->second[ii];
// 			}
// 			int label;
// 			if(count)
// 				label = WhiteMatter;
// 			else
// 				label = ZAxis;
// 			Vertex * newVert1 = new Vertex(endPoint, label);
// 			Vertex * newVert2 = new Vertex(barrelCOM, label);
// 			workingSG->addVertex(newVert1);
// 			workingSG->addVertex(newVert2);
// 			int connectionIndex[2];
// 			if(!workingSG->getNumberOfVertices())
// 			{
// 				connectionIndex[0] = 0;
// 				connectionIndex[1] = 1;
// 			}
// 			else
// 			{
// 				connectionIndex[0] = workingSG->getNumberOfVertices() - 2;
// 				connectionIndex[1] = workingSG->getNumberOfVertices() - 1;
// 			}
// // 			int noOfAxisPoints = 2;
// 			std::list< double * > axisCoords;
// 			axisCoords.push_back(endPoint);
// 			axisCoords.push_back(barrelCOM);
// 			Edge * newAxis = new Edge(connectionIndex, nrOfAxisPoints, label, axisCoords);
// 			workingSG->addEdge(newAxis);
// 		}
	/*** END MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
}

std::multimap< double, double* > CortexReconstruction::barrelAxisScores ( PolyDataPointerType piaSurface, double* barrelCentroid, double radius, double alpha )
{
	if(barrelCentroid != NULL)
	{
		PolyDataPointerType piaROI = selectSurfaceRoi(piaSurface, barrelCentroid, radius);
		if(piaROI->GetNumberOfCells())
		{
			std::vector< double > dotP;
			std::vector< double > dist;
			std::vector< double * > axis;
			std::multimap< double, double * > score;
			for(int ii = 0; ii < piaROI->GetNumberOfCells(); ++ii)
			{
				CellPointerType currentCell = piaROI->GetCell(ii);
				if(currentCell->GetNumberOfPoints() != 3)
				{
					std::cout << "Warning! vtkExtractPolyDataGeometry ouput is not triangles!" << std::endl;
				}
				PointsPointerType tmpPoints = currentCell->GetPoints();
				double centerPoint[] = {0, 0, 0};
				double axisNorm = 0;
				double * tmpAxis = new double[3];
				int nrOfPoints = tmpPoints->GetNumberOfPoints();
				for(int jj = 0; jj < nrOfPoints; ++jj)
				{
					double * tmpCoords = new double[3];
					tmpPoints->GetPoint(jj, tmpCoords);
					for(int kk = 0; kk < 3; ++kk)
						centerPoint[kk] += tmpCoords[kk];
					delete [] tmpCoords;
				}
				for(int jj = 0; jj < 3; ++jj)
					centerPoint[jj] = centerPoint[jj]/(double)nrOfPoints;
				
				for(int jj = 0; jj < 3; ++jj)
				{
					tmpAxis[jj] = centerPoint[jj] - barrelCentroid[jj];
					axisNorm += tmpAxis[jj]*tmpAxis[jj];
				}
				axisNorm = sqrt(axisNorm);
				
				PointsPointerType cellPoints = currentCell->GetPoints();
				double p1[3], p2[3], p3[3], normal[3];
				cellPoints->GetPoint(0, p1), cellPoints->GetPoint(1, p2), cellPoints->GetPoint(2, p3);
				for(int jj = 0; jj < 3; ++jj)
				{
					p3[jj] -= p1[jj];
					p2[jj] -= p1[jj];
				}
				normal[0] = p2[1]*p3[2] - p2[2]*p3[1];
				normal[1] = p2[2]*p3[0] - p2[0]*p3[2];
				normal[2] = p2[0]*p3[1] - p2[1]*p3[0];
				
				double tmpDotP = 0, normN = 0;
				for(int jj = 0; jj < 3; ++jj)
				{
					tmpDotP += tmpAxis[jj]*normal[jj];
					normN += normal[jj]*normal[jj];
				}
				normN = sqrt(normN);
				tmpDotP /= normN;
				tmpDotP /= axisNorm;
				tmpDotP = std::abs(tmpDotP);
				dotP.push_back(tmpDotP);
				dist.push_back(axisNorm);
				axis.push_back(tmpAxis);
			}
			
			if(dotP.size() == dist.size() && dotP.size() == axis.size())
			{
				if(score.size())
					score.clear();
				double maxDist = 0;
				for(int ii = 0; ii < dist.size(); ++ii)
					if(dist[ii] > maxDist)
						maxDist = dist[ii];
				for(int ii = 0; ii < dist.size(); ++ii)
				{
					double tmpScore = alpha*dotP[ii] + (1 - alpha)*(1 - dist[ii]/maxDist);
					score.insert(std::pair< double, double * >(tmpScore, axis[ii]));
				}
				return score;
			}
			else
			{
				std::cout << "Error! dotP.size() != dist.size() . Cannot calculate best axis." << std::endl;
				std::multimap< double, double * > emptyScore;
				return emptyScore;
			}
		}
	}
	std::multimap< double, double * > emptyScore;
	return emptyScore;
}

PolyDataPointerType CortexReconstruction::selectSurfaceRoi ( PolyDataPointerType surface, double* center, double radius )
{
	CylinderPointerType cylinder = CylinderPointerType::New();
	TransformPointerType cylinderTransform = TransformPointerType::New();
	TransformPointerType cylinderTranslate = TransformPointerType::New();
	ExtractPolyDataGeometryPointerType roiClipper = ExtractPolyDataGeometryPointerType::New();
	CellLocatorPointerType locator = CellLocatorPointerType::New();
	//Cylinder is centered at Center and axes of rotation is along the y-axis
	//The transformation transforms a point into the space of the implicit function (i.e., the model space).
	//=> rotation about x-axis by 90 degrees! has to be in homogeneous coordinates
	double offset[] = {0, 0, 0};
	offset[0] = -center[0];
	offset[1] = -center[1];
	offset[2] = -center[2];
	cylinderTranslate->Translate(offset[0], offset[1], offset[2]);
	cylinderTransform->RotateX(90);
	cylinderTransform->Concatenate(cylinderTranslate);
	
	cylinder->SetRadius(radius);
	cylinder->SetCenter(0, 0, 0);
	cylinder->SetTransform(cylinderTransform);
	
	roiClipper->SetImplicitFunction(cylinder);
	roiClipper->ExtractInsideOn();
	roiClipper->ExtractBoundaryCellsOn();
	roiClipper->SetInput(surface);
	roiClipper->Update();
	return roiClipper->GetOutput();
}

void CortexReconstruction::enforceAxisDivergence ( std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters )
{
	std::map< int, double * >::iterator barrelAxesIt;
	std::map< int, double * >::reverse_iterator barrelAxesRevIt;
// 	std::map< int, std::list< int > >::iterator barrelGridIt;
// 	std::map< int, std::list< int > >::reverse_iterator barrelGridRevIt;
	for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		vtkMath::Normalize(barrelAxesIt->second);
	std::map< int, std::list< int > > barrelGrid = createBarrelGrid(barrelAxes);
	
	int maxLoops = 0;
	bool change = 0;
	do
	{
// 		std::cout << "loop nr " << maxLoops << std::endl;
		change = 0;
		for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		{
			int thisBarrel = barrelAxesIt->first;
			double thisAxis[3], cumulatedChange[3];
			thisAxis[0] = barrelAxesIt->second[0], thisAxis[1] = barrelAxesIt->second[1], thisAxis[2] = barrelAxesIt->second[2];
			cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
			std::list< int >::iterator neighborIt;
			for(neighborIt = barrelGrid[barrelAxesIt->first].begin(); neighborIt != barrelGrid[barrelAxesIt->first].end(); ++neighborIt)
			{
				double neighborAxis[3];
				neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
				double direction[3];
				for(int ii = 0; ii < 3; ++ii)
					direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
				vtkMath::Normalize(direction);
				//eliminate part of neighboraxis that is parallel to this axis
				double prod = 0;
				for(int ii = 0; ii < 3; ++ii)
					prod += thisAxis[ii]*neighborAxis[ii];
				for(int ii = 0; ii < 3; ++ii)
					neighborAxis[ii] -= prod*thisAxis[ii];
				//check if other axis is pointing out of our local coordinate system
				double localDirection = 0;
				for(int ii = 0; ii < 3; ++ii)
					localDirection += direction[ii]*neighborAxis[ii];
				// < 0 means other axis is pointing towards us
				if(localDirection < 0)
				{
					change = 1;
					for(int ii = 0; ii < 3; ++ii)
					{
						cumulatedChange[ii] += 0.5*localDirection*direction[ii];
						barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
					}
					vtkMath::Normalize(barrelAxes[*neighborIt]);
				}
			}
			for(int ii = 0; ii < 3; ++ii)
				barrelAxesIt->second[ii] += cumulatedChange[ii];
			vtkMath::Normalize(barrelAxesIt->second);
		}
		if(change)
		{
			for(barrelAxesRevIt = barrelAxes.rbegin(); barrelAxesRevIt != barrelAxes.rend(); ++barrelAxesRevIt)
			{
				int thisBarrel = barrelAxesRevIt->first;
				double thisAxis[3], cumulatedChange[3];
				thisAxis[0] = barrelAxesRevIt->second[0], thisAxis[1] = barrelAxesRevIt->second[1], thisAxis[2] = barrelAxesRevIt->second[2];
				cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
				std::list< int >::iterator neighborIt;
				for(neighborIt = barrelGrid[barrelAxesRevIt->first].begin(); neighborIt != barrelGrid[barrelAxesRevIt->first].end(); ++neighborIt)
				{
					double neighborAxis[3];
					neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
					double direction[3];
					for(int ii = 0; ii < 3; ++ii)
						direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
					vtkMath::Normalize(direction);
					//eliminate part of neighboraxis that is parallel to this axis
					double prod = 0;
					for(int ii = 0; ii < 3; ++ii)
						prod += thisAxis[ii]*neighborAxis[ii];
					for(int ii = 0; ii < 3; ++ii)
						neighborAxis[ii] -= prod*thisAxis[ii];
					//check if other axis is pointing out of our local coordinate system
					double localDirection = 0;
					for(int ii = 0; ii < 3; ++ii)
						localDirection += direction[ii]*neighborAxis[ii];
					// < 0 means other axis is pointing towards us
					if(localDirection < 0)
					{
						change = 1;
						for(int ii = 0; ii < 3; ++ii)
						{
							cumulatedChange[ii] += 0.5*localDirection*direction[ii];
							barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
						}
						vtkMath::Normalize(barrelAxes[*neighborIt]);
					}
				}
				for(int ii = 0; ii < 3; ++ii)
					barrelAxesRevIt->second[ii] += cumulatedChange[ii];
				vtkMath::Normalize(barrelAxesRevIt->second);
			}
		}
		++maxLoops;
	} while(change && maxLoops < 10);
}

std::map< int, std::list< int > > CortexReconstruction::createBarrelGrid ( std::map< int, double* > barrelAxes )
{
	std::map< int, std::list< int > > grid;
	std::list< int > barrelIDs;
	std::map< int, double * >::const_iterator barrelAxesIt;
	for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		barrelIDs.push_back(barrelAxesIt->first);
	barrelIDs.sort();
	
	std::list< int > alpha;
	alpha.push_back(A1), alpha.push_back(B1), alpha.push_back(Beta);
	std::list< int > beta;
	beta.push_back(Alpha), beta.push_back(B1), beta.push_back(C1), beta.push_back(Gamma);
	std::list< int > gamma;
	gamma.push_back(Beta), gamma.push_back(C1), gamma.push_back(D1), gamma.push_back(Delta);
	std::list< int > delta;
	delta.push_back(Gamma), delta.push_back(D1), delta.push_back(E1);
	std::list< int > a1;
	a1.push_back(Alpha), a1.push_back(B1), a1.push_back(B2), a1.push_back(A2);
	std::list< int > a2;
	a2.push_back(A1), a2.push_back(B1), a2.push_back(B2), a2.push_back(B3), a2.push_back(A3);
	std::list< int > a3;
	a3.push_back(A2), a3.push_back(B2), a3.push_back(B3), a3.push_back(B4), a3.push_back(A4);
	std::list< int > a4;
	a4.push_back(A3), a4.push_back(B3), a4.push_back(B4);
	std::list< int > b1;
	b1.push_back(Beta), b1.push_back(C1), b1.push_back(C2), b1.push_back(B2), b1.push_back(A2), b1.push_back(A1), b1.push_back(Alpha);
	std::list< int > b2;
	b2.push_back(B1), b2.push_back(C1), b2.push_back(C2), b2.push_back(C3), b2.push_back(B3), b2.push_back(A3), b2.push_back(A2), b2.push_back(A1);
	std::list< int > b3;
	b3.push_back(B2), b3.push_back(C2), b3.push_back(C3), b3.push_back(C4), b3.push_back(B4), b3.push_back(A4), b3.push_back(A3), b3.push_back(A2);
	std::list< int > b4;
	b4.push_back(B3), b4.push_back(C3), b4.push_back(C4), b4.push_back(C5), b4.push_back(A4), b4.push_back(A3);
	std::list< int > c1;
	c1.push_back(Gamma), c1.push_back(D1), c1.push_back(D2), c1.push_back(C2), c1.push_back(B2), c1.push_back(B1), c1.push_back(Beta);
	std::list< int > c2;
	c2.push_back(C1), c2.push_back(D1), c2.push_back(D2), c2.push_back(D3), c2.push_back(C3), c2.push_back(B3), c2.push_back(B2), c2.push_back(B1);
	std::list< int > c3;
	c3.push_back(C2), c3.push_back(D2), c3.push_back(D3), c3.push_back(D4), c3.push_back(C4), c3.push_back(B4), c3.push_back(B3), c3.push_back(B2);
	std::list< int > c4;
	c4.push_back(C3), c4.push_back(D3), c4.push_back(D4), c4.push_back(D5), c4.push_back(C5), c4.push_back(B4), c4.push_back(B3);
	std::list< int > c5;
	c5.push_back(C4), c5.push_back(D4), c5.push_back(D5), c5.push_back(D6), c5.push_back(C6), c5.push_back(B4);
	std::list< int > c6;
	c6.push_back(C5), c6.push_back(D5), c6.push_back(D6);
	std::list< int > d1;
	d1.push_back(Delta), d1.push_back(E1), d1.push_back(E2), d1.push_back(D2), d1.push_back(C2), d1.push_back(C1), d1.push_back(Gamma);
	std::list< int > d2;
	d2.push_back(D1), d2.push_back(E1), d2.push_back(E2), d2.push_back(E3), d2.push_back(D3), d2.push_back(C3), d2.push_back(C2), d2.push_back(C1);
	std::list< int > d3;
	d3.push_back(D2), d3.push_back(E2), d3.push_back(E3), d3.push_back(E4), d3.push_back(D4), d3.push_back(C4), d3.push_back(C3), d3.push_back(C2);
	std::list< int > d4;
	d4.push_back(D3), d4.push_back(E3), d4.push_back(E4), d4.push_back(E5), d4.push_back(D5), d4.push_back(C5), d4.push_back(C4), d4.push_back(C3);
	std::list< int > d5;
	d5.push_back(D4), d5.push_back(E4), d5.push_back(E5), d5.push_back(E6), d5.push_back(D6), d5.push_back(C6), d5.push_back(C5), d5.push_back(C4);
	std::list< int > d6;
	d6.push_back(D5), d6.push_back(E5), d6.push_back(E6), d6.push_back(C6), d6.push_back(C5);
	std::list< int > e1;
	e1.push_back(Delta), e1.push_back(E2), e1.push_back(D2), e1.push_back(D1);
	std::list< int > e2;
	e2.push_back(E1), e2.push_back(E3), e2.push_back(D3), e2.push_back(D2), e2.push_back(D1);
	std::list< int > e3;
	e3.push_back(E2), e3.push_back(E4), e3.push_back(D4), e3.push_back(D3), e3.push_back(D2);
	std::list< int > e4;
	e4.push_back(E3), e4.push_back(E5), e4.push_back(D5), e4.push_back(D4), e4.push_back(D3);
	std::list< int > e5;
	e5.push_back(E4), e5.push_back(E6), e5.push_back(D6), e5.push_back(D5), e5.push_back(D4);
	std::list< int > e6;
	e6.push_back(E5), e6.push_back(D6), e6.push_back(D5);
	
	grid.insert(std::pair< int, std::list< int > >(Alpha, alpha));
	grid.insert(std::pair< int, std::list< int > >(Beta, beta));
	grid.insert(std::pair< int, std::list< int > >(Gamma, gamma));
	grid.insert(std::pair< int, std::list< int > >(Delta, delta));
	grid.insert(std::pair< int, std::list< int > >(A1, a1));
	grid.insert(std::pair< int, std::list< int > >(A2, a2));
	grid.insert(std::pair< int, std::list< int > >(A3, a3));
	grid.insert(std::pair< int, std::list< int > >(A4, a4));
	grid.insert(std::pair< int, std::list< int > >(B1, b1));
	grid.insert(std::pair< int, std::list< int > >(B2, b2));
	grid.insert(std::pair< int, std::list< int > >(B3, b3));
	grid.insert(std::pair< int, std::list< int > >(B4, b4));
	grid.insert(std::pair< int, std::list< int > >(C1, c1));
	grid.insert(std::pair< int, std::list< int > >(C2, c2));
	grid.insert(std::pair< int, std::list< int > >(C3, c3));
	grid.insert(std::pair< int, std::list< int > >(C4, c4));
	grid.insert(std::pair< int, std::list< int > >(C5, c5));
	grid.insert(std::pair< int, std::list< int > >(C6, c6));
	grid.insert(std::pair< int, std::list< int > >(D1, d1));
	grid.insert(std::pair< int, std::list< int > >(D2, d2));
	grid.insert(std::pair< int, std::list< int > >(D3, d3));
	grid.insert(std::pair< int, std::list< int > >(D4, d4));
	grid.insert(std::pair< int, std::list< int > >(D5, d5));
	grid.insert(std::pair< int, std::list< int > >(D6, d6));
	grid.insert(std::pair< int, std::list< int > >(E1, e1));
	grid.insert(std::pair< int, std::list< int > >(E2, e2));
	grid.insert(std::pair< int, std::list< int > >(E3, e3));
	grid.insert(std::pair< int, std::list< int > >(E4, e4));
	grid.insert(std::pair< int, std::list< int > >(E5, e5));
	grid.insert(std::pair< int, std::list< int > >(E6, e6));
	
	std::map< int, std::list< int > >::iterator gridIt;
	std::list< int >::iterator neighborIt;
	for(gridIt = grid.begin(); gridIt != grid.end(); )
	{
		if(std::find(barrelIDs.begin(), barrelIDs.end(), gridIt->first) != barrelIDs.end())
		{
			for(neighborIt = gridIt->second.begin(); neighborIt != gridIt->second.end(); )
			{
				if(std::find(barrelIDs.begin(), barrelIDs.end(), *neighborIt) != barrelIDs.end())
					++neighborIt;
				else
					neighborIt = gridIt->second.erase(neighborIt);
			}
			++gridIt;
		}
		else
		{
			std::map< int, std::list< int > >::iterator eraseIt = gridIt;
			++gridIt;
			grid.erase(eraseIt);
		}
	}
	
// 	for(gridIt = grid.begin(); gridIt != grid.end(); ++gridIt)
// 	{
// 		std::cout << "Barrel " << int2Labels[gridIt->first] << " has the following neighbors:" << std::endl;
// 		for(neighborIt = gridIt->second.begin(); neighborIt != gridIt->second.end(); ++neighborIt)
// 			std::cout << int2Labels[*neighborIt] << "\t";
// 		std::cout << std::endl;
// 	}
	
	return grid;
}

void CortexReconstruction::calculateBarrelCentroid ( PolyDataPointerType barrel, double centroid[3] )
{
	for(int ii = 0; ii < 3; ++ii)
		centroid[ii] = 0;
	int nrOfCells = barrel->GetNumberOfCells();
	for(int ii = 0; ii < nrOfCells; ++ii)
	{
		double centerPoint[3], paramCenter[3];
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[barrel->GetCell(ii)->GetNumberOfPoints()];
		barrel->GetCell(ii)->GetParametricCenter(pCoords);
		barrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		for(int jj = 0; jj < 3; ++jj)
			centroid[jj] += centerPoint[jj];
		delete [] weights;
	}
	for(int jj = 0; jj < 3; ++jj)
		centroid[jj] = centroid[jj]/(double)nrOfCells;
}

void CortexReconstruction::closeBarrelAlongNewAxis ( double newAxis[3], double barrelCentroid[3], PolyDataPointerType barrel, std::vector< double* >& endPoints)
{
	if(barrel->GetNumberOfCells())
	{
		PolyDataPointerType barrelTop = PolyDataPointerType::New();
		PolyDataPointerType barrelBottom = PolyDataPointerType::New();
		IdListPointerType barrelTopCellIds = IdListPointerType::New();
		IdListPointerType barrelBottomCellIds = IdListPointerType::New();
		barrelTop->Allocate();
		barrelBottom->Allocate();
		
		int minID, maxID;
		getLandmarkMinMaxIDs(barrel, minID, maxID);
		barrelTopCellIds->InsertId(0, maxID);
		barrelBottomCellIds->InsertId(0, minID);
		barrelTop->CopyCells(barrel, barrelTopCellIds);
		barrelBottom->CopyCells(barrel, barrelBottomCellIds);
		
		// case first reconstruction:
		// find extreme bottom/top point by transforming the coordinate system
		// for every point in lowest/highest polygon and looking for extreme z values
		double normA = 0;
		double alpha = 0;
		double newZAxis[3];
		newZAxis[0] = newAxis[0], newZAxis[1] = newAxis[1], newZAxis[2] = newAxis[2];
		vtkMath::Normalize(newZAxis);
		alpha = acos(std::abs(newZAxis[2]));
		TransformPointerType translate = TransformPointerType::New();
		TransformPointerType rotate = TransformPointerType::New();
		TransformPointerType inverseTranslate = TransformPointerType::New();
		
		//old axis == z axis
		//rotation axis = newAxus x oldAxis
		//then, angle of rotation = acos(oldAxis*newAxis)
		double sign = parameters.zReversed ? -1 : 1;
		double rotationAxis[3];
		rotationAxis[0] = -sign*newZAxis[1];
		rotationAxis[1] = sign*newZAxis[0];
		rotationAxis[2] = 0;
		
		double angle = sign*alpha*180/PI;	// -1 b/c we want to rotate coordinate system, not points
		
		translate->Translate(-barrelCentroid[0], -barrelCentroid[1], -barrelCentroid[2]);
		rotate->RotateWXYZ(-angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
		inverseTranslate->Translate(barrelCentroid[0], barrelCentroid[1], barrelCentroid[2]);
		
		rotate->Concatenate(translate);
		inverseTranslate->Concatenate(rotate);
		inverseTranslate->Update();
		
		// case reconstruction of registered home barrel along standard axis:
		// simply project onto new axis and look at extreme z-values
		// (all handled by bool flag HBRecon)
		
		int extremeBottomID = 0, extremeTopID = 0;
		double bottomMin = sign*1E06, topMax = -sign*1E06;
		PointsPointerType topPts = barrelTop->GetCell(0)->GetPoints();
		PointsPointerType bottomPts = barrelBottom->GetCell(0)->GetPoints();
		for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
		{
			double thisZ;
			double * tmpPt = new double[3];
			topPts->GetPoint(ii, tmpPt);
			double linePt1[3], linePt2[3], projPt[3], t;
			for(int jj = 0; jj < 3; ++jj)
			{
				linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
				linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
			}
			vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
			thisZ = projPt[2];
// 			if(HBRecon)
// 			{
// 				double linePt1[3], linePt2[3], projPt[3], t;
// 				for(int jj = 0; jj < 3; ++jj)
// 				{
// 					linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
// 					linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
// 				}
// 				vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
// 				thisZ = projPt[2];
// 			}
// 			else
// 			{
// 				double homPt[4], transPt[4];
// 				homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
// 				inverseTranslate->MultiplyPoint(homPt, transPt);
// 				thisZ = transPt[2];
// 			}
			if(parameters.zReversed ? thisZ < topMax : thisZ > topMax)
			{
				topMax = thisZ;
				extremeTopID = ii;
			}
			delete [] tmpPt;
		}
		for(int ii = 0; ii < bottomPts->GetNumberOfPoints(); ++ii)
		{
			double thisZ;
			double * tmpPt = new double[3];
			bottomPts->GetPoint(ii, tmpPt);
			double linePt1[3], linePt2[3], projPt[3], t;
			for(int jj = 0; jj < 3; ++jj)
			{
				linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
				linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
			}
			vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
			thisZ = projPt[2];
// 			if(HBRecon)
// 			{
// 				double linePt1[3], linePt2[3], projPt[3], t;
// 				for(int jj = 0; jj < 3; ++jj)
// 				{
// 					linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
// 					linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
// 				}
// 				vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
// 				thisZ = projPt[2];
// 			}
// 			else
// 			{
// 				double homPt[4], transPt[4];
// 				homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
// 				inverseTranslate->MultiplyPoint(homPt, transPt);
// 				thisZ = transPt[2];
// 			}
			if(parameters.zReversed ? thisZ > bottomMin : thisZ < bottomMin)
			{
				bottomMin = thisZ;
				extremeBottomID = ii;
			}
			delete [] tmpPt;
		}
		
		LinePointerType newZAxisLine = LinePointerType::New();
		double * extremeTopPt = new double[3];
		double * extremeBottomPt = new double[3];
		double * finalTopPt = new double[3];
		double * finalBottomPt = new double[3];
		topPts->GetPoint(extremeTopID, extremeTopPt);
		bottomPts->GetPoint(extremeBottomID, extremeBottomPt);
		double linePt1[3], linePt2[3];
		linePt1[0] = barrelCentroid[0] + 1000*newZAxis[0], linePt1[1] = barrelCentroid[1] + 1000*newZAxis[1], linePt1[2] = barrelCentroid[2] + 1000*newZAxis[2];
		linePt2[0] = barrelCentroid[0] - 1000*newZAxis[0], linePt2[1] = barrelCentroid[1] - 1000*newZAxis[1], linePt2[2] = barrelCentroid[2] - 1000*newZAxis[2];
		
		//check whether new axis actually passes through top contour!!!
		double * topBounds = barrelTop->GetCell(0)->GetBounds();
		double * bottomBounds = barrelBottom->GetCell(0)->GetBounds();
		double * axisTopCoord = new double[3];
		{axisTopCoord[0] = barrelCentroid[0] + std::abs(barrelCentroid[2]-topBounds[4])*newZAxis[0]/cos(alpha), axisTopCoord[1] = barrelCentroid[1] + std::abs(barrelCentroid[2]-topBounds[4])*newZAxis[1]/cos(alpha), axisTopCoord[2] = topBounds[4];}
		double * axisBottomCoord = new double[3];
		{axisBottomCoord[0] = barrelCentroid[0] - std::abs(bottomBounds[5] - barrelCentroid[2])*newZAxis[0]/cos(alpha), axisBottomCoord[1] = barrelCentroid[1] - std::abs(bottomBounds[5] - barrelCentroid[2])*newZAxis[1]/cos(alpha), axisBottomCoord[2] = bottomBounds[5];}
		int insideTop = 0, insideBottom = 0;
		double closestPoint1[3], closestPoint2[3];
		int subId1, subId2;
		double pCoords1[3], pCoords2[3], * weights1, * weights2, dist1, dist2, tmp1[3], tmp2[3];
		weights1 = new double[barrelTop->GetCell(0)->GetNumberOfPoints()];
		weights2 = new double[barrelBottom->GetCell(0)->GetNumberOfPoints()];
		if(axisTopCoord[0] >= topBounds[0] && axisTopCoord[0] <= topBounds[1] 
			&& axisTopCoord[1] >= topBounds[2] && axisTopCoord[1] <= topBounds[3] 
			&& axisTopCoord[2] >= topBounds[4] && axisTopCoord[2] <= topBounds[5])
		{
			insideTop = barrelTop->GetCell(0)->EvaluatePosition(axisTopCoord, closestPoint1, subId1, pCoords1, dist1, weights1);
		}
		if(axisBottomCoord[0] >= bottomBounds[0] && axisBottomCoord[0] <= bottomBounds[1] 
			&& axisBottomCoord[1] >= bottomBounds[2] && axisBottomCoord[1] <= bottomBounds[3] 
			&& axisBottomCoord[2] >= bottomBounds[4] && axisBottomCoord[2] <= bottomBounds[5])
		{
			insideBottom = barrelBottom->GetCell(0)->EvaluatePosition(axisBottomCoord, closestPoint2, subId2, pCoords2, dist2, weights2);
		}
		if(insideTop)
		{
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			endPoints.push_back(finalTopPt);
		}
		else if(!insideTop)
		{
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			double * finalTopPt2 = new double[3];
			finalTopPt2[0] = finalTopPt[0], finalTopPt2[1] = finalTopPt[1], finalTopPt2[2] = finalTopPt[2];
			endPoints.push_back(finalTopPt2);
		}
		if(insideBottom)
		{
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			endPoints.push_back(finalBottomPt);
		}
		else if(!insideBottom)
		{
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			double * finalBottomPt2 = new double[3];
			finalBottomPt2[0] = finalBottomPt[0], finalBottomPt2[1] = finalBottomPt[1], finalBottomPt2[2] = finalBottomPt[2];
			endPoints.push_back(finalBottomPt2);
		}
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel caps." << std::endl;
		return;
	}
}

void CortexReconstruction::getLandmarkMinMaxIDs ( PolyDataPointerType landmark, int& minID, int& maxID )
{
	double minZ, maxZ = landmark->GetCell(0)->GetBounds()[4];
	minZ = maxZ;
	minID = 0, maxID = 0;
	for(int ii = 1; ii < landmark->GetNumberOfCells(); ++ii)
	{
		double tmpZ = landmark->GetCell(ii)->GetBounds()[4];
		if(tmpZ > maxZ)
		{
			maxZ = tmpZ;
			maxID = ii;
		}
		if(tmpZ < minZ)
		{
			minZ = tmpZ;
			minID = ii;
		}
	}
	if(parameters.zReversed)
	{
		int tmp = minID;
		minID = maxID;
		maxID = tmp;
	}
}

int CortexReconstruction::intersectConvexCellsInPlane ( CellPointerType cell1, CellPointerType cell2, double tol, double p0[3], double p1[3] )
{
	double *x[2], pcoords[3], t, x0[3], x1[3], x2[3], x3[3];
	x[0] = p0; x[1] = p1;
	int subId, idx=0;
	double t2 = tol*tol;
	
	PointsPointerType cell1Pts = cell1->GetPoints();
	PointsPointerType cell2Pts = cell2->GetPoints();
	for(int ii = 0; ii < cell1Pts->GetNumberOfPoints(); ++ii)
	{
		cell1Pts->GetPoint(ii, x0);
		cell1Pts->GetPoint((ii+1)%cell1Pts->GetNumberOfPoints(), x1);
		PolyDataPointerType currLine = PolyDataPointerType::New();
		PointsPointerType currPts = PointsPointerType::New();
		LinePointerType line = LinePointerType::New();
		currLine->Allocate(1);
		currPts->SetDataTypeToFloat();
		line->GetPointIds()->SetNumberOfIds(2);
		line->GetPointIds()->SetId(0, 0), line->GetPointIds()->SetId(1, 1);
		currPts->InsertNextPoint(x0), currPts->InsertNextPoint(x1);
		currLine->InsertNextCell(line->GetCellType(), line->GetPointIds());
		currLine->SetPoints(currPts);
		// Loop over edges of second polygon and intersect against current edge
		for(int jj = 0; jj < cell2Pts->GetNumberOfPoints(); ++jj)
		{
			cell2Pts->GetPoint(jj, x2);
			cell2Pts->GetPoint((jj+1)%cell2Pts->GetNumberOfPoints(), x3);
			if(currLine->GetCell(0)->IntersectWithLine(x2, x3, tol, t, x[idx], pcoords, subId))
			{
				if(idx == 0)
				{
					++idx;
				}
				else if(((x[1][0]-x[0][0])*(x[1][0]-x[0][0]) + (x[1][1]-x[0][1])*(x[1][1]-x[0][1]) + 
					(x[1][2]-x[0][2])*(x[1][2]-x[0][2])) > t2)
				{
					return 2;
				}
			}//if edge intersection
		}//over all edges
	}
	
	// Evaluate what we got
	if (idx == 1)
	{
		return 1; //everything intersecting at single point
	}
	else
	{
		return 0;
	}
}

/******************************************************************************/
/*computes 3d blood vessels from vessel obects in spatialGraph by starting    */
/*from all vessels in the top plane and looking for vessels in a certain      */
/*radius in the next plane. If found, look for more vessels in deeper planes  */
/*in the direction specified by the first two vessels within tolerance        */
/******************************************************************************/
void CortexReconstruction::computeBloodVessels ( bool constrained, double minRadius )
{
	double piaCenter[3];
	if(constrained)
	{
		PolyDataPointerType pia = PolyDataPointerType::New();
		spatialGraph->extractLandmark(Pia, pia);
		double piaPCenter[3], * piaWeights =  new double[pia->GetCell(0)->GetNumberOfPoints()];
		int piaSubID;
		pia->GetCell(0)->GetParametricCenter(piaPCenter);
		pia->GetCell(0)->EvaluateLocation(piaSubID, piaPCenter, piaCenter, piaWeights);
		delete [] piaWeights;
	}
	// does not change anything significantly within the barrel field,
	// so leave at default values for now
	double lineDensity = 0.005, badAngle = 60, searchRadius = 30;	//default 0.005, 60, 30
// 	std::cout << "line density = ";
// 	std::cin >> lineDensity;
// 	std::cout << "bad angle = ";
// 	std::cin >> badAngle;
// 	std::cout << "search radius = ";
// 	std::cin >> searchRadius;
	
	spatialGraph->vesselsToPoints();
	#ifdef DEBUG
	std::flush(std::cout << "Sorting vessels by their z-position..." << std::endl);
	#endif
	std::multimap< double, Edge* > vessels;	// sort vessels by their z-position
	std::list< double > zPlanes;
	std::list< std::vector< Edge * > > vessels3D;
	std::vector< Edge * >::iterator graphIt;
	for(graphIt = spatialGraph->edgesBegin(); graphIt != spatialGraph->edgesEnd(); ++graphIt)
	{
		if((*graphIt)->label == Vessel && (*graphIt)->radius >= minRadius)
		{
			vessels.insert(std::pair< double, Edge* >(round((*graphIt)->edgePointCoordinates.back()[2]), *graphIt));
			zPlanes.push_back(round((*graphIt)->edgePointCoordinates.back()[2]));
		}
	}
	zPlanes.sort();
	std::list< double >::iterator zIter = zPlanes.begin();
	double last = *zIter;
	double curr = 0;
	++zIter;
	while(zIter != zPlanes.end())
	{
		curr = *zIter;
		if(std::abs(curr - last) < 0.5)
			zIter = zPlanes.erase(zIter);
		else
		{
			last = curr;
			++zIter;
		}
	}
	
	for(zIter = zPlanes.begin(); zIter != zPlanes.end(); ++zIter)
	{
		std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > planeVessels;
		std::multimap< double, Edge* >::iterator planeVesselIter;
		planeVessels = vessels.equal_range(*zIter);
		
		for(planeVesselIter = planeVessels.first; planeVesselIter != planeVessels.second; /*++planeVesselIter*/)
		{
			bool notInVessel = 1;
			std::vector< Edge * > tmp3DVessel;
			std::list< double >::iterator zIter2 = zIter;
			++zIter2;
			if(zIter2 != zPlanes.end())
			{
				tmp3DVessel.push_back(planeVesselIter->second);
				std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > nextPlaneVessels;
				nextPlaneVessels = vessels.equal_range(*zIter2);
				double planeDist = std::abs(*zIter - *zIter2);
				if(planeDist > 1.0)
				{
					double * vesselCenter = new double[3];
					vesselCenter = planeVesselIter->second->edgePointCoordinates.back();
					std::multimap< double, Edge* > nearestVessels = vesselDistances3D(vesselCenter, nextPlaneVessels);	// Vessels in next plane sorted by distance to current vessel
					if(!nearestVessels.size())
					{
						++planeVesselIter;
						continue;
					}
					std::multimap< double, Edge* >::iterator nearestVesselsIt;
					if(!constrained)
					{
						nearestVesselsIt = nearestVessels.begin();
						tmp3DVessel.push_back(nearestVesselsIt->second);
					}
					else
					{
						//new version: check direction before accepting
						bool goodDirection = 0;
						for(nearestVesselsIt = nearestVessels.begin(); nearestVesselsIt != nearestVessels.end(); ++nearestVesselsIt)
						{
							Edge * tmpEdge = nearestVesselsIt->second;
							double tmpVesselDirection[3], tmpCenterDirection[3];
							for(int ii = 0; ii < 2; ++ii)
							{
								tmpVesselDirection[ii] = tmpEdge->edgePointCoordinates.back()[ii] - vesselCenter[ii];
								tmpCenterDirection[ii] = vesselCenter[ii] - piaCenter[ii];
							}
							tmpVesselDirection[2] = 0, tmpCenterDirection[2] = 0;
							vtkMath::Normalize(tmpVesselDirection), vtkMath::Normalize(tmpCenterDirection);
							double angle = 0;
							for(int ii = 0; ii < 2; ++ii)
								angle += tmpVesselDirection[ii]*tmpCenterDirection[ii];
							angle = acos(angle)*180/PI;
							if(angle > badAngle)
							{
								goodDirection = 1;
								break;
							}
						}
						if(!goodDirection)
						{
							++planeVesselIter;
							continue;
						}
					}
					
					double vesselDirection[3];
					for(int ii = 0; ii < 3; ++ii)
						vesselDirection[ii] = nearestVesselsIt->second->edgePointCoordinates.back()[ii] - vesselCenter[ii];
					
					double norm = sqrt(vesselDirection[0]*vesselDirection[0] + vesselDirection[1]*vesselDirection[1] + vesselDirection[2]*vesselDirection[2]);
					for(int ii = 0; ii < 3; ++ii)
						vesselDirection[ii] /= norm;
					
					double lastPlane = *zIter2;
					double oldVesselcenter[3];
					for(int ii = 0; ii < 3; ++ii)
						oldVesselcenter[ii] = nearestVesselsIt->second->edgePointCoordinates.back()[ii];
					std::list< double >::iterator nextIt = zIter2;
					++nextIt;
					int gapCount = 1;
					while(nextIt != zPlanes.end() && gapCount < 5)
					{
						//implement search in other planes in direction of 'vesselDirection'
						nextPlaneVessels = vessels.equal_range(*nextIt);
						planeDist = std::abs(*nextIt - lastPlane);
						double cosAngle = std::abs(vesselDirection[2]);
						double nextVesselCenter[3];
						for(int ii = 0; ii < 3; ++ii)
							nextVesselCenter[ii] = oldVesselcenter[ii] + vesselDirection[ii]*planeDist/cosAngle;
						nearestVessels = vesselDistances2D(nextVesselCenter, nextPlaneVessels);
						if(!nearestVessels.size())
						{
							++nextIt;
							continue;
						}
						nearestVesselsIt = nearestVessels.begin();
						
// 						if(nearestVesselsIt->first < searchRadius*sqrt(gapCount))
						if(nearestVesselsIt->first < searchRadius*gapCount)
						{
							gapCount = 1;
							lastPlane = *nextIt;
							tmp3DVessel.push_back(nearestVesselsIt->second);
							
							for(int ii = 0; ii < 3; ++ii)
							{
								vesselDirection[ii] += (nearestVesselsIt->second->edgePointCoordinates.back()[ii] - oldVesselcenter[ii])/planeDist;
								oldVesselcenter[ii] = nearestVesselsIt->second->edgePointCoordinates.back()[ii];
							}
							
							double norm = sqrt(vesselDirection[0]*vesselDirection[0] + vesselDirection[1]*vesselDirection[1] + vesselDirection[2]*vesselDirection[2]);
							for(int ii = 0; ii < 3; ++ii)
								vesselDirection[ii] /= norm;
						}
						else
							++gapCount;
						
						++nextIt;
					}
					
					double tmpVesselVec[3];
					for(int ii = 0; ii < 3; ++ii)
						tmpVesselVec[ii] = tmp3DVessel.front()->edgePointCoordinates.back()[ii] - tmp3DVessel.back()->edgePointCoordinates.back()[ii];
					double tmpVesselLength = sqrt(tmpVesselVec[0]*tmpVesselVec[0] + tmpVesselVec[1]*tmpVesselVec[1] + tmpVesselVec[2]*tmpVesselVec[2]);
// 					normalize(tmpVesselVec);
					double lineThreshold = lineDensity;
// 					if(tmpVesselVec[2])
// 					{
// 						tmpVesselVec[2] = std::max(tmpVesselVec[2], 1/sqrt(2));	// vessels usually have max angle of ~45 deg
// 						lineThreshold *= tmpVesselVec[2];
// 					}
					if(tmp3DVessel.size() < 5 || (double)tmp3DVessel.size()/tmpVesselLength < lineThreshold)
						tmp3DVessel.clear();
					else
					{
						notInVessel = 0;
						std::vector< Edge * >::iterator tmp3DVesselIt = tmp3DVessel.begin();
						++tmp3DVesselIt;
						std::multimap< double, Edge* >::iterator eraseVesselIt;
						std::multimap< double, Edge* >::iterator eraseVesselIt2;
						while(tmp3DVesselIt != tmp3DVessel.end())
						{
							//implement search for vesselpoints in 'vessel' and remove them from the hash table
							for(eraseVesselIt = vessels.begin(); eraseVesselIt != vessels.end(); /*++eraseVesselIt*/)
							{
								double * tmp1 =  new double[3];
								double * tmp2 =  new double[3];
								tmp1 = eraseVesselIt->second->edgePointCoordinates.back();
								tmp2 = (*tmp3DVesselIt)->edgePointCoordinates.back();
								if(tmp1[0] == tmp2[0] && tmp1[1] == tmp2[1] && tmp1[2] == tmp2[2])
								{
									eraseVesselIt2 = eraseVesselIt;
									++eraseVesselIt;
									vessels.erase(eraseVesselIt2);
								}
								else
									++eraseVesselIt;
							}
							++tmp3DVesselIt;
						}
						vessels3D.push_back(tmp3DVessel);
						eraseVesselIt = planeVesselIter;
						++planeVesselIter;
						vessels.erase(eraseVesselIt);
						planeVessels = vessels.equal_range(*zIter);
					}
				}
			}
			
			if(notInVessel)
				++planeVesselIter;
		}
	}
	
	#ifdef DEBUG
	std::flush(std::cout << "converting vessels into 3D lines..." << std::endl);
	std::flush(std::cout << "nr. of vessels found in 3D: " << vessels3D.size() << std::endl);
	#endif
	std::list< std::vector< Edge * > >::iterator vessels3DIt;
	for(vessels3DIt = vessels3D.begin(); vessels3DIt != vessels3D.end(); ++vessels3DIt)
	{
		std::vector< Edge * >::iterator thisVesselIt = vessels3DIt->begin();
		std::vector< Edge * >::iterator lastVesselIt = thisVesselIt;
		double meanVessel[] = {0, 0, 0};
		++thisVesselIt;
		int zWeight = 1;	// deeper vessel points are weighted less/not at all
		double avgVesselRadius = 0;
		while(thisVesselIt != vessels3DIt->end() && zWeight <=12)
		{
			double thisVessel[3], lastVessel[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				thisVessel[ii] = (*thisVesselIt)->edgePointCoordinates.back()[ii];
				lastVessel[ii] = (*lastVesselIt)->edgePointCoordinates.back()[ii];
				meanVessel[ii] += thisVessel[ii] - lastVessel[ii];
			}
			avgVesselRadius += (*thisVesselIt)->radius;
			lastVesselIt = thisVesselIt;
			++thisVesselIt;
			++zWeight;
		}
		avgVesselRadius /= (double)zWeight;
		double norm = sqrt(meanVessel[0]*meanVessel[0] + meanVessel[1]*meanVessel[1] + meanVessel[2]*meanVessel[2]);
		for(int ii = 0; ii < 3; ++ii)
			meanVessel[ii] /= norm;
		double zDist = std::abs((*vessels3DIt->begin())->edgePointCoordinates.back()[2] - (*vessels3DIt->rbegin())->edgePointCoordinates.back()[2]);
		double * topCoords = new double[3];
		double * bottomCoords = new double[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			topCoords[ii] = (*vessels3DIt->begin())->edgePointCoordinates.back()[ii]/* - (std::abs((*vessels3DIt->begin())->edgePointCoordinates.back()[2])/meanVessel[2] + 250.0)*meanVessel[ii]*/;
			bottomCoords[ii] = (*vessels3DIt->begin())->edgePointCoordinates.back()[ii] + zDist/meanVessel[2]*meanVessel[ii];
		}
		
		Vertex *  top = new Vertex(topCoords, Vessel);
		Vertex * bottom = new Vertex(bottomCoords, Vessel);
		spatialGraph->addVertex(top);
		spatialGraph->addVertex(bottom);
		
		int connectionIndex[2] = {spatialGraph->getNumberOfVertices()-2, spatialGraph->getNumberOfVertices()-1};
		int noOfVesselPoints = 2;
		std::list< double * > vessel3DCoords;
		vessel3DCoords.push_back(topCoords);
		vessel3DCoords.push_back(bottomCoords);
		Edge * new3DVessel = new Edge(connectionIndex, noOfVesselPoints, Vessel, vessel3DCoords, avgVesselRadius);
		spatialGraph->addEdge(new3DVessel);
	}
}

/****************************************************************************/
/*vesselDistances2D() computes the 2D euclidian distance for all vessels in */
/*'vessels' from 'origin' and returns them as a hash table with the         */
/*(in ascending order) sorted distances as keys                             */
/****************************************************************************/
std::multimap< double, Edge* > CortexReconstruction::vesselDistances2D ( double origin[3], std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels )
{
	std::multimap< double, Edge* > nearestVessels;
	std::multimap< double, Edge* >::iterator nextPlaneVesselIter;
	for(nextPlaneVesselIter = vessels.first; nextPlaneVesselIter != vessels.second; ++nextPlaneVesselIter)
	{
		double tmp[3];
		for(int ii = 0; ii < 3; ++ii)
			tmp[ii] = nextPlaneVesselIter->second->edgePointCoordinates.back()[ii];
		double dist = sqrt((origin[0] - tmp[0])*(origin[0] - tmp[0]) + (origin[1] - tmp[1])*(origin[1] - tmp[1]));
		
		nearestVessels.insert(std::pair< double, Edge* >(dist, nextPlaneVesselIter->second));
	}
	
	return nearestVessels;
}

/****************************************************************************/
/*vesselDistances3D() computes the 3D euclidian distance for all vessels in */
/*'vessels' from 'origin' and returns them as a hash table with the         */
/*(in ascending order) sorted distances as keys                             */
/****************************************************************************/
std::multimap< double, Edge* > CortexReconstruction::vesselDistances3D ( double origin[3], std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels )
{
	std::multimap< double, Edge* > nearestVessels;
	std::multimap< double, Edge* >::iterator nextPlaneVesselIter;
// 	if(vessels.first == vessels.second)
// 		std::cout << "Error! no vessels in this plane!" << std::endl;
	for(nextPlaneVesselIter = vessels.first; nextPlaneVesselIter != vessels.second; ++nextPlaneVesselIter)
	{
		double tmp[3];
		for(int ii = 0; ii < 3; ++ii)
			tmp[ii] = nextPlaneVesselIter->second->edgePointCoordinates.back()[ii];
		double dist = sqrt((origin[0] - tmp[0])*(origin[0] - tmp[0]) + (origin[1] - tmp[1])*(origin[1] - tmp[1]) + (origin[2] - tmp[2])*(origin[2] - tmp[2]));
		nearestVessels.insert(std::pair< double, Edge* >(dist, nextPlaneVesselIter->second));
	}
// 	std::cout << "nearestVessels.size() = " << nearestVessels.size() << std::endl;
	return nearestVessels;
}

std::list< unsigned int > CortexReconstruction::computeConstrainingVessels ( PolyDataPointerType piaSurface )
{
	// ckeck for manually reconstructed 
	// blood vessels (only 2 pts per edge)
	bool manualVessels = 0;
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == Vessel && (*edgeIt)->edgePointCoordinates.size() == 2)
		{
			manualVessels = 1;
			#ifdef DEBUG
			std::cout << "Using manually reconstructed blood vessels for orientation!" << std::endl;
			#endif
			break;
		}
	}
	if(!manualVessels)
	{
		#ifdef DEBUG
		std::cout << "Using automatically reconstructed blood vessels for orientation!" << std::endl;
		#endif
		bool constrained = 1;
		double minRadius = 10;
		computeBloodVessels(constrained, minRadius);
	}
	
	std::list< unsigned int > goodVessels;
// 	std::list< int > vertexRemoveList;
	CellLocatorPointerType locator = CellLocatorPointerType::New();
	locator->AutomaticOn();
	locator->SetDataSet(piaSurface);
	locator->BuildLocator();
// 	locator->Print(std::cout);
	#ifdef DEBUG
	std::flush(std::cout << "Selecting good vessels" << std::endl);
	#endif
	if(spatialGraph->edgesPointer()->size())
	{
		for(unsigned int ii = 0; ii < spatialGraph->edgesPointer()->size(); ++ii)
		{
			if((*(spatialGraph->edgesPointer()))[ii]->label == Vessel && (*(spatialGraph->edgesPointer()))[ii]->numEdgePoints == 2)
			{
	// 			std::cout << "Calculating intersection of vessel with surface..." << std::endl;
				double * point1, * point2;
				point1 = (*(spatialGraph->edgesPointer()))[ii]->edgePointCoordinates.front();
				point2 = (*(spatialGraph->edgesPointer()))[ii]->edgePointCoordinates.back();
				double a0[3], a1[3], tol = 0.1, t, x[3], pcoords[3];
				int subId;
				vtkIdType cellID;
				GenericCellPointerType intersectCell = GenericCellPointerType::New();
				for(int jj = 0; jj < 3; ++jj)
				{
					a0[jj] = point1[jj];
					a1[jj] = point2[jj];
				}
				if(a0[2] < a1[2])
				{
					double direction[3];
					double dNorm = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						direction[jj] = a0[jj] - a1[jj];
						dNorm += direction[jj]*direction[jj];
					}
					dNorm = sqrt(dNorm);
					for(int jj = 0; jj < 3; ++jj)
					{
						direction[jj] /= dNorm;
						a0[jj] += direction[jj]*1000;
					}
				}
				else
				{
					double direction[3];
					double dNorm = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						direction[jj] = a1[jj] - a0[jj];
						dNorm += direction[jj]*direction[jj];
					}
					dNorm = sqrt(dNorm);
					for(int jj = 0; jj < 3; ++jj)
					{
						direction[jj] /= dNorm;
						a1[jj] += direction[jj]*1000;
					}
				}
// 				std::cout << "Vessel = [" << a0[0] << "," << a0[1] << "," << a0[2] << "] - [" << a1[0] << "," << a1[1] << "," << a1[2] << "]" << std::endl;
				int intersection = locator->IntersectWithLine(a0, a1, tol, t, x, pcoords, subId, cellID, intersectCell);
// 				std::cout << "intersection = " << intersection << std::endl;
				if(intersection)
				{
// 					#ifdef DEBUG
// 					std::flush(std::cout << "Intersection found!" << std::endl);
// 					#endif
					PointsPointerType cellPoints = intersectCell->GetPoints();
					double p1[3], p2[3], p3[3], normal[3];
					cellPoints->GetPoint(0, p1), cellPoints->GetPoint(1, p2), cellPoints->GetPoint(2, p3);
					for(int jj = 0; jj < 3; ++jj)
					{
						p3[jj] -= p1[jj];
						p2[jj] -= p1[jj];
					}
					normal[0] = p2[1]*p3[2] - p2[2]*p3[1];
					normal[1] = p2[2]*p3[0] - p2[0]*p3[2];
					normal[2] = p2[0]*p3[1] - p2[1]*p3[0];
					
					double angle = 0, normA = 0, normN = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						angle += (a0[jj] - a1[jj])*normal[jj];
						normA += (a0[jj] - a1[jj])*(a0[jj] - a1[jj]);
						normN += normal[jj]*normal[jj];
					}
					normA = sqrt(normA);
					normN = sqrt(normN);
					angle = acos(std::abs(angle)/(normA*normN))*180.0/PI;
					if(angle <= 10)
					{
						goodVessels.push_back(ii);
// 						spatialGraph->addEdge((*(spatialGraph->edgesPointer()))[ii]);
// 						(*(spatialGraph->edgesPointer()))[ii]->label = Axon;
// 						(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[ii]->edgeConnectivity[0]]->label = Axon;
// 						(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[ii]->edgeConnectivity[1]]->label = Axon;
					}
	// 				std::cout << "Surface normal @ intersection = [" << normal[0]/normN << "," << normal[1]/normN << "," << normal[2]/normN << "]" << std::endl;
				}
			}
		}
	}
	return goodVessels;
}

unsigned int CortexReconstruction::vesselsAroundBarrel ( double barrelCentroid[3], std::list< unsigned int > vessels )
{
	std::multimap< double, double * > vesselDistances;
	std::map< unsigned int, double * > neighborhoodVessels;
	std::list< unsigned int >::iterator vesselIt;
	double neighborhoodRadius = 350;
	while(!neighborhoodVessels.size() && neighborhoodRadius < 500)
	{
		for(vesselIt = vessels.begin(); vesselIt != vessels.end(); ++vesselIt)
		{
			double * point1, * point2;
			if(*vesselIt < 0 || *vesselIt >= spatialGraph->edgesPointer()->size())
			{
				std::cout << "Error! vessel list has invalid entries!" << std::endl;
				return 0;
// 				std::multimap< double, double * > emptyMap;
// 				return emptyMap;
			}
			if((*(spatialGraph->edgesPointer()))[*vesselIt]->edgePointCoordinates.size() != 2)
			{
				std::cout << "Error! vessel is not 3D vessel!" << std::endl;
				continue;
			}
			point1 = (*(spatialGraph->edgesPointer()))[*vesselIt]->edgePointCoordinates.front();
			point2 = (*(spatialGraph->edgesPointer()))[*vesselIt]->edgePointCoordinates.back();
			//make sure point 1 is closer to pia than point 2; if not, switch them
			if(point1[2] > point2[2])
			{
				double tmp[3];
				tmp[0] = point2[0], tmp[1] = point2[1], tmp[2] = point2[2];
				point2[0] = point1[0], point2[1] = point1[1], point2[2] = point1[2];
				point1[0] = tmp[0], point1[1] = tmp[1], point1[2] = tmp[2];
			}
			double * direction = new double[3];
			double dNorm = 0;
			for(int ii = 0; ii < 3; ++ii)
			{
				direction[ii] = point1[ii] - point2[ii];
				dNorm += direction[ii]*direction[ii];
			}
			dNorm = sqrt(dNorm);
			double barrelPlaneCoord[3];
			double zDist = point2[2] - barrelCentroid[2];	//order important!! in case point2[2] < barrelCentroid[2]
			for(int ii = 0; ii < 3; ++ii)
				barrelPlaneCoord[ii] = point2[ii] + direction[ii]/direction[2]*zDist;
			
			double planeDist = sqrt((barrelPlaneCoord[0] - barrelCentroid[0])*(barrelPlaneCoord[0] - barrelCentroid[0]) + (barrelPlaneCoord[1] - barrelCentroid[1])*(barrelPlaneCoord[1] - barrelCentroid[1]));
			vtkMath::Normalize(direction);
			vesselDistances.insert(std::pair< double, double * >(planeDist, direction));
			if(planeDist < neighborhoodRadius)
			{
// 				(*(spatialGraph->edgesPointer()))[*vesselIt]->label += 1;
// 				(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[*vesselIt]->edgeConnectivity[0]]->label += 1;
// 				(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[*vesselIt]->edgeConnectivity[1]]->label += 1;
				neighborhoodVessels.insert(std::pair< unsigned int, double * >(*vesselIt, direction));
			}
// 			delete [] direction;
		}
		neighborhoodRadius += 50;
	}
	
	unsigned int maxIndex = 0;
	if(neighborhoodVessels.size())
	{
		//normalization of vessel directions assumed here!!!
		double neighborhoodDirection[] = {0, 0, 0};
		std::map< unsigned int, double * >::const_iterator neighborhoodVesselIt;
		for(neighborhoodVesselIt = neighborhoodVessels.begin(); neighborhoodVesselIt != neighborhoodVessels.end(); ++neighborhoodVesselIt)
		{
			for(int ii = 0; ii < 3; ++ii)
				neighborhoodDirection[ii] += neighborhoodVesselIt->second[ii];
		}
// 		for(int ii = 0; ii < 3; ++ii)
// 			neighborhoodDirection[ii] = neighborhoodDirection[ii]/double(neighborhoodVessels.size());
		
		double maxScore = 0;
		for(neighborhoodVesselIt = neighborhoodVessels.begin(); neighborhoodVesselIt != neighborhoodVessels.end(); ++neighborhoodVesselIt)
		{
			double tmpScore = 0;
			for(int ii = 0; ii < 3; ++ii)
				tmpScore += neighborhoodVesselIt->second[ii]*neighborhoodDirection[ii];
// 			std::cout << "tmpScore = " << tmpScore << std::endl;
			if(tmpScore > maxScore)
			{
				maxScore = tmpScore;
				maxIndex = neighborhoodVesselIt->first;
			}
		}
		
// 		double * meanVessel = new double[3];
// 		meanVessel[0] = neighborhoodDirection[0], meanVessel[1] = neighborhoodDirection[1], meanVessel[2] = neighborhoodDirection[2];
// 		normalize(meanVessel);
// 		return meanVessel;
		
// 		double * endPoint = new double[3];
// 		for(int ii = 0; ii < 3; ++ii)
// 			endPoint[ii] = barrelCentroid[ii] + 100*neighborhoodDirection[ii];
// 		Vertex * newVert1 = new Vertex(endPoint, Barrel);
// 		Vertex * newVert2 = new Vertex(barrelCentroid, Barrel);
// 		spatialGraph->addVertex(newVert1);
// 		spatialGraph->addVertex(newVert2);
// 		int connectionIndex[2];
// 		if(!spatialGraph->getNumberOfVertices())
// 		{
// 			connectionIndex[0] = 0;
// 			connectionIndex[1] = 1;
// 		}
// 		else
// 		{
// 			connectionIndex[0] = spatialGraph->getNumberOfVertices() - 2;
// 			connectionIndex[1] = spatialGraph->getNumberOfVertices() - 1;
// 		}
// 		int noOfAxisPoints = 2;
// 		std::list< double * > axisCoords;
// 		axisCoords.push_back(endPoint);
// 		axisCoords.push_back(barrelCentroid);
// 		Edge * newAxis = new Edge(connectionIndex, noOfAxisPoints, Barrel, axisCoords/*maybe score*/);
// 		spatialGraph->addEdge(newAxis);
// 		std::cout << "maxIndex = " << maxIndex << std::endl;
// 		std::cout << "maxScore = " << maxScore << std::endl;
// 		(*(spatialGraph->edgesPointer()))[maxIndex]->label = Pia;
// 		(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[maxIndex]->edgeConnectivity[0]]->label = Pia;
// 		(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[maxIndex]->edgeConnectivity[1]]->label = Pia;
	}
	
	return maxIndex;
}

// PolyDataPointerType CortexReconstruction::surfaceReconstruction ( int label )
// {
// 	#ifdef DEBUG
// 	std::cout << "Reconstructing surface (label " << label << ")" << std::endl;
// 	#endif
// 	double spacing = 0;
// 	if(label == Pia)
// 		spacing = parameters.piaSpacing;
// 	else if(label == WhiteMatter)
// 		spacing = parameters.wmSpacing;
// 	ImageDataPointerType completeSurface = addTop2(label, 20, parameters.zReversed, spacing);
// 	MarchingCubesPointerType mcSurfaceFilter1 = MarchingCubesPointerType::New();
// 	mcSurfaceFilter1->SetInput(completeSurface);
// 	mcSurfaceFilter1->SetValue(0, 0);
// 	mcSurfaceFilter1->ComputeScalarsOff();
// 	mcSurfaceFilter1->ComputeGradientsOff();
// 	mcSurfaceFilter1->ComputeNormalsOff();
// 	mcSurfaceFilter1->Update();
// 	PolyDataPointerType smoothedSurface = smoothSurface(mcSurfaceFilter1->GetOutput());
// 	return smoothedSurface;
// }
// 
// ImageDataPointerType CortexReconstruction::addTop2 ( int label, int additionalSections, bool zReversed, double zSpacing )
// {
// 	ImageDataPointerType pia = piaVolume(label, additionalSections, zReversed, zSpacing);
// 	int extent[6] = {0, 0, 0, 0, 0, 0};	//minX, maxX, minY, maxY, minZ, maxZ
// 	pia->GetExtent(extent);
// 	int startZ;	//first section w/o data
// 	int stopZ = zReversed ? extent[4] : extent[5];
// // 	std::cout << "stopZ = " << stopZ << std::endl;
// 	if(zReversed)
// 	{
// 		startZ = extent[4] + additionalSections;
// // 		std::cout << "startZ = " << startZ << std::endl;
// 		if(startZ > extent[5] - 2)
// 		{
// 			std::cout << "Error! Pia volume incorrect: vertical extent smaller than data extent" << std::endl;
// 			return pia;
// 		}
// 	}
// 	else
// 	{
// 		startZ = extent[5] - additionalSections;
// // 		std::cout << "startZ = " << startZ << std::endl;
// 		if(startZ < extent[4] + 2)
// 		{
// 			std::cout << "Error! Pia volume incorrect: vertical extent smaller than data extent" << std::endl;
// 			return pia;
// 		}
// 	}
// 	int zCount = 1;
// 	int sign = zReversed ? -1 : 1;
// 	double curvature = zSpacing ? zSpacing/50.0 : 1;
// 	
// 	#ifdef REG_ACCURACY
// 	curvature *= var_gamma;
// 	#endif
// 	
// // 	if(L1flag)
// // 	{
// // 		// flat pia top
// // 		curvature *= 10000;
// // 		zCount += 1;
// // 		std::cout << ">>> Computing Pia top for L1 neurons..." << std::endl;
// // 	}
// 	
// // 	std::flush(std::cout << "sign = " << sign << std::endl);
// // 	std::flush(std::cout << "curvature = " << curvature << std::endl);
// // 	std::flush(std::cout << "startZ = " << startZ << std::endl);
// // 	std::flush(std::cout << "stopZ = " << stopZ << std::endl);
// 	for(int z = startZ; z != stopZ; z += sign)
// 	{
// 		for(int y = extent[2]; y <= extent[3]; ++y)
// 			for(int x = extent[0]; x <= extent[1]; ++x)
// 			{
// 				float * dist0 = static_cast< float * >(pia->GetScalarPointer(x, y, z));
// 				float * dist1 = static_cast< float * >(pia->GetScalarPointer(x, y, z-sign*1));
// 				float * dist2 = static_cast< float * >(pia->GetScalarPointer(x, y, z-sign*2));
// 				*dist0 = *dist1 + (*dist1 - *dist2) + (curvature*zCount)*(curvature*zCount);	// effectively calculates gradient of the distance image plus an artificial curvature
// 			}
// 		++zCount;
// 	}
// 	
// 	return pia;
// }
// 
// ImageDataPointerType CortexReconstruction::piaVolume ( int label, int additionalSections, bool zReversed, double zSpacing )
// {
// 	PolyDataPointerType polyData = PolyDataPointerType::New();
// 	if(spatialGraph->extractLandmark(label, polyData))
// 	{
// // 		polyData->Print(std::cout);
// // 		std::flush(std::cout << "Calculating isosurface..." << std::endl);
// 		double * bounds = polyData->GetBounds();
// 		int xMin = round(bounds[0]), xMax = round(bounds[1]);
// 		int yMin = round(bounds[2]), yMax = round(bounds[3]);
// 		int zMin = round(bounds[4]), zMax = round(bounds[5]);
// 		double zOffset = zSpacing ? zSpacing : 50;
// 		if(zReversed)
// 		{
// 			zMax -= zOffset;	//open @ bottom for pia & WM
// 			zMin -= zOffset*additionalSections;
// 		}
// 		else
// 		{
// 			zMax += zOffset*additionalSections;
// 			zMin += zOffset;	//open @ bottom for pia & WM
// 		}
// 		ImageDataPointerType volume = createImageVolumeFromPolyData(polyData, label, xMin, xMax, yMin, yMax, zMin, zMax, zSpacing);
// 		ImageDataPointerType distVolume = distanceTransform(volume);
// 		distVolume->SetSpacing(volume->GetSpacing());
// 		return distVolume;
// 	}
// 	
// 	else
// 	{
// 		std::cout << "Error! Empty SpatialGraph!" << std::endl;
// 		return NULL;
// 	}
// }
// 
// ImageDataPointerType CortexReconstruction::createImageVolumeFromPolyData ( PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, double zSpacing )
// {
// 	if(poly->GetNumberOfCells())
// 	{
// 		ImageDataPointerType volume = ImageDataPointerType::New();
// 		double spacing[3];
// 		if(!zSpacing)
// 			switch(label)
// 			{
// 				case Pia:
// 					spacing[0] = spacing[1] = spacing[2] = 50;
// 					break;
// 					
// 				case WhiteMatter:
// 					spacing[0] = spacing[1] = spacing[2] = 50;
// 					break;
// 					
// 				default:
// 					spacing[0] = spacing[1] = spacing[2] = 1;
// 					break;
// 			}
// 		else
// 			spacing[0] = spacing[1] = spacing[2] = zSpacing;
// 		
// 		volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
// 		int * dims = calculateExtent(xMin, xMax, yMin, yMax, zMin, zMax, label, spacing);
// // 		std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
// // 		std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
// 		volume->SetExtent(dims);
// 		volume->SetNumberOfScalarComponents(1);
// 		volume->SetScalarTypeToUnsignedChar();
// 		
// // 		poly->Print(std::cout);
// // 		volume->Print(std::cout);
// 		
// 		volume->AllocateScalars();
// 		for(int z = dims[4]; z <= dims[5]; ++z)
// 			for(int y = dims[2]; y <= dims[3]; ++y)
// 				for(int x = dims[0]; x <= dims[1]; ++x)
// 				{
// 					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
// 					*px = 0;
// 				}
// 		volume->Update();
// // 		std::flush(std::cout << "Calculating pixels inside polygons!" << std::endl);
// 		unsigned long insidePoints = 0;
// 		unsigned long outsidePoints = 0;
// 		for(int currPlane = 0; currPlane < poly->GetNumberOfCells(); ++currPlane)
// 		{
// 			double closestPoint[3];
// 			int subId;
// 			double pCoords[3];
// 			double * weights = new double[poly->GetCell(currPlane)->GetNumberOfPoints()];
// 			double dist2;
// 			double tmp[3];
// 			poly->GetCell(currPlane)->GetPoints()->GetPoint(0, tmp);
// 			double * bounds = poly->GetCell(currPlane)->GetBounds();
// // 			std::cout << "point tmp @ [" << tmp[0] << "," << tmp[1] << "," << tmp[2] << "]" << std::endl;
// 			int z = lround(tmp[2]/spacing[2]);
// // 			std::flush(std::cout << "Determining inside/outside polygon for all points plane z = " << z << std::endl);
// 			
// 			// Force all points of polygon to lie in the same plane
// 			// sometimes necessary b/c of some Amira imprecision artifacts...
// 			PolyDataPointerType thisPlanePoly = PolyDataPointerType::New();
// 			PolygonPointerType thisPolygon = PolygonPointerType::New();
// 			PointsPointerType thisPolyPoints = PointsPointerType::New();
// 			thisPlanePoly->Allocate(1);
// 			thisPolygon->GetPointIds()->SetNumberOfIds(poly->GetCell(currPlane)->GetNumberOfPoints());
// 			thisPolyPoints->SetDataTypeToFloat();
// 			thisPolyPoints->SetNumberOfPoints(poly->GetCell(currPlane)->GetNumberOfPoints());
// 			for(int ii = 0; ii < poly->GetCell(currPlane)->GetNumberOfPoints(); ++ii)
// 			{
// 				double tmp2[3];
// 				poly->GetCell(currPlane)->GetPoints()->GetPoint(ii, tmp2);
// 				tmp2[2] = 0;	// all in plane z = 0
// 				thisPolygon->GetPointIds()->SetId(ii, ii);
// 				thisPolyPoints->InsertPoint(ii, tmp2);
// 			}
// 			thisPlanePoly->InsertNextCell(thisPolygon->GetCellType(), thisPolygon->GetPointIds());
// 			thisPlanePoly->SetPoints(thisPolyPoints);
// 			thisPlanePoly->Update();
// 			
// // 			std::flush(std::cout << "physical z = " << tmp[2] << std::endl);
// // 			std::flush(std::cout << "z spacing = " << spacing[2] << std::endl);
// // 			std::flush(std::cout << "bounds = [" << bounds[0] << "," << bounds[1] << "], [" << bounds[2] << "," << bounds[3] << "],[" << bounds[4] << "," << bounds[5] << "]" << std::endl);
// // 			#pragma omp parallel for
// 			for(int y = dims[2]; y <= dims[3]; ++y)
// 			{
// 				for(int x = dims[0]; x <= dims[1]; ++x)
// 				{
// 					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
// // 					double tmpCoord[] = {x*spacing[0], y*spacing[1], z*spacing[2]};
// 					// brute force method in case pia/WM have different spacing/offset
// // 					double tmpCoord[] = {x*spacing[0], y*spacing[1], tmp[2]};
// 					double tmpCoord[] = {x*spacing[0], y*spacing[1], 0};	// check only in x-y direction!!!
// 					if(tmpCoord[0] < bounds[0] || tmpCoord[0] > bounds[1] 
// 						|| tmpCoord[1] < bounds[2] || tmpCoord[1] > bounds[3] 
// 						/*|| tmpCoord[2] < bounds[4] || tmpCoord[2] > bounds[5]*/)
// 					{
// // 						std::flush(std::cout << "out of bounds." << std::endl);
// 						++outsidePoints;
// 						*px = 0;
// 						continue;
// 					}
// // 					int insidePolygon = poly->GetCell(currPlane)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
// 					int insidePolygon = thisPlanePoly->GetCell(0)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
// 					if(insidePolygon == 1)
// 					{
// // 						std::flush(std::cout << "hit!" << std::endl);
// 						*px = 255;
// 						++insidePoints;
// 					}
// 					else
// 					{
// 						++outsidePoints;
// 						*px = 0;
// 					}
// 				}
// 			}
// // 			std::flush(std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl);
// 		}
// // 		std::flush(std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl);
// 		volume->Update();
// 		return volume;
// 	}
// 	else
// 	{
// 		std::cout << "Error! Empty PolyData! Cannot convert to vtkImageData!" << std::endl;
// 		return NULL;
// 	}
// }
// 
// ImageDataPointerType CortexReconstruction::distanceTransform ( ImageDataPointerType volume )
// {
// 	VTK2ITKImageExportPointerType v2iExport = VTK2ITKImageExportPointerType::New();
// 	ITK2VTKImageImportPointerType i2vImport = ITK2VTKImageImportPointerType::New();
// 	VTK2ITKImageImportType::Pointer v2iImport = VTK2ITKImageImportType::New();
// 	ITK2VTKCalcImageExportType::Pointer i2vExport = ITK2VTKCalcImageExportType::New();
// 	
// 	ImageType::Pointer volumeImage = ImageType::New();
// 	ImageType::Pointer contourImage = ImageType::New();
// 	CalcImageType::Pointer distanceMap = CalcImageType::New();
// 	DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
// 	
// 	v2iExport->SetInput(volume);
// 	v2iImport->SetCallbackUserData(v2iExport->GetCallbackUserData());
// 	v2iImport->SetDataExtentCallback(v2iExport->GetDataExtentCallback());
// 	v2iImport->SetNumberOfComponentsCallback(v2iExport->GetNumberOfComponentsCallback());
// 	v2iImport->SetOriginCallback(v2iExport->GetOriginCallback());
// 	v2iImport->SetPipelineModifiedCallback(v2iExport->GetPipelineModifiedCallback());
// 	v2iImport->SetPropagateUpdateExtentCallback(v2iExport->GetPropagateUpdateExtentCallback());
// 	v2iImport->SetScalarTypeCallback(v2iExport->GetScalarTypeCallback());
// 	v2iImport->SetSpacingCallback(v2iExport->GetSpacingCallback());
// 	v2iImport->SetUpdateDataCallback(v2iExport->GetUpdateDataCallback());
// 	v2iImport->SetUpdateInformationCallback(v2iExport->GetUpdateInformationCallback());
// 	v2iImport->SetWholeExtentCallback(v2iExport->GetWholeExtentCallback());
// 	v2iImport->SetBufferPointerCallback(v2iExport->GetBufferPointerCallback());
// 	
// 	volumeImage = v2iImport->GetOutput();
// 	volumeImage->Update();
// 	contourImage->SetRegions(volumeImage->GetLargestPossibleRegion());
// 	contourImage->Allocate();
// 	contourImage->FillBuffer(0);
// 	distanceMap->SetRegions(volumeImage->GetLargestPossibleRegion());
// 	distanceMap->Allocate();
// 	ImageType::SizeType radius;
// 	radius.Fill(1);
// 	SegNeighborhoodIteratorType volumeIter(radius, volumeImage, volumeImage->GetLargestPossibleRegion());
// 	SegNeighborhoodIteratorType contourIter(radius, contourImage, contourImage->GetLargestPossibleRegion());
// 	NeighborhoodOffsetVectorType offset = CreateLookUpTable();
// 	int neighborhood4[] = {10, 12, 13, 15};
// 	for(volumeIter.GoToBegin(), contourIter.GoToBegin(); !volumeIter.IsAtEnd() && !contourIter.IsAtEnd(); ++volumeIter, ++contourIter)
// 		if(volumeIter.GetCenterPixel())
// 		{
// 			contourIter.SetCenterPixel(volumeIter.GetCenterPixel());
// 			bool border[] = {0, 0, 0, 0};
// 			for(int ii = 0; ii < 4; ++ii)
// 				if(!volumeIter.GetPixel(offset[neighborhood4[ii]]))
// 					contourIter.SetPixel(offset[neighborhood4[ii]], volumeIter.GetCenterPixel());
// 		}
// 	contourImage->Update();
// 	
// 	int startZ = contourImage->GetLargestPossibleRegion().GetIndex()[2];
// 	int deltaZ = contourImage->GetLargestPossibleRegion().GetSize()[2];
// 	for(int z = startZ; z < (startZ + deltaZ); ++z)
// 	{
// 		ImageType::Pointer planeImage = ImageType::New();
// 		CalcImageType::Pointer planeDistanceMap = CalcImageType::New();
// 		ImageType::RegionType planeRegion;
// 		ImageType::RegionType planeImageRegion;
// 		ImageType::IndexType planeIndex;
// 		ImageType::SizeType planeSize;
// 		
// 		planeIndex = contourImage->GetLargestPossibleRegion().GetIndex();
// 		planeIndex[2] = z;
// 		planeSize = contourImage->GetLargestPossibleRegion().GetSize();
// 		planeSize[2] = 1;
// 		planeRegion.SetIndex(planeIndex);
// 		planeRegion.SetSize(planeSize);
// 		planeImageRegion.SetIndex(contourImage->GetLargestPossibleRegion().GetIndex());
// 		planeImageRegion.SetSize(planeSize);
// 		planeImage->SetRegions(planeImageRegion);
// 		planeImage->Allocate();
// 		
// 		ConstIteratorType copyIter(contourImage, planeRegion);
// 		IteratorType2 pasteIter(planeImage, planeImageRegion);
// 		for(copyIter.GoToBegin(), pasteIter.GoToBegin(); !copyIter.IsAtEnd() && !pasteIter.IsAtEnd(); ++copyIter, ++pasteIter)
// 			pasteIter.Set(copyIter.Get());
// 		planeImage->Update();
// 		
// 		distanceMapFilter->UseImageSpacingOff();
// 		distanceMapFilter->InsideIsPositiveOff();
// 		distanceMapFilter->SetInput(planeImage);
// 		planeDistanceMap = distanceMapFilter->GetDistanceMap();
// 		planeDistanceMap->Update();
// 		
// 		ConstCalcIteratorType copyIter2(planeDistanceMap, planeImageRegion);
// 		CalcIteratorType pasteIter2(distanceMap, planeRegion);
// 		for(copyIter2.GoToBegin(), pasteIter2.GoToBegin(); !copyIter2.IsAtEnd() && !pasteIter2.IsAtEnd(); ++copyIter2, ++pasteIter2)
// 			pasteIter2.Set(copyIter2.Get());
// 		distanceMap->Update();
// 	}
// 	
// 	i2vExport->SetInput(distanceMap);
// 	i2vImport->SetCallbackUserData(i2vExport->GetCallbackUserData());
// 	i2vImport->SetDataExtentCallback(i2vExport->GetDataExtentCallback());
// 	i2vImport->SetNumberOfComponentsCallback(i2vExport->GetNumberOfComponentsCallback());
// 	i2vImport->SetOriginCallback(i2vExport->GetOriginCallback());
// 	i2vImport->SetPipelineModifiedCallback(i2vExport->GetPipelineModifiedCallback());
// 	i2vImport->SetPropagateUpdateExtentCallback(i2vExport->GetPropagateUpdateExtentCallback());
// 	i2vImport->SetScalarTypeCallback(i2vExport->GetScalarTypeCallback());
// 	i2vImport->SetSpacingCallback(i2vExport->GetSpacingCallback());
// 	i2vImport->SetUpdateDataCallback(i2vExport->GetUpdateDataCallback());
// 	i2vImport->SetUpdateInformationCallback(i2vExport->GetUpdateInformationCallback());
// 	i2vImport->SetWholeExtentCallback(i2vExport->GetWholeExtentCallback());
// 	i2vImport->SetBufferPointerCallback(i2vExport->GetBufferPointerCallback());
// 	i2vImport->Update();
// 	
// 	ImageDataPointerType distanceVolume = ImageDataPointerType::New();
// 	distanceVolume->DeepCopy(i2vImport->GetOutput());
// 	distanceVolume->Update();
// 	return distanceVolume;
// }
// 
// int* CortexReconstruction::calculateExtent ( int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label )
// {
// 	int * extent = new int[6];
// 	double spacing[3];
// 	switch(label)
// 	{
// 		case Pia:
// 			spacing[0] = spacing[1] = spacing[2] = 50;
// 			break;
// 			
// 		case WhiteMatter:
// 			spacing[0] = spacing[1] = spacing[2] = 50;
// 			break;
// 			
// 		case Barrel:
// 			spacing[0] = spacing[1] = spacing[2] = 10;
// 			break;
// 			
// 		default:
// 			spacing[0] = spacing[1] = spacing[2] = 1;
// 			break;
// 	}
// 	
// 	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
// 	extent[0] = ((double)minX - spacing[0])/spacing[0]/* - 0.5*/;
// 	extent[1] = ((double)maxX + spacing[0])/spacing[0]/* + 0.5*/;
// 	extent[2] = ((double)minY - spacing[1])/spacing[1]/* - 0.5*/;
// 	extent[3] = ((double)maxY + spacing[1])/spacing[1]/* + 0.5*/;
// 	extent[4] = ((double)minZ - spacing[2])/spacing[2]/* - 0.5*/;
// 	extent[5] = ((double)maxZ + spacing[2])/spacing[2]/* + 0.5*/;
// 	
// 	return extent;
// }
// 
// int* CortexReconstruction::calculateExtent ( int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label, double spacing[3] )
// {
// 	int * extent = new int[6];
// 	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
// 	extent[0] = ((double)minX - spacing[0])/spacing[0]/* - 0.5*/;
// 	extent[1] = ((double)maxX + spacing[0])/spacing[0]/* + 0.5*/;
// 	extent[2] = ((double)minY - spacing[1])/spacing[1]/* - 0.5*/;
// 	extent[3] = ((double)maxY + spacing[1])/spacing[1]/* + 0.5*/;
// 	extent[4] = ((double)minZ - spacing[2])/spacing[2]/* - 0.5*/;
// 	extent[5] = ((double)maxZ + spacing[2])/spacing[2]/* + 0.5*/;
// 	// ugly fix in case minZ < 0... not sure if it's 100% bullet proof...
// 	if(minZ < 0 && abs(minZ)%lround(spacing[2]))
// 		extent[4] -= 1;
// 	
// 	return extent;
// }
// 
// NeighborhoodOffsetVectorType CortexReconstruction::CreateLookUpTable()
// {
// 	SegNeighborhoodIteratorType::OffsetType offset;
// 	NeighborhoodOffsetVectorType look_up_table(0);
// 	
// 	for(int z = -1; z <= 1 ; z++)
// 	{
// 		for(int x=-1; x<=1; x++)
// 		{
// 			for(int y=-1; y<=1; y++)
// 			{
// 				if( (x == 0) && (y == 0) && (z == 0) )
// 				{}
// 				else
// 				{
// 					offset[0]=x;
// 					offset[1]=y;
// 					offset[2]=z;
// 					look_up_table.push_back(offset);
// 				}
// 			}
// 		}
// 	}
// 	
// 	return look_up_table;
// }
// 
// PolyDataPointerType CortexReconstruction::smoothSurface ( PolyDataPointerType surface )
// {
// // 	AveragePolyDataFilterType smoothingFilter = AveragePolyDataFilterType::New();
// // 	smoothingFilter->BoundarySmoothingOff();
// // 	smoothingFilter->FeatureEdgeSmoothingOff();
// // // 	smoothingFilter->FeatureEdgeSmoothingOn();
// // // 	smoothingFilter->SetFeatureAngle(90);
// // 	smoothingFilter->SetConvergence(0.001);
// // 	smoothingFilter->SetNumberOfIterations(200);
// // 	smoothingFilter->SetInput(surface);
// // 	smoothingFilter->Update();
// // 	return smoothingFilter->GetOutput();
// 	
// 	LowpassPolyDataFilterType smoothingFilter = LowpassPolyDataFilterType::New();
// 	smoothingFilter->BoundarySmoothingOff();
// 	smoothingFilter->FeatureEdgeSmoothingOff();
// 	smoothingFilter->NormalizeCoordinatesOn();
// 	smoothingFilter->SetNumberOfIterations(20);
// 	smoothingFilter->SetPassBand(0.01);
// 	smoothingFilter->SetInput(surface);
// 	smoothingFilter->Update();
// 	return smoothingFilter->GetOutput();
// }
