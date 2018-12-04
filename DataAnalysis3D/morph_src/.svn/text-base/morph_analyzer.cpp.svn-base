#include "morph_analyzer.h"

// #define DEBUG
// #define BRANCHPOINTS

Analyzer::Analyzer(AmiraSpatialGraph * inputSpatialGraph, InputParameters parameters)
{
	SBF = new BarrelField(false);
	neuronMorphology = inputSpatialGraph;
	initializeConstants();
};

Analyzer::~Analyzer()
{
	if(neuronMorphology) delete neuronMorphology;
	if(SBF) delete SBF;
};

std::vector< double >* Analyzer::compute1DProfileLocally(const char * labelStr, double binSize)
{
	int label;
	std::string tmpStr(labelStr);
	std::map< std::string, int >::const_iterator compareIt;
	for(compareIt = neuronLabels2Int.begin(); compareIt != neuronLabels2Int.end(); ++compareIt)
	{
		if(tmpStr.compare(compareIt->first) == 0)
		{
			label = compareIt->second;
			break;
		}
	}
	if(compareIt == neuronLabels2Int.end())
	{
		std::cout << "Error! Structure label " << labelStr << " not recognized! Aborting..." << std::endl;
		std::vector< double >* emptyVec = new std::vector< double >;
		return emptyVec;
	}
	
	Profile * zProfile = new Profile(binSize);
	
	// start iteration through all edges
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != label)
			continue;
		
		// iterate through all points on this edge
		// calculate pairwise distance
		// subdivide if two points are in separate bins
		std::list< double * >::const_iterator ptIt1, ptIt2;
		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
		ptIt2 = ptIt1;
		++ptIt1;
		while(ptIt1 != (*edgeIt)->edgePointCoordinates.end())
		{
			double tmp1[] = {(*ptIt1)[0], (*ptIt1)[1], (*ptIt1)[2]};
			double tmp2[] = {(*ptIt2)[0], (*ptIt2)[1], (*ptIt2)[2]};
			double topPt[3], bottomPt[3], axis[3];
			SBF->localZAxis(tmp1, axis);
			SBF->avgPiaSurface->intersectLine(axis, tmp1);
			SBF->avgWMSurface->intersectLine(axis, tmp1);
			if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
			{
				SBF->avgPiaSurface->getLastIntersectPoint(topPt);
				SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
				assignDistanceToBins(tmp1, tmp2, topPt, bottomPt, axis, zProfile);
			}
			
			++ptIt1, ++ptIt2;
		}
	}
	
	return zProfile->getProfile();
};

std::vector< double >* Analyzer::compute1DProfileInsideS1(const char * labelStr, double binSize)
{
	int label;
	std::string tmpStr(labelStr);
	std::map< std::string, int >::const_iterator compareIt;
	for(compareIt = neuronLabels2Int.begin(); compareIt != neuronLabels2Int.end(); ++compareIt)
	{
		if(tmpStr.compare(compareIt->first) == 0)
		{
			label = compareIt->second;
			break;
		}
	}
	if(compareIt == neuronLabels2Int.end())
	{
		std::cout << "Error! Structure label " << labelStr << " not recognized! Aborting..." << std::endl;
		std::vector< double >* emptyVec = new std::vector< double >;
		return emptyVec;
	}

	Profile * zProfile = new Profile(binSize);

	// start iteration through all edges
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != label)
			continue;

		// iterate through all points on this edge
		// calculate pairwise distance
		// subdivide if two points are in separate bins
		std::list< double * >::const_iterator ptIt1, ptIt2;
		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
		ptIt2 = ptIt1;
		++ptIt1;
		while(ptIt1 != (*edgeIt)->edgePointCoordinates.end())
		{
			double tmp1[] = {(*ptIt1)[0], (*ptIt1)[1], (*ptIt1)[2]};
			double tmp2[] = {(*ptIt2)[0], (*ptIt2)[1], (*ptIt2)[2]};

			if (SBF->isInsideS1(tmp1) || SBF->isInsideS1(tmp2))
			{
				double topPt[3], bottomPt[3], axis[3];
				SBF->localZAxis(tmp1, axis);
				SBF->avgPiaSurface->intersectLine(axis, tmp1);
				SBF->avgWMSurface->intersectLine(axis, tmp1);
				if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt);
					SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
					assignDistanceToBins(tmp1, tmp2, topPt, bottomPt, axis, zProfile);
				}
			}
			++ptIt1, ++ptIt2;
		}
	}

	return zProfile->getProfile();
};

std::vector< double >* Analyzer::compute1DProfileGlobally(const char * labelStr, double binSize)
{
	int label;
	std::string tmpStr(labelStr);
	std::map< std::string, int >::const_iterator compareIt;
	for(compareIt = neuronLabels2Int.begin(); compareIt != neuronLabels2Int.end(); ++compareIt)
	{
		if(tmpStr.compare(compareIt->first) == 0)
		{
			label = compareIt->second;
			break;
		}
	}
	if(compareIt == neuronLabels2Int.end())
	{
		std::cout << "Error! Structure label " << labelStr << " not recognized! Aborting..." << std::endl;
		std::vector< double >* emptyVec = new std::vector< double >;
		return emptyVec;
	}
	
	bool hasHB, hasSoma;
	hasHB = neuronMorphology->getHomeBarrel();
	hasSoma = parameters.somaFlag;

	if(!hasHB && !hasSoma)
	{
		std::cout << "Error! No home barrel found in input SpatialGraph! Aborting..." << std::endl;
		std::vector< double > * emptyVec = new std::vector< double >;
		return emptyVec;
	}
	else if(hasSoma)
	{
		double somaPt[3];
		getPCenterOfStructure(neuronMorphology, Soma, somaPt);
		neuronMorphology->setHomeBarrel(SBF->closestBarrel(somaPt));
	}
	int HBID = neuronMorphology->getHomeBarrel();
	std::cout << "Computing z profile with respect to " << SBF->int2Labels[HBID] << std::endl;
	
	double topPt[3], bottomPt[3], axis[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		topPt[ii] = SBF->avgColumns[HBID]->top[ii];
		bottomPt[ii] = SBF->avgColumns[HBID]->bottom[ii];
	}
	SBF->localZAxis(HBID, axis);
	Profile * zProfile = new Profile(binSize);
	
	// start iteration through all edges
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != label)
			continue;
		
		// iterate through all points on this edge
		// calculate pairwise distance
		// subdivide if two points are in separate bins
		std::list< double * >::const_iterator ptIt1, ptIt2;
		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
		ptIt2 = ptIt1;
		++ptIt1;
		while(ptIt1 != (*edgeIt)->edgePointCoordinates.end())
		{
			double tmp1[] = {(*ptIt1)[0], (*ptIt1)[1], (*ptIt1)[2]};
			double tmp2[] = {(*ptIt2)[0], (*ptIt2)[1], (*ptIt2)[2]};
			assignDistanceToBins(tmp1, tmp2, topPt, bottomPt, axis, zProfile);
			
			++ptIt1, ++ptIt2;
		}
	}
	
	return zProfile->getProfile();
};

void Analyzer::computeAxonClusterParameters(const char* outputFilename, double binSize)
{
	// all lengths supra/gran/infra
	// same for branch points
	// for all columns and septum
	#ifdef DEBUG
	std::map< int, Profile * > onlyLayerProfiles;
	std::map< int, Profile * > onlyColumnProfiles;
	#endif
	std::map< int, std::map< int, Profile * > > columnProfiles;
	std::map< int, std::map< int, Profile * > > columnBPProfiles;
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
		std::map< int, Profile * > layerProfiles;
		std::map< int, Profile * > layerBPProfiles;
		for(int layer = 0; layer <= INFRA; ++layer)
		{
			Profile * newProfile = new Profile(binSize);
			Profile * newBPProfile = new Profile(binSize);
			layerProfiles.insert(std::pair< int, Profile * >(layer, newProfile));
			layerBPProfiles.insert(std::pair< int, Profile * >(layer, newBPProfile));
		}
		columnProfiles.insert(std::pair< int, std::map< int, Profile * > >(*labelIt, layerProfiles));
		columnBPProfiles.insert(std::pair< int, std::map< int, Profile * > >(*labelIt, layerBPProfiles));
		#ifdef DEBUG
		Profile * columnProfile = new Profile(binSize);
		onlyColumnProfiles.insert(std::pair< int, Profile * >(*labelIt, columnProfile));
		#endif
	}
	std::map< int, Profile * > septumLayerProfiles;
	std::map< int, Profile * > septumLayerBPProfiles;
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		Profile * newProfile = new Profile(binSize);
		Profile * newBPProfile = new Profile(binSize);
		septumLayerProfiles.insert(std::pair< int, Profile * >(layer, newProfile));
		septumLayerBPProfiles.insert(std::pair< int, Profile * >(layer, newBPProfile));
		#ifdef DEBUG
		Profile * layerProfile = new Profile(binSize);
		onlyLayerProfiles.insert(std::pair< int, Profile * >(layer, layerProfile));
		#endif
	}
	columnProfiles.insert(std::pair< int, std::map< int, Profile * > >(Septum, septumLayerProfiles));
	columnBPProfiles.insert(std::pair< int, std::map< int, Profile * > >(Septum, septumLayerBPProfiles));
	#ifdef DEBUG
	Profile * septumProfile = new Profile(binSize);
	onlyColumnProfiles.insert(std::pair< int, Profile * >(Septum, septumProfile));
	#endif
	Profile * outS1Profile = new Profile(binSize);
	Profile * outS1BPProfile = new Profile(binSize);
	
// 	#ifdef DEBUG
// 	PointsPointerType S1Pts = PointsPointerType::New();
// 	S1Pts->SetDataTypeToFloat();
// 	#endif
	#ifdef DEBUG
	PointsPointerType centerPts = PointsPointerType::New();
	centerPts->SetDataTypeToFloat();
	#endif
	
	// main loop
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != Axon)
			continue;
		
		bool oldPtS1;
		int oldPtLayer, oldPtColumn;
		std::list< double * >::const_iterator ptIt1, ptIt2;
		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
		ptIt2 = (*edgeIt)->edgePointCoordinates.begin();
		++ptIt1;
		oldPtS1 = SBF->isInsideS1(*ptIt2);
		oldPtLayer = SBF->laminarPosition(*ptIt2);
		oldPtColumn = SBF->insideColumn(*ptIt2);
		while(ptIt1 != (*edgeIt)->edgePointCoordinates.end())
		{
#ifdef DEBUG
			std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
#endif
			double * pt1, * pt2;
			bool pt1S1, pt2S1;
			int pt1Layer, pt2Layer;
			int pt1Column, pt2Column;
			
			pt1 = *ptIt1;
			pt2 = *ptIt2;
			pt1S1 = SBF->isInsideS1(pt1);
			pt2S1 = oldPtS1;
			pt1Layer = SBF->laminarPosition(pt1);
			pt2Layer = oldPtLayer;
			pt1Column = SBF->insideColumn(pt1);
			pt2Column = oldPtColumn;
			
			#ifdef DEBUG
			assignDistanceToLayers(pt1, pt2, pt1Layer, pt2Layer, onlyLayerProfiles, onlyLayerProfiles);
			#endif
			
			// one pt inside S1 is approximately enough
			// (i.e., relative length error is small)
			if(pt1S1 || pt2S1)
			{
// 				#ifdef DEBUG
// 				S1Pts->InsertNextPoint(pt1);
// 				S1Pts->InsertNextPoint(pt2);
// 				#endif
				// case 1 : both in same column/septum
				if(pt1Column == pt2Column)
				{
					assignDistanceToLayers(pt1, pt2, pt1Layer, pt2Layer, columnProfiles[pt1Column], columnProfiles[pt2Column]);
					#ifdef DEBUG
					double topPt[3], bottomPt[3], axis[3];
					SBF->localZAxis(pt1, axis);
					SBF->avgPiaSurface->intersectLine(axis, pt1);
					SBF->avgWMSurface->intersectLine(axis, pt1);
					if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
					{
						SBF->avgPiaSurface->getLastIntersectPoint(topPt);
						SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
						assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, onlyColumnProfiles[pt1Column]);
					}
					#endif
				}
				// case 2: points in different columns/septum
				else
				{
					// one point in septum, other in column:
					// find pt on column border on line
					// between pt1 and pt2
					if(!pt1Column || !pt2Column)
					{
						if(pt1Column)
						{
							double line12[3], rVec[3], r, R, angle, rNorm, lineDist;
							double t, closestPt[3], centerPt[3];
							int centerLayer;
							r = sqrt(vtkLine::DistanceToLine(pt1, SBF->avgColumns[pt1Column]->top, SBF->avgColumns[pt1Column]->bottom, t, closestPt));
							R = sqrt(SBF->avgBarrelArea[pt1Column]/PI);
							vtkMath::Subtract(pt1, closestPt, rVec);
							vtkMath::Subtract(pt2, pt1, line12);
							rNorm = vtkMath::Normalize(rVec);
							lineDist = vtkMath::Normalize(line12);
							angle = vtkMath::Dot(rVec, line12);
							
							for(int ii = 0; ii < 3; ++ii)
							{
								centerPt[ii] = pt1[ii] + (R-r)/angle*line12[ii];
							}
							// in case of long straight lines perpendicular to the radial position:
							double distance1, distance2;
							distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
							distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
							if(distance1+distance2 > lineDist+1E-06)
							{
								for(int ii = 0; ii < 3; ++ii)
								{
									centerPt[ii] = pt1[ii] + 0.5*lineDist*line12[ii];
								}
								#ifdef DEBUG
								distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
								distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
								#endif
							}
							
							centerLayer = SBF->laminarPosition(centerPt);
							assignDistanceToLayers(pt1, centerPt, pt1Layer, centerLayer, columnProfiles[pt1Column], columnProfiles[pt1Column]);
							assignDistanceToLayers(pt2, centerPt, pt2Layer, centerLayer, columnProfiles[pt2Column], columnProfiles[pt2Column]);
							
							#ifdef DEBUG
							double topPt[3], bottomPt[3], axis[3];
							SBF->localZAxis(pt1, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt1);
							SBF->avgWMSurface->intersectLine(axis, pt1);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt1, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt1Column]);
							}
							SBF->localZAxis(pt2, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt2);
							SBF->avgWMSurface->intersectLine(axis, pt2);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt2, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt2Column]);
							}
							double dist12;
							dist12 = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
							if(distance1+distance2 > dist12+1E-06)
							{
								centerPts->InsertNextPoint(pt1);
								centerPts->InsertNextPoint(pt2);
								centerPts->InsertNextPoint(centerPt);
								std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
								std::cout << "distance1 = " << distance1 << std::endl;
								std::cout << "distance2 = " << distance2 << std::endl;
								std::cout << "dist12 = " << dist12 << std::endl;
								std::cout << "pt1 @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
								std::cout << "pt2 @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
								std::cout << "rVec = [" << rVec[0] << "," << rVec[1] << "," << rVec[2] << "]" << std::endl;
								std::cout << "line12 = [" << line12[0] << "," << line12[1] << "," << line12[2] << "]" << std::endl;
								std::cout << "r = " << r << std::endl;
								std::cout << "R = " << R << std::endl;
								std::cout << "angle = " << angle << std::endl;
								std::cout << "(R-r)/angle = " << (R-r)/angle << std::endl;
								std::cout << "centerPt @ [" << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << "]" << std::endl;
								std::cout << "pt1Column = " << pt1Column << std::endl;
								std::cout << "pt2Column = " << pt2Column << std::endl;
								std::cout << "pt1Layer = " << pt1Layer << std::endl;
								std::cout << "pt2Layer = " << pt2Layer << std::endl;
								std::cout << "centerLayer = " << centerLayer << std::endl;
								std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
							}
// 							else
// 							{
// 								std::cout << "-----------------------------" << std::endl;
// 								std::cout << "-----------------------------" << std::endl;
// 							}
							#endif
						}
						else if(pt2Column)
						{
							double line12[3], rVec[3], r, R, angle, rNorm, lineDist;
							double t, closestPt[3], centerPt[3];
							int centerLayer;
							r = sqrt(vtkLine::DistanceToLine(pt2, SBF->avgColumns[pt2Column]->top, SBF->avgColumns[pt2Column]->bottom, t, closestPt));
							R = sqrt(SBF->avgBarrelArea[pt2Column]/PI);
							vtkMath::Subtract(pt2, closestPt, rVec);
							vtkMath::Subtract(pt1, pt2, line12);
							rNorm = vtkMath::Normalize(rVec);
							lineDist = vtkMath::Normalize(line12);
							angle = vtkMath::Dot(rVec, line12);
							
							for(int ii = 0; ii < 3; ++ii)
							{
								centerPt[ii] = pt2[ii] + (R-r)/angle*line12[ii];
							}
							// in case of long straight lines perpendicular to the radial position:
							double distance1, distance2;
							distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
							distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
							if(distance1+distance2 > lineDist+1E-06)
							{
								for(int ii = 0; ii < 3; ++ii)
								{
									centerPt[ii] = pt2[ii] + 0.5*lineDist*line12[ii];
								}
								#ifdef DEBUG
								distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
								distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
								#endif
							}
							
							centerLayer = SBF->laminarPosition(centerPt);
							assignDistanceToLayers(pt1, centerPt, pt1Layer, centerLayer, columnProfiles[pt1Column], columnProfiles[pt1Column]);
							assignDistanceToLayers(pt2, centerPt, pt2Layer, centerLayer, columnProfiles[pt2Column], columnProfiles[pt2Column]);
							
							#ifdef DEBUG
							double topPt[3], bottomPt[3], axis[3];
							SBF->localZAxis(pt1, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt1);
							SBF->avgWMSurface->intersectLine(axis, pt1);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt1, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt1Column]);
							}
							SBF->localZAxis(pt2, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt2);
							SBF->avgWMSurface->intersectLine(axis, pt2);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt2, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt2Column]);
							}
							double dist12;
							dist12 = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
							if(distance1+distance2 > dist12+1E-06)
							{
								centerPts->InsertNextPoint(pt1);
								centerPts->InsertNextPoint(pt2);
								centerPts->InsertNextPoint(centerPt);
								std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
								std::cout << "distance1 = " << distance1 << std::endl;
								std::cout << "distance2 = " << distance2 << std::endl;
								std::cout << "dist12 = " << dist12 << std::endl;
								std::cout << "pt1 @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
								std::cout << "pt2 @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
								std::cout << "rVec = [" << rVec[0] << "," << rVec[1] << "," << rVec[2] << "]" << std::endl;
								std::cout << "line12 = [" << line12[0] << "," << line12[1] << "," << line12[2] << "]" << std::endl;
								std::cout << "r = " << r << std::endl;
								std::cout << "R = " << R << std::endl;
								std::cout << "angle = " << angle << std::endl;
								std::cout << "(R-r)/angle = " << (R-r)/angle << std::endl;
								std::cout << "centerPt @ [" << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << "]" << std::endl;
								std::cout << "pt1Column = " << pt1Column << std::endl;
								std::cout << "pt2Column = " << pt2Column << std::endl;
								std::cout << "pt1Layer = " << pt1Layer << std::endl;
								std::cout << "pt2Layer = " << pt2Layer << std::endl;
								std::cout << "centerLayer = " << centerLayer << std::endl;
								std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
							}
							#endif
						}
					}
					// both points in columns; check if septum inbetween columns
					// or columns already overlapping
					else
					{
						double r1, R1, r2, R2, dist12, angle1, angle2;
						double t, closestPt1[3], closestPt2[3];
						double line12[3], rVec1[3], rVec2[3];
						double centerPt1[3], centerPt2[3];
						int centerLayer1, centerLayer2;
						r1 = sqrt(vtkLine::DistanceToLine(pt1, SBF->avgColumns[pt1Column]->top, SBF->avgColumns[pt1Column]->bottom, t, closestPt1));
						R1 = sqrt(SBF->avgBarrelArea[pt1Column]/PI);
						r2 = sqrt(vtkLine::DistanceToLine(pt2, SBF->avgColumns[pt2Column]->top, SBF->avgColumns[pt2Column]->bottom, t, closestPt2));
						R2 = sqrt(SBF->avgBarrelArea[pt2Column]/PI);
						vtkMath::Subtract(pt2, pt1, line12);
						dist12 = vtkMath::Normalize(line12);
						vtkMath::Subtract(pt1, closestPt1, rVec1);
						vtkMath::Subtract(pt2, closestPt2, rVec2);
						vtkMath::Normalize(rVec1);
						vtkMath::Normalize(rVec2);
						angle1 = vtkMath::Dot(rVec1, line12);
						angle2 = vtkMath::Dot(rVec2, line12);
						for(int ii = 0; ii < 3; ++ii)
						{
							centerPt1[ii] = pt1[ii] + (R1-r1)/angle1*line12[ii];
							centerPt2[ii] = pt2[ii] - (R2-r2)/angle2*line12[ii];
						}
						double distCenter1, distCenter2;
						distCenter1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt1));
						distCenter2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt2));
						// overlapping
						// split up between both columns for simplicity
						// no large systematic error
// 						if(distCenter1 > dist12 || distCenter2 > dist12)
// 						{
							double centerPt[3];
							int centerLayer;
							for(int ii = 0; ii < 3; ++ii)
								centerPt[ii] = pt1[ii] + 0.5*dist12*line12[ii];
							
// 							#ifdef DEBUG
// 							centerPts->InsertNextPoint(centerPt);
// 							#endif
							
							centerLayer = SBF->laminarPosition(centerPt);
							assignDistanceToLayers(pt1, centerPt, pt1Layer, centerLayer, columnProfiles[pt1Column], columnProfiles[pt1Column]);
							assignDistanceToLayers(pt2, centerPt, pt2Layer, centerLayer, columnProfiles[pt2Column], columnProfiles[pt2Column]);
							#ifdef DEBUG
							double topPt[3], bottomPt[3], axis[3];
							SBF->localZAxis(pt1, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt1);
							SBF->avgWMSurface->intersectLine(axis, pt1);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt1, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt1Column]);
							}
							SBF->localZAxis(pt2, axis);
							SBF->avgPiaSurface->intersectLine(axis, pt2);
							SBF->avgWMSurface->intersectLine(axis, pt2);
							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
							{
								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
								assignDistanceToBins(pt2, centerPt, topPt, bottomPt, axis, onlyColumnProfiles[pt2Column]);
							}
							#endif
// 						}
						// septum between columns
// 						else
// 						{
// 							
// // 							#ifdef DEBUG
// // 							centerPts->InsertNextPoint(centerPt1);
// // 							centerPts->InsertNextPoint(centerPt2);
// // 							#endif
// 							
// 							centerLayer1 = SBF->laminarPosition(centerPt1);
// 							centerLayer2 = SBF->laminarPosition(centerPt2);
// 							assignDistanceToLayers(pt1, centerPt1, pt1Layer, centerLayer1, columnProfiles[pt1Column], columnProfiles[pt1Column]);
// 							assignDistanceToLayers(pt2, centerPt2, pt2Layer, centerLayer2, columnProfiles[pt2Column], columnProfiles[pt2Column]);
// 							assignDistanceToLayers(centerPt1, centerPt2, centerLayer1, centerLayer2, columnProfiles[Septum], columnProfiles[Septum]);
// 							#ifdef DEBUG
// 							double topPt[3], bottomPt[3], axis[3];
// 							SBF->localZAxis(pt1, axis);
// 							SBF->avgPiaSurface->intersectLine(axis, pt1);
// 							SBF->avgWMSurface->intersectLine(axis, pt1);
// 							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
// 							{
// 								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
// 								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
// 								assignDistanceToBins(pt1, centerPt1, topPt, bottomPt, axis, onlyColumnProfiles[pt1Column]);
// 							}
// 							SBF->localZAxis(pt2, axis);
// 							SBF->avgPiaSurface->intersectLine(axis, pt2);
// 							SBF->avgWMSurface->intersectLine(axis, pt2);
// 							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
// 							{
// 								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
// 								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
// 								assignDistanceToBins(pt2, centerPt2, topPt, bottomPt, axis, onlyColumnProfiles[pt2Column]);
// 							}
// 							SBF->localZAxis(centerPt1, axis);
// 							SBF->avgPiaSurface->intersectLine(axis, centerPt1);
// 							SBF->avgWMSurface->intersectLine(axis, centerPt1);
// 							if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
// 							{
// 								SBF->avgPiaSurface->getLastIntersectPoint(topPt);
// 								SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
// 								assignDistanceToBins(centerPt1, centerPt2, topPt, bottomPt, axis, onlyColumnProfiles[Septum]);
// 							}
// 							#endif
// 						}
					}
				}
			}
			// both points outside S1 convex hull
			else
			{
				// old version: no checking whether points below standard surfaces
				double topPt[3], bottomPt[3], axis[3];
// 				SBF->localZAxis(pt1, axis);
// 				SBF->avgPiaSurface->intersectLine(axis, pt1);
// 				SBF->avgWMSurface->intersectLine(axis, pt1);
// 				if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
// 				{
// 					SBF->avgPiaSurface->getLastIntersectPoint(topPt);
// 					SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
// 					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, outS1Profile);
// 					#ifdef DEBUG
// 					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, onlyColumnProfiles[Septum]);
// 					#endif
// 				}
				
				// new version: includes checking whether points below standard surfaces
				SBF->localZAxis(pt1, axis);
				SBF->avgPiaSurface->intersectLine(axis, pt1);
				if(SBF->avgPiaSurface->isIntersectionFound())
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt);
					for(int ii = 0; ii < 3; ++ii)
					{
						bottomPt[ii] = topPt[ii] - 4000*axis[ii];
					}
					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, outS1Profile);
					#ifdef DEBUG
					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, onlyColumnProfiles[Septum]);
					#endif
				}
				//heuristic safeguard in case points are outside of interpolated surfaces
				//maybe fix better by interpolating surfaces further
				else
				{
				//in case there is structure outside of standard surfaces -> use closest pt
					double closestPt1[3], closestPt2[3], correction[3] = {0,0,0};
					//use closestPt for all calculations and simply shift everything back with correction vector at the end
					closestPt1[0] = pt1[0], closestPt1[1] = pt1[1], closestPt1[2] = pt1[2];
					SBF->localZAxis(pt1, closestPt1, axis);
					for(int ii = 0; ii < 3; ++ii)
					{
						correction[ii] = pt1[ii] - closestPt1[ii];
					}
					SBF->avgPiaSurface->intersectLine(axis, closestPt1);
					if(SBF->avgPiaSurface->isIntersectionFound())
					{
						SBF->avgPiaSurface->getLastIntersectPoint(topPt);
						for(int ii = 0; ii < 3; ++ii)
						{
							bottomPt[ii] = topPt[ii] - 4000*axis[ii];
							closestPt2[ii] = pt2[ii] - correction[ii];
						}
						assignDistanceToBins(closestPt1, closestPt2, topPt, bottomPt, axis, outS1Profile);
#ifdef DEBUG
						std::cout << "---------------------------------" << std::endl;
						std::cout << "Found segment outside of standard surfaces!" << std::endl;
						std::cout << "Original length:\t" << sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2)) << std::endl;
						std::cout << "Projected length:\t" << sqrt(vtkMath::Distance2BetweenPoints(closestPt1, closestPt2)) << std::endl;
						std::cout << "---------------------------------" << std::endl;
#endif
					}
					else
					{
						std::cout << "Warning! Could not determine depth location of axon segment!" << std::endl;
						std::cout << "Adding segment to first bin below pia..." << std::endl;
						double dist2 = vtkMath::Distance2BetweenPoints(pt1, pt2);
						outS1Profile->addSegment(sqrt(dist2), 0);
					}
				}
			}
			
			oldPtS1 = pt1S1;
			oldPtLayer = pt1Layer;
			oldPtColumn = pt1Column;
			++ptIt1, ++ptIt2;
		}
	}
#ifdef DEBUG
		std::cout << "\nDone processing axon length!" << std::endl;
#endif
	
#ifdef BRANCHPOINTS
	//loop through branch points
	std::list< int > branchPtList = getBranchPointIDs(Axon);
#ifdef DEBUG
		std::flush(std::cout << "Nr. of branch points = " << branchPtList.size() << std::endl);
#endif
	std::list< int >::const_iterator bpListIt;
	for(bpListIt = branchPtList.begin(); bpListIt != branchPtList.end(); ++bpListIt)
	{
		double * bp = (*(neuronMorphology->edgesPointer()))[*bpListIt]->edgePointCoordinates.back();
		int layer, column;
		bool inS1;
		Profile * locationProfile;
		inS1 = SBF->isInsideS1(bp);
		if(inS1)
		{
			layer = SBF->laminarPosition(bp);
			column = SBF->insideColumn(bp);
			locationProfile = columnBPProfiles[column][layer];
		}
		else
		{
			locationProfile = outS1BPProfile;
		}
		double topPt[3], axis[3], binSize;
		binSize = locationProfile->getBinSize();
		SBF->localZAxis(bp, axis);
		SBF->avgPiaSurface->intersectLine(axis, bp);
		if(SBF->avgPiaSurface->isValid())
		{
			SBF->avgPiaSurface->getLastIntersectPoint(topPt);
			double depth = sqrt(vtkMath::Distance2BetweenPoints(bp, topPt));
			unsigned int bin = (unsigned int)(depth/binSize);	// round down
			locationProfile->incrementBin(bin);
		}
	}
#endif
	
	// all parameters to be computed:
	int HC = -1;
	if(neuronMorphology->isLabelInSpatialGraph(Soma))
	{
		PolyDataPointerType soma = PolyDataPointerType::New();
		if(!neuronMorphology->extractLandmark(Soma, soma))
		{
			std::cout << "Error! Could not find soma in SpatialGraph!" << std::endl;
			return;
		}
		int subID;
		double somaPos[3], pCoords[3], * weights;
		weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
		soma->GetCell(0)->GetParametricCenter(pCoords);
		soma->GetCell(0)->EvaluateLocation(subID, pCoords, somaPos, weights);
		delete [] weights;
		
		HC = SBF->closestBarrel(somaPos);
#ifdef DEBUG
		std::flush(std::cout << "Home column: " << SBF->int2Labels[HC] << std::endl);
#endif
	}
	
	double totalLength = 0, S1Length = 0, HCLength = 0, otherColLength = 0, septumLength = 0;
	double onlyLayerTotal = 0, onlyColumnTotal = 0;
	double totalVolume = 0;
#ifdef BRANCHPOINTS
	unsigned int nrBranchPoints = 0;
	nrBranchPoints = branchPtList.size();
#endif
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		#ifdef DEBUG
		onlyLayerTotal += onlyLayerProfiles[layer]->getIntegral();
		#endif
		if(HC > 0)
			HCLength += columnProfiles[HC][layer]->getIntegral();
		septumLength += columnProfiles[Septum][layer]->getIntegral();
		S1Length += columnProfiles[Septum][layer]->getIntegral();
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			if(*labelIt == C5 || *labelIt == C6
				|| *labelIt == D5 || *labelIt == D6
				|| *labelIt == E5 || *labelIt == E6)
			{
				continue;
			}
			int column = *labelIt;
			S1Length += columnProfiles[column][layer]->getIntegral();
			if(column != HC)
				otherColLength += columnProfiles[column][layer]->getIntegral();
			#ifdef DEBUG
			if(!layer)
				onlyColumnTotal += onlyColumnProfiles[column]->getIntegral();
			#endif
		}
	}
	double axonBounds[6];
	neuronMorphology->getBoundingBox(Axon, axonBounds);
	totalVolume = (axonBounds[1] - axonBounds[0])*(axonBounds[3] - axonBounds[2])*(axonBounds[5] - axonBounds[4]);
	
	#ifdef DEBUG
	onlyColumnTotal += onlyColumnProfiles[Septum]->getIntegral();
	#endif
	totalLength = S1Length + outS1Profile->getIntegral();
	std::string summaryOutName(outputFilename);
	summaryOutName += "_axon_length_summary.csv";
	std::ofstream summaryFile;
	summaryFile.open(summaryOutName.c_str());
	std::cout << "*************************************************" << std::endl;
	#ifdef DEBUG
	std::cout << "Output filename:        " << outputFilename << std::endl;
	#endif
#ifdef BRANCHPOINTS
	std::cout << "Nr of branch points:    " << nrBranchPoints << std::endl;
	summaryFile << "Nr of branch points:\t" << nrBranchPoints << std::endl;
#endif
	#ifdef DEBUG
	std::cout << "Only layer axon length: " << onlyLayerTotal << std::endl;
	std::cout << "Only col axon length:   " << onlyColumnTotal << std::endl;
	#endif
	std::cout << "Total axon length:      " << totalLength << std::endl;
	summaryFile << "Total axon length:\t" << totalLength << std::endl;
	std::cout << "S1 axon length:         " << S1Length << std::endl;
	summaryFile << "S1 axon length:\t" << S1Length << std::endl;
	if(HC > 0)
	{
		std::cout << "PC (" << SBF->int2Labels[HC] << ") length:         " << HCLength << std::endl;
		summaryFile << "PC (" << SBF->int2Labels[HC] << ") length:\t" << HCLength << std::endl;
	}
	std::cout << "Other columns length:   " << otherColLength << std::endl;
	summaryFile << "Other columns length:\t" << otherColLength << std::endl;
	std::cout << "Septum axon length:     " << septumLength << std::endl;
	summaryFile << "Septum axon length:\t" << septumLength << std::endl;
	std::cout << "Out S1 axon length:     " << outS1Profile->getIntegral() << std::endl;
	summaryFile << "Out S1 axon length:\t" << outS1Profile->getIntegral() << std::endl;
#ifdef DEBUG
	std::cout << "Bounding box: " << axonBounds[0] << " " << axonBounds[1] << " " << axonBounds[2] << " " << axonBounds[3] << " " << axonBounds[4] << " " << axonBounds[5] << std::endl;
#endif
	std::cout << "Total volume(BBox):     " << totalVolume << std::endl;
	summaryFile << "Total volume(BBox):\t" << totalVolume << std::endl;
	std::cout << "*************************************************" << std::endl;
	summaryFile.close();
	#ifdef DEBUG
	std::string centerPtsOutName(outputFilename);
	centerPtsOutName += "_center_pts";
	Reader * centerLandmarkWriter = new Reader(centerPtsOutName.c_str(), centerPtsOutName.c_str());
	centerLandmarkWriter->writeLandmarkFile(centerPts);
	delete centerLandmarkWriter;
// 	std::string S1PtsOutName(outputFilename);
// 	S1PtsOutName += "_S1_pts";
// 	Reader * S1LandmarkWriter = new Reader(S1PtsOutName.c_str(), S1PtsOutName.c_str());
// 	S1LandmarkWriter->writeLandmarkFile(S1Pts);
// 	delete S1LandmarkWriter;
	#endif
	
	#ifndef DEBUG
	// write output files
	
	Profile * inS1Profile = new Profile(binSize);
	Profile * inHCProfile = new Profile(binSize);
	Profile * totalProfile = new Profile(binSize);
	
	std::string axonOutName(outputFilename);
	axonOutName += "_axon_length.csv";
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
			
			inS1Profile->addProfile(columnProfiles[*labelIt][layer]);
			totalProfile->addProfile(columnProfiles[*labelIt][layer]);
			
			AxonFile << columnProfiles[*labelIt][layer]->getIntegral() << "\t";
		}
		
		if(HC > 0)
			inHCProfile->addProfile(columnProfiles[HC][layer]);
		inS1Profile->addProfile(columnProfiles[Septum][layer]);
		totalProfile->addProfile(columnProfiles[Septum][layer]);
		
		AxonFile << columnProfiles[Septum][layer]->getIntegral() << "\t";
		
		if(!layer)
		{
			totalProfile->addProfile(outS1Profile);
			AxonFile << outS1Profile->getIntegral();
			delete outS1Profile;
		}
		AxonFile << std::endl;
	}
	AxonFile.close();
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			if(*labelIt == C5 || *labelIt == C6
				|| *labelIt == D5 || *labelIt == D6
				|| *labelIt == E5 || *labelIt == E6)
			{
				continue;
			}
			delete columnProfiles[*labelIt][layer];
		}
		delete columnProfiles[Septum][layer];
	}
	
	std::string axonS1OutName(outputFilename);
	axonS1OutName += "_axon_S1_profile.csv";
	inS1Profile->writeProfile(axonS1OutName.c_str(), 0.5*binSize);
	delete inS1Profile;
	
	std::string axonHCOutName(outputFilename);
	axonHCOutName += "_axon_PC_profile.csv";
	inHCProfile->writeProfile(axonHCOutName.c_str(), 0.5*binSize);
	delete inHCProfile;
	
	std::string axonTotalOutName(outputFilename);
	axonTotalOutName += "_axon_total_profile.csv";
	totalProfile->writeProfile(axonTotalOutName.c_str(), 0.5*binSize);
	delete totalProfile;

#ifdef BRANCHPOINTS
	std::string bpOutName(outputFilename);
	bpOutName += "_branchpoints.csv";
	std::ofstream BPFile;
	BPFile.open(bpOutName.c_str());
	BPFile << "Layer\t";
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
		BPFile << SBF->int2Labels[*labelIt] << "\t";
	}
	BPFile << "Septum\t";
	BPFile << "Outside S1";
	BPFile << std::endl;
	for(int layer = 0; layer <= INFRA; ++layer)
	{
		if(!layer)
			BPFile << "Other\t";
		if(layer == SUPRA)
			BPFile << "Supra\t";
		if(layer == GRAN)
			BPFile << "Gran\t";
		if(layer == INFRA)
			BPFile << "Infra\t";
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
			BPFile << columnBPProfiles[*labelIt][layer]->getIntegral() << "\t";
			delete columnBPProfiles[*labelIt][layer];
		}
		
		BPFile << columnBPProfiles[Septum][layer]->getIntegral() << "\t";
		delete columnBPProfiles[Septum][layer];
		
		if(!layer)
		{
			BPFile << outS1BPProfile->getIntegral();
			delete outS1BPProfile;
		}
		BPFile << std::endl;
	}
	BPFile.close();
#endif
	#endif
}

void Analyzer::assignDistanceToLayers ( double pt1[3], double pt2[3], int pt1Layer, int pt2Layer, std::map< int, Profile* >& pt1Profile, std::map< int, Profile* >& pt2Profile )
{
	if(pt1Layer == pt2Layer)
	{
		double topPt[3], bottomPt[3], axis[3];
		SBF->localZAxis(pt1, axis);
		SBF->avgPiaSurface->intersectLine(axis, pt1);
// 		SBF->avgWMSurface->intersectLine(axis, pt1);
		if(SBF->avgPiaSurface->isValid()/* && SBF->avgWMSurface->isValid()*/)
		{
			SBF->avgPiaSurface->getLastIntersectPoint(topPt);
// 			SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
			for(int ii = 0; ii < 3; ++ii)
			{
				bottomPt[ii] = topPt[ii] - 4000*axis[ii];
			}
			assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, pt1Profile[pt1Layer]);
		}
	}
	// assume they're in adjacent layers
	else
	{
		if(pt1Layer && pt2Layer)
		{
			double line12[3], centerPt[3];
			double topPt1[3], topPt2[3], bottomPt1[3], bottomPt2[3], axis1[3], axis2[3];
			Surface * boundary;
			if(pt1Layer == SUPRA && pt2Layer == GRAN || pt1Layer == GRAN && pt2Layer == SUPRA)
				boundary = SBF->avgL4UpperSurface;
			else if(pt1Layer == GRAN && pt2Layer == INFRA || pt1Layer == INFRA && pt2Layer == GRAN)
				boundary = SBF->avgL4LowerSurface;
			else if(pt1Layer == SUPRA && pt2Layer == INFRA || pt1Layer == INFRA && pt2Layer == SUPRA)
			{
				
				vtkMath::Subtract(pt1, pt2, line12);
				vtkMath::Normalize(line12);
				SBF->avgL4UpperSurface->intersectLine(line12, pt2);
				SBF->avgL4LowerSurface->intersectLine(line12, pt2);
				if(SBF->avgL4UpperSurface->isIntersectionFound() && SBF->avgL4LowerSurface->isIntersectionFound())
				{
					double newPt1[3], newPt2[3];
					SBF->avgL4UpperSurface->getLastIntersectPoint(newPt1);
					SBF->avgL4LowerSurface->getLastIntersectPoint(newPt2);
					
					if(pt1Layer == SUPRA && pt2Layer == INFRA)
					{
						assignDistanceToLayers(pt1, newPt1, pt1Layer, GRAN, pt1Profile, pt1Profile);
						assignDistanceToLayers(newPt1, newPt2, GRAN, GRAN, pt1Profile, pt1Profile);
						assignDistanceToLayers(newPt2, pt2, GRAN, pt2Layer, pt2Profile, pt2Profile);
					}
					else if(pt2Layer == SUPRA && pt1Layer == INFRA)
					{
						assignDistanceToLayers(pt2, newPt1, pt2Layer, GRAN, pt2Profile, pt2Profile);
						assignDistanceToLayers(newPt1, newPt2, GRAN, GRAN, pt2Profile, pt2Profile);
						assignDistanceToLayers(newPt2, pt1, GRAN, pt1Layer, pt1Profile, pt1Profile);
					}
					return;
				}
				// in the unlikely case this fails, just define the boundary to be L4Lower
				else
				{
					std::cout << "Warning! Could not split up axon segment between all three layers." << std::endl;
					std::cout << "Splitting segment at granular/infragranular boundary instead..." << std::endl;
					boundary = SBF->avgL4LowerSurface;
				}
			}
			vtkMath::Subtract(pt1, pt2, line12);
			vtkMath::Normalize(line12);
			SBF->localZAxis(pt1, axis1);
			SBF->localZAxis(pt2, axis2);
			boundary->intersectLine(line12, pt2);
			if(boundary->isValid())
			{
				boundary->getLastIntersectPoint(centerPt);
				
				double distance1, distance2, dist12;
				distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
				distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
				dist12 = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
				if(distance1+distance2 > dist12+1E-06)
				{
					for(int ii = 0; ii < 3; ++ii)
					{
						centerPt[ii] = pt2[ii] + 0.5*dist12*line12[ii];
					}
					#ifdef DEBUG
					distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, centerPt));
					distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, centerPt));
					#endif
				}
				#ifdef DEBUG
				if(distance1+distance2 > dist12+1E-06)
				{
// 					centerPts->InsertNextPoint(pt1);
// 					centerPts->InsertNextPoint(pt2);
// 					centerPts->InsertNextPoint(centerPt);
					std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
					std::cout << "distance1 = " << distance1 << std::endl;
					std::cout << "distance2 = " << distance2 << std::endl;
					std::cout << "dist12 = " << dist12 << std::endl;
					std::cout << "pt1 @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
					std::cout << "pt2 @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
					std::cout << "line12 = [" << line12[0] << "," << line12[1] << "," << line12[2] << "]" << std::endl;
					std::cout << "centerPt @ [" << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << "]" << std::endl;
					std::cout << "pt1Layer = " << pt1Layer << std::endl;
					std::cout << "pt2Layer = " << pt2Layer << std::endl;
					std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
				}
				#endif
				
				SBF->avgPiaSurface->intersectLine(axis1, pt1);
// 				SBF->avgWMSurface->intersectLine(axis1, pt1);
				if(SBF->avgPiaSurface->isValid()/* && SBF->avgWMSurface->isValid()*/)
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt1);
// 					SBF->avgWMSurface->getLastIntersectPoint(bottomPt1);
					for(int ii = 0; ii < 3; ++ii)
					{
						bottomPt1[ii] = topPt1[ii] - 4000*axis1[ii];
					}
					assignDistanceToBins(pt1, centerPt, topPt1, bottomPt1, axis1, pt1Profile[pt1Layer]);
				}
				SBF->avgPiaSurface->intersectLine(axis2, pt2);
// 				SBF->avgWMSurface->intersectLine(axis2, pt2);
				if(SBF->avgPiaSurface->isValid()/* && SBF->avgWMSurface->isValid()*/)
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt2);
// 					SBF->avgWMSurface->getLastIntersectPoint(bottomPt2);
					for(int ii = 0; ii < 3; ++ii)
					{
						bottomPt2[ii] = topPt2[ii] - 4000*axis2[ii];
					}
					assignDistanceToBins(pt2, centerPt, topPt2, bottomPt2, axis2, pt2Profile[pt2Layer]);
				}
			}
		}
		// if one of them is outside, just assign
		// everything to the layer of the other point
		else if(!pt1Layer || !pt2Layer)
		{
			if(!pt2Layer)
			{
				double topPt[3], bottomPt[3], axis[3];
				SBF->localZAxis(pt1, axis);
				SBF->avgPiaSurface->intersectLine(axis, pt1);
				SBF->avgWMSurface->intersectLine(axis, pt1);
				if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt);
					SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, pt1Profile[pt1Layer]);
				}
			}
			if(!pt1Layer)
			{
				double topPt[3], bottomPt[3], axis[3];
				SBF->localZAxis(pt2, axis);
				SBF->avgPiaSurface->intersectLine(axis, pt2);
				SBF->avgWMSurface->intersectLine(axis, pt2);
				if(SBF->avgPiaSurface->isValid() && SBF->avgWMSurface->isValid())
				{
					SBF->avgPiaSurface->getLastIntersectPoint(topPt);
					SBF->avgWMSurface->getLastIntersectPoint(bottomPt);
					assignDistanceToBins(pt1, pt2, topPt, bottomPt, axis, pt2Profile[pt2Layer]);
				}
			}
		}
	}
}

void Analyzer::assignDistanceToBins(double pt1[3], double pt2[3], double topPt[3], double bottomPt[3], double axis[3], Profile * zProfile)
{
	double binSize = zProfile->getBinSize();
	double tmp1Proj[3], tmp2Proj[3], t1, t2;
	vtkLine::DistanceToLine(pt1, topPt, bottomPt, t1, tmp1Proj);
	vtkLine::DistanceToLine(pt2, topPt, bottomPt, t2, tmp2Proj);
	double depth1 = sqrt(vtkMath::Distance2BetweenPoints(tmp1Proj, topPt));
	double depth2 = sqrt(vtkMath::Distance2BetweenPoints(tmp2Proj, topPt));
	unsigned int bin1 = (unsigned int)(depth1/binSize);	// round down
	unsigned int bin2 = (unsigned int)(depth2/binSize);
	
	if(bin1 == bin2)
	{
		zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2)), bin1);
	}
	else if(bin1 > bin2)
	{
		double direction[3], angle;
		for(int ii = 0; ii < 3; ++ii)
			direction[ii] = pt2[ii] - pt1[ii];
		vtkMath::Normalize(direction);
		angle = vtkMath::Dot(direction, axis);
		
		// beware of overestimating nearly parallel sections!
		// angle < 1E-7 together with a really small distance to bin border
		// may be numerically unfavorable; better be safe than sorry here!
		if(angle < 1E-7)
			zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2)), bin1);
		else
		{
			int binDiff = bin1 - bin2;
			double lastDepth, lastPt[3], borderPt[3];
			lastDepth = depth1;
			for(int ii = 0; ii < 3; ++ii)
				lastPt[ii] = pt1[ii];
			for(int n = 0; n < binDiff; ++n)
			{
				double binDist = 0;
				if(n == 0)
					binDist = lastDepth - bin1*binSize;
				else
					binDist = binSize;
				for(int ii = 0; ii < 3; ++ii)
					borderPt[ii] = lastPt[ii] + binDist/angle*direction[ii];
				
				zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(lastPt, borderPt)), (bin1 - n));
				
				for(int ii = 0; ii < 3; ++ii)
					lastPt[ii] = borderPt[ii];
				double tmp3Proj[3], t3;
				vtkLine::DistanceToLine(lastPt, topPt, bottomPt, t3, tmp3Proj);
				lastDepth = sqrt(vtkMath::Distance2BetweenPoints(tmp3Proj, topPt));
			}
			zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(pt2, lastPt)), bin2);
		}
	}
	else if(bin2 > bin1)
	{
		double direction[3], angle;
		for(int ii = 0; ii < 3; ++ii)
			direction[ii] = pt1[ii] - pt2[ii];
		vtkMath::Normalize(direction);
		angle = vtkMath::Dot(direction, axis);
		
		// beware of overestimating nearly parallel sections!
		// angle < 1E-7 together with a really small distance to bin border
		// may be numerically unfavorable; better be safe than sorry here!
		if(angle < 1E-7)
			zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2)), bin1);
		else
		{
			int binDiff = bin2 - bin1;
			double lastDepth, lastPt[3], borderPt[3];
			lastDepth = depth2;
			for(int ii = 0; ii < 3; ++ii)
				lastPt[ii] = pt2[ii];
			for(int n = 0; n < binDiff; ++n)
			{
				double binDist = 0;
				if(n == 0)
					binDist = lastDepth - bin2*binSize;
				else
					binDist = binSize;
				for(int ii = 0; ii < 3; ++ii)
					borderPt[ii] = lastPt[ii] + binDist/angle*direction[ii];
				
				zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(lastPt, borderPt)), (bin2 - n));
				
				for(int ii = 0; ii < 3; ++ii)
					lastPt[ii] = borderPt[ii];
				double tmp3Proj[3], t3;
				vtkLine::DistanceToLine(lastPt, topPt, bottomPt, t3, tmp3Proj);
				lastDepth = sqrt(vtkMath::Distance2BetweenPoints(tmp3Proj, topPt));
			}
			zProfile->addSegment(sqrt(vtkMath::Distance2BetweenPoints(pt1, lastPt)), bin1);
		}
	}
};

std::list< int > Analyzer::getBranchPointIDs ( int label )
{
	std::map< int, int > nrOfChildBranches;
	std::list< int > branchPtIDs;
	
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == label)
		{
			int fatherID = (*edgeIt)->fatherID;
			if(nrOfChildBranches.find(fatherID) != nrOfChildBranches.end())
			{
				nrOfChildBranches[fatherID]++;
			}
			else
			{
				nrOfChildBranches.insert(std::pair< int, int >(fatherID, 1));
			}
		}
	}
	
	std::map< int, int >::const_iterator nrBranchIt;
	for(nrBranchIt = nrOfChildBranches.begin(); nrBranchIt != nrOfChildBranches.end(); ++nrBranchIt)
	{
		if(nrBranchIt->second > 1)
		{
			branchPtIDs.push_back(nrBranchIt->first);
		}
	}
	
	return branchPtIDs;
}

void Analyzer::getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3])
{
	PolyDataPointerType structure = PolyDataPointerType::New();
	if(!sg->extractLandmark(ID, structure))
	{
		std::cout << "Error! Could not find structure with ID " << ID <<" in SpatialGraph!" << std::endl;
		return;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
	structure->GetCell(0)->GetParametricCenter(pCoords);
	structure->GetCell(0)->EvaluateLocation(subID, pCoords, centerPt, weights);
	delete [] weights;
};

void Analyzer::initializeConstants()
{
	if(neuronLabels2Int.size())
		neuronLabels2Int.clear();
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("Neuron"), Neuron));
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("Dendrite"), Dendrite));
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("ApicalDendrite"), ApicalDendrite));
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("BasalDendrite"), BasalDendrite));
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("Axon"), Axon));
	neuronLabels2Int.insert(std::pair< std::string, int >(std::string("Soma"), Soma));
};
