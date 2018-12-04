/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      geometry.cpp                                                  */
/*                                                                          */
/* Purpose:   class providing all methods for computing geometric           */
/*            parameters and structures: 3D vessels, surfaces, barrel       */
/*            parameters                                                    */
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
#include "geometry.h"

Geometry::Geometry(AmiraSpatialGraph * inputSpatialGraph)
{
	spatialGraph = inputSpatialGraph;
	barrelLabels.push_back(Alpha);
	barrelLabels.push_back(A1);
	barrelLabels.push_back(A2);
	barrelLabels.push_back(A3);
	barrelLabels.push_back(A4);
	barrelLabels.push_back(Beta);
	barrelLabels.push_back(B1);
	barrelLabels.push_back(B2);
	barrelLabels.push_back(B3);
	barrelLabels.push_back(B4);
	barrelLabels.push_back(Gamma);
	barrelLabels.push_back(C1);
	barrelLabels.push_back(C2);
	barrelLabels.push_back(C3);
	barrelLabels.push_back(C4);
	barrelLabels.push_back(C5);
	barrelLabels.push_back(C6);
	barrelLabels.push_back(Delta);
	barrelLabels.push_back(D1);
	barrelLabels.push_back(D2);
	barrelLabels.push_back(D3);
	barrelLabels.push_back(D4);
	barrelLabels.push_back(D5);
	barrelLabels.push_back(D6);
	barrelLabels.push_back(E1);
	barrelLabels.push_back(E2);
	barrelLabels.push_back(E3);
	barrelLabels.push_back(E4);
	barrelLabels.push_back(E5);
	barrelLabels.push_back(E6);
	int2Labels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2Labels.insert(std::pair< int, const char * >(A1, "A1"));
	int2Labels.insert(std::pair< int, const char * >(A2, "A2"));
	int2Labels.insert(std::pair< int, const char * >(A3, "A3"));
	int2Labels.insert(std::pair< int, const char * >(A4, "A4"));
	int2Labels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2Labels.insert(std::pair< int, const char * >(B1, "B1"));
	int2Labels.insert(std::pair< int, const char * >(B2, "B2"));
	int2Labels.insert(std::pair< int, const char * >(B3, "B3"));
	int2Labels.insert(std::pair< int, const char * >(B4, "B4"));
	int2Labels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2Labels.insert(std::pair< int, const char * >(C1, "C1"));
	int2Labels.insert(std::pair< int, const char * >(C2, "C2"));
	int2Labels.insert(std::pair< int, const char * >(C3, "C3"));
	int2Labels.insert(std::pair< int, const char * >(C4, "C4"));
	int2Labels.insert(std::pair< int, const char * >(C5, "C5"));
	int2Labels.insert(std::pair< int, const char * >(C6, "C6"));
	int2Labels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2Labels.insert(std::pair< int, const char * >(D1, "D1"));
	int2Labels.insert(std::pair< int, const char * >(D2, "D2"));
	int2Labels.insert(std::pair< int, const char * >(D3, "D3"));
	int2Labels.insert(std::pair< int, const char * >(D4, "D4"));
	int2Labels.insert(std::pair< int, const char * >(D5, "D5"));
	int2Labels.insert(std::pair< int, const char * >(D6, "D6"));
	int2Labels.insert(std::pair< int, const char * >(E1, "E1"));
	int2Labels.insert(std::pair< int, const char * >(E2, "E2"));
	int2Labels.insert(std::pair< int, const char * >(E3, "E3"));
	int2Labels.insert(std::pair< int, const char * >(E4, "E4"));
	int2Labels.insert(std::pair< int, const char * >(E5, "E5"));
	int2Labels.insert(std::pair< int, const char * >(E6, "E6"));
};

Geometry::Geometry()
{
	barrelLabels.push_back(Alpha);
	barrelLabels.push_back(A1);
	barrelLabels.push_back(A2);
	barrelLabels.push_back(A3);
	barrelLabels.push_back(A4);
	barrelLabels.push_back(Beta);
	barrelLabels.push_back(B1);
	barrelLabels.push_back(B2);
	barrelLabels.push_back(B3);
	barrelLabels.push_back(B4);
	barrelLabels.push_back(Gamma);
	barrelLabels.push_back(C1);
	barrelLabels.push_back(C2);
	barrelLabels.push_back(C3);
	barrelLabels.push_back(C4);
	barrelLabels.push_back(C5);
	barrelLabels.push_back(C6);
	barrelLabels.push_back(Delta);
	barrelLabels.push_back(D1);
	barrelLabels.push_back(D2);
	barrelLabels.push_back(D3);
	barrelLabels.push_back(D4);
	barrelLabels.push_back(D5);
	barrelLabels.push_back(D6);
	barrelLabels.push_back(E1);
	barrelLabels.push_back(E2);
	barrelLabels.push_back(E3);
	barrelLabels.push_back(E4);
	barrelLabels.push_back(E5);
	barrelLabels.push_back(E6);
	int2Labels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2Labels.insert(std::pair< int, const char * >(A1, "A1"));
	int2Labels.insert(std::pair< int, const char * >(A2, "A2"));
	int2Labels.insert(std::pair< int, const char * >(A3, "A3"));
	int2Labels.insert(std::pair< int, const char * >(A4, "A4"));
	int2Labels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2Labels.insert(std::pair< int, const char * >(B1, "B1"));
	int2Labels.insert(std::pair< int, const char * >(B2, "B2"));
	int2Labels.insert(std::pair< int, const char * >(B3, "B3"));
	int2Labels.insert(std::pair< int, const char * >(B4, "B4"));
	int2Labels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2Labels.insert(std::pair< int, const char * >(C1, "C1"));
	int2Labels.insert(std::pair< int, const char * >(C2, "C2"));
	int2Labels.insert(std::pair< int, const char * >(C3, "C3"));
	int2Labels.insert(std::pair< int, const char * >(C4, "C4"));
	int2Labels.insert(std::pair< int, const char * >(C5, "C5"));
	int2Labels.insert(std::pair< int, const char * >(C6, "C6"));
	int2Labels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2Labels.insert(std::pair< int, const char * >(D1, "D1"));
	int2Labels.insert(std::pair< int, const char * >(D2, "D2"));
	int2Labels.insert(std::pair< int, const char * >(D3, "D3"));
	int2Labels.insert(std::pair< int, const char * >(D4, "D4"));
	int2Labels.insert(std::pair< int, const char * >(D5, "D5"));
	int2Labels.insert(std::pair< int, const char * >(D6, "D6"));
	int2Labels.insert(std::pair< int, const char * >(E1, "E1"));
	int2Labels.insert(std::pair< int, const char * >(E2, "E2"));
	int2Labels.insert(std::pair< int, const char * >(E3, "E3"));
	int2Labels.insert(std::pair< int, const char * >(E4, "E4"));
	int2Labels.insert(std::pair< int, const char * >(E5, "E5"));
	int2Labels.insert(std::pair< int, const char * >(E6, "E6"));
};

Geometry::~Geometry()
{
	barrelLabels.clear();
	int2Labels.clear();
};

/******************************************************************************/
/*computes 3d blood vessels from vessel obects in spatialGraph by starting    */
/*from all vessels in the top plane and looking for vessels in a certain      */
/*radius in the next plane. If found, look for more vessels in deeper planes  */
/*in the direction specified by the first two vessels within tolerance        */
/******************************************************************************/
void Geometry::computeBloodVessels(bool constrained)
{
	double piaCenter[3];
	if(constrained)
	{
		PolyDataPointerType pia = PolyDataPointerType::New();
		std::list< int > zIndexList;
		spatialGraph->extractLandmark(Pia, pia, zIndexList);
		double piaPCenter[3], * piaWeights =  new double[pia->GetCell(0)->GetNumberOfPoints()];
		int piaSubID;
		pia->GetCell(0)->GetParametricCenter(piaPCenter);
		pia->GetCell(0)->EvaluateLocation(piaSubID, piaPCenter, piaCenter, piaWeights);
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
	std::flush(std::cout << "Sorting vessels by their z-position..." << std::endl);
	std::multimap< double, Edge* > vessels;	// sort vessels by their z-position
	std::list< double > zPlanes;
	std::list< std::vector< Edge * > > vessels3D;
	std::vector< Edge * >::iterator graphIt;
	for(graphIt = spatialGraph->edgesBegin(); graphIt != spatialGraph->edgesEnd(); ++graphIt)
	{
		if((*graphIt)->label == Vessel)
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
							normalize(tmpVesselDirection), normalize(tmpCenterDirection);
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
	std::flush(std::cout << "converting vessels into 3D lines..." << std::endl);
	std::list< std::vector< Edge * > >::iterator vessels3DIt;
	std::flush(std::cout << "nr. of vessels found in 3D: " << vessels3D.size() << std::endl);
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
};

/****************************************************************************/
/*vesselDistances3D() computes the 3D euclidian distance for all vessels in */
/*'vessels' from 'origin' and returns them as a hash table with the         */
/*(in ascending order) sorted distances as keys                             */
/****************************************************************************/
std::multimap< double, Edge* > Geometry::vesselDistances3D(double * origin, std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels)
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
};

/****************************************************************************/
/*vesselDistances2D() computes the 2D euclidian distance for all vessels in */
/*'vessels' from 'origin' and returns them as a hash table with the         */
/*(in ascending order) sorted distances as keys                             */
/****************************************************************************/
std::multimap< double, Edge* > Geometry::vesselDistances2D(double * origin, std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels)
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
};

/****************************************************************************/
/*computes a histogram of angles between the reconstructed blood vessels    */
/*and the surface normals at the intersection points between pia and vessel */
/****************************************************************************/
void Geometry::computeVesselToNormalAngleHistogram(const char * graphFName, const char * surfFName, const char * oFName, double binSize)
{
	Reader * graphReader = new Reader(graphFName);
	graphReader->readSpatialGraphFile(0);
	spatialGraph = graphReader->getSpatialGraph();
	PolyDataPointerType pia = PolyDataPointerType::New();
	std::list< int > zIndexList;
	if(spatialGraph->extractLandmark(Pia, pia, zIndexList))
	{
		double piaPCenter[3], piaCenter[3], * piaWeights =  new double[pia->GetCell(0)->GetNumberOfPoints()];
		int piaSubID;
		pia->GetCell(0)->GetParametricCenter(piaPCenter);
		pia->GetCell(0)->EvaluateLocation(piaSubID, piaPCenter, piaCenter, piaWeights);
		
		computeBloodVessels(1);
		
		Reader * surfaceReader = new Reader(surfFName);
		PolyDataPointerType surface = surfaceReader->readAmiraSurfaceFile();
		CellLocatorPointerType locator = CellLocatorPointerType::New();
		locator->AutomaticOn();
		locator->SetDataSet(surface);
		locator->BuildLocator();
	// 	locator->Print(std::cout);
		
		unsigned int noOfBins = (unsigned int)(90/binSize + 0.5) + 1;
		unsigned int * histogram = new unsigned int[noOfBins];
		for(int ii = 0; ii < noOfBins; ++ii)
			histogram[ii] = 0;
		
		std::string ofString(oFName);
		ofString += "_vessels.csv";
		std::ofstream VesselData( ofString.c_str() );
		VesselData << "vessel-surface intersection\tvessel direction\tsurface normal\tradial coord of intersection pt" << std::endl;
		
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label == Vessel && (*edgeIt)->numEdgePoints == 2)
			{
	// 			std::cout << "Calculating intersection of vessel with surface..." << std::endl;
				double * point1, * point2;
				point1 = (*edgeIt)->edgePointCoordinates.front();
				point2 = (*edgeIt)->edgePointCoordinates.back();
				double a0[3], a1[3], tol = 0.1, t, x[3], pcoords[3], vesselDirection[3];
				int subId;
				vtkIdType cellID;
				GenericCellPointerType intersectCell = GenericCellPointerType::New();
				for(int ii = 0; ii < 3; ++ii)
				{
					a0[ii] = point1[ii];
					a1[ii] = point2[ii];
				}
				if(a0[2] < a1[2])
				{
	// 				double direction[3];
					double dNorm = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						vesselDirection[jj] = a0[jj] - a1[jj];
						dNorm += vesselDirection[jj]*vesselDirection[jj];
					}
					dNorm = sqrt(dNorm);
					for(int jj = 0; jj < 3; ++jj)
					{
						vesselDirection[jj] /= dNorm;
						a0[jj] += vesselDirection[jj]*1000;
					}
				}
				else
				{
	// 				double direction[3];
					double dNorm = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						vesselDirection[jj] = a1[jj] - a0[jj];
						dNorm += vesselDirection[jj]*vesselDirection[jj];
					}
					dNorm = sqrt(dNorm);
					for(int jj = 0; jj < 3; ++jj)
					{
						vesselDirection[jj] /= dNorm;
						a1[jj] += vesselDirection[jj]*1000;
					}
				}
	// 			std::cout << "Vessel = [" << a0[0] << "," << a0[1] << "," << a0[2] << "] - [" << a1[0] << "," << a1[1] << "," << a1[2] << "]" << std::endl;
				int intersection = locator->IntersectWithLine(a0, a1, tol, t, x, pcoords, subId, cellID, intersectCell);
	// 			std::cout << "intersection = " << intersection << std::endl;
				if(intersection)
				{
	// 				std::cout << "Intersection found!" << std::endl;
					PointsPointerType cellPoints = intersectCell->GetPoints();
					double p1[3], p2[3], p3[3], normal[3];
					cellPoints->GetPoint(0, p1), cellPoints->GetPoint(1, p2), cellPoints->GetPoint(2, p3);
					for(int ii = 0; ii < 3; ++ii)
					{
						p3[ii] -= p1[ii];
						p2[ii] -= p1[ii];
					}
					normal[0] = p2[1]*p3[2] - p2[2]*p3[1];
					normal[1] = p2[2]*p3[0] - p2[0]*p3[2];
					normal[2] = p2[0]*p3[1] - p2[1]*p3[0];
					
					double angle = 0, normA = 0, normN = 0;
					for(int ii = 0; ii < 3; ++ii)
					{
						angle += (a0[ii] - a1[ii])*normal[ii];
						normA += (a0[ii] - a1[ii])*(a0[ii] - a1[ii]);
						normN += normal[ii]*normal[ii];
					}
					normA = sqrt(normA);
					normN = sqrt(normN);
					if(normN)
						normal[0] /= normN, normal[1] /= normN, normal[2] /= normN;
					angle = acos(std::abs(angle)/(normA*normN))*180.0/PI;
					double radialDirection[3];
					radialDirection[0] = x[0] - piaCenter[0], radialDirection[1] = x[1] - piaCenter[1], radialDirection[2] = 0;
					
					++histogram[(unsigned int)(angle/binSize)];
	// 				std::cout << "Surface normal @ intersection = [" << normal[0]/normN << "," << normal[1]/normN << "," << normal[2]/normN << "]" << std::endl;
					
					VesselData << x[0] << "," << x[1] << "," << x[2] << "\t";
					VesselData << vesselDirection[0] << "," << vesselDirection[1] << "," << vesselDirection[2] << "\t";
					VesselData << normal[0] << "," << normal[1] << "," << normal[2] << "\t";
					VesselData << radialDirection[0] << "," << radialDirection[1] << "," << radialDirection[2] << std::endl;
				}
			}
		}
		VesselData.close();
		
		std::string format = oFName;
		format += "_vessel_histo.csv";
		
		std::ofstream VesselHistoData( format.c_str() );
		for(int ii = 0; ii < noOfBins; ++ii)
			VesselHistoData << std::fixed << ii*binSize << "\t" << histogram[ii] << std::endl;
		VesselHistoData.close();
		
		delete graphReader;
		delete surfaceReader;
	}
	else
		std::cout << "Error! Could not find landmark with label Pia in SpatialGraph " << graphFName << std::endl;
};

/****************************************************************************/
/*computes angles between the reconstructed blood vessels                   */
/*and the surface normals at the intersection points between pia and vessel.*/
/*only returns vessels that have an angle to the surface normal of less     */
/*than 10 degrees                                                           */
/****************************************************************************/
std::list< unsigned int > Geometry::computeConstrainingVessels(PolyDataPointerType piaSurface)
{
	bool constrained = 1;
	computeBloodVessels(constrained);
	std::list< unsigned int > goodVessels;
// 	std::list< int > vertexRemoveList;
	CellLocatorPointerType locator = CellLocatorPointerType::New();
	locator->AutomaticOn();
	locator->SetDataSet(piaSurface);
	locator->BuildLocator();
// 	locator->Print(std::cout);
	if(spatialGraph->edgesPointer()->size())
	{
		for(unsigned int ii = 0; ii < spatialGraph->edgesPointer()->size();/* ++ii*/)
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
				// 			std::cout << "Vessel = [" << a0[0] << "," << a0[1] << "," << a0[2] << "] - [" << a1[0] << "," << a1[1] << "," << a1[2] << "]" << std::endl;
				int intersection = locator->IntersectWithLine(a0, a1, tol, t, x, pcoords, subId, cellID, intersectCell);
				// 			std::cout << "intersection = " << intersection << std::endl;
				if(intersection)
				{
					// 				std::cout << "Intersection found!" << std::endl;
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
						spatialGraph2->addEdge((*(spatialGraph->edgesPointer()))[ii]);
// 						(*(spatialGraph->edgesPointer()))[ii]->label = Axon;
// 						(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[ii]->edgeConnectivity[0]]->label = Axon;
// 						(*(spatialGraph->verticesPointer()))[(*(spatialGraph->edgesPointer()))[ii]->edgeConnectivity[1]]->label = Axon;
					}
	// 				std::cout << "Surface normal @ intersection = [" << normal[0]/normN << "," << normal[1]/normN << "," << normal[2]/normN << "]" << std::endl;
				}
			}
			++ii;
		}
	}
	return goodVessels;
};

/******************************************************************************/
/*pipeline for extracting contours with a certain label, converting them into */
/*the VTK image format (regular spacing) and running a marching cubes contour */
/*filter. The MC filter output is then passed through a low-pass filter for   */
/*smoothing of the surface. Optinally, a (heuristic) top can be added         */
/******************************************************************************/
void Geometry::computeSurfaces(const char * filename, int interpolParam)
{
// 	ImageDataPointerType piaTop = addTop(Pia, 3);
// 	ImageDataPointerType pia3D = piaVolume(Pia, filename);
// 	ImageDataPointerType completePia = mergeStructures(piaTop, pia3D);
	ImageDataPointerType completePia = addTop2(Pia, 20);
// 	ImageDataPointerType wmTop = addTop(WhiteMatter, 3);
// 	ImageDataPointerType wm3D = piaVolume(WhiteMatter, filename);
// 	ImageDataPointerType completeWM = mergeStructures(wmTop, wm3D);
// 	ImageDataPointerType completeWM = addTop2(WhiteMatter, 10);
	MarchingCubesPointerType mcSurfaceFilter = MarchingCubesPointerType::New();
	mcSurfaceFilter->SetInput(completePia);
// 	mcSurfaceFilter->SetInput(completePia);
// 	mcSurfaceFilter->SetInput(completeWM);
// 	mcSurfaceFilter->SetInput(piaTop);
	std::flush(std::cout << "runnning vtkMarchingCubes..." << std::endl);
	mcSurfaceFilter->SetValue(0, 0);
	mcSurfaceFilter->ComputeScalarsOff();
	mcSurfaceFilter->ComputeGradientsOff();
	mcSurfaceFilter->ComputeNormalsOff();
	mcSurfaceFilter->Update();
// 	mcSurfaceFilter->GetOutput()->Print(std::cout);
	std::string output(filename);
	output += "_mc_vtk";
	std::flush(std::cout << "writing smoothed vtkMarchingCubes output to ascii file..." << std::endl);
	Reader * amiraSurfaceWriter = new Reader(output.c_str(), output.c_str());
	amiraSurfaceWriter->writeAmiraSurfaceFile(smoothSurface(mcSurfaceFilter->GetOutput()));
// 	amiraSurfaceWriter->writeAmiraSurfaceFile(mcSurfaceFilter->GetOutput());
	delete amiraSurfaceWriter;
};

void Geometry::computeBarrelSurfaces(const char * filename, double alpha)
{
// 	spatialGraph->applyTransformation();
	ImageDataPointerType completePia = addTop2(Pia, 20);
	ImageDataPointerType completeWM = addTop2(WhiteMatter, 10);
	MarchingCubesPointerType mcSurfaceFilter1 = MarchingCubesPointerType::New();
	MarchingCubesPointerType mcSurfaceFilter2 = MarchingCubesPointerType::New();
	mcSurfaceFilter1->SetInput(completePia);
	mcSurfaceFilter1->SetValue(0, 0);
	mcSurfaceFilter1->ComputeScalarsOff();
	mcSurfaceFilter1->ComputeGradientsOff();
	mcSurfaceFilter1->ComputeNormalsOff();
	mcSurfaceFilter1->Update();
	mcSurfaceFilter2->SetInput(completeWM);
	mcSurfaceFilter2->SetValue(0, 0);
	mcSurfaceFilter2->ComputeScalarsOff();
	mcSurfaceFilter2->ComputeGradientsOff();
	mcSurfaceFilter2->ComputeNormalsOff();
	mcSurfaceFilter2->Update();
	PolyDataPointerType piaSmooth = smoothSurface(mcSurfaceFilter1->GetOutput());
	PolyDataPointerType wmSmooth = smoothSurface(mcSurfaceFilter2->GetOutput());
	
	spatialGraph2 = new AmiraSpatialGraph;
	std::map< int, PolyDataPointerType > barrels;
	std::map< int, PolyDataPointerType > avgBarrels;
	std::map< int, Column * > columns;
	std::map< int, double * > barrelAxes;
	std::map< int, double * > barrelCenters;
	std::list< unsigned int > goodVessels = computeConstrainingVessels(piaSmooth);
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			std::cout << "Computing parameters of Barrel " << int2Labels[*labelIt] << std::endl;
			PolyDataPointerType barrel = smoothBarrelInZ2(*labelIt);
			double * newAxis = newBarrelAxis(barrel, goodVessels, piaSmooth, filename, alpha);
			double * barrelCenter = calculateBarrelCentroid(barrel);
			barrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, barrel));
			barrelCenters.insert(std::pair< int, double * >(*labelIt, barrelCenter));
			barrelAxes.insert(std::pair< int, double * >(*labelIt, newAxis));
		}
		else
			std::cout << "Barrel " << int2Labels[*labelIt] << " not found in SpatialGraph!" << std::endl;
	}
	// apply divergence constraint on barrel axis vector field
	enforceAxisDivergence(barrelAxes, barrelCenters);
	
	// calculate avg barrel contours and ensure mutual non-overlap
	// of those in barrel field (they can overlap in deeper layers...)
	std::map< int, std::vector< double * > > endPointMap;
	calculateAvgContours(barrels, barrelAxes, barrelCenters, avgBarrels, endPointMap);
	
	std::string outFile(filename);
	outFile += "_parameter.csv";
	std::ofstream ParameterFile(outFile.c_str());
	ParameterFile << "Barrel\tbottom\ttop\theight\tWM\tarea" << std::endl;
	std::map< int, std::vector< double > > barrelParameters;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			PolyDataPointerType barrel = barrels[*labelIt];
			double * finalAxis = barrelAxes[*labelIt];
			double * barrelCenter = barrelCenters[*labelIt];
			std::vector< double * > endPoints = endPointMap[*labelIt];
// 			std::vector< double * > endPoints;
// 			closeBarrelAlongNewAxis(finalAxis, barrelCenter, barrel, endPoints);
			std::vector< double > parameters = computeBarrelParameters(barrel, piaSmooth, wmSmooth, finalAxis, barrelCenter, endPoints, *labelIt, columns, avgBarrels);
			
			if(parameters.size() == 5)
				ParameterFile << int2Labels[*labelIt] << "\t" << parameters[0] << "\t" << parameters[1] << "\t" << parameters[2] << "\t" << parameters[3] << "\t" << parameters[4] << std::endl;
			
			barrelParameters.insert(std::pair< int, std::vector< double > >(*labelIt, parameters));
			spatialGraph2->addPolyDataObject(barrel, *labelIt);
			spatialGraph2->addPolyDataObject(columns[*labelIt]->contours, *labelIt);
		}
	}
	ParameterFile.close();
	
	std::string outFile2(filename);
	outFile2 += "_column_overlap.csv";
	std::ofstream ColumnFile(outFile2.c_str());
	ColumnFile << "column overlap" << std::endl;
	ColumnFile << "center barrel\tneighbor barrel\tdepth";
	for(int ii = 0; ii < 45; ++ii)
		ColumnFile << "\t" << ii*50;
	ColumnFile << std::endl;
	std::string outFile3(filename);
	outFile3 += "_septal_distances.csv";
	std::ofstream SeptaFile(outFile3.c_str());
	SeptaFile << "septal distances" << std::endl;
	SeptaFile << "center barrel\tneighbor barrel\tdepth";
	for(int ii = 0; ii < 45; ++ii)
		SeptaFile << "\t" << ii*50;
	SeptaFile << std::endl;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			std::cout << "Computing overlap and septa of column " << int2Labels[*labelIt] << std::endl;
			std::map< int, std::vector< double > > overlapProfiles = calculateColumnOverlap(columns, *labelIt);
			std::map< int, std::vector< double > >::const_iterator neighborColIt;
			for(neighborColIt = overlapProfiles.begin(); neighborColIt != overlapProfiles.end(); ++neighborColIt)
			{
				ColumnFile << int2Labels[*labelIt] << "\t" << int2Labels[neighborColIt->first] << "\toverlap ratio";
				int ii;
				for(ii = 0; ii < neighborColIt->second.size(); ++ii)
					ColumnFile << "\t" << neighborColIt->second[ii];
				while(ii < 45)
				{
					ColumnFile << "\t" << 0;
					++ii;
				}
				ColumnFile << std::endl;
			}
			
			std::map< int, std::vector< double > > septalProfiles = calculateSeptalDistances(columns, *labelIt);
			for(neighborColIt = septalProfiles.begin(); neighborColIt != septalProfiles.end(); ++neighborColIt)
			{
				SeptaFile << int2Labels[*labelIt] << "\t" << int2Labels[neighborColIt->first] << "\tseptal distance";
				int ii;
				for(ii = 0; ii < neighborColIt->second.size(); ++ii)
					SeptaFile << "\t" << neighborColIt->second[ii];
				while(ii < 45)
				{
					SeptaFile << "\t" << 0;
					++ii;
				}
				SeptaFile << std::endl;
			}
		}
	}
	ColumnFile.close();
	SeptaFile.close();
	
// 	std::vector< Edge * > vesselVec;
// 	std::vector< Edge * >::const_iterator sgVesselIt;
// 	for(sgVesselIt = spatialGraph->edgesBegin(); sgVesselIt != spatialGraph->edgesEnd(); ++sgVesselIt)
// 		if((*sgVesselIt)->label == Vessel && (*sgVesselIt)->edgePointCoordinates.size() == 2)
// 			vesselVec.push_back(*sgVesselIt);
// 	std::map< int, std::list< int > > vesselRadii = computeBarrelVesselCorrelations(avgBarrels, barrelCenters, vesselVec);
// 	std::string outFile3(filename);
// 	outFile3 += "_vessel_correlation.csv";
// 	std::ofstream VesselFile(outFile3.c_str());
// 	VesselFile << "Barrel\theight\tarea\tclosest vessel radii" << std::endl;
// 	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 	{
// 		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
// 		{
// 			VesselFile << int2Labels[*labelIt] << "\t" << barrelParameters[*labelIt][2] << "\t" << barrelParameters[*labelIt][4];
// 			std::list< int >::const_iterator vesselIDListIt;
// 			for(vesselIDListIt = vesselRadii[*labelIt].begin(); vesselIDListIt != vesselRadii[*labelIt].end(); ++vesselIDListIt)
// 				VesselFile << "\t" << *vesselIDListIt;
// 			VesselFile << std::endl;
// 		}
// 	}
// 	VesselFile << std::endl;
// 	VesselFile << "Vessel ID\tradius" << std::endl;
// 	for(int ii = 0; ii < vesselVec.size(); ++ii)
// 		VesselFile << ii << "\t" << vesselVec[ii]->radius << std::endl;
// 	VesselFile.close();
	
	
	std::string piaSurfFile(filename);
	piaSurfFile += "_pia";
	Reader * piaSurfaceWriter = new Reader(piaSurfFile.c_str(), piaSurfFile.c_str());
	piaSurfaceWriter->writeAmiraSurfaceFile(piaSmooth);
	delete piaSurfaceWriter;
	
	std::string wmSurfFile(filename);
	wmSurfFile += "_wm";
	Reader * wmSurfaceWriter = new Reader(wmSurfFile.c_str(), wmSurfFile.c_str());
	wmSurfaceWriter->writeAmiraSurfaceFile(wmSmooth);
	delete wmSurfaceWriter;
};

/******************************************************************************/
/*pipeline for calcuting the total volume of the region populated by the      */
/*columns in the input file & the total volume inside the columns             */
/******************************************************************************/
void Geometry::computeTotalVolumes(const char * filename)
{
	double totalBounds[] = {1E8, -1E8, 1E8, -1E8, 1E8, -1E8};
	std::map< int, Column * > columns;
	AppendPolyDataPointerType mergeFilter = AppendPolyDataPointerType::New();
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		PolyDataPointerType contour = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(*labelIt, contour))
		{
			double * top = new double[3], * bottom = new double[3];
			int subId, subId2;
			double pCoords[3], pCoords2[3], * weights = new double[contour->GetCell(0)->GetNumberOfPoints()], * weights2 = new double[contour->GetCell(1)->GetNumberOfPoints()];
			contour->GetCell(0)->GetParametricCenter(pCoords);
			contour->GetCell(0)->EvaluateLocation(subId, pCoords, top, weights);
			contour->GetCell(1)->GetParametricCenter(pCoords2);
			contour->GetCell(1)->EvaluateLocation(subId2, pCoords2, bottom, weights2);
// 			contour->GetCell(0)->Print(std::cout);
// 			contour->GetCell(1)->Print(std::cout);
			double colBounds[6];
			contour->GetBounds(colBounds);
			for(int ii = 0; ii < 3; ++ii)
			{
				if(colBounds[2*ii] < totalBounds[2*ii])
					totalBounds[2*ii] = colBounds[2*ii];
				if(colBounds[2*ii+1] > totalBounds[2*ii+1])
					totalBounds[2*ii+1] = colBounds[2*ii+1];
			}
			Column * thisCol = new Column(contour, top, bottom);
			columns.insert(std::pair< int, Column * >(*labelIt, thisCol));
			mergeFilter->AddInput(contour);
		}
	}
	mergeFilter->Update();
	PolyDataPointerType allColumns = mergeFilter->GetOutput();
	ConvexHullFilterPointerType convexHullFilter = ConvexHullFilterPointerType::New();
	convexHullFilter->SetInput(allColumns);
// 	convexHullFilter->AddCubeFacePlanes();
	convexHullFilter->AddRecursiveSpherePlanes(2);
	convexHullFilter->Update();
	PolyDataPointerType barrelFieldHull = convexHullFilter->GetOutput();
	PolyDataNormalsPointerType normalsFilter = PolyDataNormalsPointerType::New();
	normalsFilter->SetInput(barrelFieldHull);
	normalsFilter->ConsistencyOn();
	normalsFilter->ComputePointNormalsOff();
	normalsFilter->ComputeCellNormalsOn();
	normalsFilter->SplittingOff();
	normalsFilter->Update();
	barrelFieldHull = normalsFilter->GetOutput();
	spatialGraph->addPolyDataObject(barrelFieldHull, Barrel);
	
	std::flush(std::cout << "Creating voxel grid..." << std::endl);
	ImageDataPointerType totalVolume = createImageVolume(Barrel, totalBounds[0], totalBounds[1], totalBounds[2], totalBounds[3], totalBounds[4], totalBounds[5]);
// 	totalVolume->Print(std::cout);
	int * extent = totalVolume->GetExtent();
	std::flush(std::cout << "Marking voxels inside columns..." << std::endl);
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(spatialGraph->isLabelInSpatialGraph(ID))
		{
			std::flush(std::cout << "Checking for voxels in column " << int2Labels[ID] << std::endl);
			Column * thisCol = columns[ID];
			double minDist = 1E06, maxDist = 0;
			std::vector< double * > contourPts;
			for(int ii = 0; ii < thisCol->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
			{
				double pt[3];
				double dist = 0, t;
				double * closestPt = new double[3];
				thisCol->contours->GetCell(0)->GetPoints()->GetPoint(ii, pt);
				dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, closestPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				for(int jj = 0; jj < 3; ++jj)
					closestPt[jj] = pt[jj] - closestPt[jj];
				contourPts.push_back(closestPt);
			}
	// 		std::flush(std::cout << "minDist = " << minDist << " - maxDist = " << maxDist << std::endl);
			for(int z = extent[4]; z <= extent[5]; ++z)
				for(int y = extent[2]; y <= extent[3]; ++y)
					for(int x = extent[0]; x <= extent[1]; ++x)
					{
						unsigned char * px = static_cast< unsigned char * >(totalVolume->GetScalarPointer(x, y, z));
						if(*px)
							continue;	// no need to waste computing time if this px is already inside a column
						double pt[3];
						pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
						double t, projectedPt[3];
						double dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, projectedPt);
						dist = sqrt(dist);
						if(dist > maxDist || t < 0 || t > 1)
							continue;
						if(dist < minDist)
						{
							*px = Barrel;
							continue;
						}
						
						PolyDataPointerType polyData = PolyDataPointerType::New();
						polyData->Allocate(1);
						PointsPointerType points = PointsPointerType::New();
						points->SetDataTypeToFloat();
						PolygonPointerType poly = PolygonPointerType::New();
						poly->GetPointIds()->SetNumberOfIds(contourPts.size());
						for(int ii = 0; ii < contourPts.size(); ++ii)
						{
							double tmp[3];
							for(int jj = 0; jj < 3; ++jj)
								tmp[jj] = projectedPt[jj] + contourPts[ii][jj];
							points->InsertNextPoint(tmp);
							poly->GetPointIds()->SetId(ii, ii);
						}
						polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
						polyData->SetPoints(points);
						polyData->Update();
						double closestPoint[3], pCoords[3], dist2;
						int subId;
						double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
						int insidePolygon = polyData->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pCoords, dist2, weights);
						if(insidePolygon == 1)
							*px = Barrel;
					}
			totalVolume->Update();
		}
	}
	std::flush(std::cout << "Counting inside voxels..." << std::endl);
	std::list< double * > colPtList;
	std::list< double * > sepPtList;
	unsigned long insideCol = 0;
	unsigned long insideSep = 0;
	PlanePointerType * hullPlanes = new PlanePointerType[barrelFieldHull->GetNumberOfCells()];
	for(int ii = 0; ii < barrelFieldHull->GetNumberOfCells(); ++ii)
	{
		hullPlanes[ii] = PlanePointerType::New();
		double normal[3], origin[3], pCoords[3], * weights = new double[barrelFieldHull->GetCell(ii)->GetNumberOfPoints()];
		int subId;
		barrelFieldHull->GetCell(ii)->GetParametricCenter(pCoords);
		barrelFieldHull->GetCell(ii)->EvaluateLocation(subId, pCoords, origin, weights);
		vtkPolygon::ComputeNormal(barrelFieldHull->GetCell(ii)->GetPoints(), normal);
		hullPlanes[ii]->SetNormal(normal);
		hullPlanes[ii]->SetOrigin(origin);
	}
	for(int z = extent[4]; z <= extent[5]; ++z)
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				bool inside = 1;
				double pt[3];
				pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
				for(int ii = 0; ii < barrelFieldHull->GetNumberOfCells(); ++ii)
				{
					if(hullPlanes[ii]->EvaluateFunction(pt) > 0)
					{
						inside = 0;
						break;
					}
				}
				if(inside)
				{
					unsigned char * px = static_cast< unsigned char * >(totalVolume->GetScalarPointer(x, y, z));
					if(*px == Barrel)
					{
						++insideCol;
						double * coords = new double[3];
						coords[0] = 10*x, coords[1] = 10*y, coords[2] = 10*z;
						colPtList.push_back(coords);
					}
					else
					{
						++insideSep;
						double * coords = new double[3];
						coords[0] = 10*x, coords[1] = 10*y, coords[2] = 10*z;
						sepPtList.push_back(coords);
					}
				}
			}
	
	std::string outFileName1(filename);
	outFileName1 += "_voxels.csv";
	std::ofstream VoxelSizeFile(outFileName1.c_str());
	VoxelSizeFile << "# Nr of Voxels in columns and septa" << std::endl;
	VoxelSizeFile << "# Voxel size = 10x10x10 um^3" << std::endl;
	VoxelSizeFile << "Columns\t" << insideCol << std::endl;
	VoxelSizeFile << "Septa\t" << insideSep << std::endl;
	
	std::flush(std::cout << "Writing landmark files..." << std::endl);
	std::string landmarkFileName1(filename);
	landmarkFileName1 += "_col_voxels.landmarkAscii";
	std::ofstream LandmarkFile1(landmarkFileName1.c_str());
	LandmarkFile1 << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	LandmarkFile1 << std::endl;
	LandmarkFile1 << "# Volume of S1 bounding box:" << std::endl;
	LandmarkFile1 << "# " << extent[1]-extent[0]+1 << " x " << extent[3]-extent[2]+1 << " x " << extent[5]-extent[4]+1  << std::endl;
	LandmarkFile1 << "# Voxel size = 10x10x10 um^3" << std::endl;
	LandmarkFile1 << std::endl;
	LandmarkFile1 << "define Markers " << colPtList.size() << std::endl;
	LandmarkFile1 << std::endl;
	LandmarkFile1 << "Parameters {" << std::endl;
	LandmarkFile1 << "\tNumSets 1," << std::endl;
	LandmarkFile1 << "\tContentType \"LandmarkSet\"" << std::endl;
	LandmarkFile1 << "}" << std::endl;
	LandmarkFile1 << std::endl;
	LandmarkFile1 << "Markers { float[3] Coordinates } @1" << std::endl;
	LandmarkFile1 << std::endl;
	LandmarkFile1 << "# Data section follows" << std::endl;
	LandmarkFile1 << "@1" << std::endl;
	std::list< double * >::const_iterator ptListIt;
	for(ptListIt = colPtList.begin(); ptListIt != colPtList.end(); ++ptListIt)
	{
		LandmarkFile1 << (*ptListIt)[0] << " " << (*ptListIt)[1] << " " << (*ptListIt)[2] << std::endl;
		delete *ptListIt;
	}
	colPtList.clear();
	LandmarkFile1.close();
	
	std::string landmarkFileName2(filename);
	landmarkFileName2 += "_sep_voxels.landmarkAscii";
	std::ofstream LandmarkFile2(landmarkFileName2.c_str());
	LandmarkFile2 << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	LandmarkFile2 << std::endl;
	LandmarkFile2 << "# Volume of S1 bounding box:" << std::endl;
	LandmarkFile2 << "# " << extent[1]-extent[0]+1 << " x " << extent[3]-extent[2]+1 << " x " << extent[5]-extent[4]+1  << std::endl;
	LandmarkFile2 << "# Voxel size = 10x10x10 um^3" << std::endl;
	LandmarkFile2 << std::endl;
	LandmarkFile2 << "define Markers " << sepPtList.size() << std::endl;
	LandmarkFile2 << std::endl;
	LandmarkFile2 << "Parameters {" << std::endl;
	LandmarkFile2 << "\tNumSets 1," << std::endl;
	LandmarkFile2 << "\tContentType \"LandmarkSet\"" << std::endl;
	LandmarkFile2 << "}" << std::endl;
	LandmarkFile2 << std::endl;
	LandmarkFile2 << "Markers { float[3] Coordinates } @1" << std::endl;
	LandmarkFile2 << std::endl;
	LandmarkFile2 << "# Data section follows" << std::endl;
	LandmarkFile2 << "@1" << std::endl;
	std::list< double * >::const_iterator sepListIt;
	for(sepListIt = sepPtList.begin(); sepListIt != sepPtList.end(); ++sepListIt)
	{
		LandmarkFile2 << (*sepListIt)[0] << " " << (*sepListIt)[1] << " " << (*sepListIt)[2] << std::endl;
		delete *sepListIt;
	}
	sepPtList.clear();
	LandmarkFile2.close();
};

PolyDataPointerType Geometry::smoothSurface(PolyDataPointerType surface)
{
// 	AveragePolyDataFilterType smoothingFilter = AveragePolyDataFilterType::New();
// 	smoothingFilter->BoundarySmoothingOff();
// 	smoothingFilter->FeatureEdgeSmoothingOn();
// 	smoothingFilter->SetFeatureAngle(90);
// 	smoothingFilter->SetConvergence(0.01);
// 	smoothingFilter->SetNumberOfIterations(20);
// 	smoothingFilter->SetInput(surface);
// 	smoothingFilter->Update();
// 	return smoothingFilter->GetOutput();
	
	LowpassPolyDataFilterType smoothingFilter = LowpassPolyDataFilterType::New();
	smoothingFilter->BoundarySmoothingOff();
	smoothingFilter->FeatureEdgeSmoothingOff();
// 	smoothingFilter->FeatureEdgeSmoothingOn();
// 	smoothingFilter->SetFeatureAngle(90);
	smoothingFilter->NormalizeCoordinatesOn();
	smoothingFilter->SetNumberOfIterations(20);
	smoothingFilter->SetPassBand(0.1);
	smoothingFilter->SetInput(surface);
	smoothingFilter->Update();
	return smoothingFilter->GetOutput();
};

/******************************************************************************/
/*convert pia contours into vtk polygon (one for each plane),                 */
/*set voxel grid with appropriate resolution                                  */
/******************************************************************************/
ImageDataPointerType Geometry::piaVolume(int label, const char * filename)
{
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	if(spatialGraph->extractLandmark(label, planeEdgePointList, zIndexList))
	{
		// 	std::flush(std::cout << "Calculating bounding box for image..." << std::endl);
		int xMin = 1E06, xMax = -1E06, yMin = 1E06, yMax = -1E06, zMin = 1E06, zMax = -1E06;
		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planeEdgePointList.begin(); planeEdgePointListIt != planeEdgePointList.end(); ++planeEdgePointListIt)
		{
			std::list< double * >::iterator planeVertexIt;
			for(planeVertexIt = planeEdgePointListIt->begin(); planeVertexIt != planeEdgePointListIt->end(); ++planeVertexIt)
			{
				double * tmp = *planeVertexIt;
				if(tmp[0] < xMin)
					xMin = lround(tmp[0]);
				if(tmp[0] > xMax)
					xMax = lround(tmp[0]);
				if(tmp[1] < yMin)
					yMin = lround(tmp[1]);
				if(tmp[1] > yMax)
					yMax = lround(tmp[1]);
				if(tmp[2] < zMin)
					zMin = lround(tmp[2]);
				if(tmp[2] > zMax)
					zMax = lround(tmp[2]);
			}
		}
		zMax -= 50;	//open @ bottom for pia & WM
		PolyDataPointerType polyData = createPolyDataFromPlanePointList(planeEdgePointList);
		ImageDataPointerType volume = createImageVolumeFromPolyData(polyData, label, xMin, xMax, yMin, yMax, zMin, zMax);
		ImageDataPointerType distVolume = distanceTransform(volume);
		distVolume->SetSpacing(volume->GetSpacing());
		return distVolume;
	}
	
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		return NULL;
	}
};

ImageDataPointerType Geometry::piaVolume(int label, int additionalSections)
{
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	if(spatialGraph->extractLandmark(label, planeEdgePointList, zIndexList))
	{
		std::flush(std::cout << "Calculating isosurface..." << std::endl);
		int xMin = 1E06, xMax = -1E06, yMin = 1E06, yMax = -1E06, zMin = 1E06, zMax = -1E06;
		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planeEdgePointList.begin(); planeEdgePointListIt != planeEdgePointList.end(); ++planeEdgePointListIt)
		{
			std::list< double * >::iterator planeVertexIt;
			for(planeVertexIt = planeEdgePointListIt->begin(); planeVertexIt != planeEdgePointListIt->end(); ++planeVertexIt)
			{
				double * tmp = *planeVertexIt;
				if(tmp[0] < xMin)
					xMin = lround(tmp[0]);
				if(tmp[0] > xMax)
					xMax = lround(tmp[0]);
				if(tmp[1] < yMin)
					yMin = lround(tmp[1]);
				if(tmp[1] > yMax)
					yMax = lround(tmp[1]);
				if(tmp[2] < zMin)
					zMin = lround(tmp[2]);
				if(tmp[2] > zMax)
					zMax = lround(tmp[2]);
			}
		}
		zMax -= 50;	//open @ bottom for pia & WM
		zMin -= 50*additionalSections;
		PolyDataPointerType polyData = createPolyDataFromPlanePointList(planeEdgePointList);
		ImageDataPointerType volume = createImageVolumeFromPolyData(polyData, label, xMin, xMax, yMin, yMax, zMin, zMax);
		ImageDataPointerType distVolume = distanceTransform(volume);
		distVolume->SetSpacing(volume->GetSpacing());
		return distVolume;
	}
	
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		return NULL;
	}
};

int * Geometry::calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label)
{
	int * extent = new int[6];
	double spacing[3];
	switch(label)
	{
		case Pia:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case WhiteMatter:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case Barrel:
			spacing[0] = spacing[1] = spacing[2] = 10;
			break;

		default:
			spacing[0] = spacing[1] = spacing[2] = 1;
			break;
	}
	
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = ((double)minX - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = ((double)maxX + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = ((double)minY - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = ((double)maxY + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = ((double)minZ - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = ((double)maxZ + spacing[2])/spacing[2]/* + 0.5*/;
	
	return extent;
};

/******************************************************************************/
/*determine contour pixels for each plane, then compute plane-wise            */
/*distance-transform of image. Implemented in ITK.                            */
/******************************************************************************/
ImageDataPointerType Geometry::distanceTransform(ImageDataPointerType volume)
{
// 	std::flush(std::cout << "Importing vtkImage into ITK..." << std::endl);
	VTK2ITKImageExportPointerType v2iExport = VTK2ITKImageExportPointerType::New();
	ITK2VTKImageImportPointerType i2vImport = ITK2VTKImageImportPointerType::New();
	VTK2ITKImageImportType::Pointer v2iImport = VTK2ITKImageImportType::New();
	ITK2VTKCalcImageExportType::Pointer i2vExport = ITK2VTKCalcImageExportType::New();
	
	ImageType::Pointer volumeImage = ImageType::New();
	ImageType::Pointer contourImage = ImageType::New();
	CalcImageType::Pointer distanceMap = CalcImageType::New();
	DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
	
	v2iExport->SetInput(volume);
	v2iImport->SetCallbackUserData(v2iExport->GetCallbackUserData());
	v2iImport->SetDataExtentCallback(v2iExport->GetDataExtentCallback());
	v2iImport->SetNumberOfComponentsCallback(v2iExport->GetNumberOfComponentsCallback());
	v2iImport->SetOriginCallback(v2iExport->GetOriginCallback());
	v2iImport->SetPipelineModifiedCallback(v2iExport->GetPipelineModifiedCallback());
	v2iImport->SetPropagateUpdateExtentCallback(v2iExport->GetPropagateUpdateExtentCallback());
	v2iImport->SetScalarTypeCallback(v2iExport->GetScalarTypeCallback());
	v2iImport->SetSpacingCallback(v2iExport->GetSpacingCallback());
	v2iImport->SetUpdateDataCallback(v2iExport->GetUpdateDataCallback());
	v2iImport->SetUpdateInformationCallback(v2iExport->GetUpdateInformationCallback());
	v2iImport->SetWholeExtentCallback(v2iExport->GetWholeExtentCallback());
	v2iImport->SetBufferPointerCallback(v2iExport->GetBufferPointerCallback());
	
	volumeImage = v2iImport->GetOutput();
	volumeImage->Update();
	contourImage->SetRegions(volumeImage->GetLargestPossibleRegion());
	contourImage->Allocate();
	contourImage->FillBuffer(0);
	distanceMap->SetRegions(volumeImage->GetLargestPossibleRegion());
	distanceMap->Allocate();
	ImageType::SizeType radius;
	radius.Fill(1);
	SegNeighborhoodIteratorType volumeIter(radius, volumeImage, volumeImage->GetLargestPossibleRegion());
	SegNeighborhoodIteratorType contourIter(radius, contourImage, contourImage->GetLargestPossibleRegion());
	NeighborhoodOffsetVectorType offset = CreateLookUpTable();
	int neighborhood4[] = {10, 12, 13, 15};
	for(volumeIter.GoToBegin(), contourIter.GoToBegin(); !volumeIter.IsAtEnd() && !contourIter.IsAtEnd(); ++volumeIter, ++contourIter)
		if(volumeIter.GetCenterPixel())
		{
			contourIter.SetCenterPixel(volumeIter.GetCenterPixel());
			bool border[] = {0, 0, 0, 0};
			for(int ii = 0; ii < 4; ++ii)
				if(!volumeIter.GetPixel(offset[neighborhood4[ii]]))
					contourIter.SetPixel(offset[neighborhood4[ii]], volumeIter.GetCenterPixel());
		}
	contourImage->Update();
	
// 	std::flush(std::cout << "Calculating 2D contours and distance transform..." << std::endl);
	int startZ = contourImage->GetLargestPossibleRegion().GetIndex()[2];
	int deltaZ = contourImage->GetLargestPossibleRegion().GetSize()[2];
	for(int z = startZ; z < (startZ + deltaZ); ++z)
	{
// 		std::cout << "calculating distance transform in plane " << z << std::endl;
		ImageType::Pointer planeImage = ImageType::New();
		CalcImageType::Pointer planeDistanceMap = CalcImageType::New();
		ImageType::RegionType planeRegion;
		ImageType::RegionType planeImageRegion;
		ImageType::IndexType planeIndex;
		ImageType::SizeType planeSize;
		
		planeIndex = contourImage->GetLargestPossibleRegion().GetIndex();
		planeIndex[2] = z;
		planeSize = contourImage->GetLargestPossibleRegion().GetSize();
		planeSize[2] = 1;
		planeRegion.SetIndex(planeIndex);
		planeRegion.SetSize(planeSize);
		planeImageRegion.SetIndex(contourImage->GetLargestPossibleRegion().GetIndex());
		planeImageRegion.SetSize(planeSize);
		planeImage->SetRegions(planeImageRegion);
		planeImage->Allocate();
		
		ConstIteratorType copyIter(contourImage, planeRegion);
		IteratorType2 pasteIter(planeImage, planeImageRegion);
		for(copyIter.GoToBegin(), pasteIter.GoToBegin(); !copyIter.IsAtEnd() && !pasteIter.IsAtEnd(); ++copyIter, ++pasteIter)
			pasteIter.Set(copyIter.Get());
		planeImage->Update();
		
		distanceMapFilter->UseImageSpacingOff();
		distanceMapFilter->InsideIsPositiveOff();
		distanceMapFilter->SetInput(planeImage);
		planeDistanceMap = distanceMapFilter->GetDistanceMap();
		planeDistanceMap->Update();
		
		ConstCalcIteratorType copyIter2(planeDistanceMap, planeImageRegion);
		CalcIteratorType pasteIter2(distanceMap, planeRegion);
		for(copyIter2.GoToBegin(), pasteIter2.GoToBegin(); !copyIter2.IsAtEnd() && !pasteIter2.IsAtEnd(); ++copyIter2, ++pasteIter2)
			pasteIter2.Set(copyIter2.Get());
		distanceMap->Update();
	}
	
// 	std::flush(std::cout << "Importing itkImage into VTK..." << std::endl);
	i2vExport->SetInput(distanceMap);
	i2vImport->SetCallbackUserData(i2vExport->GetCallbackUserData());
	i2vImport->SetDataExtentCallback(i2vExport->GetDataExtentCallback());
	i2vImport->SetNumberOfComponentsCallback(i2vExport->GetNumberOfComponentsCallback());
	i2vImport->SetOriginCallback(i2vExport->GetOriginCallback());
	i2vImport->SetPipelineModifiedCallback(i2vExport->GetPipelineModifiedCallback());
	i2vImport->SetPropagateUpdateExtentCallback(i2vExport->GetPropagateUpdateExtentCallback());
	i2vImport->SetScalarTypeCallback(i2vExport->GetScalarTypeCallback());
	i2vImport->SetSpacingCallback(i2vExport->GetSpacingCallback());
	i2vImport->SetUpdateDataCallback(i2vExport->GetUpdateDataCallback());
	i2vImport->SetUpdateInformationCallback(i2vExport->GetUpdateInformationCallback());
	i2vImport->SetWholeExtentCallback(i2vExport->GetWholeExtentCallback());
	i2vImport->SetBufferPointerCallback(i2vExport->GetBufferPointerCallback());
	i2vImport->Update();
	
	ImageDataPointerType distanceVolume = ImageDataPointerType::New();
	distanceVolume->DeepCopy(i2vImport->GetOutput());
	distanceVolume->Update();
// 	distanceVolume->Print(std::cout);
// 	int * dims = distanceVolume->GetExtent();
// 	unsigned int contourPxCount = 0, insidePxCount = 0, outsidePxCount = 0;
// 	for(int z = dims[4]; z <= dims[5]; ++z)
// 		for(int y = dims[2]; y <= dims[3]; ++y)
// 			for(int x = dims[0]; x <= dims[1]; ++x)
// 			{
// 				double * px = static_cast< double * >(distanceVolume->GetScalarPointer(x, y, z));
// 				if(*px == 0)
// 				{
// 					++contourPxCount;
// // 					double * thisCoord = new double[3];
// // 					thisCoord[0] = x*50, thisCoord[1] = y*50, thisCoord[2] = z*50;
// // 					Vertex newVertex(thisCoord, 2);
// // 					spatialGraph->addVertex(newVertex);
// 				}
// 				if(*px > 0)
// 					++outsidePxCount;
// 				if(*px < 0)
// 					++insidePxCount;
// 			}
// 	std::cout << "#contour pixels = " << contourPxCount << std::endl;
// 	std::cout << "#inside pixels = " << insidePxCount << std::endl;
// 	std::cout << "#outside pixels = " << outsidePxCount << std::endl;
	return distanceVolume;
};

//calculate deviation of automatic barrel detection from manual barrel outlines
//absolute deviation of barrel centers/centroids
void Geometry::reconstructionError(const char * filename)
{
	ImageDataPointerType completePia = addTop2(Pia, 20);
	MarchingCubesPointerType mcSurfaceFilter1 = MarchingCubesPointerType::New();
	mcSurfaceFilter1->SetInput(completePia);
	mcSurfaceFilter1->SetValue(0, 0);
	mcSurfaceFilter1->ComputeScalarsOff();
	mcSurfaceFilter1->ComputeGradientsOff();
	mcSurfaceFilter1->ComputeNormalsOff();
	mcSurfaceFilter1->Update();
	PolyDataPointerType piaSmooth = smoothSurface(mcSurfaceFilter1->GetOutput());
	ImageDataPointerType completeWM;
	PolyDataPointerType wmSmooth;
	if(spatialGraph->isLabelInSpatialGraph(WhiteMatter))
	{
		completeWM = addTop2(WhiteMatter, 10);
		MarchingCubesPointerType mcSurfaceFilter2 = MarchingCubesPointerType::New();
		mcSurfaceFilter2->SetInput(completeWM);
		mcSurfaceFilter2->SetValue(0, 0);
		mcSurfaceFilter2->ComputeScalarsOff();
		mcSurfaceFilter2->ComputeGradientsOff();
		mcSurfaceFilter2->ComputeNormalsOff();
		mcSurfaceFilter2->Update();
		wmSmooth = smoothSurface(mcSurfaceFilter2->GetOutput());
	}
	
	spatialGraph2 = new AmiraSpatialGraph;
	std::map< int, PolyDataPointerType > barrels;
	std::map< int, PolyDataPointerType > avgBarrels;
	std::map< int, double * > barrelAxes;
	std::map< int, double * > barrelCenters;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		PolyDataPointerType barrel = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(*labelIt, barrel))
		{
			std::cout << "Computing parameters of Barrel " << int2Labels[*labelIt] << std::endl;
			double alpha = 0.5;
			double * newAxis = simpleBarrelAxis(barrel, piaSmooth, alpha);
			double * barrelCenter = calculateBarrelCentroid(barrel);
			barrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, barrel));
			barrelCenters.insert(std::pair< int, double * >(*labelIt, barrelCenter));
			barrelAxes.insert(std::pair< int, double * >(*labelIt, newAxis));
		}
		else
			std::cout << "Barrel " << int2Labels[*labelIt] << " not found in SpatialGraph!" << std::endl;
	}
	// apply divergence constraint on barrel axis vector field
// 	enforceAxisDivergence(barrelAxes, barrelCenters);
	
	// calculate avg barrel contours and ensure mutual non-overlap
	// of those in barrel field (they can overlap in deeper layers...)
	std::map< int, std::vector< double * > > endPointMap;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			std::vector< double * > endPoints;
			closeBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPoints);
			endPointMap.insert(std::pair< int, std::vector< double * > >(*labelIt, endPoints));
			computeAverageHomeBarrel(barrels[*labelIt], barrelCenters[*labelIt], barrelAxes[*labelIt]);
// 			avgBarrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, smoothBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPointMap[*labelIt])));
		}
	
	std::string outFile(filename);
	outFile += "_parameter.csv";
	std::ofstream ParameterFile(outFile.c_str());
	ParameterFile << "Barrel\tbottom\ttop\theight\tWM" << std::endl;
	std::map< int, std::vector< double > > barrelParameters;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		ParameterFile << int2Labels[*labelIt];
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			PolyDataPointerType barrel = barrels[*labelIt];
			double * finalAxis = barrelAxes[*labelIt];
			double * barrelCenter = barrelCenters[*labelIt];
			std::vector< double * > endPoints = endPointMap[*labelIt];
// 			std::vector< double * > endPoints;
// 			closeBarrelAlongNewAxis(finalAxis, barrelCenter, barrel, endPoints);
			std::vector< double > parameters = computeManualBarrelParameters(barrel, piaSmooth, wmSmooth, finalAxis, barrelCenter, endPoints);
			for(int ii = 0; ii < parameters.size(); ++ii)
				ParameterFile << "\t" << parameters[ii];
			
			barrelParameters.insert(std::pair< int, std::vector< double > >(*labelIt, parameters));
			spatialGraph2->addPolyDataObject(barrel, *labelIt);
		}
		ParameterFile << std::endl;
	}
	ParameterFile.close();
	
	std::string piaSurfFile(filename);
	piaSurfFile += "_pia";
	Reader * piaSurfaceWriter = new Reader(piaSurfFile.c_str(), piaSurfFile.c_str());
	piaSurfaceWriter->writeAmiraSurfaceFile(piaSmooth);
	delete piaSurfaceWriter;
	
	if(wmSmooth)
	{
		std::string wmSurfFile(filename);
		wmSurfFile += "_wm";
		Reader * wmSurfaceWriter = new Reader(wmSurfFile.c_str(), wmSurfFile.c_str());
		wmSurfaceWriter->writeAmiraSurfaceFile(wmSmooth);
		delete wmSurfaceWriter;
	}
};

NeighborhoodOffsetVectorType Geometry::CreateLookUpTable()
{
	SegNeighborhoodIteratorType::OffsetType offset;
	NeighborhoodOffsetVectorType look_up_table(0);
	
	for(int z = -1; z <= 1 ; z++)
	{
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				if( (x == 0) && (y == 0) && (z == 0) )
				{}
				else
				{
					offset[0]=x;
					offset[1]=y;
					offset[2]=z;
					look_up_table.push_back(offset);
				}
			}
		}
	}
	
	return look_up_table;
};

/******************************************************************************/
/*interpolate barrel/pia/WM top by radial sampling of #interpolationSlices    */
/*slices below top and linear interpolation through sampled points.           */
/******************************************************************************/
ImageDataPointerType Geometry::addTop(int label, int interpolationSlices)
{
	if(interpolationSlices < 2)
	{
		std::cout << "Error! At least 2 slices required for top interpolation!" << std::endl;
		return NULL;
	}
	
	//first, extract pia contours
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	if(spatialGraph->extractLandmark(label, planeEdgePointList, zIndexList))
	{
		double * topBounds = new double[6];
		topBounds[0] = 1E06, topBounds[1] = -1E06, topBounds[2] = 1E06, topBounds[3] = -1E06;
		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planeEdgePointList.begin(); planeEdgePointListIt != planeEdgePointList.end(); ++planeEdgePointListIt)
		{
			std::list< double * >::iterator planeVertexIt;
			for(planeVertexIt = planeEdgePointListIt->begin(); planeVertexIt != planeEdgePointListIt->end(); ++planeVertexIt)
			{
				double * tmp = *planeVertexIt;
				if(tmp[0] < topBounds[0])
					topBounds[0] = round(tmp[0]);
				if(tmp[0] > topBounds[1])
					topBounds[1] = round(tmp[0]);
				if(tmp[1] < topBounds[2])
					topBounds[2] = round(tmp[1]);
				if(tmp[1] > topBounds[3])
					topBounds[3] = round(tmp[1]);
			}
		}
		
// 		zIndexList.sort();	//already returned sorted
		int * zInterpolPlanes = new int[interpolationSlices];
		std::list< int >::iterator zIndexListIter = zIndexList.begin();
		for(int ii = 0; ii < interpolationSlices && zIndexListIter != zIndexList.end(); ++ii, ++zIndexListIter)
		{
			zInterpolPlanes[ii] = *zIndexListIter;
// 			std::cout << "Plane " << ii << " @ z = " << zInterpolPlanes[ii] << std::endl;
		}
		
		// only work with top slices from now on
		PolyDataPointerType polyData = PolyDataPointerType::New();
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
// 		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planeEdgePointList.begin(); planeEdgePointListIt != planeEdgePointList.end(); ++planeEdgePointListIt)
		{
			bool interpolPlane = 0;
			for(int ii = 0; ii < interpolationSlices; ++ii)
				if(planeEdgePointListIt->back()[2] == zInterpolPlanes[ii])
				{
					interpolPlane = 1;
					break;
				}
			if(!interpolPlane)
				continue;
// 			std::flush(std::cout << "Allocating memory for " << planeEdgePointListIt->size() - 1 << " points in plane " << planeEdgePointListIt->back()[2] << "..." << std::endl);
			int end = planeEdgePointListIt->size() - 1;	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			
			std::list< double * >::iterator planeVertexIt = planeEdgePointListIt->begin();
			for(int ii = 0; ii < end; ++planeVertexIt, ++ii)
			{
				double * tmp = *planeVertexIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
		}
		polyData->SetPoints(points);
		polyData->Update();
// 		polyData->Print(std::cout);

		// radial sampling of top contours
		std::cout << "starting radial sampling..." << std::endl;
		double * centerPoint = new double[3];	// calculate centerPoint as center of top plane
		for(int ii = 0; ii < polyData->GetNumberOfCells(); ++ii)
		{
			double * bounds = polyData->GetCell(ii)->GetBounds();
// 			std::cout << "bounds[4] = " << bounds[4] << " --- zInterpolPlanes[0] = " << zInterpolPlanes[0] << std::endl;
			if(bounds[4] == zInterpolPlanes[0])
			{
				int subID;
				double pCoords[3], * weights;
				weights = new double[polyData->GetCell(ii)->GetNumberOfPoints()];
				polyData->GetCell(ii)->GetParametricCenter(pCoords);
				polyData->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
// 				std::cout << "centerPoint @ [" << centerPoint[0] << "," << centerPoint[1] << "," << centerPoint[2] << "]" << std::endl;
			}
		}
// 		double diag = polyData->GetLength();
		int noAngles = 20;
		std::vector< PointsPointerType > sampledPoints;
		std::vector< double * > sampledDerivatives;
		for(int ii = 0; ii < noAngles/2; ++ii)
		{
			PointsPointerType thisAnglePoints = PointsPointerType::New();
			thisAnglePoints->SetDataTypeToFloat();
			int pointID = 0;
			std::map< int, int > planeIDs;
			for(int jj = 0; jj < interpolationSlices; ++jj)
				planeIDs.insert(std::pair< int, int >(zInterpolPlanes[jj], jj));
// 			std::cout << "Sampling ray @ angle " << ii*360/noAngles << std::endl;
			centerPoint[2] = zInterpolPlanes[0];
			
			double angle = ii*2*PI/noAngles;
			double direction[3];
			direction[0] = sin(angle);	// plane normal!!!
			direction[1] = cos(angle);
			direction[2] = 0;
			
			PlanePointerType thisPlane = PlanePointerType::New();
			CutterPointerType planeCutter = CutterPointerType::New();
			thisPlane->SetOrigin(centerPoint);
			thisPlane->SetNormal(direction);
			planeCutter->SetCutFunction(thisPlane);
			planeCutter->SetInput(polyData);
			planeCutter->Update();
			PolyDataPointerType cutPoints = planeCutter->GetOutput();
			thisAnglePoints = cutPoints->GetPoints();
			
			//get rid of extra points created by vtkCutter:
			//calculate pairwise distance between all points; the 2 points with max pairwise distance
			//are the 2 outermost points that we need
			PointsPointerType extremePoints = PointsPointerType::New();
			extremePoints->SetNumberOfPoints(2*interpolationSlices);
			std::map< double, std::pair< int, int > > * distanceIDs = new std::map< double, std::pair< int, int > >[interpolationSlices];	//map automatically sorts pairs by their labels (==distance)
			for(int kk = 0; kk < thisAnglePoints->GetNumberOfPoints(); ++kk)
			{
				double * tmpCoords1 = new double[3];
				thisAnglePoints->GetPoint(kk, tmpCoords1);
				int currPlane = lround(tmpCoords1[2]);
				for(int ll = kk; ll < thisAnglePoints->GetNumberOfPoints(); ++ll)
				{
					double * tmpCoords2 = new double[3];
					thisAnglePoints->GetPoint(ll, tmpCoords2);
					if(lround(tmpCoords2[2]) != currPlane)
						continue;
					double distance = sqrt((tmpCoords1[0]-tmpCoords2[0])*(tmpCoords1[0]-tmpCoords2[0]) + (tmpCoords1[1]-tmpCoords2[1])*(tmpCoords1[1]-tmpCoords2[1]));
					distanceIDs[planeIDs[currPlane]].insert(std::pair< double, std::pair< int, int> >(distance, std::pair< int, int >(kk, ll)));
					delete [] tmpCoords2;
				}
				delete [] tmpCoords1;
			}
			for(int kk = 0; kk < interpolationSlices; ++kk)
			{
				double * tmpCoords1 = new double[3];
				double * tmpCoords2 = new double[3];
				int tmpID1 = distanceIDs[kk].rbegin()->second.first;	//maps sort their keys in ascending order
				int tmpID2 = distanceIDs[kk].rbegin()->second.second;
				thisAnglePoints->GetPoint(tmpID1, tmpCoords1);
				thisAnglePoints->GetPoint(tmpID2, tmpCoords2);
				extremePoints->InsertPoint(2*kk, tmpCoords1);
				extremePoints->InsertPoint(2*kk+1, tmpCoords2);
// 				Vertex newPoint1(tmpCoords1, Pia);
// 				Vertex newPoint2(tmpCoords2, Pia);
// 				spatialGraph->addVertex(newPoint1);
// 				spatialGraph->addVertex(newPoint2);
// 				std::list< double * > edgePts;
// 				edgePts.push_back(tmpCoords1);
// 				edgePts.push_back(tmpCoords2);
// 				int connectivity[2];
// 				if(!spatialGraph->getNumberOfVertices())
// 				{
// 					connectivity[0] = 0;
// 					connectivity[1] = 1;
// 				}
// 				else
// 				{
// 					connectivity[0] = spatialGraph->getNumberOfVertices() - 2;
// 					connectivity[1] = spatialGraph->getNumberOfVertices() - 1;
// 				}
// 				Edge pt1(connectivity, edgePts.size(), Pia, edgePts);
// 				spatialGraph->addEdge(pt1);
			}
			sampledPoints.push_back(extremePoints);
// 			extremePoints->Print(std::cout);
		}
		
		//calculate finite differences for each ray
		for(int ii = 0; ii < sampledPoints.size(); ++ii)
		{
// 			sampledPoints[ii]->Print(std::cout);
			//sort sampled points -> important for correct finite difference calculation!!!
			// IDs 0,1,2: first 3 points out->in
			// IDs 3,4,5: 3 points out->in on other end of ray
			int maxID = 0;
			double maxDist = 0;
			double * maxCoords;
			std::map< double, int > distSortedIDs;
			for(int jj = 0; jj < sampledPoints[ii]->GetNumberOfPoints(); ++jj)
			{
				double * tmpCoords = new double[3];
				sampledPoints[ii]->GetPoint(jj, tmpCoords);
				double currDist = sqrt((tmpCoords[0] - centerPoint[0])*(tmpCoords[0] - centerPoint[0]) + (tmpCoords[1] - centerPoint[1])*(tmpCoords[1] - centerPoint[1]));
				if(currDist > maxDist)
				{
					maxID = jj;
					maxCoords = tmpCoords;
				}
				else
					delete [] tmpCoords;
			}
			for(int jj = 0; jj < sampledPoints[ii]->GetNumberOfPoints(); ++jj)
			{
				double * tmpCoords = new double[3];
				sampledPoints[ii]->GetPoint(jj, tmpCoords);
				double currDist = sqrt((tmpCoords[0] - maxCoords[0])*(tmpCoords[0] - maxCoords[0]) + (tmpCoords[1] - maxCoords[1])*(tmpCoords[1] - maxCoords[1]));
				distSortedIDs.insert(std::pair< double, int >(currDist, jj));	// maxID is going to be @ first position b/c currDist == 0 then
				delete [] tmpCoords;
			}
			
			std::map< double, int >::iterator distSortedIDsIter;
			int sortedIDs[] = {0, 1, 2, 5, 4, 3};
			int kk = 0;
			PointsPointerType orderedPoints = PointsPointerType::New();
			orderedPoints->SetNumberOfPoints(6);
			orderedPoints->SetDataTypeToDouble();
			for(distSortedIDsIter = distSortedIDs.begin(); distSortedIDsIter != distSortedIDs.end(); ++distSortedIDsIter, ++kk)
			{
				// IDs 0,1,2: first 3 points out->in
				// IDs 3,4,5: 3 points out->in on other end of ray
				double * tmpCoords = new double[3];
				sampledPoints[ii]->GetPoint(distSortedIDsIter->second, tmpCoords);
// 				std::cout << "[" << tmpCoords[0] << "," << tmpCoords[1] << "," << tmpCoords[2] << "]" << std::endl;
				orderedPoints->InsertPoint(sortedIDs[kk], tmpCoords);
			}
			sampledPoints[ii] = orderedPoints;	//order needed again later for calculateTop()
			
			sampledDerivatives.push_back(samplingRayDerivatives(orderedPoints));
		}
		
		PolyDataPointerType top = calculateTop(sampledDerivatives, sampledPoints);
		//write top contours to spatial graph file
// 		delete spatialGraph;
// 		spatialGraph = new AmiraSpatialGraph();
// 		PointsPointerType tmpPoints = top->GetPoints();
// 		double * begin = new double[3];
// 		double * end = new double[3];
// 		tmpPoints->GetPoint(0, begin);
// 		tmpPoints->GetPoint(0, end);
// 		Vertex * newVert1 = new Vertex(begin, label);
// 		Vertex * newVert2 = new Vertex(end, label);
// 		spatialGraph->addVertex(newVert1);
// 		spatialGraph->addVertex(newVert2);
// 		std::list< double * > edgePtList;
// 		for(int ii = 0; ii < tmpPoints->GetNumberOfPoints(); ++ii)
// 		{
// 			double * tmpPoint = new double[3];
// 			tmpPoints->GetPoint(ii, tmpPoint);
// 			edgePtList.push_back(tmpPoint);
// 		}
// 		edgePtList.push_back(end);
// 		int connectivity[2];
// 		if(!spatialGraph->getNumberOfVertices())
// 		{
// 			connectivity[0] = 0;
// 			connectivity[1] = 1;
// 		}
// 		else
// 		{
// 			connectivity[0] = spatialGraph->getNumberOfVertices() - 2;
// 			connectivity[1] = spatialGraph->getNumberOfVertices() - 1;
// 		}
// 		Edge * newEdge = new Edge(connectivity, edgePtList.size(), label, edgePtList);
// 		spatialGraph->addEdge(newEdge);
		
		topBounds[4] = top->GetBounds()[4];
		topBounds[5] = top->GetBounds()[5];
		int xMin = lround(topBounds[0]), xMax = lround(topBounds[1]);
		int yMin = lround(topBounds[2]), yMax = lround(topBounds[3]);
		int zMin = lround(topBounds[4]), zMax = lround(topBounds[5]);
		ImageDataPointerType volume = createImageVolumeFromPolyData(top, label, xMin, xMax, yMin, yMax, zMin, zMax);
		ImageDataPointerType distVolume = distanceTransform(volume);
		distVolume->SetSpacing(volume->GetSpacing());
		return distVolume;
	}
	
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		return NULL;
	}
};

/******************************************************************************/
/*interpolate barrel/pia/WM top by following the gradient of the distance     */
/*transform of the existing landmark                                          */
/******************************************************************************/
ImageDataPointerType Geometry::addTop2(int label, int additionalSections)
{
	ImageDataPointerType pia = piaVolume(label, additionalSections);
	int extent[6] = {0, 0, 0, 0, 0, 0};	//minX, maxX, minY, maxY, minZ, maxZ
	pia->GetExtent(extent);
	int startZ = extent[4] + additionalSections;	//first section w/o data
	if(startZ > extent[5] - 2)
	{
		std::cout << "Error! Pia volume incorrect: vertical extent smaller than data extent" << std::endl;
		return pia;
	}
	int zCount = 0;
	for(int z = startZ; z >= extent[4]; --z)
	{
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				float * dist0 = static_cast< float * >(pia->GetScalarPointer(x, y, z));
				float * dist1 = static_cast< float * >(pia->GetScalarPointer(x, y, z+1));
				float * dist2 = static_cast< float * >(pia->GetScalarPointer(x, y, z+2));
				*dist0 = *dist1 + (*dist1 - *dist2) + zCount*zCount;	// effectively calculates gradient of the distance image plus an artificial curvature
			}
			++zCount;
	}
	
	return pia;
};

//calculate centered differences from 3 given points
//ias an approximation to slope and curvature at top of pia
//point order:
// IDs 0,1,2: first 3 points out->in
// IDs 3,4,5: 3 points out->in on other end of ray
double * Geometry::samplingRayDerivatives(PointsPointerType rayPoints)
{
	if(rayPoints->GetNumberOfPoints() != 6)
	{
		std::cout << "Error! 6 points per sampling ray required to compute derivatives!" << std::endl;
		return NULL;
	}
	
	double * derivs = new double[4];
	
	double out1[3], center1[3], in1[3], out2[3], center2[3], in2[3];
	rayPoints->GetPoint(0, out1);
	rayPoints->GetPoint(1, center1);
	rayPoints->GetPoint(2, in1);
	rayPoints->GetPoint(3, out2);
	rayPoints->GetPoint(4, center2);
	rayPoints->GetPoint(5, in2);
	
	double h1in = sqrt((in1[0]-center1[0])*(in1[0]-center1[0]) + (in1[1]-center1[1])*(in1[1]-center1[1]));
	double h1out = sqrt((out1[0]-center1[0])*(out1[0]-center1[0]) + (out1[1]-center1[1])*(out1[1]-center1[1]));
	double h2in = sqrt((in2[0]-center2[0])*(in2[0]-center2[0]) + (in2[1]-center2[1])*(in2[1]-center2[1]));
	double h2out = sqrt((out2[0]-center2[0])*(out2[0]-center2[0]) + (out2[1]-center2[1])*(out2[1]-center2[1]));
	
	derivs[0] = (in1[2] - out1[2])/(h1in + h1out);
	derivs[1] = (in1[2] + out1[2] - 2*(derivs[0]*(h1in + h1out)/2))/(h1in*h1out);
	derivs[2] = (in2[2] - out2[2])/(h2in + h2out);
	derivs[3] = (in2[2] + out2[2] - 2*(derivs[2]*(h1in + h1out)/2))/(h2in*h2out);
// 	std::cout << "in1[2] = " << in1[2] << std::endl;
// 	std::cout << "out1[2] = " << out1[2] << std::endl;
// 	std::cout << "center1[2] = " << center1[2] << std::endl;
// 	std::cout << "in2[2] = " << in2[2] << std::endl;
// 	std::cout << "out2[2] = " << out2[2] << std::endl;
// 	std::cout << "center2[2] = " << center2[2] << std::endl;
// 	std::cout << "derivs[0] = " << derivs[0] << std::endl;
// 	std::cout << "derivs[1] = " << derivs[1] << std::endl;
// 	std::cout << "derivs[2] = " << derivs[2] << std::endl;
// 	std::cout << "derivs[3] = " << derivs[3] << std::endl;
	
	return derivs;
};

//point order:
// IDs 0,1,2: first 3 points out->in
// IDs 3,4,5: 3 points out->in on other end of ray
PolyDataPointerType Geometry::calculateTop(std::vector< double * > derivatives, std::vector< PointsPointerType > existingPoints)
{
	double averageCurvature = 0;
	double * averageSlope = new double[derivatives.size()*2];
	for(int ii = 0; ii < derivatives.size(); ++ii)
	{
		averageSlope[2*ii] = derivatives[ii][0];
		averageSlope[2*ii+1] = derivatives[ii][2];
		if(ii > 0)
		{
			averageSlope[2*ii] += derivatives[ii-1][0];
			averageSlope[2*ii+1] += derivatives[ii-1][2];
		}
		else
		{
			averageSlope[2*ii] += derivatives[derivatives.size()-1][0];
			averageSlope[2*ii+1] += derivatives[derivatives.size()-1][2];
		}
		if(ii < derivatives.size()-1)
		{
			averageSlope[2*ii] += derivatives[ii+1][0];
			averageSlope[2*ii+1] += derivatives[ii+1][2];
		}
		else
		{
			averageSlope[2*ii] += derivatives[0][0];
			averageSlope[2*ii+1] += derivatives[0][2];
		}
		averageSlope[2*ii] /= 3;
		averageSlope[2*ii+1] /= 3;
		averageCurvature += derivatives[ii][1];
		averageCurvature += derivatives[ii][3];
	}
	averageCurvature = averageCurvature/(2*derivatives.size());
	
	double zDirection;	// +/- 1; indicates whether top should be in positive or negative z-Direction
	double startZ;
	std::vector< double > samplingRadii;
	std::vector< double > totalInwardStepDist;
	std::vector< double * > samplingDirections;	// [2*existingPointsIndex] is inward direction of ray w/ pointID 2, [2*existingPointsIndex+1] -> pointID 5
	for(int ii = 0; ii < existingPoints.size(); ++ii)
	{
		double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
		existingPoints[ii]->GetPoint(1, tmp3);
		existingPoints[ii]->GetPoint(2, tmp1);
		existingPoints[ii]->GetPoint(4, tmp4);
		existingPoints[ii]->GetPoint(5, tmp2);
		
		if(ii == 0)
		{
			zDirection = (tmp1[2] - tmp3[2]) < 0 ? -1 : 1;
			startZ = tmp1[2];
		}
		
		double diam = sqrt((tmp1[0]-tmp2[0])*(tmp1[0]-tmp2[0]) + (tmp1[1]-tmp2[1])*(tmp1[1]-tmp2[1]));
		samplingRadii.push_back(diam/2);
// 		std::cout << "sampling radius " << ii << " = " << samplingRadii[ii] << std::endl;
		totalInwardStepDist.push_back(0);
		totalInwardStepDist.push_back(0);
		double * direction1 = new double[2];
		double * direction2 = new double[2];
		direction1[0] = tmp1[0] - tmp3[0], direction1[1] = tmp1[1] - tmp3[1];
		direction2[0] = tmp2[0] - tmp4[0], direction2[1] = tmp2[1] - tmp4[1];
		double direction1Norm = sqrt(direction1[0]*direction1[0] + direction1[1]*direction1[1]);
		double direction2Norm = sqrt(direction2[0]*direction2[0] + direction2[1]*direction2[1]);
		direction1[0] /= direction1Norm, direction1[1] /= direction1Norm;
		direction2[0] /= direction2Norm, direction2[1] /= direction2Norm;
		samplingDirections.push_back(direction1);
		samplingDirections.push_back(direction2);
	}
// 	std::cout << "startZ = " << startZ << std::endl;
	PolyDataPointerType top = PolyDataPointerType::New();
	PointsPointerType topPoints = PointsPointerType::New();
	top->Allocate(1);
	topPoints->SetDataTypeToFloat();
	std::list< std::vector< double * > > allPlanePoints;
	int lastID = 0;
	bool validContour = 1;
	int zStep = 1;
	double zSamplingSize = 50;
// 	std::cout << "zSamplingSize = " << zSamplingSize << std::endl;
	while(validContour /*&& zStep < 2*/)
	{
		std::vector< double * > thisPlanePoints;
// 		thisPlanePoints.assign(derivatives.size()*2, NULL);
		for(int ii = 0; ii < derivatives.size(); ++ii)
		{
			double root1 = 0, root2 = 0, stepSize = 0;
			double a = averageCurvature*0.5;
			double b = averageSlope[2*ii];
			
			stepSize = std::abs(zSamplingSize/b);
			totalInwardStepDist[2*ii] += std::abs(stepSize);
// 			std::cout << "stepSize = " << stepSize << std::endl;
// 			std::cout << "totalInwardStepDist[2*ii] = " << totalInwardStepDist[2*ii] << std::endl;
			if(totalInwardStepDist[2*ii] >= /*0.75**/samplingRadii[ii])
			{
				validContour = 0;
				break;
			}
			
			double edgePoint[3];
			existingPoints[ii]->GetPoint(2, edgePoint);
			double * newPoint1 = new double[3];
			newPoint1[0] = edgePoint[0] + totalInwardStepDist[2*ii]*samplingDirections[2*ii][0];
			newPoint1[1] = edgePoint[1] + totalInwardStepDist[2*ii]*samplingDirections[2*ii][1];
			newPoint1[2] = round(startZ + zDirection*zStep*zSamplingSize);
// 			std::cout << "edgePoint = [" << edgePoint[0] << "," << edgePoint[1] << "," << edgePoint[2] << "]" << std::endl;
// 			std::cout << "newPoint1 = [" << newPoint1[0] << "," << newPoint1[1] << "," << newPoint1[2] << "]" << std::endl;
// 			topPoints->InsertPoint(2*ii, newPoint1);
// 			thisPlanePoints[2*ii] = newPoint1;
			thisPlanePoints.push_back(newPoint1);
			
			// other point on this sampling ray
			b = averageSlope[2*ii+1];
			stepSize = std::abs(zSamplingSize/b);
			totalInwardStepDist[2*ii+1] += std::abs(stepSize);
// 			std::cout << "stepSize = " << stepSize << std::endl;
// 			std::cout << "totalInwardStepDist[2*ii+1] = " << totalInwardStepDist[2*ii+1] << std::endl;
			if(totalInwardStepDist[2*ii+1] >= /*0.75**/samplingRadii[ii])
			{
				validContour = 0;
				break;
			}
			
			existingPoints[ii]->GetPoint(5, edgePoint);
			double * newPoint2 = new double[3];
			newPoint2[0] = edgePoint[0] + totalInwardStepDist[2*ii+1]*samplingDirections[2*ii+1][0];
			newPoint2[1] = edgePoint[1] + totalInwardStepDist[2*ii+1]*samplingDirections[2*ii+1][1];
			newPoint2[2] = round(startZ + zDirection*zStep*zSamplingSize);
// 			std::cout << "edgePoint = [" << edgePoint[0] << "," << edgePoint[1] << "," << edgePoint[2] << "]" << std::endl;
// 			std::cout << "newPoint2 = [" << newPoint2[0] << "," << newPoint2[1] << "," << newPoint2[2] << "]" << std::endl;
// 			topPoints->InsertPoint(2*ii+1, newPoint2);
// 			thisPlanePoints[2*ii+1] = newPoint2;
			thisPlanePoints.push_back(newPoint2);
		}
		if(!validContour)
			break;
		
		allPlanePoints.push_back(thisPlanePoints);
		lastID += 2*derivatives.size();
		++zStep;
	}
	unsigned int totalNoPoints = 0;
	for(int ii = 0; ii < allPlanePoints.size(); ++ii)
		totalNoPoints += existingPoints.size()*2;
	topPoints->SetNumberOfPoints(totalNoPoints);
	//we want each polygon to be non-intersecting (topologically, like a circle)
	//look for closest point of the two points on neighboring sampling ray
	//CAUTION: may give wrong results when sampling from very few rays
	//but still better than brute-force solving the traveling salesman...
	std::list< std::vector< double * > >::iterator allPlanePointsIter;
	int currID = 0;
	for(allPlanePointsIter = allPlanePoints.begin(); allPlanePointsIter != allPlanePoints.end(); ++allPlanePointsIter)
	{
		PolygonPointerType thisPlanePoly = PolygonPointerType::New();
		thisPlanePoly->GetPointIds()->SetNumberOfIds(derivatives.size()*2);
		int polyPtID = 0;
// 		std::cout << "no. of points in this plane: " << allPlanePointsIter->size() << std::endl;
		for(int ii = 0; ii < allPlanePointsIter->size()/2; ++ii)
		{
			if(ii == 0)
			{
				topPoints->InsertPoint(currID, (*allPlanePointsIter)[0]);
				thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
				topPoints->InsertPoint(currID + derivatives.size(), (*allPlanePointsIter)[1]);
				thisPlanePoly->GetPointIds()->SetId(polyPtID + derivatives.size(), currID + derivatives.size());
				++currID;
				++polyPtID;
// 				std::cout << "pt1 = [" << (*allPlanePointsIter)[0][0] << "," << (*allPlanePointsIter)[0][1] << "," << (*allPlanePointsIter)[0][2] << "]" << std::endl;
// 				std::cout << "pt2 = [" << (*allPlanePointsIter)[1][0] << "," << (*allPlanePointsIter)[1][1] << "," << (*allPlanePointsIter)[1][2] << "]" << std::endl;
			}
			else
			{
				double * previousPt = new double[3];
				double * pt1/* = new double[3]*/;
				double * pt2/* = new double[3]*/;
				pt1 = (*allPlanePointsIter)[2*ii];
				pt2 = (*allPlanePointsIter)[2*ii+1];
				topPoints->GetPoint(currID - 1, previousPt);
				double dist1 = sqrt((previousPt[0] - pt1[0])*(previousPt[0] - pt1[0]) + (previousPt[1] - pt1[1])*(previousPt[1] - pt1[1]));
				double dist2 = sqrt((previousPt[0] - pt2[0])*(previousPt[0] - pt2[0]) + (previousPt[1] - pt2[1])*(previousPt[1] - pt2[1]));
				if(dist1 < dist2)
				{
					topPoints->InsertPoint(currID, pt1);
					thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
					topPoints->InsertPoint(currID + derivatives.size(), pt2);
					thisPlanePoly->GetPointIds()->SetId(polyPtID + derivatives.size(), currID + derivatives.size());
				}
				else
				{
					topPoints->InsertPoint(currID, pt2);
					thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
					topPoints->InsertPoint(currID + derivatives.size(), pt1);
					thisPlanePoly->GetPointIds()->SetId(polyPtID + derivatives.size(), currID + derivatives.size());
				}
				++currID;
				++polyPtID;
// 				std::cout << "pt1 = [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
// 				std::cout << "pt2 = [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
			}
		}
		currID += derivatives.size();
// 		thisPlanePoly->Print(std::cout);
		top->InsertNextCell(thisPlanePoly->GetCellType(), thisPlanePoly->GetPointIds());
	}
// 	topPoints->Print(std::cout);
	top->SetPoints(topPoints);
// 	top->Print(std::cout);
	return top;
};

// smooth Barrel contours in z-Direction by averaging them
// with all points in +- 25(?) um in z
PolyDataPointerType Geometry::smoothBarrelInZ()
{
	int label = Barrel;
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	PolyDataPointerType barrel = PolyDataPointerType::New();
	if(spatialGraph->extractLandmark(label, barrel, zIndexList))
	{
		std::vector< std::vector< double * > > allPlaneAvgPoints;
		double * barrelTopBounds = barrel->GetCell(0)->GetBounds();
		int zDirection = 0;	//+/-1 depending on whether zIndexList parallel/antiparallel to orientation of barrel
		std::vector< int > zIndexVector;
		std::cout << "barrel bounds z = " << barrelTopBounds[4] << " --- zIndexList.front() = " << zIndexList.front() << " --- zIndexList.back() = " << zIndexList.back() << std::endl;
		if(lround(barrelTopBounds[4]) == zIndexList.front())
		{
			zDirection = 1;
			std::list< int >::iterator zListIt;
			for(zListIt = zIndexList.begin(); zListIt != zIndexList.end(); ++zListIt)
				zIndexVector.push_back(*zListIt);
			std::cout << "Barrel and zIndexList parallel!" << std::endl;
		}
		else if(lround(barrelTopBounds[4]) == zIndexList.back())
		{
			zDirection = -1;
			std::list< int >::reverse_iterator zListIt;
			for(zListIt = zIndexList.rbegin(); zListIt != zIndexList.rend(); ++zListIt)
				zIndexVector.push_back(*zListIt);
			std::cout << "Barrel and zIndexList antiparallel!" << std::endl;
		}
		else
		{
			std::cout << "Error! Could not determine z-Direction of barrel PolyData!" << std::endl;
			return NULL;
		}
		std::cout << "Calculating z-smoothing for " << barrel->GetNumberOfCells() << " slices..." << std::endl;
		for(int currentCell = 0/*85*/; currentCell < /*86*/barrel->GetNumberOfCells(); ++currentCell)
		{
			PolyDataPointerType barrelSubset = PolyDataPointerType::New();
			IdListPointerType barrelSubsetCellIds = IdListPointerType::New();
			barrelSubset->Allocate();
			int startID = std::max(currentCell - 25, 0);
// 			int stopID = std::min(currentCell + 25, int(barrel->GetNumberOfCells() - 1));
			int stopID = startID;
			int maxZ = lround(barrel->GetCell(currentCell)->GetBounds()[4]) + 25;
			int currZ = lround(barrel->GetCell(currentCell)->GetBounds()[4]);
			int currentCellID = currentCell - startID;
			std::map< int, int > planeIDs;
			int tmpIDCount = 0;
			unsigned int nrOfPlanes = 0;
			while(nrOfPlanes <= 50 && currZ < maxZ && stopID < barrel->GetNumberOfCells())
			{
				barrelSubsetCellIds->InsertId(stopID-startID, stopID);
				// right now there is no correction for several cells lying in the same plane.
				// current solution: insert all cells in one plane and count nr of planes instead of nr of cells
				currZ = lround(barrel->GetCell(stopID)->GetBounds()[4]);
				if(!planeIDs.count(currZ))
				{
					planeIDs.insert(std::pair< int, int >(currZ, tmpIDCount));	//makes all cells with same z point to one index
					++tmpIDCount;
					++nrOfPlanes;
				}
				++stopID;
			}
// 			for(int ii = startID; ii <= /*currentCell + 1*/stopID; ++ii)
// 			{
// 				barrelSubsetCellIds->InsertId(ii-startID, ii);
// 				// right now there is no correction for several cells lying in the same plane
// 				// maybe just insert all cells in one plane and count nr of planes instead of nr of cells
// // 				planeIDs.insert(std::pair< int, int >(zIndexVector[ii-startID], ii-startID));
// 				int currZ = lround(barrel->GetCell(ii)->GetBounds()[4]);
// 				if(!planeIDs.count(currZ))
// 				{
// 					planeIDs.insert(std::pair< int, int >(currZ, tmpIDCount));	//makes all cells with same z point to one index
// 					++tmpIDCount;
// 				}
// 			}
			barrelSubset->CopyCells(barrel, barrelSubsetCellIds);
// 			barrelSubset->Print(std::cout);
			unsigned int nrOfCells = barrelSubset->GetNumberOfCells();
// // 			unsigned int nrOfPlanes = planeIDs.size();
// 			std::cout << "nrOfCells = " << nrOfCells << std::endl;
// 			std::cout << "nrOfPlanes = " << nrOfPlanes << std::endl;
// 			std::cout << "planeIDs:" << std::endl;
			std::map< int, int >::const_iterator planeIDIt;
// 			for(planeIDIt = planeIDs.begin(); planeIDIt != planeIDs.end(); ++planeIDIt)
// 				std::cout << " z = " << planeIDIt->first << "\tid = " << planeIDIt->second << std::endl;
			// radial sampling of barrel contours
			std::cout << "starting radial sampling in plane " << currentCell << "..." << std::endl;
			double * centerPoint = new double[3];	// calculate centerPoint as center of top plane
			int subID;
			double pCoords[3], * weights;
			weights = new double[barrelSubset->GetCell(currentCellID)->GetNumberOfPoints()];
			barrelSubset->GetCell(currentCellID)->GetParametricCenter(pCoords);
			barrelSubset->GetCell(currentCellID)->EvaluateLocation(subID, pCoords, centerPoint, weights);
// 			std::cout << "centerPoint of cell w/ ID " << currentCellID << " @ [" << centerPoint[0] << "," << centerPoint[1] << "," << centerPoint[2] << "]" << std::endl;
			
			int noAngles = 36;
			std::vector< double * > thisPlaneAvgPoints;
			PlanePointerType thisPlane = PlanePointerType::New();
			CutterPointerType planeCutter = CutterPointerType::New();
			PointsPointerType thisAnglePoints = PointsPointerType::New();
			for(int ii = 0; ii < /*1*/noAngles/2; ++ii)
			{
// 				std::cout << "Sampling ray @ angle " << ii*360/noAngles << std::endl;
				
				double angle = ii*2*PI/noAngles;
				double direction[3];
				direction[0] = sin(angle);	// plane normal!!!
				direction[1] = cos(angle);
				direction[2] = 0;
				
				thisPlane->SetOrigin(centerPoint);
				thisPlane->SetNormal(direction);
				planeCutter->SetNumberOfContours(1);
				planeCutter->SetCutFunction(thisPlane);
				planeCutter->SetInput(barrelSubset);
				planeCutter->Update();
				PolyDataPointerType cutLines = planeCutter->GetOutput();	// cutter creates contoured lines from cutting planes w/ planes
				thisAnglePoints = cutLines->GetPoints();
// 				cutLines->Print(std::cout);
// 				thisAnglePoints->Print(std::cout);
				
				//getting rid of extra points created by vtkCutter:
				//calculate length of each line that lies in a plane; 
				//longest line has the two points of interest as endpoints
				std::map< double, std::pair< int, int > > * distanceIDs = new std::map< double, std::pair< int, int > >[nrOfPlanes];	//map automatically sorts pairs by their labels (==distance)
				std::vector< int > zAveragePoints1; // two vectors containing the IDs for all points belonging to one of the two respective z-averages
				std::vector< int > zAveragePoints2;
				for(int jj = 0; jj < thisAnglePoints->GetNumberOfPoints(); ++jj)
				{
					double * tmpCoords1 = new double[3];
					thisAnglePoints->GetPoint(jj, tmpCoords1);
					int currPlane = lround(tmpCoords1[2]);
					for(int ll = jj + 1; ll < thisAnglePoints->GetNumberOfPoints(); ++ll)
					{
						double * tmpCoords2 = new double[3];
						thisAnglePoints->GetPoint(ll, tmpCoords2);
						if(lround(tmpCoords2[2]) != currPlane)
						{
							delete [] tmpCoords2;
							continue;
						}
						double distance = sqrt((tmpCoords1[0]-tmpCoords2[0])*(tmpCoords1[0]-tmpCoords2[0]) + (tmpCoords1[1]-tmpCoords2[1])*(tmpCoords1[1]-tmpCoords2[1]));
						distanceIDs[planeIDs[currPlane]].insert(std::pair< double, std::pair< int, int> >(distance, std::pair< int, int >(jj, ll)));
						delete [] tmpCoords2;
					}
					delete [] tmpCoords1;
				}
// 				std::cout << "distanceIDs[0].size() = " << distanceIDs[0].size() << std::endl;
				//assign points to their respective z-average by comparing the x-y-distance from point from the previous plane.
				//because they are sampled from opposite sides of the barrel (which is cut approx. horizontally),
				//the distances should be significantly different
				if(!distanceIDs[0].size())
					continue;
				zAveragePoints1.push_back(distanceIDs[0].rbegin()->second.first);
				zAveragePoints2.push_back(distanceIDs[0].rbegin()->second.second);
				for(int jj = 1; jj < nrOfPlanes; ++jj)
				{
					if(distanceIDs[jj].size())
					{
						double * thisPlaneCoords1 = new double[3];
						double * thisPlaneCoords2 = new double[3];
						double * lastPlaneCoords = new double[3];
						int thisPlaneID1 = distanceIDs[jj].rbegin()->second.first;	//maps sort their keys in ascending order
						int thisPlaneID2 = distanceIDs[jj].rbegin()->second.second;
						int lastPlaneID = zAveragePoints1[zAveragePoints1.size()-1];	//in case last plane did not push back any ID
						thisAnglePoints->GetPoint(thisPlaneID1, thisPlaneCoords1);
						thisAnglePoints->GetPoint(thisPlaneID2, thisPlaneCoords2);
						thisAnglePoints->GetPoint(lastPlaneID, lastPlaneCoords);
						
						double dist1 = sqrt((thisPlaneCoords1[0] - lastPlaneCoords[0])*(thisPlaneCoords1[0] - lastPlaneCoords[0]) + (thisPlaneCoords1[1] - lastPlaneCoords[1])*(thisPlaneCoords1[1] - lastPlaneCoords[1]));
						double dist2 = sqrt((thisPlaneCoords2[0] - lastPlaneCoords[0])*(thisPlaneCoords2[0] - lastPlaneCoords[0]) + (thisPlaneCoords2[1] - lastPlaneCoords[1])*(thisPlaneCoords2[1] - lastPlaneCoords[1]));
						if(dist1 < dist2)
						{
							zAveragePoints1.push_back(thisPlaneID1);
							zAveragePoints2.push_back(thisPlaneID2);
						}
						else
						{
							zAveragePoints1.push_back(thisPlaneID2);
							zAveragePoints2.push_back(thisPlaneID1);
						}
						
						delete [] thisPlaneCoords1, delete[] thisPlaneCoords2, delete [] lastPlaneCoords;
					}
					else
					{
						std::flush(std::cout << "Warning! No z-average points found in plane " << jj + startID << " @ angle " << ii*360/noAngles << std::endl);
// 						std::flush(std::cout << "currentCellID = " << currentCellID <<std::endl);
						std::flush(std::cout << "currentCell = " << currentCell <<std::endl);
					}
				}
// 				std::flush(std::cout << "zAveragePoints1.size() = " << zAveragePoints1.size() << std::endl);
// 				std::flush(std::cout << "zAveragePoints2.size() = " << zAveragePoints2.size() << std::endl);
				if(zAveragePoints1.size() != zAveragePoints2.size())
					std::flush(std::cout << "Warning! z-average point lists not consistent!" << std::endl);
				
				double * avgCoords1 = new double[3];
				double * avgCoords2 = new double[3];
				for(int jj = 0; jj < 3; ++jj)
				{
					avgCoords1[jj] = 0;
					avgCoords2[jj] = 0;
				}
				int weight = 0;
				for(int jj = 0; jj < zAveragePoints1.size() && jj < zAveragePoints2.size(); ++jj, ++weight)
				{
					double * tmp1 = new double[3];
					double * tmp2 = new double[3];
					thisAnglePoints->GetPoint(zAveragePoints1[jj], tmp1);
					thisAnglePoints->GetPoint(zAveragePoints2[jj], tmp2);
					avgCoords1[0] += tmp1[0], avgCoords1[1] += tmp1[1];
					avgCoords2[0] += tmp2[0], avgCoords2[1] += tmp2[1];
					delete [] tmp1, delete [] tmp2;
				}
				if(weight)
					for(int jj = 0; jj < 2; ++jj)
					{
						avgCoords1[jj] /= weight;
						avgCoords2[jj] /= weight;
					}
				avgCoords1[2] = centerPoint[2], avgCoords2[2] = centerPoint[2];
				thisPlaneAvgPoints.push_back(avgCoords1);
				thisPlaneAvgPoints.push_back(avgCoords2);
				zAveragePoints1.clear();
				zAveragePoints2.clear();
			}
			allPlaneAvgPoints.push_back(thisPlaneAvgPoints);
		}
		PolyDataPointerType avgBarrel = PolyDataPointerType::New();
		PointsPointerType avgBarrelPoints = PointsPointerType::New();
		avgBarrel->Allocate(1);
		avgBarrelPoints->SetDataTypeToFloat();
		unsigned int totalNoPoints = 0;
		for(int ii = 0; ii < allPlaneAvgPoints.size(); ++ii)
			totalNoPoints += allPlaneAvgPoints[ii].size();
		avgBarrelPoints->SetNumberOfPoints(totalNoPoints);
		//we want each polygon to be non-intersecting (topologically, like a circle)
		//look for closest point of the two points on neighboring sampling ray
		//CAUTION: may give wrong results when sampling from very few rays
		//but still better than brute-force solving the traveling salesman...
		std::vector< std::vector< double * > >::iterator allPlanePointsIter;
		int currID = 0;
		for(allPlanePointsIter = allPlaneAvgPoints.begin(); allPlanePointsIter != allPlaneAvgPoints.end(); ++allPlanePointsIter)
		{
			PolygonPointerType thisPlanePoly = PolygonPointerType::New();
			thisPlanePoly->GetPointIds()->SetNumberOfIds(allPlanePointsIter->size());
			int polyPtID = 0;
			// 		std::cout << "no. of points in this plane: " << allPlanePointsIter->size() << std::endl;
			for(int ii = 0; ii < allPlanePointsIter->size()/2; ++ii)
			{
				if(ii == 0)
				{
					avgBarrelPoints->InsertPoint(currID, (*allPlanePointsIter)[0]);
					thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
					avgBarrelPoints->InsertPoint(currID + allPlanePointsIter->size()/2, (*allPlanePointsIter)[1]);
					thisPlanePoly->GetPointIds()->SetId(polyPtID + allPlanePointsIter->size()/2, currID + allPlanePointsIter->size()/2);
					++currID;
					++polyPtID;
				}
				else
				{
					double * previousPt = new double[3];
					double * pt1/* = new double[3]*/;
					double * pt2/* = new double[3]*/;
					pt1 = (*allPlanePointsIter)[2*ii];
					pt2 = (*allPlanePointsIter)[2*ii+1];
					avgBarrelPoints->GetPoint(currID - 1, previousPt);
					double dist1 = sqrt((previousPt[0] - pt1[0])*(previousPt[0] - pt1[0]) + (previousPt[1] - pt1[1])*(previousPt[1] - pt1[1]));
					double dist2 = sqrt((previousPt[0] - pt2[0])*(previousPt[0] - pt2[0]) + (previousPt[1] - pt2[1])*(previousPt[1] - pt2[1]));
					if(dist1 < dist2)
					{
						avgBarrelPoints->InsertPoint(currID, pt1);
						thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
						avgBarrelPoints->InsertPoint(currID + allPlanePointsIter->size()/2, pt2);
						thisPlanePoly->GetPointIds()->SetId(polyPtID + allPlanePointsIter->size()/2, currID + allPlanePointsIter->size()/2);
					}
					else
					{
						avgBarrelPoints->InsertPoint(currID, pt2);
						thisPlanePoly->GetPointIds()->SetId(polyPtID, currID);
						avgBarrelPoints->InsertPoint(currID + allPlanePointsIter->size()/2, pt1);
						thisPlanePoly->GetPointIds()->SetId(polyPtID + allPlanePointsIter->size()/2, currID + allPlanePointsIter->size()/2);
					}
					++currID;
					++polyPtID;
				}
			}
			currID += allPlanePointsIter->size()/2;
// 			thisPlanePoly->Print(std::cout);
			avgBarrel->InsertNextCell(thisPlanePoly->GetCellType(), thisPlanePoly->GetPointIds());
		}
// 		topPoints->Print(std::cout);
		avgBarrel->SetPoints(avgBarrelPoints);
		avgBarrel->Update();
// 		avgBarrel->Print(std::cout);
		delete spatialGraph;
		spatialGraph = new AmiraSpatialGraph();
		for(int ii = 0; ii < avgBarrel->GetNumberOfCells(); ++ii)
		{
			IdListPointerType ptIDs = avgBarrel->GetCell(ii)->GetPointIds();
			double * ptCoords = new double[3];
			avgBarrel->GetPoint(ptIDs->GetId(0), ptCoords);
			Vertex * pt1 = new Vertex(ptCoords, Pia);
			Vertex * pt2 = new Vertex(ptCoords, Pia);
			spatialGraph->addVertex(pt1);
			spatialGraph->addVertex(pt2);
			std::list< double * > edgePts;
			for(int jj = 0; jj < ptIDs->GetNumberOfIds(); ++jj)
			{
				double * tmpCoords = new double[3];
				avgBarrel->GetPoint(ptIDs->GetId(jj), tmpCoords);
				edgePts.push_back(tmpCoords);
			}
			edgePts.push_back(ptCoords);
			int connectivity[2];
			if(!spatialGraph->getNumberOfVertices())
			{
				connectivity[0] = 0;
				connectivity[1] = 1;
			}
			else
			{
				connectivity[0] = spatialGraph->getNumberOfVertices() - 2;
				connectivity[1] = spatialGraph->getNumberOfVertices() - 1;
			}
			Edge * edge = new Edge(connectivity, edgePts.size(), Pia, edgePts);
			spatialGraph->addEdge(edge);
		}
		
// 		ImageDataPointerType volume = createImageVolumeFromPolyData(polyData, label, xMin, xMax, yMin, yMax, zMin, zMax);
// 		ImageDataPointerType distVolume = distanceTransform(volume);
// 		distVolume->SetSpacing(volume->GetSpacing());
// 		return distVolume;
		return avgBarrel;
	}
	
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		return NULL;
	}
};

/****************************************************************************/
/*smoothes barrel contours along old z axis: compute points on every        */
/*contour in 10 degree steps and average all points in one direction within */
/*+- 25 micron in z                                                         */
/****************************************************************************/
PolyDataPointerType Geometry::smoothBarrelInZ2(int label)
{
// 	int label = Barrel;
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	PolyDataPointerType barrel = PolyDataPointerType::New();
	if(spatialGraph->extractLandmark(label, barrel, zIndexList))
	{
		barrel = correctOutlierContours(barrel);
		std::vector< std::vector< double * > > allPlaneAvgPoints;
		double * barrelTopBounds = barrel->GetCell(0)->GetBounds();
		int zDirection = 0;	//+/-1 depending on whether zIndexList parallel/antiparallel to orientation of barrel
		std::vector< int > zIndexVector;
// 		std::cout << "barrel bounds z = " << barrelTopBounds[4] << " --- zIndexList.front() = " << zIndexList.front() << " --- zIndexList.back() = " << zIndexList.back() << std::endl;
		if(lround(barrelTopBounds[4]) == zIndexList.front())
		{
			zDirection = 1;
			std::list< int >::iterator zListIt;
			for(zListIt = zIndexList.begin(); zListIt != zIndexList.end(); ++zListIt)
				zIndexVector.push_back(*zListIt);
// 			std::cout << "Barrel and zIndexList parallel!" << std::endl;
		}
		else if(lround(barrelTopBounds[4]) == zIndexList.back())
		{
			zDirection = -1;
			std::list< int >::reverse_iterator zListIt;
			for(zListIt = zIndexList.rbegin(); zListIt != zIndexList.rend(); ++zListIt)
				zIndexVector.push_back(*zListIt);
// 			std::cout << "Barrel and zIndexList antiparallel!" << std::endl;
		}
		else
		{
			std::cout << "Error! Could not determine z-Direction of barrel PolyData!" << std::endl;
			return NULL;
		}
		
		std::multimap< int, int > zToIDMap;
		std::vector< int > IDList;
		for(int ii = 0; ii < barrel->GetNumberOfCells(); ++ii)
			zToIDMap.insert(std::pair< int, int >(lround(barrel->GetCell(ii)->GetBounds()[4]), ii));
		for(int z = barrel->GetBounds()[4]; z <= barrel->GetBounds()[5]; ++z)
		{
			std::multimap< int, int >::iterator zToIDIt = zToIDMap.find(z);
			if(zToIDIt != zToIDMap.end())
				IDList.push_back(zToIDIt->second);
		}
		std::sort(IDList.begin(), IDList.end());
// 		std::cout << "Calculating z-smoothing for " << IDList.size() << " slices..." << std::endl;
// 		for(int currentCell = 0; currentCell < barrel->GetNumberOfCells(); ++currentCell)
		for(int curr = 0; curr < IDList.size(); ++curr)
		{
			std::vector< int >::iterator IDListIt;
			int currentCell = IDList[curr];
			PolyDataPointerType barrelSubset = PolyDataPointerType::New();
			IdListPointerType barrelSubsetCellIds = IdListPointerType::New();
			barrelSubset->Allocate();
			// automatic version
			int startID = std::max(curr - 25, 0);
// 			int stopID = std::min(currentCell + 25, int(barrel->GetNumberOfCells() - 1));
			int stopID = curr + 25;
			// BEGIN parameters for manual contours (delta z = 50)
// 			int startID = std::max(curr - 1, 0);
// // 			int stopID = std::min(currentCell + 25, int(barrel->GetNumberOfCells() - 1));
// 			int stopID = curr + 1;
			// END parameters for manual contours (delta z = 50)
			if(stopID >= IDList.size())
				stopID = IDList.size() - 1;
			for(int ii = startID; ii <= stopID; ++ii)
				barrelSubsetCellIds->InsertId(ii-startID, IDList[ii]);
			barrelSubset->CopyCells(barrel, barrelSubsetCellIds);
// 			barrelSubset->Print(std::cout);
// 			unsigned int nrOfCells = barrelSubset->GetNumberOfCells();
// 			std::cout << "nrOfCells = " << nrOfCells << std::endl;
			
			// radial sampling of barrel contours
// 			std::cout << "starting radial sampling in plane " << currentCell << "..." << std::endl;
			
			int nrAngles = 36;
			std::vector< std::vector< double * > > thisPlaneSamplingPoints;	// one vector for each sampling angle
			for(int ii = 0; ii < nrAngles; ++ii)
			{
				std::vector< double * > angleSamplingVec;
				thisPlaneSamplingPoints.push_back(angleSamplingVec);
			}
			for(int ii = 0; ii < barrelSubset->GetNumberOfCells(); ++ii)
			{
				double * centerPoint = new double[3];	// calculate centerPoint as center of individual planes
				double paramCenter[3];
				//polygon centroid
// 				IdTypeArrayPointerType idArray = IdTypeArrayPointerType::New();
// 				idArray->SetArray(barrelSubset->GetCell(ii)->GetPointIds()->GetPointer(0), barrelSubset->GetCell(ii)->GetNumberOfPoints(), 1);
// 				vtkPolygon::ComputeCentroid(idArray, barrelSubset->GetCell(ii)->GetPoints(), centerPoint);
				
				//parametric center
				int subID;
				double pCoords[3], * weights;
				weights = new double[barrelSubset->GetCell(ii)->GetNumberOfPoints()];
				barrelSubset->GetCell(ii)->GetParametricCenter(pCoords);
				barrelSubset->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
				
// 				double dist = sqrt((paramCenter[0] - centerPoint[0])*(paramCenter[0] - centerPoint[0]) + (paramCenter[1] - centerPoint[1])*(paramCenter[1] - centerPoint[1]) + (paramCenter[2] - centerPoint[2])*(paramCenter[2] - centerPoint[2]));
// 				if(curr == IDList.size()/2 && ii == barrelSubset->GetNumberOfCells()/2)
// 					std::cout << "Dist between parametric center and centroid = " << dist << " um" << std::endl;
				
				// identify polygon normals so each contour can be traversed in same direction
				double normal[3];
				vtkPolygon::ComputeNormal(barrelSubset->GetCell(ii)->GetPoints(), normal);
				bool reverseDirection;
				if(normal[2] < 0)
					reverseDirection = 1;
				else
					reverseDirection = 0;
				
				PlanePointerType thisPlane = PlanePointerType::New();
				PointsPointerType thisCellPoints = barrelSubset->GetCell(ii)->GetPoints();
				for(int jj = 0; jj < nrAngles; ++jj)	// has to be unique b/c order of the points is not clear
				{
	// 				std::cout << "Sampling ray @ angle " << ii*360/noAngles << std::endl;
					double angle = jj*2*PI/nrAngles;
					double planeNormal[3];
					if(reverseDirection)
					{
						planeNormal[0] = -sin(angle);	// plane normal!!!
						planeNormal[1] = -cos(angle);
					}
					else
					{
						planeNormal[0] = sin(angle);	// plane normal!!!
						planeNormal[1] = cos(angle);
					}
					planeNormal[2] = 0;
					thisPlane->SetOrigin(centerPoint);
					thisPlane->SetNormal(planeNormal);
					
					double * firstPt = new double[3];
					thisCellPoints->GetPoint(0, firstPt);
					double lastVal = thisPlane->EvaluateFunction(firstPt);
					delete [] firstPt;
					for(int kk = 1; kk <= thisCellPoints->GetNumberOfPoints(); ++kk)
					{
						double * thisPt = new double[3];
						thisCellPoints->GetPoint(kk%thisCellPoints->GetNumberOfPoints(), thisPt);
						double ptVal = thisPlane->EvaluateFunction(thisPt);
// 						std::cout << "value of plane fct @ point " << kk << " = " << ptVal << std::endl;
						if((ptVal < 0 && lastVal > 0)/* || (ptVal > 0 && lastVal < 0)*/)	// has to be unique b/c order of the points is not clear
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
						}
						else if(ptVal == 0 && lastVal > 0)
						{
							// case: this point is directly on plane
							double * samplePt =  new double[3];
							for(int ll = 0; ll < 3; ++ll)
								samplePt[ll] = thisPt[ll];
							thisPlaneSamplingPoints[jj].push_back(samplePt);
						}
// 						else if(ptVal != 0 && lastVal == 0)
// 						{
// 							;// do nothing
// 						}
// 						else if(ptVal == 0 && lastVal == 0)
// 						{
// 							;// either an error or both points lie directly on the plane (extremely unlikely)
// 						}
						lastVal = ptVal;
						delete [] thisPt;
					} // all cell points
				} // all angles
			} // all barrel subset cells
			
			std::vector< double * > thisPlaneAvgPoints;
			for(int ii = 0; ii < nrAngles; ++ii)
			{
				if(thisPlaneSamplingPoints[ii].size())
				{
					double * avgPt = new double[3];
					avgPt[0] = 0, avgPt[1] = 0, avgPt[2] = barrel->GetCell(currentCell)->GetBounds()[4];
					for(int jj = 0; jj < thisPlaneSamplingPoints[ii].size(); ++jj)
					{
						double * tmpPt = thisPlaneSamplingPoints[ii][jj];
						for(int kk = 0; kk < 2; ++kk)
							avgPt[kk] += tmpPt[kk];
					}
					for(int jj = 0; jj < 2; ++jj)
						avgPt[jj] = avgPt[jj]/double(thisPlaneSamplingPoints[ii].size());
					thisPlaneAvgPoints.push_back(avgPt);
				}
			}
			allPlaneAvgPoints.push_back(thisPlaneAvgPoints);
		} // all barrel cells
		
		PolyDataPointerType avgBarrel = PolyDataPointerType::New();
		PointsPointerType avgBarrelPoints = PointsPointerType::New();
		avgBarrel->Allocate();
		avgBarrelPoints->SetDataTypeToFloat();
		unsigned int totalNoPoints = 0;
		for(int ii = 0; ii < allPlaneAvgPoints.size(); ++ii)
			totalNoPoints += allPlaneAvgPoints[ii].size();
		avgBarrelPoints->SetNumberOfPoints(totalNoPoints);
		std::vector< std::vector< double * > >::iterator allPlanePointsIter;
		int currID = 0;
		for(allPlanePointsIter = allPlaneAvgPoints.begin(); allPlanePointsIter != allPlaneAvgPoints.end(); ++allPlanePointsIter)
		{
			PolygonPointerType thisPlanePoly = PolygonPointerType::New();
			thisPlanePoly->GetPointIds()->SetNumberOfIds(allPlanePointsIter->size());
	// 		std::cout << "no. of points in this plane: " << allPlanePointsIter->size() << std::endl;
			for(int ii = 0; ii < allPlanePointsIter->size(); ++ii)
			{
				avgBarrelPoints->InsertPoint(currID, (*allPlanePointsIter)[ii]);
				thisPlanePoly->GetPointIds()->SetId(ii, currID);
				++currID;
			}
			currID += allPlanePointsIter->size();
			avgBarrel->InsertNextCell(thisPlanePoly->GetCellType(), thisPlanePoly->GetPointIds());
		}
		avgBarrel->SetPoints(avgBarrelPoints);
		avgBarrel->Update();
// 		delete spatialGraph;
// 		spatialGraph = new AmiraSpatialGraph;
// 		spatialGraph->addPolyDataObject(avgBarrel, Soma);
		return avgBarrel;
	} // if(spatialGraph->extractLandmark())
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		PolyDataPointerType avgBarrel = PolyDataPointerType::New();
		avgBarrel->Allocate();
		return avgBarrel;
	}
};

/******************************************************************************/
/*simple avg contour by sampling barrel contours at regular angular intervals */
/*and then averaging in z                                                     */
/******************************************************************************/
void Geometry::computeAverageHomeBarrel(PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3])
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
		double * centerPoint = new double[3];
		double paramCenter[3];
		
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[completeBarrel->GetCell(ii)->GetNumberOfPoints()];
		completeBarrel->GetCell(ii)->GetParametricCenter(pCoords);
		completeBarrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		
		// identify polygon normals so each contour can be traversed in same direction
		double normal[3], rotAxis[3], zUnitVec[3] = {0,0,-1};
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
		
		// rotate sampling plane so it is perpendicular to 
		// the plane defined by the contour
		vtkMath::Cross(normal, zUnitVec, rotAxis);
		normalize(rotAxis);
		rotAngle = std::acos(vtkMath::Dot(normal, zUnitVec))*180/PI;
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
	
	if(spatialGraph2) spatialGraph2->addPolyDataObject(avgHBContour, Barrel);
};

/****************************************************************************/
/*shrink large outlying contours b/c they are segmentation artifacts,       */
/*not real barrel outlines.                                                 */
/****************************************************************************/
PolyDataPointerType Geometry::correctOutlierContours(PolyDataPointerType barrel)
{
	double * barrelCOM = calculateBarrelCentroid(barrel);
	double * barrelRadius = calculateBarrelRadius(barrel, barrelCOM);
// 	std::cout << "radius = " << barrelRadius[0] << " +- " << barrelRadius[1] << std::endl;
// 	std::cout << "innerThresh = " << innerThresh << std::endl;
	PolyDataPointerType filteredBarrel = PolyDataPointerType::New();
	filteredBarrel->Allocate();
	PointsPointerType filteredPoints = PointsPointerType::New();
	filteredPoints->SetDataTypeToFloat();
	int lastID = 0;
	int nrOfCells = barrel->GetNumberOfCells();
	for(int ii = 0; ii < nrOfCells; ++ii)
	{
		PointsPointerType cellPts = barrel->GetCell(ii)->GetPoints();
		PolygonPointerType poly = PolygonPointerType::New();
		poly->GetPointIds()->SetNumberOfIds(cellPts->GetNumberOfPoints());
		std::map< int, double * > orderedPoints;
		
		std::list< int > outlierList;
		std::list< int > regularList;
		bool outliers = 0;
		double meanDist = 0;
		for(int jj = 0; jj < cellPts->GetNumberOfPoints(); ++jj)
		{
			double * pt = new double[3];
			cellPts->GetPoint(jj, pt);
			double dist = sqrt((pt[0] - barrelCOM[0])*(pt[0] - barrelCOM[0]) + (pt[1] - barrelCOM[1])*(pt[1] - barrelCOM[1]));
			if(dist > barrelRadius[0] + 1.7*barrelRadius[1])
			{
				outliers = 1;
				outlierList.push_back(jj);
			}
			else
			{
				meanDist += dist;
				regularList.push_back(jj);
				orderedPoints.insert(std::pair< int, double * >(jj, pt));
			}
		}
		if(outliers)
		{
			if(regularList.size())
				meanDist = meanDist/double(regularList.size());
			else
				meanDist = barrelRadius[0];
			
			std::list< int > outlierListCopy(outlierList);
			std::list< int >::iterator outlierListCopyIt;
			for(outlierListCopyIt = outlierListCopy.begin(); outlierListCopyIt != outlierListCopy.end(); ++outlierListCopyIt)
			{
				for(int offset = -3; offset <= 3; ++offset)
				{
					int tmp = *outlierListCopyIt + offset;
					if(tmp >= 0 && tmp < cellPts->GetNumberOfPoints())
					{
						std::list< int >::iterator tmpIt = std::find(regularList.begin(), regularList.end(), tmp);
						if(tmpIt != regularList.end())
						{
							outlierList.push_back(tmp);
							regularList.erase(tmpIt);
							orderedPoints.erase(tmp);
						}
					}
				}
			}
			
			std::list< int >::iterator outlierListIt;
			for(outlierListIt = outlierList.begin(); outlierListIt != outlierList.end(); ++outlierListIt)
			{
				double pt[3], direction[3];
				double normD = 0;
				cellPts->GetPoint(*outlierListIt, pt);
				for(int kk = 0; kk < 2; ++kk)
				{
					direction[kk] = pt[kk] - barrelCOM[kk];
					normD += direction[kk]*direction[kk];
				}
				normD = sqrt(normD);
				if(normD)
					direction[0] /= normD, direction[1] /= normD;
				direction[2] = 0;
				double * newPt = new double[3];
				newPt[0] = barrelCOM[0] + meanDist*direction[0];
				newPt[1] = barrelCOM[1] + meanDist*direction[1];
				newPt[2] = pt[2];
				orderedPoints.insert(std::pair< int, double * >(*outlierListIt, newPt));
			}
		}
		
		std::map< int, double * >::iterator orderedPointsIt;
		for(orderedPointsIt = orderedPoints.begin(); orderedPointsIt != orderedPoints.end(); ++orderedPointsIt)
		{
			filteredPoints->InsertNextPoint(orderedPointsIt->second);
			poly->GetPointIds()->SetId(orderedPointsIt->first, orderedPointsIt->first + lastID);
		}
		filteredBarrel->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
		lastID += cellPts->GetNumberOfPoints();
	}
	filteredBarrel->SetPoints(filteredPoints);
	filteredBarrel->Update();
	
	return filteredBarrel;
};

/****************************************************************************/
/*computes new barrel axis from possible candidate axes in a x-y radius of  */
/*2000micron. rates candidate axes based on distance to Pia, angle to Pia   */
/*normal at intersection point and orientation wrt neighboring vessels      */
/****************************************************************************/
// std::vector< double > Geometry::newBarrelAxis(int label, std::list< unsigned int > vessels, PolyDataPointerType piaSurface, const char * outputFilename, double alpha)
double * Geometry::newBarrelAxis(PolyDataPointerType barrel, std::list< unsigned int > vessels, PolyDataPointerType piaSurface, const char * outputFilename, double alpha)
{
// 	std::vector< double > barrelParameters;
// 	PolyDataPointerType barrel = PolyDataPointerType::New();
// 	std::list< int > zIndexList;
// 	if(spatialGraph->extractLandmark(label, barrel, zIndexList))
// 	{
		double * barrelCOM = calculateBarrelCentroid(barrel);
		std::multimap< double, double * > scores = barrelAxisScores(piaSurface, barrelCOM, 2000, alpha);
		// following line for best vessel
		unsigned int directionVessel = vesselsAroundBarrel(barrelCOM, vessels);
		// following line for mean vessel
// 		double * vessel = vesselsAroundBarrel(barrelCOM, vessels);
		
		std::multimap< double, double * >::reverse_iterator scoreIt;
		std::multimap< double, double * >::iterator vesselIt;
		for(scoreIt = scores.rbegin(); scoreIt != scores.rend(); ++scoreIt)
		{
			bool keepAxis = 1;
			double * axis = scoreIt->second;
			/****** begin best vessel ******/
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
			/****** end best vessel ******/
			
			/****** begin mean vessel ******/
// 			double aNorm = 0, vNorm = 0, angle = 0;
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				aNorm += axis[ii]*axis[ii];
// 				vNorm += vessel[ii]*vessel[ii];
// 				angle += axis[ii]*vessel[ii];
// 			}
// 			aNorm = sqrt(aNorm);
// 			vNorm = sqrt(vNorm);
// 			angle = std::acos(std::abs(angle)/(aNorm*vNorm))*180.0/PI;
// 			if(angle > 5)
// 				keepAxis = 0;
			/****** end mean vessel ******/
			
			if(keepAxis)
				break;
		}
		if(scoreIt == scores.rend())
		{
			std::cout << "Warning! Could not find barrel axis satisfying the vessel constraints. Selecting axis with highest score." << std::endl;
			scoreIt = scores.rbegin();
		}
		
// 		std::vector< double * > endPoints;
// 		closeBarrelAlongNewAxis(newAxis, barrelCOM, barrel, endPoints);
// 		std::vector< double > parameters = computeBarrelParameters(barrel, newAxis, barrelCOM, endPoints, label);
// 		
// 		spatialGraph2->addPolyDataObject(barrel, label);
		
		
		/*** BEGIN MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
	// 	std::cout << "Writing 5 highest scoring axes out of " << scores.size() << std::endl;
// 		std::multimap< double, double * >::const_reverse_iterator axesIt;
// 		int count = 0;
// 		int nrOfAxisPoints = 2;
// 		for(axesIt = scores.rbegin(); axesIt != scores.rend() && count < 5; ++axesIt, ++count)
// 		{
// 			if(axesIt == scoreIt)
// 				continue;
// 			
// 			double * endPoint = new double[3];
// 			double * bottomPoint = new double[3];
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				endPoint[ii] = barrelCOM[ii] + axesIt->second[ii];
// 				bottomPoint[ii] = barrelCOM[ii] - 3*axesIt->second[ii];
// 			}
// // 			Vertex * newVert1 = new Vertex(endPoint, label);
// // 			Vertex * newVert2 = new Vertex(barrelCOM, label);
// // 			spatialGraph->addVertex(newVert1);
// // 			spatialGraph->addVertex(newVert2);
// 			int connectionIndex[2];
// 			if(!spatialGraph2->getNumberOfVertices())
// 			{
// 				connectionIndex[0] = 0;
// 				connectionIndex[1] = 1;
// 			}
// 			else
// 			{
// 				connectionIndex[0] = spatialGraph2->getNumberOfVertices() - 2;
// 				connectionIndex[1] = spatialGraph2->getNumberOfVertices() - 1;
// 			}
// // 			int noOfAxisPoints = 2;
// 			std::list< double * > axisCoords;
// 			axisCoords.push_back(endPoint);
// 			axisCoords.push_back(barrelCOM);
// 			Edge * newAxis = new Edge(connectionIndex, nrOfAxisPoints, ApicalDendrite, axisCoords);
// 			spatialGraph2->addEdge(newAxis);
// 		}
// 		//selected axis
// 		double * endPoint = new double[3];
// 		double * bottomPoint = new double[3];
// 		for(int ii = 0; ii < 3; ++ii)
// 		{
// 			endPoint[ii] = barrelCOM[ii] + scoreIt->second[ii];
// 			bottomPoint[ii] = barrelCOM[ii] - 3*scoreIt->second[ii];
// 		}
// // 		Vertex * newVert1 = new Vertex(endPoint, 0);
// // 		Vertex * newVert2 = new Vertex(bottomPoint, 0);
// // 		spatialGraph2->addVertex(newVert1);
// // 		spatialGraph2->addVertex(newVert2);
// 		int connectionIndex[2];
// 		if(!spatialGraph2->getNumberOfVertices())
// 		{
// 			connectionIndex[0] = 0;
// 			connectionIndex[1] = 1;
// 		}
// 		else
// 		{
// 			connectionIndex[0] = spatialGraph2->getNumberOfVertices() - 2;
// 			connectionIndex[1] = spatialGraph2->getNumberOfVertices() - 1;
// 		}
// // 		int noOfAxisPoints = 2;
// 		std::list< double * > axisCoords;
// 		axisCoords.push_back(endPoint);
// 		axisCoords.push_back(bottomPoint);
// 		Edge * newEdge = new Edge(connectionIndex, nrOfAxisPoints, Axon, axisCoords);
// 		spatialGraph2->addEdge(newEdge);
// 		
// 		double oldBottomPt[3], oldTopPt[3];
// 		oldBottomPt[0] = barrelCOM[0], oldBottomPt[1] = barrelCOM[1], oldBottomPt[2] = barrelCOM[2] - 3*scoreIt->second[2];
// 		oldTopPt[0] = barrelCOM[0], oldTopPt[1] = barrelCOM[1], oldTopPt[2] = barrelCOM[2] + scoreIt->second[2];
// // 		int connectionIndex[2];
// 		if(!spatialGraph2->getNumberOfVertices())
// 		{
// 			connectionIndex[0] = 0;
// 			connectionIndex[1] = 1;
// 		}
// 		else
// 		{
// 			connectionIndex[0] = spatialGraph2->getNumberOfVertices() - 2;
// 			connectionIndex[1] = spatialGraph2->getNumberOfVertices() - 1;
// 		}
// // 		int nrOfAxisPoints = 2;
// 		std::list< double * > oldAxisCoords;
// 		oldAxisCoords.push_back(oldBottomPt);
// 		oldAxisCoords.push_back(oldTopPt);
// 		Edge * oldAxis = new Edge(connectionIndex, nrOfAxisPoints, Neuron, oldAxisCoords);
// 		spatialGraph->addEdge(oldAxis);
		/*** END MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
		
		double * newAxis = new double[3];
		newAxis[0] = scoreIt->second[0], newAxis[1] = scoreIt->second[1], newAxis[2] = scoreIt->second[2];
		return newAxis;
// 		
// 		double newAxisTopZ, newAxisBottomZ;
// 		newAxisTopZ = sqrt((endPoint[0] - endPoints[0][0])*(endPoint[0] - endPoints[0][0]) + (endPoint[1] - endPoints[0][1])*(endPoint[1] - endPoints[0][1]) + (endPoint[2] - endPoints[0][2])*(endPoint[2] - endPoints[0][2]));
// 		newAxisBottomZ = sqrt((endPoint[0] - endPoints[1][0])*(endPoint[0] - endPoints[1][0]) + (endPoint[1] - endPoints[1][1])*(endPoint[1] - endPoints[1][1]) + (endPoint[2] - endPoints[1][2])*(endPoint[2] - endPoints[1][2]));
// 		barrelParameters.push_back(newAxisBottomZ);
// 		barrelParameters.push_back(newAxisTopZ);
// 		barrelParameters.push_back(parameters[0]);
// 		barrelParameters.push_back(parameters[1]);
// 	}
// 	return barrelParameters;
};

double * Geometry::simpleBarrelAxis(PolyDataPointerType barrel, PolyDataPointerType piaSurface, double alpha)
{
	double * barrelCOM = calculateBarrelCentroid(barrel);
	std::multimap< double, double * > scores = barrelAxisScores(piaSurface, barrelCOM, 2000, alpha);
	std::multimap< double, double * >::reverse_iterator scoreIt = scores.rbegin();
	double * newAxis = new double[3];
	newAxis[0] = scoreIt->second[0], newAxis[1] = scoreIt->second[1], newAxis[2] = scoreIt->second[2];
	return newAxis;
};

/******************************************************************************/
/*returns pairs of (double score, double * vec) where score is the score of   */
/*the potential barrel axis in direction of vec from barrelCentroid           */
/*to the Pia                                                                  */
/*radius defines the radius of the cylinder defining the ROI                  */
/*of Pia surface patches                                                      */
/******************************************************************************/
// PolyDataPointerType Geometry::barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha)
std::multimap< double, double * > Geometry::barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha)
{
	if(barrelCentroid != NULL)
	{
// 		std::cout << "Calculating pia ROI..." << std::endl;
		PolyDataPointerType piaROI = selectSurfaceRoi(piaSurface, barrelCentroid, radius);
		if(piaROI->GetNumberOfCells())
		{
// 			std::cout << "Calculating score values for " << piaROI->GetNumberOfCells() << " cells" << std::endl;
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
// 				double * centerPoint = new double[3];
// 				double * tmpAxis = new double[3];
// 				double axisNorm = 0;
// 				int subID;
// 				double pCoords[3], * weights;
// 				currentCell->GetParametricCenter(pCoords);
// 				currentCell->EvaluateLocation(subID, pCoords, centerPoint, weights);
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
// 				delete [] centerPoint;
			}
			
// 			std::flush(std::cout << "Done calculating score values" << std::endl);
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
};

PolyDataPointerType Geometry::selectSurfaceRoi(PolyDataPointerType surface, double * center, double radius)
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
// 	roiClipper->GetOutput()->Print(std::cout);
	return roiClipper->GetOutput();
};

/******************************************************************************/
/*returns pairs of (double dist, double * vec) where dist is the 2D distance  */
/*in the plane of the barrel center of all "good" vessels passed as IDs       */
/*in std::list< unsigned int > vessels. vec is direction of vessel            */
/******************************************************************************/
unsigned int Geometry::vesselsAroundBarrel(double * barrelCentroid, std::list< unsigned int > vessels)
// double * Geometry::vesselsAroundBarrel(double * barrelCentroid, std::list< unsigned int > vessels)
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
			for(int ii = 0; ii < 3; ++ii)
				direction[ii] /= dNorm;
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
};

void Geometry::enforceAxisDivergence(std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters)
{
	std::map< int, double * >::iterator barrelAxesIt;
	std::map< int, double * >::reverse_iterator barrelAxesRevIt;
// 	std::map< int, std::list< int > >::iterator barrelGridIt;
// 	std::map< int, std::list< int > >::reverse_iterator barrelGridRevIt;
	for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		normalize(barrelAxesIt->second);
	std::map< int, std::list< int > > barrelGrid = createBarrelGrid(barrelAxes);
	
	int maxLoops = 0;
	bool change = 0;
	do
	{
// 		std::cout << "loop nr " << maxLoops << std::endl;
		change = 0;
// 		for(barrelAxesRevIt = barrelAxes.rbegin(); barrelAxesRevIt != barrelAxes.rend(); ++barrelAxesRevIt)
// 		{
// 			int thisBarrel = barrelAxesRevIt->first;
// 			double thisAxis[3], cumulatedChange[3];
// 			thisAxis[0] = barrelAxesRevIt->second[0], thisAxis[1] = barrelAxesRevIt->second[1], thisAxis[2] = barrelAxesRevIt->second[2];
// 			cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
// 			std::list< int >::iterator neighborIt;
// 			for(neighborIt = barrelGrid[barrelAxesRevIt->first].begin(); neighborIt != barrelGrid[barrelAxesRevIt->first].end(); ++neighborIt)
// 			{
// 				double neighborAxis[3];
// 				neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
// 				double direction[3];
// 				for(int ii = 0; ii < 3; ++ii)
// 					direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
// 				normalize(direction);
// 				//eliminate part of neighboraxis that is parallel to this axis
// 				double prod = 0;
// 				for(int ii = 0; ii < 3; ++ii)
// 					prod += thisAxis[ii]*neighborAxis[ii];
// 				for(int ii = 0; ii < 3; ++ii)
// 					neighborAxis[ii] -= prod*thisAxis[ii];
// 				//check if other axis is pointing out of our local coordinate system
// 				double localDirection = 0;
// 				for(int ii = 0; ii < 3; ++ii)
// 					localDirection += direction[ii]*neighborAxis[ii];
// 				// < 0 means other axis is pointing towards us
// 				if(localDirection < 0)
// 				{
// 					change = 1;
// 					for(int ii = 0; ii < 3; ++ii)
// 					{
// 						cumulatedChange[ii] += 0.5*localDirection*direction[ii];
// 						barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
// 					}
// 					normalize(barrelAxes[*neighborIt]);
// 				}
// 			}
// 			for(int ii = 0; ii < 3; ++ii)
// 				barrelAxesRevIt->second[ii] += cumulatedChange[ii];
// 			normalize(barrelAxesRevIt->second);
// 		}
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
				normalize(direction);
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
					normalize(barrelAxes[*neighborIt]);
				}
			}
			for(int ii = 0; ii < 3; ++ii)
				barrelAxesIt->second[ii] += cumulatedChange[ii];
			normalize(barrelAxesIt->second);
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
					normalize(direction);
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
						normalize(barrelAxes[*neighborIt]);
					}
				}
				for(int ii = 0; ii < 3; ++ii)
					barrelAxesRevIt->second[ii] += cumulatedChange[ii];
				normalize(barrelAxesRevIt->second);
			}
// 			for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
// 			{
// 				int thisBarrel = barrelAxesIt->first;
// 				double thisAxis[3], cumulatedChange[3];
// 				thisAxis[0] = barrelAxesIt->second[0], thisAxis[1] = barrelAxesIt->second[1], thisAxis[2] = barrelAxesIt->second[2];
// 				cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
// 				std::list< int >::iterator neighborIt;
// 				for(neighborIt = barrelGrid[barrelAxesIt->first].begin(); neighborIt != barrelGrid[barrelAxesIt->first].end(); ++neighborIt)
// 				{
// 					double neighborAxis[3];
// 					neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
// 					double direction[3];
// 					for(int ii = 0; ii < 3; ++ii)
// 						direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
// 					normalize(direction);
// 					//eliminate part of neighboraxis that is parallel to this axis
// 					double prod = 0;
// 					for(int ii = 0; ii < 3; ++ii)
// 						prod += thisAxis[ii]*neighborAxis[ii];
// 					for(int ii = 0; ii < 3; ++ii)
// 						neighborAxis[ii] -= prod*thisAxis[ii];
// 					//check if other axis is pointing out of our local coordinate system
// 					double localDirection = 0;
// 					for(int ii = 0; ii < 3; ++ii)
// 						localDirection += direction[ii]*neighborAxis[ii];
// 					// < 0 means other axis is pointing towards us
// 					if(localDirection < 0)
// 					{
// 						change = 1;
// 						for(int ii = 0; ii < 3; ++ii)
// 						{
// 							cumulatedChange[ii] += 0.5*localDirection*direction[ii];
// 							barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
// 						}
// 						normalize(barrelAxes[*neighborIt]);
// 					}
// 				}
// 				for(int ii = 0; ii < 3; ++ii)
// 					barrelAxesIt->second[ii] += cumulatedChange[ii];
// 				normalize(barrelAxesIt->second);
// 			}
		}
		++maxLoops;
	} while(change && maxLoops < 10);
	
};

/******************************************************************************/
/*create grid of barrels actually present in data: each barrel label has a    */
/*LUT of neighbors associated with it (i.e., create a graph)                  */
/******************************************************************************/
std::map< int, std::list< int > > Geometry::createBarrelGrid(std::map< int, double * > barrelAxes)
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
};

/******************************************************************************/
/*create grid of barrels actually present in data: each barrel label has a    */
/*LUT of neighbors associated with it (i.e., create a graph)                  */
/******************************************************************************/
std::map< int, std::list< int > > Geometry::createBarrelGrid(std::map< int, Column * > barrelColumns)
{
	std::map< int, std::list< int > > grid;
	std::list< int > barrelIDs;
	std::map< int, Column * >::const_iterator barrelColsIt;
	for(barrelColsIt = barrelColumns.begin(); barrelColsIt != barrelColumns.end(); ++barrelColsIt)
		barrelIDs.push_back(barrelColsIt->first);
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
};

/******************************************************************************/
/*create grid of barrels actually present in data: each barrel label has a    */
/*LUT of neighbors associated with it (i.e., create a graph)                  */
/******************************************************************************/
std::map< int, std::list< int > > Geometry::createNearestNeighborBarrelGrid(std::map< int, double * > barrelCenters)
{
	std::map< int, std::list< int > > grid;
	std::list< int > barrelIDs;
	std::map< int, double * >::const_iterator barrelColsIt;
	for(barrelColsIt = barrelCenters.begin(); barrelColsIt != barrelCenters.end(); ++barrelColsIt)
		barrelIDs.push_back(barrelColsIt->first);
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
	a1.push_back(Alpha), a1.push_back(B1), a1.push_back(A2);
	std::list< int > a2;
	a2.push_back(A1), a2.push_back(B2), a2.push_back(A3);
	std::list< int > a3;
	a3.push_back(A2), a3.push_back(B3), a3.push_back(A4);
	std::list< int > a4;
	a4.push_back(A3), a4.push_back(B4);
	std::list< int > b1;
	b1.push_back(Beta), b1.push_back(C1), b1.push_back(B2), b1.push_back(A1), b1.push_back(Alpha);
	std::list< int > b2;
	b2.push_back(B1), b2.push_back(C2), b2.push_back(B3), b2.push_back(A2);
	std::list< int > b3;
	b3.push_back(B2), b3.push_back(C3), b3.push_back(B4), b3.push_back(A3);
	std::list< int > b4;
	b4.push_back(B3), b4.push_back(C4), b4.push_back(A4);
	std::list< int > c1;
	c1.push_back(Gamma), c1.push_back(D1), c1.push_back(C2), c1.push_back(B1), c1.push_back(Beta);
	std::list< int > c2;
	c2.push_back(C1), c2.push_back(D2), c2.push_back(C3), c2.push_back(B2);
	std::list< int > c3;
	c3.push_back(C2), c3.push_back(D3), c3.push_back(C4), c3.push_back(B3);
	std::list< int > c4;
	c4.push_back(C3), c4.push_back(D4), c4.push_back(C5), c4.push_back(B4);
	std::list< int > c5;
	c5.push_back(C4), c5.push_back(D5), c5.push_back(C6);
	std::list< int > c6;
	c6.push_back(C5), c6.push_back(D6);
	std::list< int > d1;
	d1.push_back(Delta), d1.push_back(E1), d1.push_back(D2), d1.push_back(C1), d1.push_back(Gamma);
	std::list< int > d2;
	d2.push_back(D1), d2.push_back(E2), d2.push_back(D3), d2.push_back(C2);
	std::list< int > d3;
	d3.push_back(D2), d3.push_back(E3), d3.push_back(D4), d3.push_back(C3);
	std::list< int > d4;
	d4.push_back(D3), d4.push_back(E4), d4.push_back(D5), d4.push_back(C4);
	std::list< int > d5;
	d5.push_back(D4), d5.push_back(E5), d5.push_back(D6), d5.push_back(C5);
	std::list< int > d6;
	d6.push_back(D5), d6.push_back(E6), d6.push_back(C6);
	std::list< int > e1;
	e1.push_back(Delta), e1.push_back(E2), e1.push_back(D1);
	std::list< int > e2;
	e2.push_back(E1), e2.push_back(E3), e2.push_back(D2);
	std::list< int > e3;
	e3.push_back(E2), e3.push_back(E4), e3.push_back(D3);
	std::list< int > e4;
	e4.push_back(E3), e4.push_back(E5), e4.push_back(D4);
	std::list< int > e5;
	e5.push_back(E4), e5.push_back(E6), e5.push_back(D5);
	std::list< int > e6;
	e6.push_back(E5), e6.push_back(D6);
	
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
};

double * Geometry::calculateBarrelCentroid(PolyDataPointerType barrel)
{
// 	std::cout << "Calculating barrel centroid" << std::endl;
	double * centroid = new double[3];
	for(int ii = 0; ii < 3; ++ii)
		centroid[ii] = 0;
	int nrOfCells = barrel->GetNumberOfCells();
// 		std::cout << "nrOfCells = " << nrOfCells << std::endl;
	for(int ii = 0; ii < nrOfCells; ++ii)
	{
// 			std::cout << "cellID = " << ii << std::endl;
		PointsPointerType tmpPoints = barrel->GetCell(ii)->GetPoints();
		double centerPoint[] = {0, 0, 0};
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
		{
			centerPoint[jj] = centerPoint[jj]/(double)nrOfPoints;
			centroid[jj] += centerPoint[jj];
		}
	}
	for(int jj = 0; jj < 3; ++jj)
		centroid[jj] = centroid[jj]/(double)nrOfCells;
// 		std::flush(std::cout << "centroid = [" << centroid[0] << "," << centroid[1] << "," << centroid[2] << "]" << std::endl);
	return centroid;
};

double * Geometry::calculateBarrelRadius(PolyDataPointerType barrel, double * barrelCenter)
{
// 	std::flush(std::cout << "Calculating mean barrel radius" << std::endl);
	std::list< double > radiusList;
	int nrOfCells = barrel->GetNumberOfCells();
	for(int ii = 0; ii < nrOfCells; ++ii)
	{
// 		double * centerPoint = new double[3];	// calculate centerPoint as center of individual planes
// 		int subID;
// 		double pCoords[3], * weights;
// 		weights = new double[barrel->GetCell(ii)->GetNumberOfPoints()];
// 		barrel->GetCell(ii)->GetParametricCenter(pCoords);
// 		barrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		PointsPointerType cellPts = barrel->GetCell(ii)->GetPoints();
		for(int jj = 0; jj < cellPts->GetNumberOfPoints(); ++jj)
		{
			double * pt = new double[3];
			cellPts->GetPoint(jj, pt);
// 			double radius = sqrt((pt[0] - centerPoint[0])*(pt[0] - centerPoint[0]) + (pt[1] - centerPoint[1])*(pt[1] - centerPoint[1]));
			double radius = sqrt((pt[0] - barrelCenter[0])*(pt[0] - barrelCenter[0]) + (pt[1] - barrelCenter[1])*(pt[1] - barrelCenter[1]));
			radiusList.push_back(radius);
			delete [] pt;
		}
	}
	
	double mean = 0, stddev = 0;
	std::list< double >::iterator radiusListIt;
	for(radiusListIt = radiusList.begin(); radiusListIt != radiusList.end(); ++radiusListIt)
		mean += *radiusListIt;
	if(radiusList.size())
		mean = mean/double(radiusList.size());
	for(radiusListIt = radiusList.begin(); radiusListIt != radiusList.end(); ++radiusListIt)
		stddev += (*radiusListIt - mean)*(*radiusListIt - mean);
	if(radiusList.size())
		stddev = stddev/double(radiusList.size());
	stddev = sqrt(stddev);
	
	double * stats = new double[2];
	stats[0] = mean, stats[1] = stddev;
	
	return stats;
};

void Geometry::transformBarrelZAxis(double * newAxis, double * barrelCentroid, PolyDataPointerType barrel)
{
	TransformPointerType translate = TransformPointerType::New();
	TransformPointerType rotate = TransformPointerType::New();
	TransformPointerType inverseTranslate = TransformPointerType::New();
	TransformFilterType barrelTransform = TransformFilterType::New();
	
	//old axis == z axis
	//rotation axis = newAxus x oldAxis
	//then, angle of rotation = acos(oldAxis*newAxis)
	double rotationAxis[3];
	rotationAxis[0] = newAxis[1];
	rotationAxis[1] = -newAxis[0];
	rotationAxis[2] = 0;
	
	double axisNorm = sqrt(newAxis[0]*newAxis[0] + newAxis[1]*newAxis[1]+ newAxis[2]*newAxis[2]);
	double angle = 0;
	if(axisNorm)
		angle = acos(-newAxis[2]/axisNorm)*180/PI;
	
	translate->Translate(-barrelCentroid[0], -barrelCentroid[1], -barrelCentroid[2]);
	rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
	inverseTranslate->Translate(barrelCentroid[0], barrelCentroid[1], barrelCentroid[2]);
	
	rotate->Concatenate(translate);
	inverseTranslate->Concatenate(rotate);
	inverseTranslate->Update();
// 	inverseTranslate->Print(std::cout);
	
	barrelTransform->SetTransform(inverseTranslate);
	barrelTransform->SetInput(barrel);
	barrelTransform->Update();
	barrelTransform->GetOutput()->Print(std::cout);
// 	barrel = barrelTransform->GetOutput();
	barrel->DeepCopy(barrelTransform->GetOutput());
	barrel->Update();
	barrel->Print(std::cout);
};

/******************************************************************************/
/*calculates avg contour of center contours (all points in middle third of    */
/*vertical extent of barrel along new axis are included)                      */
/******************************************************************************/
PolyDataPointerType Geometry::smoothBarrelAlongNewAxis(double* newAxis, double* barrelCentroid, PolyDataPointerType barrel, std::vector< double* > endPoints)
{
	if(barrel->GetNumberOfPoints())
	{
// 		barrel->Print(std::cout);
// 		std::flush(std::cout << "Calculate avg barrel contour" << std::endl);
		if(barrel->GetCell(0)->GetNumberOfPoints() != 36)
		{
			std::cout << "Warning! Barrel has not been sampled correctly! max barrel contour may be corrupted..." << std::endl;
		}
		double zAxis[3];
		zAxis[0] = newAxis[0], zAxis[1] = newAxis[1], zAxis[2] = newAxis[2];
		normalize(zAxis);
		double axisPt1[3], axisPt2[3];
		axisPt1[0] = barrelCentroid[0] + 1000*zAxis[0], axisPt1[1] = barrelCentroid[1] + 1000*zAxis[1], axisPt1[2] = barrelCentroid[2] + 1000*zAxis[2];
		axisPt2[0] = barrelCentroid[0] - 1000*zAxis[0], axisPt2[1] = barrelCentroid[1] - 1000*zAxis[1], axisPt2[2] = barrelCentroid[2] - 1000*zAxis[2];
		double zExtent = sqrt((endPoints[0][0] - endPoints[1][0])*(endPoints[0][0] - endPoints[1][0]) + (endPoints[0][1] - endPoints[1][1])*(endPoints[0][1] - endPoints[1][1]) + (endPoints[0][2] - endPoints[1][2])*(endPoints[0][2] - endPoints[1][2]));
// 		double * barrelRadius = calculateBarrelRadius(barrel, barrelCentroid);
// 		std::flush(std::cout << "barrelRadius = " << barrelRadius[0] << " +- " << barrelRadius[1] << std::endl);
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
				double dist1 = sqrt((endPoints[0][0] - closestPt[0])*(endPoints[0][0] - closestPt[0]) + (endPoints[0][1] - closestPt[1])*(endPoints[0][1] - closestPt[1]) + (endPoints[0][2] - closestPt[2])*(endPoints[0][2] - closestPt[2]));
				double dist2 = sqrt((endPoints[1][0] - closestPt[0])*(endPoints[1][0] - closestPt[0]) + (endPoints[1][1] - closestPt[1])*(endPoints[1][1] - closestPt[1]) + (endPoints[1][2] - closestPt[2])*(endPoints[1][2] - closestPt[2]));
// 				if(dist1 > zExtent*0.33 && dist2 > zExtent*0.33)
				{
					for(int kk = 0; kk < 3; ++kk)
						closestPt[kk] = pt[kk] - closestPt[kk];
					avgPoints[jj].push_back(closestPt);
					avgDistances[jj].push_back(dist);
				}
// 				else
// 					delete [] closestPt;
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
				normalize(avgVec);
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
				normalize(avgVec);
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
		std::cout << "Error! PolyData barrel is empty! Could not calculate max barrel contour." << std::endl;
		return NULL;
	}
};

/******************************************************************************/
/*adds polygons on top and bottom to make barrel extend up to extreme points  */
/*of existing contours measured along the new z axis.                         */
/*Stores the extreme points in vector endpoints in the order top, bottom      */
/******************************************************************************/
void Geometry::closeBarrelAlongNewAxis(double* newAxis, double* barrelCentroid, PolyDataPointerType barrel, std::vector< double* >& endPoints)
{
	if(barrel->GetNumberOfCells())
	{
		PolyDataPointerType barrelTop = PolyDataPointerType::New();
		PolyDataPointerType barrelBottom = PolyDataPointerType::New();
		IdListPointerType barrelTopCellIds = IdListPointerType::New();
		IdListPointerType barrelBottomCellIds = IdListPointerType::New();
		barrelTop->Allocate();
		barrelBottom->Allocate();
		
		// estimate how many planes should be used
		// to be rotated for barrel caps: ~ base * sin (a)
		// for base use ~ diagonal of bounding box / sqrt(2) (maybe *factor???)
		double normA = 0;
		double alpha = 0;
		double newZAxis[3];
		newZAxis[0] = newAxis[0], newZAxis[1] = newAxis[1], newZAxis[2] = newAxis[2];
		normalize(newZAxis);
		alpha = acos(std::abs(newZAxis[2]));
		
		int nrOfTopCapPlanes = 1;
		int nrOfBottomCapPlanes = 1;
		nrOfTopCapPlanes = std::min(nrOfTopCapPlanes, (int)barrel->GetNumberOfCells());
		nrOfBottomCapPlanes = std::min(nrOfBottomCapPlanes, (int)barrel->GetNumberOfCells());
		
		// insert IDs in order out->in so that ID 0 is either the global top or bottom polygon
		for(int ii = 0; ii < nrOfTopCapPlanes; ++ii)
			barrelTopCellIds->InsertId(ii, ii);
		for(int ii = 0, jj = barrel->GetNumberOfCells() - 1; ii < nrOfBottomCapPlanes && jj >= 0; ++ii, --jj)
			barrelBottomCellIds->InsertId(ii, jj);
		
		barrelTop->CopyCells(barrel, barrelTopCellIds);
		barrelBottom->CopyCells(barrel, barrelBottomCellIds);
		
		// now, careful! barrel top has globally lowest z coordinates and vice versa
		// find extreme bottom/top point by transforming the coordinate system
		// for every point in lowest/highest polygon and looking for extreme z values
		TransformPointerType translate = TransformPointerType::New();
		TransformPointerType rotate = TransformPointerType::New();
		TransformPointerType inverseTranslate = TransformPointerType::New();
		
		//old axis == z axis
		//rotation axis = newAxus x oldAxis
		//then, angle of rotation = acos(oldAxis*newAxis)
		double rotationAxis[3];
		rotationAxis[0] = newZAxis[1];
		rotationAxis[1] = -newZAxis[0];
		rotationAxis[2] = 0;
		
		double angle = -alpha*180/PI;	// -1 b/c we want to rotate coordinate system, not points
		
		translate->Translate(-barrelCentroid[0], -barrelCentroid[1], -barrelCentroid[2]);
		rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
		inverseTranslate->Translate(barrelCentroid[0], barrelCentroid[1], barrelCentroid[2]);
		
		rotate->Concatenate(translate);
		inverseTranslate->Concatenate(rotate);
		inverseTranslate->Update();
		
		int extremeBottomID = 0, extremeTopID = 0;
		double bottomMax = -1E06, topMin = 1E06;
		PointsPointerType topPts = barrelTop->GetCell(0)->GetPoints();
		PointsPointerType bottomPts = barrelBottom->GetCell(0)->GetPoints();
		for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			double homPt[4], transPt[4];
			topPts->GetPoint(ii, tmpPt);
			homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
			inverseTranslate->MultiplyPoint(homPt, transPt);
			if(transPt[2] < topMin)
			{
				topMin = transPt[2];
				extremeTopID = ii;
			}
			delete [] tmpPt;
		}
		for(int ii = 0; ii < bottomPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			double homPt[4], transPt[4];
			bottomPts->GetPoint(ii, tmpPt);
			homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
			inverseTranslate->MultiplyPoint(homPt, transPt);
			if(transPt[2] > bottomMax)
			{
				bottomMax = transPt[2];
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
// 			std::cout << "inside top!!!" << std::endl;
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			endPoints.push_back(finalTopPt);
			//could do something here to make it look nice...not crucial...
// 			for(int ii = 0; ii < 3; ++ii)
// 				finalTopPt[ii] += finalTopPt[ii] - extremeTopPt[ii];
		}
		else if(!insideTop)
		{
// 			std::cout << "not inside top!!!" << std::endl;
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			double * finalTopPt2 = new double[3];
			finalTopPt2[0] = finalTopPt[0], finalTopPt2[1] = finalTopPt[1], finalTopPt2[2] = finalTopPt[2];
			endPoints.push_back(finalTopPt2);
			double * centerPoint = new double[3];
			int subID;
			double pCoords[3], * weights;
			weights = new double[barrelTop->GetCell(0)->GetNumberOfPoints()];
			barrelTop->GetCell(0)->GetParametricCenter(pCoords);
			barrelTop->GetCell(0)->EvaluateLocation(subID, pCoords, centerPoint, weights);
			PlanePointerType topPlane = PlanePointerType::New();
			topPlane->SetOrigin(finalTopPt);
			topPlane->SetNormal(newZAxis);
			double planeDist = std::abs(topPlane->DistanceToPlane(centerPoint));
			finalTopPt[0] = centerPoint[0], finalTopPt[1] = centerPoint[1], finalTopPt[2] = centerPoint[2];
			if(cos(alpha))
				finalTopPt[2] -= planeDist/cos(alpha);
		}
		if(insideBottom)
		{
// 			std::cout << "inside bottom!!!" << std::endl;
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			endPoints.push_back(finalBottomPt);
			//could do something here to make it look nice...not crucial...
// 			for(int ii = 0; ii < 3; ++ii)
// 				finalBottomPt[ii] += finalBottomPt[ii] - extremeBottomPt[ii];
		}
		else if(!insideBottom)
		{
// 			std::cout << "not inside bottom!!!" << std::endl;
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			double * finalBottomPt2 = new double[3];
			finalBottomPt2[0] = finalBottomPt[0], finalBottomPt2[1] = finalBottomPt[1], finalBottomPt2[2] = finalBottomPt[2];
			endPoints.push_back(finalBottomPt2);
			double * centerPoint = new double[3];
			int subID;
			double pCoords[3], * weights;
			weights = new double[barrelBottom->GetCell(0)->GetNumberOfPoints()];
			barrelBottom->GetCell(0)->GetParametricCenter(pCoords);
			barrelBottom->GetCell(0)->EvaluateLocation(subID, pCoords, centerPoint, weights);
			PlanePointerType bottomPlane = PlanePointerType::New();
			bottomPlane->SetOrigin(finalBottomPt);
			bottomPlane->SetNormal(newZAxis);
			double planeDist = std::abs(bottomPlane->DistanceToPlane(centerPoint));
			finalBottomPt[0] = centerPoint[0], finalBottomPt[1] = centerPoint[1], finalBottomPt[2] = centerPoint[2];
			if(cos(alpha))
				finalBottomPt[2] += planeDist/cos(alpha);
		}
		// start linear interpolation of contours up to final points
// 		std::flush(std::cout << "start linear interpolation of contours up to final points" << std::endl);
		int nrOfTopPolys = 0, nrOfBottomPolys = 0;
		nrOfTopPolys = int(std::abs(finalTopPt[2] - extremeTopPt[2]));
		nrOfBottomPolys = int(std::abs(finalBottomPt[2] - extremeBottomPt[2]));
		std::vector< double * > topContourSlopes;
		std::vector< double * > bottomContourSlopes;
		for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			topPts->GetPoint(ii, tmpPt);
			for(int jj = 0; jj < 3; ++jj)
				tmpPt[jj] = (finalTopPt[jj] - tmpPt[jj])/(double)(nrOfTopPolys);
			topContourSlopes.push_back(tmpPt);
// 			std::flush(std::cout << "topContourSlopes[" << ii << "] = " << "[" << topContourSlopes.back()[0] << "," << topContourSlopes.back()[1] << "," << topContourSlopes.back()[2] << "]" << std::endl);
		}
		for(int ii = 0; ii < bottomPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			bottomPts->GetPoint(ii, tmpPt);
			for(int jj = 0; jj < 3; ++jj)
				tmpPt[jj] = (finalBottomPt[jj] - tmpPt[jj])/(double)(nrOfBottomPolys);
			bottomContourSlopes.push_back(tmpPt);
		}
// 		std::flush(std::cout << "topContourSlopes.size() = " << topContourSlopes.size() << std::endl);
// 		std::flush(std::cout << "bottomContourSlopes.size() = " << bottomContourSlopes.size() << std::endl);
		
		// first top
// 		std::flush(std::cout << "first top: " << nrOfTopPolys << " addtl polygons" << std::endl);
		PolyDataPointerType newBarrelTop = PolyDataPointerType::New();
		PointsPointerType newBarrelTopPoints = PointsPointerType::New();
		newBarrelTop->Allocate();
		newBarrelTopPoints->SetDataTypeToFloat();
		int lastID = 0;
		double startZ = 0;
		for(int ii = 1; ii < nrOfTopPolys; ++ii)	// leave out last polygon on purpose to make sure we do not have some self-intersecting mess
		{
			PolygonPointerType newPoly = PolygonPointerType::New();
			newPoly->GetPointIds()->SetNumberOfIds(topPts->GetNumberOfPoints());
			for(int jj = 0; jj < topPts->GetNumberOfPoints(); ++jj)
			{
				double * tmpPt = new double[3];
				topPts->GetPoint(jj, tmpPt);
				for(int kk = 0; kk < 3; ++kk)
					tmpPt[kk] += ii*topContourSlopes[jj][kk];
				if(ii == 1)
				{
					tmpPt[2] = round(tmpPt[2]);
					startZ = tmpPt[2];
				}
				else
					tmpPt[2] = startZ - ii + 1;
				newBarrelTopPoints->InsertNextPoint(tmpPt);
				newPoly->GetPointIds()->SetId(jj, jj + lastID);
			}
			lastID += topPts->GetNumberOfPoints();
			newBarrelTop->InsertNextCell(newPoly->GetCellType(), newPoly->GetPointIds());
		}
		newBarrelTop->SetPoints(newBarrelTopPoints);
		newBarrelTop->Update();
		
		//now bottom
// 		std::flush(std::cout << "now bottom: " << nrOfBottomPolys << " addtl polygons" << std::endl);
		PolyDataPointerType newBarrelBottom = PolyDataPointerType::New();
		PointsPointerType newBarrelBottomPoints = PointsPointerType::New();
		newBarrelBottom->Allocate();
		newBarrelBottomPoints->SetDataTypeToFloat();
		lastID = 0;
		for(int ii = 1; ii < nrOfBottomPolys; ++ii)	// leave out last polygon on purpose to make sure we do not have some self-intersecting mess
		{
			PolygonPointerType newPoly = PolygonPointerType::New();
			newPoly->GetPointIds()->SetNumberOfIds(bottomPts->GetNumberOfPoints());
			for(int jj = 0; jj < bottomPts->GetNumberOfPoints(); ++jj)
			{
				double * tmpPt = new double[3];
				bottomPts->GetPoint(jj, tmpPt);
				for(int kk = 0; kk < 3; ++kk)
					tmpPt[kk] += ii*bottomContourSlopes[jj][kk];
				if(ii == 1)
				{
					tmpPt[2] = round(tmpPt[2]);
					startZ = tmpPt[2];
				}
				else
					tmpPt[2] = startZ + ii - 1;
				newBarrelBottomPoints->InsertNextPoint(tmpPt);
				newPoly->GetPointIds()->SetId(jj, jj + lastID);
			}
			lastID += topPts->GetNumberOfPoints();
			newBarrelBottom->InsertNextCell(newPoly->GetCellType(), newPoly->GetPointIds());
		}
		newBarrelBottom->SetPoints(newBarrelBottomPoints);
		newBarrelBottom->Update();
// 		newBarrelTop->Print(std::cout);
// 		barrel->Print(std::cout);
// 		newBarrelBottom->Print(std::cout);
		/*** BEGIN MANUAL BARREL COMMENT ***/
		// fill empty contours (for whatever reason they might be empty...)
		// gaps should be small, b/c this is just a constant continuation
		std::list< int > zList;
		std::list< int > missingZList;
		int minZ = 0, maxZ = 0;
		
		if(newBarrelTop->GetNumberOfPoints())
			minZ = lround(newBarrelTop->GetCell(newBarrelTop->GetNumberOfCells()-1)->GetBounds()[4]);
		else
			minZ = lround(barrel->GetBounds()[4]);
		if(newBarrelBottom->GetNumberOfPoints())
			maxZ = lround(newBarrelBottom->GetCell(newBarrelBottom->GetNumberOfCells()-1)->GetBounds()[4]);
		else
			maxZ = lround(barrel->GetBounds()[5]);
		
		for(int ii = 0; ii < barrel->GetNumberOfCells(); ++ii)
			zList.push_back(lround(barrel->GetCell(ii)->GetBounds()[4]));
		for(int ii = 0; ii < newBarrelTop->GetNumberOfCells(); ++ii)
			zList.push_back(lround(newBarrelTop->GetCell(ii)->GetBounds()[4]));
		for(int ii = 0; ii < newBarrelBottom->GetNumberOfCells(); ++ii)
			zList.push_back(lround(newBarrelBottom->GetCell(ii)->GetBounds()[4]));
		zList.sort();
		for(int ii = minZ; ii <= maxZ; ++ii)
		{
			if(std::find(zList.begin(), zList.end(), ii) == zList.end())
			{
				missingZList.push_back(ii);
// 				std::cout << "z = " << ii << " missing..." << std::endl;
			}
		}
		PolyDataPointerType missingBarrelPlanes = PolyDataPointerType::New();
		PointsPointerType missingBarrelPoints = PointsPointerType::New();
		missingBarrelPlanes->Allocate();
		missingBarrelPoints->SetDataTypeToFloat();
		if(missingZList.size())
		{
			std::list< int > topMissingPlanes;
			std::list< int > bottomMissingPlanes;
			std::list< int > centerMissingPlanes;
			std::list< int >::iterator missingZListIt;
			for(missingZListIt = missingZList.begin(); missingZListIt != missingZList.end(); ++missingZListIt)
			{
				if(newBarrelTop->GetNumberOfPoints() && *missingZListIt < newBarrelTop->GetBounds()[5])
					topMissingPlanes.push_back(*missingZListIt);
				else if(newBarrelBottom->GetNumberOfPoints() && *missingZListIt > newBarrelBottom->GetBounds()[4])
					bottomMissingPlanes.push_back(*missingZListIt);
				else
					centerMissingPlanes.push_back(*missingZListIt);
			}
			lastID = 0;
			if(topMissingPlanes.size())
			{
				std::list< int >::iterator topMissingIt;
				for(topMissingIt = topMissingPlanes.begin(); topMissingIt != topMissingPlanes.end(); ++topMissingIt)
				{
					int minDist = 1E06;
					int closestCell = 0;
					for(int ii = 0; ii < newBarrelTop->GetNumberOfCells(); ++ii)
					{
						int tmpDist = std::abs(newBarrelTop->GetCell(ii)->GetBounds()[4] - *topMissingIt);
						if(tmpDist < minDist)
						{
							minDist = tmpDist;
							closestCell = ii;
						}
					}
					if(minDist == 0)
						continue;
					
					PointsPointerType copyPts = newBarrelTop->GetCell(closestCell)->GetPoints();
					PolygonPointerType pastePoly = PolygonPointerType::New();
					pastePoly->GetPointIds()->SetNumberOfIds(copyPts->GetNumberOfPoints());
					for(int ii = 0; ii < copyPts->GetNumberOfPoints(); ++ii)
					{
						double * pt = new double[3];
						copyPts->GetPoint(ii, pt);
						pt[2] += round(*topMissingIt - pt[2]);
						missingBarrelPoints->InsertNextPoint(pt);
						pastePoly->GetPointIds()->SetId(ii, ii + lastID);
					}
					lastID += copyPts->GetNumberOfPoints();
					missingBarrelPlanes->InsertNextCell(pastePoly->GetCellType(), pastePoly->GetPointIds());
				}
			}
			if(bottomMissingPlanes.size())
			{
				std::list< int >::iterator bottomMissingIt;
				for(bottomMissingIt = bottomMissingPlanes.begin(); bottomMissingIt != bottomMissingPlanes.end(); ++bottomMissingIt)
				{
					int minDist = 1E06;
					int closestCell = 0;
					for(int ii = 0; ii < newBarrelBottom->GetNumberOfCells(); ++ii)
					{
						int tmpDist = std::abs(newBarrelBottom->GetCell(ii)->GetBounds()[4] - *bottomMissingIt);
						if(tmpDist < minDist)
						{
							minDist = tmpDist;
							closestCell = ii;
						}
					}
					if(minDist == 0)
						continue;
					
					PointsPointerType copyPts = newBarrelBottom->GetCell(closestCell)->GetPoints();
					PolygonPointerType pastePoly = PolygonPointerType::New();
					pastePoly->GetPointIds()->SetNumberOfIds(copyPts->GetNumberOfPoints());
					for(int ii = 0; ii < copyPts->GetNumberOfPoints(); ++ii)
					{
						double * pt = new double[3];
						copyPts->GetPoint(ii, pt);
						pt[2] += round(*bottomMissingIt - pt[2]);
						missingBarrelPoints->InsertNextPoint(pt);
						pastePoly->GetPointIds()->SetId(ii, ii + lastID);
					}
					lastID += copyPts->GetNumberOfPoints();
					missingBarrelPlanes->InsertNextCell(pastePoly->GetCellType(), pastePoly->GetPointIds());
				}
			}
			if(centerMissingPlanes.size())
			{
				std::list< int >::iterator centerMissingIt;
				for(centerMissingIt = centerMissingPlanes.begin(); centerMissingIt != centerMissingPlanes.end(); ++centerMissingIt)
				{
					int minDist = 1E06;
					int closestCell = 0;
					for(int ii = 0; ii < barrel->GetNumberOfCells(); ++ii)
					{
						int tmpDist = std::abs(barrel->GetCell(ii)->GetBounds()[4] - *centerMissingIt);
						if(tmpDist < minDist)
						{
							minDist = tmpDist;
							closestCell = ii;
						}
					}
					if(minDist == 0)
						continue;
					
					PointsPointerType copyPts = barrel->GetCell(closestCell)->GetPoints();
					PolygonPointerType pastePoly = PolygonPointerType::New();
					pastePoly->GetPointIds()->SetNumberOfIds(copyPts->GetNumberOfPoints());
					for(int ii = 0; ii < copyPts->GetNumberOfPoints(); ++ii)
					{
						double * pt = new double[3];
						copyPts->GetPoint(ii, pt);
						pt[2] += round(*centerMissingIt - pt[2]);
						missingBarrelPoints->InsertNextPoint(pt);
						pastePoly->GetPointIds()->SetId(ii, ii + lastID);
					}
					lastID += copyPts->GetNumberOfPoints();
					missingBarrelPlanes->InsertNextCell(pastePoly->GetCellType(), pastePoly->GetPointIds());
				}
			}
			missingBarrelPlanes->SetPoints(missingBarrelPoints);
			missingBarrelPlanes->Update();
		}
		/*** END MANUAL BARREL COMMENT ***/
// 		std::flush(std::cout << "appending bottom and top to existing barrel..." << std::endl);
		AppendPolyDataPointerType appendToBarrel = AppendPolyDataPointerType::New();
		appendToBarrel->AddInput(barrel);
		appendToBarrel->AddInput(newBarrelTop);
		appendToBarrel->AddInput(newBarrelBottom);
		/*** BEGIN MANUAL BARREL COMMENT ***/
		if(missingZList.size())
		{
// 			std::flush(std::cout << "setting missing planes" << std::endl);
			appendToBarrel->AddInput(missingBarrelPlanes);
		}
		/*** END MANUAL BARREL COMMENT ***/
		appendToBarrel->Update();
		barrel->DeepCopy(appendToBarrel->GetOutput());
		barrel->Update();
		
// // 		Vertex * newVert1 = new Vertex(extremeTopPt, Soma);
// // 		Vertex * newVert2 = new Vertex(extremeBottomPt, Soma);
// 		Vertex * newVert1 = new Vertex(axisTopCoord, Soma);
// 		Vertex * newVert2 = new Vertex(axisBottomCoord, Soma);
// 		spatialGraph2->addVertex(newVert1);
// 		spatialGraph2->addVertex(newVert2);
// 		int connectionIndex[2];
// 		if(!spatialGraph2->getNumberOfVertices())
// 		{
// 			connectionIndex[0] = 0;
// 			connectionIndex[1] = 1;
// 		}
// 		else
// 		{
// 			connectionIndex[0] = spatialGraph2->getNumberOfVertices() - 2;
// 			connectionIndex[1] = spatialGraph2->getNumberOfVertices() - 1;
// 		}
// 		std::list< double * > axisCoords;
// // 		axisCoords.push_back(extremeTopPt);
// 		axisCoords.push_back(axisTopCoord);
// 		axisCoords.push_back(finalTopPt);
// 		axisCoords.push_back(finalBottomPt);
// 		axisCoords.push_back(axisBottomCoord);
// // 		axisCoords.push_back(extremeBottomPt);
// 		int noOfAxisPoints = axisCoords.size();
// 		Edge * newAxis = new Edge(connectionIndex, noOfAxisPoints, Soma, axisCoords/*maybe score*/);
// 		spatialGraph2->addEdge(newAxis);
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel caps." << std::endl;
		return;
	}
};

void Geometry::calculateAvgContours(std::map< int, PolyDataPointerType > barrels, std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double * > >& endPointMap)
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
		{
			std::vector< double * > endPoints;
			closeBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPoints);
			endPointMap.insert(std::pair< int, std::vector< double * > >(*labelIt, endPoints));
			avgBarrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, smoothBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPointMap[*labelIt])));
		}
	// now, enforce non-overlap constraint
	std::map< int, std::list< int > > barrelGrid = createBarrelGrid(barrelAxes);
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
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
					// not accurate enough
					// use actual center of cell instead
// 					centerBottom[ii] = endPointMap[centerBarrel][1][ii];
					neighborBottom[ii] = endPointMap[neighborBarrel][1][ii];
					projectedNeighborBottom[ii] = endPointMap[neighborBarrel][1][ii];
					centerTop[ii] = endPointMap[centerBarrel][0][ii];
					neighborTop[ii] = endPointMap[neighborBarrel][0][ii];
				}
				int cSubID;
				double cPCoords[3], * cWeights = new double[avgBarrels[centerBarrel]->GetCell(1)->GetNumberOfPoints()];
				avgBarrels[centerBarrel]->GetCell(1)->GetParametricCenter(cPCoords);
				avgBarrels[centerBarrel]->GetCell(1)->EvaluateLocation(cSubID, cPCoords, centerBottom, cWeights);
				
				// barrel axes point upwards towards Pia by construction
				// project neighbor pts on bottom plane of center barrel
				normalize(centerAxis);
				normalize(neighborAxis);
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
					normalize(intersectVec);
					pNormal[0] = centerAxis[1]*intersectVec[2] - centerAxis[2]*intersectVec[1];
					pNormal[1] = centerAxis[2]*intersectVec[0] - centerAxis[0]*intersectVec[2];
					pNormal[2] = centerAxis[0]*intersectVec[1] - centerAxis[1]*intersectVec[0];
					normalize(pNormal);
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
							newDist = std::max(newDist - intersectPlane->DistanceToPlane(pt) - 5.0, 1.0);
							for(int jj = 0; jj < 3; ++jj)
								pt[jj] = pt[jj] - centerBottom[jj];
							normalize(pt);
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
							newDist = std::max(newDist - intersectPlane->DistanceToPlane(pt) - 5.0, 1.0);
							for(int jj = 0; jj < 3; ++jj)
								pt[jj] = pt[jj] - projectedNeighborBottom[jj];
							normalize(pt);
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
	
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		if(spatialGraph->isLabelInSpatialGraph(*labelIt))
			spatialGraph2->addPolyDataObject(avgBarrels[*labelIt], *labelIt);
};

/******************************************************************************/
/*Basically same method as vtkPolygon::IntersectConvex2DCells(), but for      */
/*coplanar polygons. Finds intersection by brute-force intersecting all       */
/*pairs of edges. Returns (int)nr of intersection points in {0,1,2}           */
/******************************************************************************/
int Geometry::intersectConvexCellsInPlane(CellPointerType cell1, CellPointerType cell2, double tol, double p0[3], double p1[3])
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
};

std::vector< double > Geometry::computeBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double* newAxis, double* barrelCenter, std::vector< double* > endPoints, int label, std::map< int, Column* >& barrelColumns, std::map< int, PolyDataPointerType >& avgBarrels)
{
	std::vector< double > params;
	if(barrel->GetNumberOfPoints())
	{
		double height, area, topPiaDist, bottomPiaDist, piaWMDist;
		double * topPt = new double[3];
		double * bottomPt = new double[3];
		topPt[0] = endPoints[0][0], topPt[1] = endPoints[0][1], topPt[2] = endPoints[0][2];
		bottomPt[0] = endPoints[1][0], bottomPt[1] = endPoints[1][1], bottomPt[2] = endPoints[1][2];
		height = sqrt((topPt[0] - bottomPt[0])*(topPt[0] - bottomPt[0]) + (topPt[1] - bottomPt[1])*(topPt[1] - bottomPt[1]) + (topPt[2] - bottomPt[2])*(topPt[2] - bottomPt[2]));
		
		PolyDataPointerType avgContours = avgBarrels[label];
		double normal[3];
		// don't use GetCell(0)->GetPoints b/c points in PointIDList are not sorted in strictly ascending order for overlap-corrected contours!
		area = vtkPolygon::ComputeArea(avgContours->GetPoints(), avgContours->GetCell(0)->GetNumberOfPoints(), avgContours->GetCell(0)->GetPointIds()->GetPointer(0), normal);
		
// 		PolyDataPointerType maxContours = maxBarrelContour(barrel, newAxis, barrelCenter, endPoints);
// 		spatialGraph2->addPolyDataObject(maxContours, label);
// 		double normal[3];
// 		area = vtkPolygon::ComputeArea(maxContours->GetCell(0)->GetPoints(), maxContours->GetCell(0)->GetNumberOfPoints(), maxContours->GetCell(0)->GetPointIds()->GetPointer(0), normal);
		
		double * piaIntersection = axisSurfaceIntersection(pia, newAxis, barrelCenter);
		bottomPiaDist = sqrt((piaIntersection[0] - bottomPt[0])*(piaIntersection[0] - bottomPt[0]) + (piaIntersection[1] - bottomPt[1])*(piaIntersection[1] - bottomPt[1]) + (piaIntersection[2] - bottomPt[2])*(piaIntersection[2] - bottomPt[2]));
		topPiaDist = sqrt((piaIntersection[0] - topPt[0])*(piaIntersection[0] - topPt[0]) + (piaIntersection[1] - topPt[1])*(piaIntersection[1] - topPt[1]) + (piaIntersection[2] - topPt[2])*(piaIntersection[2] - topPt[2]));
		
		double * wmIntersection = axisSurfaceIntersection(wm, newAxis, barrelCenter);
		piaWMDist = sqrt((piaIntersection[0] - wmIntersection[0])*(piaIntersection[0] - wmIntersection[0]) + (piaIntersection[1] - wmIntersection[1])*(piaIntersection[1] - wmIntersection[1]) + (piaIntersection[2] - wmIntersection[2])*(piaIntersection[2] - wmIntersection[2]));
		
		barrelColumns.insert(std::pair< int, Column * >(label, createBarrelColumn(avgContours, piaIntersection, wmIntersection)));
		
		//selected axis
		normalize(newAxis);
		double * topAxisPt = new double[3];
		double * bottomAxisPoint = new double[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			topAxisPt[ii] = piaIntersection[ii];
			bottomAxisPoint[ii] = barrelCenter[ii] - 1500*newAxis[ii];
		}
		Vertex * newVert1 = new Vertex(topAxisPt, ZAxis);
		Vertex * newVert2 = new Vertex(bottomAxisPoint, ZAxis);
		spatialGraph2->addVertex(newVert1);
		spatialGraph2->addVertex(newVert2);
		int connectionIndex[2];
		if(!spatialGraph2->getNumberOfVertices())
		{
			connectionIndex[0] = 0;
			connectionIndex[1] = 1;
		}
		else
		{
			connectionIndex[0] = spatialGraph2->getNumberOfVertices() - 2;
			connectionIndex[1] = spatialGraph2->getNumberOfVertices() - 1;
		}
		int noOfAxisPoints = 2;
		std::list< double * > axisCoords;
		axisCoords.push_back(topAxisPt);
		axisCoords.push_back(bottomAxisPoint);
		Edge * newEdge = new Edge(connectionIndex, noOfAxisPoints, ZAxis, axisCoords/*maybe score*/);
		spatialGraph2->addEdge(newEdge);
		
		params.push_back(bottomPiaDist);
		params.push_back(topPiaDist);
		params.push_back(height);
		params.push_back(piaWMDist);
		params.push_back(area);
		
		delete [] topPt, delete [] bottomPt, delete [] piaIntersection, delete [] wmIntersection;
		return params;
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel parameters." << std::endl;
		return params;
	}
};
std::vector< double > Geometry::computeManualBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double * newAxis, double * barrelCenter, std::vector< double * > endPoints)
{
	std::vector< double > params;
	if(barrel->GetNumberOfPoints())
	{
		double height, topPiaDist, bottomPiaDist, piaWMDist = 0;
		double * topPt = new double[3];
		double * bottomPt = new double[3];
		topPt[0] = endPoints[0][0], topPt[1] = endPoints[0][1], topPt[2] = endPoints[0][2];
		bottomPt[0] = endPoints[1][0], bottomPt[1] = endPoints[1][1], bottomPt[2] = endPoints[1][2];
		height = sqrt((topPt[0] - bottomPt[0])*(topPt[0] - bottomPt[0]) + (topPt[1] - bottomPt[1])*(topPt[1] - bottomPt[1]) + (topPt[2] - bottomPt[2])*(topPt[2] - bottomPt[2]));
		
		double * piaIntersection = axisSurfaceIntersection(pia, newAxis, barrelCenter);
		bottomPiaDist = sqrt((piaIntersection[0] - bottomPt[0])*(piaIntersection[0] - bottomPt[0]) + (piaIntersection[1] - bottomPt[1])*(piaIntersection[1] - bottomPt[1]) + (piaIntersection[2] - bottomPt[2])*(piaIntersection[2] - bottomPt[2]));
		topPiaDist = sqrt((piaIntersection[0] - topPt[0])*(piaIntersection[0] - topPt[0]) + (piaIntersection[1] - topPt[1])*(piaIntersection[1] - topPt[1]) + (piaIntersection[2] - topPt[2])*(piaIntersection[2] - topPt[2]));
		
		double * wmIntersection;
		if(wm)
		{
			wmIntersection = axisSurfaceIntersection(wm, newAxis, barrelCenter);
			if(wmIntersection)
				piaWMDist = sqrt((piaIntersection[0] - wmIntersection[0])*(piaIntersection[0] - wmIntersection[0]) + (piaIntersection[1] - wmIntersection[1])*(piaIntersection[1] - wmIntersection[1]) + (piaIntersection[2] - wmIntersection[2])*(piaIntersection[2] - wmIntersection[2]));
		}
		
		
		params.push_back(bottomPiaDist);
		params.push_back(topPiaDist);
		params.push_back(height);
		params.push_back(piaWMDist);
		
		delete [] topPt, delete [] bottomPt, delete [] piaIntersection;
		if(wm) delete [] wmIntersection;
		return params;
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel parameters." << std::endl;
		return params;
	}
};

Column * Geometry::createBarrelColumn(PolyDataPointerType avgBarrel, double * top, double * bottom)
{
	double axis[3];
	axis[0] = top[0] - bottom[0], axis[1] = top[1] - bottom[1], axis[2] = top[2] - bottom[2];
	//normalize(axis);
	std::vector< double * > pts;
	PointsPointerType cellPts = avgBarrel->GetCell(0)->GetPoints();
	for(int jj = 0; jj < cellPts->GetNumberOfPoints(); ++jj)
	{
		double t;
		double * closestPt = new double[3];
		double pt[3];
		cellPts->GetPoint(jj, pt);
		vtkLine::DistanceToLine(pt, top, bottom, t, closestPt);
		for(int kk = 0; kk < 3; ++kk)
			closestPt[kk] = pt[kk] - closestPt[kk];
		pts.push_back(closestPt);
	}
	PolyDataPointerType avgContours = PolyDataPointerType::New();
	PointsPointerType avgContourPoints = PointsPointerType::New();
	PolygonPointerType topPoly = PolygonPointerType::New();
	PolygonPointerType bottomPoly = PolygonPointerType::New();
	avgContours->Allocate();
	avgContourPoints->SetDataTypeToFloat();
	topPoly->GetPointIds()->SetNumberOfIds(cellPts->GetNumberOfPoints());
	bottomPoly->GetPointIds()->SetNumberOfIds(cellPts->GetNumberOfPoints());
	for(int ii = 0; ii < cellPts->GetNumberOfPoints(); ++ii)
	{
		double * topPt = new double[3];
		topPt[0] = top[0], topPt[1] = top[1], topPt[2] = top[2];
		for(int jj = 0; jj < 3; ++jj)
			topPt[jj] += pts[ii][jj];
		avgContourPoints->InsertNextPoint(topPt);
		topPoly->GetPointIds()->SetId(ii, ii);
	}
	for(int ii = 0; ii < cellPts->GetNumberOfPoints(); ++ii)
	{
		double * bottomPt = new double[3];
		bottomPt[0] = bottom[0], bottomPt[1] = bottom[1], bottomPt[2] = bottom[2];
		for(int jj = 0; jj < 3; ++jj)
			bottomPt[jj] += pts[ii][jj];
		avgContourPoints->InsertNextPoint(bottomPt);
		bottomPoly->GetPointIds()->SetId(ii, ii + cellPts->GetNumberOfPoints());
	}
	avgContours->InsertNextCell(topPoly->GetCellType(), topPoly->GetPointIds());
	avgContours->InsertNextCell(bottomPoly->GetCellType(), bottomPoly->GetPointIds());
	avgContours->SetPoints(avgContourPoints);
	avgContours->Update();

	Column * newColumn = new Column(avgContours, top, bottom);
	return newColumn;
};

std::map< int, std::vector< double > > Geometry::calculateColumnOverlap(std::map< int, Column * > barrelColumns, int label)
{
	std::map< int, std::list< int > > grid = createBarrelGrid(barrelColumns);
	Column * thisCol = barrelColumns[label];
	
	double colHeight = thisCol->getHeight();
	double zBinSize = 50;
	int nrOfZBins = (int)(colHeight/zBinSize + 1);
	std::map< int, std::vector< double > > colOverlapProfiles;
	
	double minDist, maxDist;
	std::vector< double * > contourPts;
	ImageDataPointerType columnVoxelGrid = createColumnVoxels(thisCol, contourPts, minDist, maxDist);
//	columnVoxelGrid->Print(std::cout);
	int * extent = columnVoxelGrid->GetExtent();
	int insidePx = 0;
	for(int z = extent[4]; z <= extent[5]; ++z)
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				double pt[3];
				pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
				double t, projectedPt[3];
				double dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, projectedPt);
				dist = sqrt(dist);
				if(dist > maxDist || t < 0 || t > 1)
				{
					continue;
				}
				if(dist < minDist)
				{
					unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
					*px = label;
					++insidePx;
					continue;
				}
				
				PolyDataPointerType polyData = PolyDataPointerType::New();
				polyData->Allocate(1);
				PointsPointerType points = PointsPointerType::New();
				points->SetDataTypeToFloat();
				PolygonPointerType poly = PolygonPointerType::New();
				poly->GetPointIds()->SetNumberOfIds(contourPts.size());
				for(int ii = 0; ii < contourPts.size(); ++ii)
				{
					double tmp[3];
					for(int jj = 0; jj < 3; ++jj)
						tmp[jj] = projectedPt[jj] + contourPts[ii][jj];
					points->InsertNextPoint(tmp);
					poly->GetPointIds()->SetId(ii, ii);
				}
				polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
				polyData->SetPoints(points);
				polyData->Update();
				double closestPoint[3], pCoords[3], dist2;
				int subId;
				double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
				int insidePolygon = polyData->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pCoords, dist2, weights);
				if(insidePolygon == 1)
				{
					unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
					*px = label;
					++insidePx;
				}
			}
// 	std::cout << "nr. inside px = " << insidePx << std::endl;
// 	std::cout << "this is " << (double)insidePx/(double)((extent[1] - extent[0])*(extent[3] - extent[2])*(extent[5] - extent[4]))*100 << "% of the column bounding box" << std::endl;
	columnVoxelGrid->Update();
	
	int overlapPx = 0;
	std::list< int >::const_iterator neighborColIt;
	for(neighborColIt = grid[label].begin(); neighborColIt != grid[label].end(); ++neighborColIt)
	{
		double neighborMinDist, neighborMaxDist;
		std::vector< double * > neighborContourPts;
		ImageDataPointerType neighborColumnVoxelGrid = createColumnVoxels(barrelColumns[*neighborColIt], neighborContourPts, neighborMinDist, neighborMaxDist);
		int * neighborExtent = neighborColumnVoxelGrid->GetExtent();
		int insidePx = 0;
		for(int z = neighborExtent[4]; z <= neighborExtent[5]; ++z)
			for(int y = neighborExtent[2]; y <= neighborExtent[3]; ++y)
				for(int x = neighborExtent[0]; x <= neighborExtent[1]; ++x)
				{
					if(x < extent[0] || x > extent[1] || y < extent[2] || y > extent[3] || z < extent[4] || z > extent[5])
						continue;
					
					unsigned char * currPx = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
					if(!*currPx)
						continue;
					
					double pt[3];
					pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
					double t, projectedPt[3];
					double dist = vtkLine::DistanceToLine(pt, barrelColumns[*neighborColIt]->top, barrelColumns[*neighborColIt]->bottom, t, projectedPt);
					dist = sqrt(dist);
					if(dist > neighborMaxDist || t < 0 || t > 1)
					{
						continue;
					}
					if(dist < neighborMinDist)
					{
						unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
						*px = *neighborColIt;
						++overlapPx;
						continue;
					}
					
					PolyDataPointerType polyData = PolyDataPointerType::New();
					polyData->Allocate(1);
					PointsPointerType points = PointsPointerType::New();
					points->SetDataTypeToFloat();
					PolygonPointerType poly = PolygonPointerType::New();
					poly->GetPointIds()->SetNumberOfIds(neighborContourPts.size());
					for(int ii = 0; ii < neighborContourPts.size(); ++ii)
					{
						double tmp[3];
						for(int jj = 0; jj < 3; ++jj)
							tmp[jj] = projectedPt[jj] + neighborContourPts[ii][jj];
						points->InsertNextPoint(tmp);
						poly->GetPointIds()->SetId(ii, ii);
					}
					polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
					polyData->SetPoints(points);
					polyData->Update();
					double closestPoint[3], pCoords[3], dist2;
					int subId;
					double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
					int insidePolygon = polyData->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pCoords, dist2, weights);
					if(insidePolygon == 1)
					{
						unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
						*px = *neighborColIt;
						++overlapPx;
					}
				}
		//calculate overlap profile for this specific neighbor column
		unsigned int * totalPxProfile = new unsigned int[nrOfZBins];
		unsigned int * sharedPxProfile = new unsigned int[nrOfZBins];
		for(int ii = 0; ii < nrOfZBins; ++ii)
		{
			totalPxProfile[ii] = 0;
			sharedPxProfile[ii] = 0;
		}
		for(int z = extent[4]; z <= extent[5]; ++z)
			for(int y = extent[2]; y <= extent[3]; ++y)
				for(int x = extent[0]; x <= extent[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
					if(*px)
					{
						double * pt = new double[3];
						pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
						double t, projectedPt[3];
						double dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, projectedPt);
						if(t >= 0 && t <= 1)
						{
							int currZBin = (int)(t*colHeight/zBinSize);
							++totalPxProfile[currZBin];
							if(*px == *neighborColIt)
								++sharedPxProfile[currZBin];
						}
					}
				}
		std::vector< double > overlapProfile;
		for(int ii = 0; ii < nrOfZBins; ++ii)
		{
			double overlapRatio = (double)sharedPxProfile[ii]/(double)totalPxProfile[ii];
			overlapProfile.push_back(overlapRatio);
		}
		colOverlapProfiles.insert(std::pair< int, std::vector< double > >(*neighborColIt, overlapProfile));
	}
// 	std::cout << "nr. overlap px = " << overlapPx << std::endl;
// 	std::cout << "this is " << (double)overlapPx/(double)insidePx*100 << "% of the column volume" << std::endl;
	
	std::vector< double > totalOverlapProfile;
	unsigned int * totalPxProfile = new unsigned int[nrOfZBins];
	unsigned int * sharedPxProfile = new unsigned int[nrOfZBins];
	for(int ii = 0; ii < nrOfZBins; ++ii)
	{
		totalPxProfile[ii] = 0;
		sharedPxProfile[ii] = 0;
	}
	for(int z = extent[4]; z <= extent[5]; ++z)
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				unsigned char * px = static_cast< unsigned char * >(columnVoxelGrid->GetScalarPointer(x, y, z));
				if(*px)
				{
					double * pt = new double[3];
					pt[0] = 10*x, pt[1] = 10*y, pt[2] = 10*z;
					double t, projectedPt[3];
					double dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, projectedPt);
					if(t >= 0 && t <= 1)
					{
						int currZBin = (int)(t*colHeight/zBinSize);
						++totalPxProfile[currZBin];
						if(*px != label)
							++sharedPxProfile[currZBin];
					}
				}
			}
	for(int ii = 0; ii < nrOfZBins; ++ii)
	{
		double overlapRatio = (double)sharedPxProfile[ii]/(double)totalPxProfile[ii];
		totalOverlapProfile.push_back(overlapRatio);
	}
	colOverlapProfiles.insert(std::pair< int, std::vector< double > >(label, totalOverlapProfile));
	
	return colOverlapProfiles;
};

ImageDataPointerType Geometry::createColumnVoxels(Column * column, std::vector< double * >& contour, double& minDist, double& maxDist)
{
	minDist = 1E06, maxDist = 0;
	for(int ii = 0; ii < column->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
	{
		double pt[3];
		double dist = 0, t;
		double * closestPt = new double[3];
		column->contours->GetCell(0)->GetPoints()->GetPoint(ii, pt);
		dist = vtkLine::DistanceToLine(pt, column->top, column->bottom, t, closestPt);
		for(int jj = 0; jj < 3; ++jj)
		{
			closestPt[jj] = pt[jj] - closestPt[jj];
		}
		dist = sqrt(dist);
		if(dist < minDist)
			minDist = dist;
		if(dist > maxDist)
			maxDist = dist;
		contour.push_back(closestPt);
	}
	//	std::cout << "minDist = " << minDist << std::endl;
	//	std::cout << "maxDist = " << maxDist << std::endl;
	double * bounds = column->contours->GetBounds();
	int * intBounds = new int[6];
	for(int ii = 0; ii < 6; ++ii)
		intBounds[ii] = lround(bounds[ii]);
	return createImageVolume(Barrel, intBounds[0], intBounds[1], intBounds[2], intBounds[3], intBounds[4], intBounds[5]);
};

std::map< int, std::vector< double > > Geometry::calculateSeptalDistances(std::map< int, Column * > barrelColumns, int label)
{
	std::map< int, std::list< int > > grid = createBarrelGrid(barrelColumns);
	Column * thisCol = barrelColumns[label];
	
	double colHeight = thisCol->getHeight();
	double zBinSize = 50;
	int nrOfZBins = (int)(colHeight/zBinSize + 1);
	std::map< int, std::vector< double > > septalDistProfiles;
	
	std::vector< double * > radialContour;
	PointsPointerType topPts = thisCol->contours->GetCell(0)->GetPoints();
	for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
	{
		double * tmpPt = new double[3];
		topPts->GetPoint(ii, tmpPt);
		for(int jj = 0; jj < 3; ++jj)
			tmpPt[jj] = tmpPt[jj] - thisCol->top[jj];
		radialContour.push_back(tmpPt);
	}
	
	std::list< int >::const_iterator neighborColIt;
	for(neighborColIt = grid[label].begin(); neighborColIt != grid[label].end(); ++neighborColIt)
	{
		std::vector< double * > neighborRadialContour;
		PointsPointerType neighborTopPts = barrelColumns[*neighborColIt]->contours->GetCell(0)->GetPoints();
		for(int ii = 0; ii < neighborTopPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			neighborTopPts->GetPoint(ii, tmpPt);
			for(int jj = 0; jj < 3; ++jj)
				tmpPt[jj] = tmpPt[jj] - barrelColumns[*neighborColIt]->top[jj];
			neighborRadialContour.push_back(tmpPt);
		}
		
		double closestLine[3];
		for(int ii = 0; ii < 3; ++ii)
			closestLine[ii] = barrelColumns[*neighborColIt]->top[ii] - thisCol->top[ii];
		normalize(closestLine);
		
		int thisClosestIndex, neighborClosestIndex;
		double max = 0;
		for(int ii = 0; ii < radialContour.size(); ++ii)
		{
			double tmp = 0, tmpNorm = 0;
			for(int jj = 0; jj < 3; ++jj)
			{
				tmp += radialContour[ii][jj]*closestLine[jj];
				tmpNorm += radialContour[ii][jj]*radialContour[ii][jj];
			}
			tmpNorm = sqrt(tmpNorm);
			if(tmpNorm)
				tmp /= tmpNorm;
			if(tmp > max)
			{
				max = tmp;
				thisClosestIndex = ii;
			}
		}
		max = 0;
		for(int ii = 0; ii < neighborRadialContour.size(); ++ii)
		{
			double tmp = 0, tmpNorm = 0;
			for(int jj = 0; jj < 3; ++jj)
			{
				tmp += -neighborRadialContour[ii][jj]*closestLine[jj];
				tmpNorm += neighborRadialContour[ii][jj]*neighborRadialContour[ii][jj];
			}
			tmpNorm = sqrt(tmpNorm);
			if(tmpNorm)
				tmp /= tmpNorm;
			if(tmp > max)
			{
				max = tmp;
				neighborClosestIndex = ii;
			}
		}
		
		double thisRadius = 0, neighborRadius = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			thisRadius += radialContour[thisClosestIndex][ii]*radialContour[thisClosestIndex][ii];
			neighborRadius += neighborRadialContour[neighborClosestIndex][ii]*neighborRadialContour[neighborClosestIndex][ii];
		}
		thisRadius = sqrt(thisRadius);
		neighborRadius = sqrt(neighborRadius);
		
		double thisZAxis[3], neighborZAxis[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			thisZAxis[ii] = thisCol->bottom[ii] - thisCol->top[ii];
			neighborZAxis[ii] = barrelColumns[*neighborColIt]->bottom[ii] - barrelColumns[*neighborColIt]->top[ii];
		}
		normalize(thisZAxis);
		normalize(neighborZAxis);
		
		//calculate part of radial contour vector parallel to line connecting the two columns
		//by subtracting the part perpendicular to the plane spanned by the connecting line
		//and the respective z axis
		double thisPerpRadius[3], neighborPerpRadius[3];
		double thisNormal[3], neighborNormal[3];
		thisNormal[0] = closestLine[1]*thisZAxis[2] - closestLine[2]*thisZAxis[1];
		thisNormal[1] = closestLine[2]*thisZAxis[0] - closestLine[0]*thisZAxis[2];
		thisNormal[2] = closestLine[0]*thisZAxis[1] - closestLine[1]*thisZAxis[0];
		neighborNormal[0] = closestLine[1]*neighborZAxis[2] - closestLine[2]*neighborZAxis[1];
		neighborNormal[1] = closestLine[2]*neighborZAxis[0] - closestLine[0]*neighborZAxis[2];
		neighborNormal[2] = closestLine[0]*neighborZAxis[1] - closestLine[1]*neighborZAxis[0];
		normalize(thisNormal);
		normalize(neighborNormal);
		double thisParallelRadius = 0, neighborParallelRadius = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			thisParallelRadius += radialContour[thisClosestIndex][ii]*thisNormal[ii];
			neighborParallelRadius += neighborRadialContour[neighborClosestIndex][ii]*neighborNormal[ii];
		}
		double thisParallelRadiusVec[3], neighborParallelRadiusVec[3];
		double thisCorrectionAngle = 0, neighborCorrectionAngle = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			thisParallelRadiusVec[ii] = radialContour[thisClosestIndex][ii] - thisParallelRadius*thisNormal[ii];
			neighborParallelRadiusVec[ii] = neighborRadialContour[neighborClosestIndex][ii] - neighborParallelRadius*neighborNormal[ii];
		}
		normalize(thisParallelRadiusVec);
		normalize(neighborParallelRadiusVec);
		for(int ii = 0; ii < 3; ++ii)
		{
			thisCorrectionAngle += thisParallelRadiusVec[ii]*closestLine[ii];
			neighborCorrectionAngle += neighborParallelRadiusVec[ii]*closestLine[ii];
		}
		thisCorrectionAngle = std::abs(thisCorrectionAngle);
		neighborCorrectionAngle = std::abs(neighborCorrectionAngle);
		if(thisCorrectionAngle)
			thisRadius /= thisCorrectionAngle;
		if(neighborCorrectionAngle)
			neighborRadius /= neighborCorrectionAngle;
		
		std::vector< double > zDistProfile;
		for(int ii = 0; ii < nrOfZBins; ++ii)
		{
			double thisCurrPt[3], neighborCurrPt[3];
			double dist = 0;
			for(int jj = 0; jj < 3; ++jj)
			{
				thisCurrPt[jj] = thisCol->top[jj] + ii*zBinSize*thisZAxis[jj];
				neighborCurrPt[jj] = barrelColumns[*neighborColIt]->top[jj] + ii*zBinSize*neighborZAxis[jj];
				dist += (neighborCurrPt[jj] - thisCurrPt[jj])*(neighborCurrPt[jj] - thisCurrPt[jj]);
			}
			dist = sqrt(dist);
			dist = dist - thisRadius - neighborRadius;
			zDistProfile.push_back(std::max(0.0, dist));
		}
		septalDistProfiles.insert(std::pair< int, std::vector< double > >(*neighborColIt, zDistProfile));
	}
	
	return septalDistProfiles;
};

// calculate avg septal distance between adjacent barrels
// avgBarrels: top = cell 0, bottom = cell 1
void Geometry::septalDistances(std::map< int, PolyDataPointerType > avgBarrels, std::map< int, double* > barrelCenters)
{
	std::map< int, std::list< int > > barrelGrid = createNearestNeighborBarrelGrid(barrelCenters);
	
	std::list< int > barrelIDList;
	std::map< int, std::list< int > >::const_iterator readBarrelIDsIt;
	for(readBarrelIDsIt = barrelGrid.begin(); readBarrelIDsIt != barrelGrid.end(); ++readBarrelIDsIt)
		barrelIDList.push_back(readBarrelIDsIt->first);
	
	std::list< int >::const_iterator barrelIDIt;
	for(barrelIDIt = barrelIDList.begin(); barrelIDIt != barrelIDList.end(); ++barrelIDIt)
	{
		int thisBarrel = *barrelIDIt;
		if(avgBarrels[thisBarrel]->GetNumberOfCells())
		{
			std::list< int >::const_iterator neighborListIt;
			for(neighborListIt = barrelGrid[thisBarrel].begin(); neighborListIt != barrelGrid[thisBarrel].end(); ++neighborListIt)
			{
				int currNeighbor = *neighborListIt;
				double centerVec[] = {0, 0, 0};
				double normCenterVec[3];
				double centerDist = 0;
				for(int ii = 0; ii < 3; ++ii)
				{
					centerVec[ii] = barrelCenters[currNeighbor][ii] - barrelCenters[thisBarrel][ii];
					centerDist += centerVec[ii]*centerVec[ii];
					normCenterVec[ii] = centerVec[ii];
				}
				centerDist = sqrt(centerDist);
				normalize(normCenterVec);
				//calculate avg distance @top and @bottom
				//use fact that only top & bottom are stored in avg contour
				for(int ii = 0; ii < 2; ++ii)
				{
					int subId;
					double pCoords[3], thisCenterPt[3], * weights = new double[avgBarrels[thisBarrel]->GetCell(ii)->GetNumberOfPoints()];
					avgBarrels[thisBarrel]->GetCell(ii)->GetParametricCenter(pCoords);
					avgBarrels[thisBarrel]->GetCell(ii)->EvaluateLocation(subId, pCoords, thisCenterPt, weights);
					double pCoords2[3], currCenterPt[3], * weights2 = new double[avgBarrels[currNeighbor]->GetCell(ii)->GetNumberOfPoints()];
					avgBarrels[currNeighbor]->GetCell(ii)->GetParametricCenter(pCoords2);
					avgBarrels[currNeighbor]->GetCell(ii)->EvaluateLocation(subId, pCoords2, currCenterPt, weights2);
					//collect IDs of all pts facing the septum
					std::list< int > thisBarrelFacePts;
					std::list< int > currBarrelFacePts;
					PointsPointerType thisBarrelPts = avgBarrels[thisBarrel]->GetCell(ii)->GetPoints();
					PointsPointerType currBarrelPts = avgBarrels[currNeighbor]->GetCell(ii)->GetPoints();
					for(int jj = 0; jj < thisBarrelPts->GetNumberOfPoints(); ++jj)
					{
						double tmpPt[3];
						thisBarrelPts->GetPoint(jj, tmpPt);
						for(int kk = 0; kk < 3; ++kk)
							tmpPt[kk] = tmpPt[kk] - thisCenterPt[kk];
						normalize(tmpPt);
						double angle = 0;
						for(int kk = 0; kk < 3; ++kk)
							angle += tmpPt[kk]*normCenterVec[kk];
						if(acos(angle) < 30)
							thisBarrelFacePts.push_back(jj);
					}
					for(int jj = 0; jj < currBarrelPts->GetNumberOfPoints(); ++jj)
					{
						double tmpPt[3];
						currBarrelPts->GetPoint(jj, tmpPt);
						for(int kk = 0; kk < 3; ++kk)
							tmpPt[kk] = tmpPt[kk] - currCenterPt[kk];
						normalize(tmpPt);
						double angle = 0;
						for(int kk = 0; kk < 3; ++kk)
							angle += tmpPt[kk]*normCenterVec[kk];
						if(acos(angle) < 30)
							currBarrelFacePts.push_back(jj);
					}
					//calculate rough avg face
					//use fact that pts are ordered in contour
					double thisFaceAvg[3], neighborFaceAvg[3];
					for(int jj = 0; jj < 3; ++jj)
					{
						thisFaceAvg[jj] = normCenterVec[jj];
						neighborFaceAvg[jj] = -normCenterVec[jj];
					}
					std::list< int >::const_iterator thisFacePtsIt;
					std::list< int >::const_iterator currFacePtsIt;
					double thisFaceAvgDist = 0, currFaceAvgDist = 0;
					for(thisFacePtsIt = thisBarrelFacePts.begin(); thisFacePtsIt != thisBarrelFacePts.end(); ++thisFacePtsIt)
					{
						double tmpPt[3];
						thisBarrelPts->GetPoint(*thisFacePtsIt, tmpPt);
						double tmpProd = 0;
						for(int jj = 0; jj < 3; ++jj)
						{
							tmpPt[jj] = tmpPt[jj] - thisCenterPt[jj];
							tmpProd += tmpPt[jj]*normCenterVec[jj];
						}
						tmpProd = std::abs(tmpProd);
						thisFaceAvgDist += tmpProd;
					}
					if(thisBarrelFacePts.size())
						thisFaceAvgDist = thisFaceAvgDist/(double)thisBarrelFacePts.size();
					for(currFacePtsIt = currBarrelFacePts.begin(); currFacePtsIt != currBarrelFacePts.end(); ++currFacePtsIt)
					{
						double tmpPt[3];
						currBarrelPts->GetPoint(*currFacePtsIt, tmpPt);
						double tmpProd = 0;
						for(int jj = 0; jj < 3; ++jj)
						{
							tmpPt[jj] = tmpPt[jj] - currCenterPt[jj];
							tmpProd += tmpPt[jj]*normCenterVec[jj];
						}
						tmpProd = std::abs(tmpProd);
						currFaceAvgDist += tmpProd;
					}
					if(currBarrelFacePts.size())
						currFaceAvgDist = currFaceAvgDist/(double)currBarrelFacePts.size();
					//store endpoints of avg faces
					std::list< double * > thisFace;
					std::list< double * > currFace;
					double * thisPt1 = new double[3], * thisPt2 = new double[3];
					double * currPt1 = new double [3], * currPt2 = new double[3];
					double thisPt1Parallel = 0, thisPt2Parallel = 0;
					double currPt1Parallel = 0, currPt2Parallel = 0;
					double thisFaceListPt1[3], thisFaceListPt2[3];
					double currFaceListPt1[3], currFaceListPt2[3];
					thisBarrelPts->GetPoint(thisBarrelFacePts.front(), thisFaceListPt1);
					thisBarrelPts->GetPoint(thisBarrelFacePts.back(), thisFaceListPt2);
					currBarrelPts->GetPoint(currBarrelFacePts.front(), currFaceListPt1);
					currBarrelPts->GetPoint(currBarrelFacePts.back(), currFaceListPt2);
					for(int jj = 0; jj < 3; ++jj)
					{
						thisFaceListPt1[jj] = thisFaceListPt1[jj] - thisCenterPt[jj];
						thisFaceListPt2[jj] = thisFaceListPt2[jj] - thisCenterPt[jj];
						currFaceListPt1[jj] = currFaceListPt1[jj] - currCenterPt[jj];
						currFaceListPt2[jj] = currFaceListPt2[jj] - currCenterPt[jj];
						thisPt1Parallel += thisFaceListPt1[jj]*normCenterVec[jj];
						thisPt2Parallel += thisFaceListPt2[jj]*normCenterVec[jj];
						currPt1Parallel += -currFaceListPt1[jj]*normCenterVec[jj];
						currPt2Parallel += -currFaceListPt2[jj]*normCenterVec[jj];
					}
					//compute actual line as avg distance line between endpoints of face list
					for(int jj = 0; jj < 3; ++jj)
					{
						thisPt1[jj] = thisCenterPt[jj] + (thisFaceAvgDist - thisPt1Parallel)*normCenterVec[jj] + thisFaceListPt1[jj];
						thisPt2[jj] = thisCenterPt[jj] + (thisFaceAvgDist - thisPt2Parallel)*normCenterVec[jj] + thisFaceListPt2[jj];
						currPt1[jj] = currCenterPt[jj] + (currPt1Parallel - currFaceAvgDist)*normCenterVec[jj] + currFaceListPt1[jj];
						currPt2[jj] = currCenterPt[jj] + (currPt2Parallel - currFaceAvgDist)*normCenterVec[jj] + currFaceListPt2[jj];
					}
					thisFace.push_back(thisPt1), thisFace.push_back(thisPt2);
					currFace.push_back(currPt1), currFace.push_back(currPt2);
					
				}
			}
		}
		else
			std::cout << "Warning! Average barrel " << int2Labels[thisBarrel] << " is empty; should contain avg contours!" << std::endl;
	}
};

// assign all vessels that are uniquely contained in one barrel
// to that barrel and all other vessels to all barrels closer than
// a certain max distance (~250-350 mu?)
std::map< int, std::list< int > > Geometry::computeBarrelVesselCorrelations(std::map< int, PolyDataPointerType > avgBarrels, std::map< int, double * > barrelCenters, std::vector< Edge * > vesselVec)
{
	std::cout << "maxDist = ";
	double maxDist;
	std::cin >> maxDist;
// 	double maxDist = 300;
	
// 	std::vector< Edge * > vesselVec;
// 	std::vector< Edge * >::const_iterator sgVesselIt;
// 	for(sgVesselIt = spatialGraph->edgesBegin(); sgVesselIt != spatialGraph->edgesEnd(); ++sgVesselIt)
// 		if((*sgVesselIt)->label == Vessel && (*sgVesselIt)->edgePointCoordinates.size() == 2)
// 			vesselVec.push_back(*sgVesselIt);
	
	std::map< int, std::list< int > > barrelVesselMap;
	std::list< int > uniqueVesselList;
	std::list< int >::const_iterator barrelIt;
	// assign unique vessels
	for(barrelIt = barrelLabels.begin(); barrelIt != barrelLabels.end(); ++barrelIt)
	{
		if(spatialGraph->isLabelInSpatialGraph(*barrelIt))
		{
			std::list< int > vesselIDList;
			if(avgBarrels[*barrelIt]->GetNumberOfCells())
			{
				CellLocatorPointerType locator = CellLocatorPointerType::New();
				locator->AutomaticOn();
				locator->SetDataSet(avgBarrels[*barrelIt]);
				locator->BuildLocator();
				double * a0, * a1, intersectPt[3], tol = 0.1, t, pcoords[3];
				int subId;
				vtkIdType cellID;
				GenericCellPointerType intersectCell = GenericCellPointerType::New();
				for(int ii = 0; ii < vesselVec.size(); ++ii)
				{
					a0 = vesselVec[ii]->edgePointCoordinates.front(), a1 = vesselVec[ii]->edgePointCoordinates.back();
					if(locator->IntersectWithLine(a0, a1, tol, t, intersectPt, pcoords, subId, cellID, intersectCell))
					{
						vesselIDList.push_back(ii);
						uniqueVesselList.push_back(ii);
						vesselVec[ii]->label = *barrelIt;
						spatialGraph2->addEdge(vesselVec[ii]);
					}
				}
			}
			else
				std::cout << "Warning! Average barrel " << int2Labels[*barrelIt] << " is empty; should contain avg contours!" << std::endl;
			
			barrelVesselMap.insert(std::pair< int, std::list< int > >(*barrelIt, vesselIDList));
		}
	}
	// now assign remaining vessels to all neighbors
	for(barrelIt = barrelLabels.begin(); barrelIt != barrelLabels.end(); ++barrelIt)
	{
		if(spatialGraph->isLabelInSpatialGraph(*barrelIt))
		{
			if(avgBarrels[*barrelIt]->GetNumberOfCells())
			{
				for(int ii = 0; ii < vesselVec.size(); ++ii)
				{
					if(std::find(uniqueVesselList.begin(), uniqueVesselList.end(), ii) != uniqueVesselList.end())
						continue;
					
					double t, closestPoint[3];
					double dist = vtkLine::DistanceToLine(barrelCenters[*barrelIt], vesselVec[ii]->edgePointCoordinates.front(), vesselVec[ii]->edgePointCoordinates.back(), t, closestPoint);
					dist = sqrt(dist);
					if(t >= 0 && t <= 1 && dist <= maxDist)
					{
						barrelVesselMap[*barrelIt].push_back(ii);
						Edge * newVessel = new Edge(vesselVec[ii]);
						newVessel->label = *barrelIt;
						spatialGraph2->addEdge(newVessel);
					}
				}
			}
			else
				std::cout << "Warning! Average barrel " << int2Labels[*barrelIt] << " is empty; should contain avg contours!" << std::endl;
		}
	}
	return barrelVesselMap;
};

PolyDataPointerType Geometry::maxBarrelContour(PolyDataPointerType barrel, double * newAxis, double * barrelCenter, std::vector< double * > endPoints)
{
	if(barrel->GetNumberOfPoints())
	{
		if(barrel->GetCell(0)->GetNumberOfPoints() != 36)
		{
			std::cout << "Warning! Barrel has not been sampled correctly! max barrel contour may be corrupted..." << std::endl;
		}
		double * zAxis = new double[3];
		zAxis[0] = newAxis[0], zAxis[1] = newAxis[1], zAxis[2] = newAxis[2];
		normalize(zAxis);
		double axisPt1[3], axisPt2[3];
		axisPt1[0] = barrelCenter[0] + 1000*zAxis[0], axisPt1[1] = barrelCenter[1] + 1000*zAxis[1], axisPt1[2] = barrelCenter[2] + 1000*zAxis[2];
		axisPt2[0] = barrelCenter[0] - 1000*zAxis[0], axisPt2[1] = barrelCenter[1] - 1000*zAxis[1], axisPt2[2] = barrelCenter[2] - 1000*zAxis[2];
		std::vector< std::map< double, double * > > maxPoints;
		for(int ii = 0; ii < 36; ++ii)
		{
			std::map< double, double * > ptMap;
			ptMap.insert(std::pair< double, double * >(0, new double[3]));
			maxPoints.push_back(ptMap);
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
				maxPoints[jj].insert(std::pair< double, double * >(dist, closestPt));
			}
		}
		
// 		for(int ii = 0; ii < 36; ++ii)
// 			std::flush(std::cout << "maxPoints[" << ii << "].size() = " << maxPoints[ii].size() << std::endl);
		// create top & bottom max contours @ extreme points 
		PolyDataPointerType maxContours = PolyDataPointerType::New();
		PointsPointerType maxContourPoints = PointsPointerType::New();
		PolygonPointerType topPoly = PolygonPointerType::New();
		PolygonPointerType bottomPoly = PolygonPointerType::New();
		maxContours->Allocate();
		maxContourPoints->SetDataTypeToFloat();
		topPoly->GetPointIds()->SetNumberOfIds(36);
		bottomPoly->GetPointIds()->SetNumberOfIds(36);
		for(int ii = 0; ii < 36; ++ii)
		{
			if(maxPoints[ii].size())
			{
				double * topPt = new double[3];
				topPt[0] = endPoints[0][0], topPt[1] = endPoints[0][1], topPt[2] = endPoints[0][2];
				for(int jj = 0; jj < 3; ++jj)
					topPt[jj] += maxPoints[ii].rbegin()->second[jj];
				maxContourPoints->InsertNextPoint(topPt);
				topPoly->GetPointIds()->SetId(ii, ii);
			}
		}
		for(int ii = 0; ii < 36; ++ii)
		{
			if(maxPoints[ii].size())
			{
				double * bottomPt = new double[3];
				bottomPt[0] = endPoints[1][0], bottomPt[1] = endPoints[1][1], bottomPt[2] = endPoints[1][2];
				for(int jj = 0; jj < 3; ++jj)
					bottomPt[jj] += maxPoints[ii].rbegin()->second[jj];
				maxContourPoints->InsertNextPoint(bottomPt);
				bottomPoly->GetPointIds()->SetId(ii, ii + 36);
			}
		}
		maxContours->InsertNextCell(topPoly->GetCellType(), topPoly->GetPointIds());
		maxContours->InsertNextCell(bottomPoly->GetCellType(), bottomPoly->GetPointIds());
		maxContours->SetPoints(maxContourPoints);
		maxContours->Update();
		
		delete [] zAxis;
		return maxContours;
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate max barrel contour." << std::endl;
		return NULL;
	}
};

double * Geometry::axisSurfaceIntersection(PolyDataPointerType surface, double * axis, double * center)
{
	PolyDataPointerType surfRoi = selectSurfaceRoi(surface, center, 2000);
	if(surfRoi->GetNumberOfPoints())
	{
		CellLocatorPointerType locator = CellLocatorPointerType::New();
		locator->AutomaticOn();
		locator->SetDataSet(surfRoi);
		locator->BuildLocator();
		double * intersectPt = new double[3];
		double a0[3], a1[3], tol = 0.1, t, pcoords[3];
		int subId;
		vtkIdType cellID;
		GenericCellPointerType intersectCell = GenericCellPointerType::New();
		for(int jj = 0; jj < 3; ++jj)
		{
			a0[jj] = center[jj];
			a1[jj] = center[jj];
		}
		double * direction = new double[3];
		for(int jj = 0; jj < 3; ++jj)
			direction[jj] = axis[jj];
		normalize(direction);
		for(int jj = 0; jj < 3; ++jj)
		{
			a0[jj] += direction[jj]*2000;
			a1[jj] -= direction[jj]*2000;
		}
		delete [] direction;
// 		std::cout << "Vessel = [" << a0[0] << "," << a0[1] << "," << a0[2] << "] - [" << a1[0] << "," << a1[1] << "," << a1[2] << "]" << std::endl;
		int intersection = locator->IntersectWithLine(a0, a1, tol, t, intersectPt, pcoords, subId, cellID, intersectCell);
// 		std::cout << "intersection = " << intersection << std::endl;
		if(intersection)
		{
// 			std::cout << "Intersection found!" << std::endl;
			return intersectPt;
		}
		else
		{
			std::cout << "Warning! Could not find intersection point of axis with surface!" << std::endl;
			intersectPt[0] = 0, intersectPt[1] = 0, intersectPt[2] = 0;
			return intersectPt;
		}
	}
	else
	{
		std::flush(std::cout << "Error! Could not select surface ROI for calculation of intersection point." << std::endl);
		return NULL;
	}
};

void Geometry::writeBestBarrelAxes(std::multimap< double, double * > axes, double * barrelCentroid)
{
	if(axes.size())
	{
		std::cout << "Writing 5 highest scoring axes out of " << axes.size() << std::endl;
		std::multimap< double, double * >::const_reverse_iterator axesIt;
		int count = 0;
		for(axesIt = axes.rbegin(); axesIt != axes.rend() && count < 5; ++axesIt, ++count)
		{
			double * endPoint = new double[3];
			for(int ii = 0; ii < 3; ++ii)
				endPoint[ii] = barrelCentroid[ii] + axesIt->second[ii];
			Vertex * newVert1 = new Vertex(endPoint, C1);
			Vertex * newVert2 = new Vertex(barrelCentroid, C1);
			spatialGraph->addVertex(newVert1);
			spatialGraph->addVertex(newVert2);
			int connectionIndex[2];
			if(!spatialGraph->getNumberOfVertices())
			{
				connectionIndex[0] = 0;
				connectionIndex[1] = 1;
			}
			else
			{
				connectionIndex[0] = spatialGraph->getNumberOfVertices() - 2;
				connectionIndex[1] = spatialGraph->getNumberOfVertices() - 1;
			}
			int noOfAxisPoints = 2;
			std::list< double * > axisCoords;
			axisCoords.push_back(endPoint);
			axisCoords.push_back(barrelCentroid);
			Edge * newAxis = new Edge(connectionIndex, noOfAxisPoints, C1, axisCoords/*maybe score*/);
			spatialGraph->addEdge(newAxis);
		}
	}
};

// merge two ImageData structures of the SAME SPACING!!!
// no information overlap assumed, boundaries are allowed to overlap
ImageDataPointerType Geometry::mergeStructures(ImageDataPointerType structure1, ImageDataPointerType structure2)
{
	double * spacing1 = structure1->GetSpacing();
	double * spacing2 = structure2->GetSpacing();
	
	bool sameSpacing = 1;
	for(int ii = 0; ii < 3; ++ii)
		if(spacing1[ii] != spacing2[ii])
			sameSpacing = 0;
	if(!sameSpacing)
	{
		std::cout << "Error! Structures must have the same spacing in order to merge them!" << std::endl;
		return NULL;
	}
	if(structure1->GetScalarType() != structure2->GetScalarType())
	{
		std::cout << "Error! Structures must have the same scalar type in order to merge them!" << std::endl;
		return NULL;
	}
	
	// xMin, xMax, yMin, yMax, zMin, zMax
	int * bounds1 = structure1->GetExtent();
	int * bounds2 = structure2->GetExtent();
	
	//determine direction in which the structures are to be merged
	int mergeDirection = -1;
	for(int ii = 0; ii < 3; ++ii)
		if(bounds1[2*ii] != bounds2[2*ii] || bounds1[2*ii+1] != bounds2[2*ii+1])
			mergeDirection = ii;
	
	//determine whether the structures overlap or not
	//if not, merging is trivial
	bool overlap = 0;
	int topIndex = 0;
	if(bounds1[2*mergeDirection] < bounds2[2*mergeDirection+1] && bounds1[2*mergeDirection+1] > bounds2[2*mergeDirection+1])
	{
		overlap = 1;
		topIndex = 1;
	}
	else if(bounds2[2*mergeDirection] < bounds1[2*mergeDirection+1] && bounds2[2*mergeDirection+1] > bounds1[2*mergeDirection+1])
	{
		overlap = 1;
		topIndex = 2;
	}
	
// 	std::cout << "overlap = " << overlap << std::endl;
// 	std::cout << "topIndex = " << topIndex << std::endl;
// 	std::cout << "mergeDirection = " << mergeDirection << std::endl;
	
	int * totalExtent = new int[6];
	for(int ii = 0; ii < 3; ++ii)
	{
		if(ii == mergeDirection)
		{
			totalExtent[2*ii] = std::min(bounds1[2*ii], bounds2[2*ii]);
			totalExtent[2*ii+1] = std::max(bounds1[2*ii+1], bounds2[2*ii+1]);
		}
		else
		{
			totalExtent[2*ii] = bounds1[2*ii];
			totalExtent[2*ii+1] = bounds1[2*ii+1];
// 			std::cout << "bounds1[" << 2*ii << "] = " << bounds1[2*ii] << std::endl;
// 			std::cout << "bounds1[" << 2*ii+1 << "] = " << bounds1[2*ii+1] << std::endl;
// 			std::cout << "totalExtent[" << 2*ii << "] = " << totalExtent[2*ii] << std::endl;
// 			std::cout << "totalExtent[" << 2*ii+1 << "] = " << totalExtent[2*ii+1] << std::endl;
		}
	}
// 	std::cout << "bounds1 = [" << bounds1[0] << "," << bounds1[1] << "] x [" << bounds1[2] << "," << bounds1[3] << "] x [" << bounds1[4] << "," << bounds1[5] << "]" << std::endl;
// 	std::cout << "bounds2 = [" << bounds2[0] << "," << bounds2[1] << "] x [" << bounds2[2] << "," << bounds2[3] << "] x [" << bounds2[4] << "," << bounds2[5] << "]" << std::endl;
// 	std::cout << "extent = [" << totalExtent[0] << "," << totalExtent[1] << "] x [" << totalExtent[2] << "," << totalExtent[3] << "] x [" << totalExtent[4] << "," << totalExtent[5] << "]" << std::endl;
	ImageDataPointerType mergeVolume = ImageDataPointerType::New();
	mergeVolume->SetExtent(totalExtent);
	mergeVolume->SetSpacing(structure1->GetSpacing());
	mergeVolume->SetScalarType(structure1->GetScalarType());
	mergeVolume->AllocateScalars();
	const char * scalarType = mergeVolume->GetScalarTypeAsString();
// 	std::cout << "scalarType = " << scalarType << std::endl;
	if(!overlap)
	{
		if(scalarType == "unsigned char")
		{
			for(int z = totalExtent[4]; z <= totalExtent[5]; ++z)
				for(int y = totalExtent[2]; y <= totalExtent[3]; ++y)
					for(int x = totalExtent[0]; x <= totalExtent[1]; ++x)
					{
						unsigned char * px = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
						*px = 0;
					}
			for(int z = bounds1[4]; z <= bounds1[5]; ++z)
				for(int y = bounds1[2]; y <= bounds1[3]; ++y)
					for(int x = bounds1[0]; x <= bounds1[1]; ++x)
					{
						unsigned char * getPx = static_cast< unsigned char * >(structure1->GetScalarPointer(x, y, z));
						unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
			for(int z = bounds2[4]; z <= bounds2[5]; ++z)
				for(int y = bounds2[2]; y <= bounds2[3]; ++y)
					for(int x = bounds2[0]; x <= bounds2[1]; ++x)
					{
						unsigned char * getPx = static_cast< unsigned char * >(structure2->GetScalarPointer(x, y, z));
						unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
		}
		else if(scalarType == "char")
		{
			for(int z = totalExtent[4]; z <= totalExtent[5]; ++z)
				for(int y = totalExtent[2]; y <= totalExtent[3]; ++y)
					for(int x = totalExtent[0]; x <= totalExtent[1]; ++x)
					{
						char * px = static_cast< char * >(mergeVolume->GetScalarPointer(x, y, z));
						*px = 0;
					}
			for(int z = bounds1[4]; z <= bounds1[5]; ++z)
				for(int y = bounds1[2]; y <= bounds1[3]; ++y)
					for(int x = bounds1[0]; x <= bounds1[1]; ++x)
					{
						char * getPx = static_cast< char * >(structure1->GetScalarPointer(x, y, z));
						char * setPx = static_cast< char * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
			for(int z = bounds2[4]; z <= bounds2[5]; ++z)
				for(int y = bounds2[2]; y <= bounds2[3]; ++y)
					for(int x = bounds2[0]; x <= bounds1[1]; ++x)
					{
						char * getPx = static_cast< char * >(structure2->GetScalarPointer(x, y, z));
						char * setPx = static_cast< char * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
		}
		else if(scalarType == "float")
		{
			for(int z = totalExtent[4]; z <= totalExtent[5]; ++z)
				for(int y = totalExtent[2]; y <= totalExtent[3]; ++y)
					for(int x = totalExtent[0]; x <= totalExtent[1]; ++x)
					{
						float * px = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
						*px = 0;
					}
			for(int z = bounds1[4]; z <= bounds1[5]; ++z)
				for(int y = bounds1[2]; y <= bounds1[3]; ++y)
					for(int x = bounds1[0]; x <= bounds1[1]; ++x)
					{
						float * getPx = static_cast< float * >(structure1->GetScalarPointer(x, y, z));
						float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
			for(int z = bounds2[4]; z <= bounds2[5]; ++z)
				for(int y = bounds2[2]; y <= bounds2[3]; ++y)
					for(int x = bounds2[0]; x <= bounds1[1]; ++x)
					{
						float * getPx = static_cast< float * >(structure2->GetScalarPointer(x, y, z));
						float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
		}
		else if(scalarType == "double")
		{
			for(int z = totalExtent[4]; z <= totalExtent[5]; ++z)
				for(int y = totalExtent[2]; y <= totalExtent[3]; ++y)
					for(int x = totalExtent[0]; x <= totalExtent[1]; ++x)
					{
						double * px = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
						*px = 0;
					}
			for(int z = bounds1[4]; z <= bounds1[5]; ++z)
				for(int y = bounds1[2]; y <= bounds1[3]; ++y)
					for(int x = bounds1[0]; x <= bounds1[1]; ++x)
					{
						double * getPx = static_cast< double * >(structure1->GetScalarPointer(x, y, z));
						double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
			for(int z = bounds2[4]; z <= bounds2[5]; ++z)
				for(int y = bounds2[2]; y <= bounds2[3]; ++y)
					for(int x = bounds2[0]; x <= bounds1[1]; ++x)
					{
						double * getPx = static_cast< double * >(structure2->GetScalarPointer(x, y, z));
						double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
						*setPx = *getPx;
					}
		}
	}
	
	else if(overlap)
	{
		if(scalarType == "unsigned char")
		{
			if(topIndex == 1)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds2[0], maxX = bounds2[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2], maxY = bounds2[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4], maxZ = bounds2[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							unsigned char * getPx = static_cast< unsigned char * >(structure2->GetScalarPointer(x, y, z));
							unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds1[0] + 1, maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2] + 1, maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4] + 1, maxZ = bounds1[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							unsigned char * getPx = static_cast< unsigned char * >(structure1->GetScalarPointer(x, y, z));
							unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
			if(topIndex == 2)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds1[0], maxX = bounds1[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							unsigned char * getPx = static_cast< unsigned char * >(structure1->GetScalarPointer(x, y, z));
							unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds2[0] + 1, maxX = bounds2[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2] + 1, maxY = bounds2[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4] + 1, maxZ = bounds2[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							unsigned char * getPx = static_cast< unsigned char * >(structure2->GetScalarPointer(x, y, z));
							unsigned char * setPx = static_cast< unsigned char * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
		}
		if(scalarType == "float")
		{
			if(topIndex == 1)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds2[0], maxX = bounds2[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2], maxY = bounds2[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4], maxZ = bounds2[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							float * getPx = static_cast< float * >(structure2->GetScalarPointer(x, y, z));
							float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds1[0] + 1, maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2] + 1, maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4] + 1, maxZ = bounds1[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							float * getPx = static_cast< float * >(structure1->GetScalarPointer(x, y, z));
							float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
			if(topIndex == 2)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds1[0], maxX = bounds1[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							float * getPx = static_cast< float * >(structure1->GetScalarPointer(x, y, z));
							float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds2[0] + 1, maxX = bounds2[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2] + 1, maxY = bounds2[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4] + 1, maxZ = bounds2[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							float * getPx = static_cast< float * >(structure2->GetScalarPointer(x, y, z));
							float * setPx = static_cast< float * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
		}
		if(scalarType == "double")
		{
			if(topIndex == 1)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds2[0], maxX = bounds2[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2], maxY = bounds2[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4], maxZ = bounds2[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							double * getPx = static_cast< double * >(structure2->GetScalarPointer(x, y, z));
							double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds1[0] + 1, maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2] + 1, maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4] + 1, maxZ = bounds1[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							double * getPx = static_cast< double * >(structure1->GetScalarPointer(x, y, z));
							double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
			if(topIndex == 2)
			{
				int maxX, maxY, maxZ, minX, minY, minZ;
				if(mergeDirection == 0)
				{
					minX = bounds1[0], maxX = bounds1[1] - 1;
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3] - 1;
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5] - 1;
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							double * getPx = static_cast< double * >(structure1->GetScalarPointer(x, y, z));
							double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
				
				if(mergeDirection == 0)
				{
					minX = bounds2[0] + 1, maxX = bounds2[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 1)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds2[2] + 1, maxY = bounds2[3];
					minZ = bounds1[4], maxZ = bounds1[5];
				}
				if(mergeDirection == 2)
				{
					minX = bounds1[0], maxX = bounds1[1];
					minY = bounds1[2], maxY = bounds1[3];
					minZ = bounds2[4] + 1, maxZ = bounds2[5];
				}
				for(int z = minZ; z <= maxZ; ++z)
					for(int y = minY; y <= maxY; ++y)
						for(int x = minX; x <= maxX; ++x)
						{
							double * getPx = static_cast< double * >(structure2->GetScalarPointer(x, y, z));
							double * setPx = static_cast< double * >(mergeVolume->GetScalarPointer(x, y, z));
							*setPx = *getPx;
						}
			}
		}
	}
	
	mergeVolume->Update();
	return mergeVolume;
};

ImageDataPointerType Geometry::createImageVolumeFromPolyData(PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
{
	if(poly->GetNumberOfCells())
	{
		ImageDataPointerType volume = ImageDataPointerType::New();
		double spacing[3];
		switch(label)
		{
			case Pia:
				spacing[0] = spacing[1] = spacing[2] = 50;
				break;
			
			case WhiteMatter:
				spacing[0] = spacing[1] = spacing[2] = 50;
				break;
				
			default:
				spacing[0] = spacing[1] = spacing[2] = 1;
				break;
		}
			
		volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
		int * dims = calculateExtent(xMin, xMax, yMin, yMax, zMin, zMax, label);
// 				std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
// 				std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
		volume->SetExtent(dims);
		volume->SetNumberOfScalarComponents(1);
		volume->SetScalarTypeToUnsignedChar();
		// 	volume->Print(std::cout);
		volume->AllocateScalars();
		for(int z = dims[4]; z <= dims[5]; ++z)
			for(int y = dims[2]; y <= dims[3]; ++y)
				for(int x = dims[0]; x <= dims[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
					*px = 0;
				}
				volume->Update();
// 		std::cout << "Calculating pixels inside polygons!" << std::endl;
		unsigned long insidePoints = 0;
		unsigned long outsidePoints = 0;
		for(int currPlane = 0; currPlane < poly->GetNumberOfCells(); ++currPlane)
		{
			double closestPoint[3];
			int subId;
			double pCoords[3];
			double * weights = new double[poly->GetCell(currPlane)->GetNumberOfPoints()];
			double dist2;
			double tmp[3];
			poly->GetCell(currPlane)->GetPoints()->GetPoint(0, tmp);
			double * bounds = poly->GetCell(currPlane)->GetBounds();
// 			std::cout << "point tmp @ [" << tmp[0] << "," << tmp[1] << "," << tmp[2] << "]" << std::endl;
			int z = lround(tmp[2]/spacing[2]);
// 			std::flush(std::cout << "Determining inside/outside polygon for all points plane z = " << z << std::endl);
// 			std::flush(std::cout << "bounds = [" << bounds[0] << "," << bounds[1] << "], [" << bounds[2] << "," << bounds[3] << "],[" << bounds[4] << "," << bounds[5] << "]" << std::endl);
// 			#pragma omp parallel for
			for(int y = dims[2]; y <= dims[3]; ++y)
			{
				for(int x = dims[0]; x <= dims[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
					double tmpCoord[] = {x*spacing[0], y*spacing[1], z*spacing[2]};
					if(tmpCoord[0] < bounds[0] || tmpCoord[0] > bounds[1] 
						|| tmpCoord[1] < bounds[2] || tmpCoord[1] > bounds[3] 
						|| tmpCoord[2] < bounds[4] || tmpCoord[2] > bounds[5])
					{
						// 						std::flush(std::cout << "out of bounds." << std::endl);
						++outsidePoints;
						*px = 0;
						continue;
					}
					int insidePolygon = poly->GetCell(currPlane)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
					if(insidePolygon == 1)
					{
						// 						std::flush(std::cout << "hit!" << std::endl);
						*px = 255;
						++insidePoints;
					}
					else
					{
						++outsidePoints;
						*px = 0;
					}
				}
			}
// 			std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl;
		}
// 				std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl;
		volume->Update();
		return volume;
	}
	else
	{
		std::cout << "Error! Empty PolyData! Cannot convert to vtkImageData!" << std::endl;
		return NULL;
	}
};

ImageDataPointerType Geometry::createImageVolume(int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
{

	ImageDataPointerType volume = ImageDataPointerType::New();
	double spacing[3];
	if(Barrel <= label && label <= E6)
		label = Barrel;
	switch(label)
	{
		case Pia:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
		
		case WhiteMatter:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;

		case Barrel:
			spacing[0] = spacing[1] = spacing[2] = 10;
			break;
			
		default:
			spacing[0] = spacing[1] = spacing[2] = 1;
			break;
	}
		
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	int * dims = calculateExtent(xMin, xMax, yMin, yMax, zMin, zMax, label);
	// 		std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
	// 		std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToUnsignedChar();
	// 	volume->Print(std::cout);
	volume->AllocateScalars();
	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
				*px = 0;
			}

	volume->Update();
	return volume;
};


//deprecated; use AmiraSpatialGraph::extractLandmark() directly instead
PolyDataPointerType Geometry::createPolyDataFromPlanePointList(std::list< std::list< double * > > planewisePointList)
{
	if(planewisePointList.size())
	{
		PolyDataPointerType polyData = PolyDataPointerType::New();
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planewisePointList.begin(); planeEdgePointListIt != planewisePointList.end(); ++planeEdgePointListIt)
		{
			// 			std::flush(std::cout << "Allocating memory for " << planeEdgePointListIt->size() - 1 << " points in plane " << planeEdgePointListIt->back()[2] << "..." << std::endl);
			int end = planeEdgePointListIt->size() - 1;	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			
			std::list< double * >::iterator planeVertexIt = planeEdgePointListIt->begin();
			for(int ii = 0; ii < end; ++planeVertexIt, ++ii)
			{
				double * tmp = *planeVertexIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
		}
		polyData->SetPoints(points);
		polyData->Update();
		return polyData;
	}
	else
	{
		std::cout << "Error! Empty PointList! Cannot create vtkPolyData!" << std::endl;
		return NULL;
	}
};

// does what you think it does
void Geometry::normalize(double * vec)
{
	double norm = 0;
	for(int ii = 0; ii < 3; ++ii)
		norm += vec[ii]*vec[ii];
	norm = sqrt(norm);
	if(norm)
		for(int ii = 0; ii < 3; ++ii)
			vec[ii] = vec[ii]/norm;
};










