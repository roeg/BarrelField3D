#include "spatialGraph.h"

// Basics::Basics()
// {
// 	inputFilename	= 0;
// 	outputFilename	= 0;
// };
// 
// Basics::~Basics()
// {
// 	
// };

Vertex::Vertex(float * _coordinates, int _label)
{
	coordinates[0] = (double)_coordinates[0];
	coordinates[1] = (double)_coordinates[1];
	coordinates[2] = (double)_coordinates[2];
	label = _label;
};

Vertex::Vertex(double * _coordinates, int _label)
{
	coordinates[0] = _coordinates[0];
	coordinates[1] = _coordinates[1];
	coordinates[2] = _coordinates[2];
	label = _label;
};

Vertex::~Vertex()
{
	//tbd
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = 0.0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, float _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = _radius;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = 0.0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = _radius;
};

Edge::~Edge()
{
	//tbd
}

AmiraSpatialGraph::AmiraSpatialGraph()
{
	
};

AmiraSpatialGraph::~AmiraSpatialGraph()
{
	std::list< Edge * >::iterator it1;
	std::list< Vertex * >::iterator it2;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
		delete *it1;
	for(it2 = vertices.begin(); it2 != vertices.end(); ++it2)
		delete *it2;
	
	vertices.clear();
	edges.clear();
};

unsigned int AmiraSpatialGraph::getNumberOfPoints()
{
	unsigned int number = 0;
	std::list< Edge * >::iterator it1;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
	{
		number += (*it1)->numEdgePoints;
	}
	
	return number;
}

void AmiraSpatialGraph::addVertex(Vertex * newVertex)
{
	vertices.push_back(newVertex);
};

void AmiraSpatialGraph::addEdge(Edge * newEdge)
{
	edges.push_back(newEdge);
};

std::list< Vertex* >::iterator AmiraSpatialGraph::verticesBegin()
{
	std::list< Vertex * >::iterator it = vertices.begin();
	return it;
};

std::list< Edge* >::iterator AmiraSpatialGraph::edgesBegin()
{
	std::list< Edge * >::iterator it = edges.begin();
	return it;
};

std::list< Vertex* >::iterator AmiraSpatialGraph::verticesEnd()
{
	std::list< Vertex * >::iterator it = vertices.end();
	return it;
};

std::list< Edge* >::iterator AmiraSpatialGraph::edgesEnd()
{
	std::list< Edge * >::iterator it = edges.end();
	return it;
};

void AmiraSpatialGraph::vesselsToPoints()
{
	std::list< Edge * >::iterator edgeIt;
	std::list< Vertex * >::iterator vertexIt;
	std::vector< Vertex * > verticesVec;
	
	std::flush(std::cout << "converting blood vessels to points..." << std::endl);
	
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
	{
		verticesVec.push_back(*vertexIt);
	}
	
	for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
	{
		if((*edgeIt)->label == Vessel)
		{
			double circ = 0;
			double * center = new double[3];
			center[0] = 0, center[1] = 0, center[2] = 0;
			double * curr;
			double * last;
			std::list< double * >::iterator it = (*edgeIt)->edgePointCoordinates.begin();
			last = *it;
			++it;
			while(it != (*edgeIt)->edgePointCoordinates.end())
			{
				curr = *it;
				circ += std::sqrt((curr[0] - last[0])*(curr[0] - last[0]) + (curr[1] - last[1])*(curr[1] - last[1]) + (curr[2] - last[2])*(curr[2] - last[2]));
				center[0] += curr[0];
				center[1] += curr[1];
				center[2] += curr[2];
				last = curr;
				++it;
			}
			center[0] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[1] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[2] /= (double)((*edgeIt)->numEdgePoints - 1);
			(*edgeIt)->numEdgePoints = 1;
			(*edgeIt)->edgePointCoordinates.clear();
			(*edgeIt)->edgePointCoordinates.push_back(center);
			(*edgeIt)->radius = circ/(2*PI);
			
			for(int ii = 0; ii < 3; ++ii)
			{
				verticesVec[(*edgeIt)->edgeConnectivity[0]]->coordinates[ii] = center[ii];
				verticesVec[(*edgeIt)->edgeConnectivity[1]]->coordinates[ii] = center[ii];
			}
		}
	}
	
	int ii = 0;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end() && ii < verticesVec.size(); ++vertexIt, ++ii)
	{
		for(int jj = 0; jj < 3; ++jj)
			(*vertexIt)->coordinates[jj] = verticesVec[ii]->coordinates[jj];
	}
	
	removeDuplicatePoints();
};

void AmiraSpatialGraph::removeDuplicatePoints()
{
	std::list< Edge* >::iterator edgeIt = edges.begin();
	std::list< Vertex* >::iterator vertexIt = vertices.begin();
	bool duplicate = 0;
	while(edgeIt != edges.end() && vertexIt != vertices.end())
	{
		bool duplicateVertex = (bool)((*edgeIt)->edgeConnectivity[0] != (*edgeIt)->edgeConnectivity[1]);
		duplicate = 0;
		
		if((*edgeIt)->label == Vessel)
		{
			std::list< Edge* >::iterator checkIt = edgeIt;
			std::list< Vertex* >::iterator checkIt2 = vertexIt;
			++checkIt;
			++checkIt2;
			if(duplicateVertex)
				++checkIt2;
			
			while(checkIt != edges.end() && checkIt2 != vertices.end())
			{
				if((*edgeIt)->label == Vessel)
				{
					double * tmp1 = (*edgeIt)->edgePointCoordinates.front();
					double * tmp2 = (*checkIt)->edgePointCoordinates.front();
					
					if(std::abs(tmp1[2] - tmp2[2]) > 1.0)
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
						continue;
					}
					
					float dist = std::sqrt((tmp1[0] - tmp2[0])*(tmp1[0] - tmp2[0]) + (tmp1[1] - tmp2[1])*(tmp1[1] - tmp2[1]) + (tmp1[2] - tmp2[2])*(tmp1[2] - tmp2[2]));
					if(dist < (*edgeIt)->radius || dist < (*checkIt)->radius)
					{
						if((*edgeIt)->radius < (*checkIt)->radius)
						{
							duplicate = 1;
							edgeIt = edges.erase(edgeIt);
							vertexIt = vertices.erase(vertexIt);
							if(duplicateVertex)
								vertexIt = vertices.erase(vertexIt);
							
							std::list< Edge* >::iterator updateIt;
							for(updateIt = edgeIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
							++checkIt;
							++checkIt2;
							if(duplicateVertex)
								++checkIt2;
						}
						else
						{
							duplicate = 0;
							checkIt = edges.erase(checkIt);
							checkIt2 = vertices.erase(checkIt2);
							if(duplicateVertex)
								checkIt2 = vertices.erase(checkIt2);
							
							std::list< Edge* >::iterator updateIt;
							for(updateIt = checkIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
						}
					}
					else
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
					}
				}
				else
				{
					++checkIt;
					++checkIt2;
					if(duplicateVertex)
						++checkIt2;
				}
			}
		}
		
		if(!duplicate)
		{
			++edgeIt;
			++vertexIt;
			if(duplicateVertex)
				++vertexIt;
		}
	}
};


//extract all points of landmark 'label' planewise; return their plane indices in zIndexList; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList)
{
	std::list< Edge* >::iterator edgeIt;
	edgeIt = this->edges.begin();
	if(edgeIt != this->edges.end())
	{
		int zIndex = 0;
		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
			++edgeIt;
		
		zIndex = lround((*edgeIt)->edgePointCoordinates.back()[2]);
		zIndexList.push_back(zIndex);
		std::list< double * > planeCoordinates;
		std::list< double * >::iterator pointListIt;
		while(edgeIt != this->edges.end())
		{
			if((*edgeIt)->label != label)
			{
				++edgeIt;
				continue;
			}
			int tmpZ = lround((*edgeIt)->edgePointCoordinates.back()[2]);
			// 			std::flush(std::cout << "in plane " << tmpZ << std::endl);
			if(tmpZ == zIndex)
			{
				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
				{
					double * tmpPoint = *pointListIt;
					// 					tmpPoint[2] = (double)(int)(tmpPoint[2] + 0.5);
					tmpPoint[2] = round(tmpPoint[2]);
					planeCoordinates.push_back(*pointListIt);
					planeCoordinates.back()[2] = tmpPoint[2];
				}
			}
			else
			{
				std::list< double * > tmpList(planeCoordinates);
				planewisePointList.push_back(tmpList);
				planeCoordinates.clear();
				
				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
				{
					double * tmpPoint = *pointListIt;
					tmpPoint[2] = round(tmpPoint[2]);
					planeCoordinates.push_back(*pointListIt);
					planeCoordinates.back()[2] = tmpPoint[2];
				}
				zIndex = tmpZ;
				zIndexList.push_back(zIndex);
			}
			++edgeIt;
		}
		zIndexList.sort();
		planewisePointList.push_back(planeCoordinates);
		std::list< std::list< double * > >::const_iterator controlit;
		for(controlit = planewisePointList.begin(); controlit != planewisePointList.end(); ++controlit)
			if(controlit->size())
				return 1;
		
		return 0;
	}
	
	return 0;
};

//extract all points of landmark 'label' as PolyData; return their plane indices in zIndexList; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList)
{
	std::list< Edge* >::iterator edgeIt;
	edgeIt = this->edges.begin();
	if(edgeIt != this->edges.end())
	{
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
			++edgeIt;
		
		std::list< double * >::iterator pointListIt;
		while(edgeIt != this->edges.end())
		{
			if((*edgeIt)->label != label)
			{
				++edgeIt;
				continue;
			}
			// 			std::flush(std::cout << "Allocating memory for " << planeEdgePointListIt->size() - 1 << " points in plane " << planeEdgePointListIt->back()[2] << "..." << std::endl);
			int end = (*edgeIt)->edgePointCoordinates.size() - 1;	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			pointListIt = (*edgeIt)->edgePointCoordinates.begin();
			for(int ii = 0; ii < end; ++pointListIt, ++ii)
			{
				double * tmp = *pointListIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			++edgeIt;
		}
		polyData->SetPoints(points);
		polyData->Update();
		
		for(int ii = 0; ii < polyData->GetNumberOfCells(); ++ii)
		{
			double * bounds = polyData->GetCell(ii)->GetBounds();
			if(lround(bounds[4]) == lround(bounds[5]))
				zIndexList.push_back(lround(bounds[4]));
		}
		zIndexList.sort();
		std::cout << "nr of z-elements incl. duplicates = " << zIndexList.size() << std::endl;
		std::list< int >::iterator removeDuplicatesIt;
		int currentZ = zIndexList.front();
		removeDuplicatesIt = zIndexList.begin();
		++removeDuplicatesIt;
		while(removeDuplicatesIt != zIndexList.end())
		{
			if(*removeDuplicatesIt == currentZ)
				removeDuplicatesIt = zIndexList.erase(removeDuplicatesIt);
			else
			{
				currentZ = *removeDuplicatesIt;
				++removeDuplicatesIt;
			}
		}
		std::cout << "nr of z-elements after removing duplicates = " << zIndexList.size() << std::endl;
		
		if(polyData->GetNumberOfPoints())
			return 1;
			
		return 0;
	}
	
	return 0;
};




















