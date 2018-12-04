#include "vpmreconstruction.h"

VPMReconstruction::VPMReconstruction ( AmiraSpatialGraph * inputSG, InputParameters parameters )
{
	this->parameters = parameters;
	this->spatialGraph = new AmiraSpatialGraph;
	this->spatialGraph->mergeSpatialGraph(inputSG);
	SBF = new BarrelField(false);
}

VPMReconstruction::~VPMReconstruction()
{
	if(spatialGraph) delete spatialGraph;
	if(SBF) delete SBF;
}

std::map< int, ClosedSurface * > VPMReconstruction::getBarreloids()
{
	return barreloidField;
}

void VPMReconstruction::startReconstruction ( const char* outputFilename )
{
	std::cout << "Starting 3D reconstruction of barreloids!" << std::endl;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		PolyDataPointerType barreloid = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(ID, barreloid))
		{
			BarreloidSurfaceReconstruct * barreloidReconstruct = new BarreloidSurfaceReconstruct(parameters, barreloid);
			ClosedSurface * barreloidSurface = new ClosedSurface(barreloidReconstruct->surfaceReconstruction(ID));
			barreloidField.insert(std::pair< int, ClosedSurface * >(ID, barreloidSurface));
			
			std::vector< double > barreloidCenter = barreloidReconstruct->getEllipsoidCenter();
			std::vector< std::vector< double > > barreloidAxes = barreloidReconstruct->getPrincipalAxes();
			barreloidCenters.insert(std::pair< int, std::vector< double > >(ID, barreloidCenter));
			barreloidEVecs.insert(std::pair< int, std::vector< std::vector< double > > >(ID, barreloidAxes));
			
			delete barreloidReconstruct;
			
			std::string barreloidName(outputFilename);
			barreloidName += "_";
			barreloidName += SBF->int2Labels[ID];
			Reader * surfWriter = new Reader(barreloidName.c_str(), barreloidName.c_str());
			surfWriter->writeAmiraSurfaceFile(barreloidSurface->ptr());
			delete surfWriter;
		}
	}
	
	computeBarreloidDimensions();
	
	writeBarreloidParameters(outputFilename);
}

void VPMReconstruction::writeBarreloidParameters ( const char* outputFilename )
{
	std::string ofName(outputFilename);
	ofName += "_parameters.csv";
	std::ofstream BarreloidParameters;
	BarreloidParameters.open(ofName.c_str());
	BarreloidParameters << "# Barreloid\tVolume\tz extent\trow extent\tarc extent" << std::endl;
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barreloidField.find(ID) != barreloidField.end())
		{
			double vol = barreloidVolume(barreloidField[ID]->ptr());
			BarreloidParameters << SBF->int2Labels[ID] << "\t" << vol << "\t";
			if(barreloidDimensions.find(ID) != barreloidDimensions.end())
			{
				BarreloidParameters << barreloidDimensions[ID][0] << "\t";
				BarreloidParameters << barreloidDimensions[ID][1] << "\t";
				BarreloidParameters << barreloidDimensions[ID][2];
			}
			BarreloidParameters << std::endl;
		}
	}
	
	BarreloidParameters.close();
}

// compute volume of each barreloid
// by summing the signed volumes of all
// tetraheda formed by all surface triangles
// with the origin (see Amira documentation).
// works since all faces are oriented consistently
double VPMReconstruction::barreloidVolume ( PolyDataPointerType barreloid )
{
	double totalVol = 0;
	double origin[] = {0,0,0};
	for(int ii = 0; ii < barreloid->GetNumberOfCells(); ++ii)
	{
		double * pts = new double[9];
		PointsPointerType cellPts = barreloid->GetCell(ii)->GetPoints();
		for(int jj = 0; jj < 3; ++jj)
			cellPts->GetPoint(jj, pts+3*jj);
		
		double vol = vtkTetra::ComputeVolume(pts, pts+3, pts+6, origin);
		totalVol += vol;
		
		delete [] pts;
	}
	
	return fabs(totalVol);
}

// compute extent of barreloid along
// eigenvectors centered on center
void VPMReconstruction::computeBarreloidDimensions()
{
	const unsigned int Z = 0;
	const unsigned int ROW = 1;
	const unsigned int ARC = 2;
	
	if(barreloidCenters.size() && barreloidEVecs.size())
	{
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(barreloidCenters.find(ID) != barreloidCenters.end() && barreloidEVecs.find(ID) != barreloidEVecs.end())
			{
				// find eigenvectors corresponding
				// to z/row/arc
				// using inner product of
				// eigenvectors with z axis
				// and row direction vector
				int eVecOrder[3];
				
				double maxScore = 0;
				int maxIndex;
				for(int ii = 0; ii < 3; ++ii)
				{
					double tmpScore = fabs(barreloidEVecs[ID][ii][2]);
					if(tmpScore > maxScore)
					{
						maxScore = tmpScore;
						maxIndex = ii;
					}
				}
				eVecOrder[Z] = maxIndex;
				
				int previousBarrel = ID-1, nextBarrel = ID+1;
				if(std::find(SBF->barrelLabels.begin(), SBF->barrelLabels.end(), previousBarrel) == SBF->barrelLabels.end()
					|| barreloidCenters.find(previousBarrel) == barreloidCenters.end())
					previousBarrel == ID;
				if(std::find(SBF->barrelLabels.begin(), SBF->barrelLabels.end(), nextBarrel) == SBF->barrelLabels.end()
					|| barreloidCenters.find(nextBarrel) == barreloidCenters.end())
					nextBarrel == ID;
				if(previousBarrel == nextBarrel)
				{
					std::cout << "Warning! Cannot determine row/arc orientation of " << SBF->int2Labels[ID];
					std::cout << ". Writing row/arc dimensions in random order." << std::endl;
					eVecOrder[ROW] = (eVecOrder[Z]+1)%3;
					eVecOrder[ARC] = (eVecOrder[Z]+2)%3;
				}
				else
				{
					double direction[3];
					for(int ii = 0; ii < 3; ++ii)
						direction[ii] = barreloidCenters[nextBarrel][ii] - barreloidCenters[previousBarrel][ii];
					vtkMath::Normalize(direction);
					
					int index1, index2;
					double score1, score2;
					double vec1[3], vec2[3];
					
					index1 = (eVecOrder[Z]+1)%3;
					index2 = (eVecOrder[Z]+2)%3;
					
					for(int ii = 0; ii < 3; ++ii)
					{
						vec1[ii] = barreloidEVecs[ID][index1][ii];
						vec2[ii] = barreloidEVecs[ID][index2][ii];
					}
					
					score1 = fabs(vtkMath::Dot(vec1, direction));
					score2 = fabs(vtkMath::Dot(vec2, direction));
					
					if(score1 > score2)
					{
						eVecOrder[ROW] = index1;
						eVecOrder[ARC] = index2;
					}
					else
					{
						eVecOrder[ARC] = index1;
						eVecOrder[ROW] = index2;
					}
				}
				
				// compute dimensions
				std::vector< double > dimensions;
				double zDirection[3], rowDirection[3], arcDirection[3];
				for(int ii = 0; ii < 3; ++ii)
				{
					zDirection[ii] = barreloidEVecs[ID][eVecOrder[Z]][ii];
					rowDirection[ii] = barreloidEVecs[ID][eVecOrder[ROW]][ii];
					arcDirection[ii] = barreloidEVecs[ID][eVecOrder[ARC]][ii];
				}
				
				// find extent along z axis
				// by selecting points with min/max z
				// and projecting them on axis through center
				double maxPt[3], minPt[3];
				double maxZ, minZ;
				PointsPointerType barreloidPts = barreloidField[ID]->ptr()->GetPoints();
				barreloidPts->GetPoint(0, maxPt);
				barreloidPts->GetPoint(0, minPt);
				maxZ = maxPt[2], minZ = minPt[2];
				
				for(int ii = 0; ii < barreloidPts->GetNumberOfPoints(); ++ii)
				{
					double pt[3];
					barreloidPts->GetPoint(ii, pt);
					if(pt[2] < minZ)
					{
						minZ = pt[2];
						minPt[0] = pt[0], minPt[1] = pt[1], minPt[2] = pt[2];
					}
					if(pt[2] > maxZ)
					{
						maxZ = pt[2];
						maxPt[0] = pt[0], maxPt[1] = pt[1], maxPt[2] = pt[2];
					}
				}
				
				double zLinePt1[3], zLinePt2[3], zExtent1[3], zExtent2[3], t;
				for(int ii = 0; ii < 3; ++ii)
				{
					zLinePt1[ii] = barreloidCenters[ID][ii] + 5*(maxZ-minZ)*zDirection[ii];
					zLinePt2[ii] = barreloidCenters[ID][ii] - 5*(maxZ-minZ)*zDirection[ii];
				}
				vtkLine::DistanceToLine(maxPt, zLinePt1, zLinePt2, t, zExtent1);
				vtkLine::DistanceToLine(minPt, zLinePt1, zLinePt2, t, zExtent2);
				
				// find extent along row/arc
				double center[3], rowExtent1[3], rowExtent2[3], arcExtent1[3], arcExtent2[3];
				center[0] = barreloidCenters[ID][0];
				center[1] = barreloidCenters[ID][1];
				center[2] = barreloidCenters[ID][2];
				
				barreloidField[ID]->intersectLineInDirection(rowDirection, center);
				barreloidField[ID]->getLastIntersectPoint(rowExtent1);
				vtkMath::MultiplyScalar(rowDirection, -1.0);
				barreloidField[ID]->intersectLineInDirection(rowDirection, center);
				barreloidField[ID]->getLastIntersectPoint(rowExtent2);
				
				barreloidField[ID]->intersectLineInDirection(arcDirection, center);
				barreloidField[ID]->getLastIntersectPoint(arcExtent1);
				vtkMath::MultiplyScalar(arcDirection, -1.0);
				barreloidField[ID]->intersectLineInDirection(arcDirection, center);
				barreloidField[ID]->getLastIntersectPoint(arcExtent2);
				
				dimensions.push_back(sqrt(vtkMath::Distance2BetweenPoints(zExtent1, zExtent2)));
				dimensions.push_back(sqrt(vtkMath::Distance2BetweenPoints(rowExtent1, rowExtent2)));
				dimensions.push_back(sqrt(vtkMath::Distance2BetweenPoints(arcExtent1, arcExtent2)));
				barreloidDimensions.insert(std::pair< int, std::vector< double > >(ID, dimensions));
			}
		}
		
	}
}

