/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
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

PolyDataPointerType firstSamplingPass(std::vector< Surface * > surfaceVec);
PolyDataPointerType firstSamplingPassWithNormals(std::vector< Surface * > surfaceVec);

void createL4FromAvgPia(PolyDataPointerType avgPia, PolyDataPointerType L4Upper, PolyDataPointerType L4Lower, AmiraSpatialGraph * zAxes);
std::multimap< double, double * > barrelAxisScores(PolyDataPointerType piaSurface, double * centroid, double alpha);
int closestBarrel(double samplePt[3], std::map< int, double * > avgCenters);
double * computeLayerThicknesses(int barrel, std::map< int, Column * > avgColumns, std::map< int, Column * > avgBarrels);
PolyDataPointerType readStandardBarrelField(std::map< int, Column * >& avgColumns,
			     std::map< int, Column * >& avgBarrels, std::map< int, double * >& avgAxes, std::map< int, double * >& avgCenters);

int main( int argc , char * argv[])
{
	if(argc > 3)
	{
		std::vector< const char * > inputNameVec;
		std::vector< Reader * > amiraReaderVec;
		std::vector< Surface * > surfaceData;
		
		const char * outputFilename = argv[argc-1];
		std::string outString(outputFilename);
		
		for(int ii = 1; ii < argc-1; ++ii)
		{
			inputNameVec.push_back(argv[ii]);
			Reader * newReader = new Reader(inputNameVec[ii-1], inputNameVec[ii-1]);
			Surface * newSurface = new Surface(newReader->readAmiraSurfaceFile());
			surfaceData.push_back(newSurface);
			delete newReader;
		}
		
		PolyDataPointerType avgSurface = firstSamplingPass(surfaceData);
// 		PolyDataPointerType avgSurface = firstSamplingPassWithNormals(surfaceData);
		LowpassPolyDataFilterType smoothingFilter = LowpassPolyDataFilterType::New();
		smoothingFilter->BoundarySmoothingOff();
		smoothingFilter->FeatureEdgeSmoothingOff();
		smoothingFilter->NormalizeCoordinatesOn();
		smoothingFilter->SetNumberOfIterations(20);
		smoothingFilter->SetPassBand(0.1);
		smoothingFilter->SetInput(avgSurface);
		smoothingFilter->Update();
		avgSurface = smoothingFilter->GetOutput();
		
		Reader * surfWriter = new Reader(outString.c_str(), outString.c_str());
		surfWriter->writeAmiraSurfaceFile(avgSurface);
		delete surfWriter;
	}
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		Reader * surfReader = new Reader(inputFilename, inputFilename);
		PolyDataPointerType avgPia = surfReader->readAmiraSurfaceFile();
		delete surfReader;
		
		PolyDataPointerType L4Upper, L4Lower;
		L4Upper = PolyDataPointerType::New();
		L4Lower = PolyDataPointerType::New();
		L4Upper->Allocate(1);
		L4Lower->Allocate(1);
		AmiraSpatialGraph * zAxesFromAvgPia = new AmiraSpatialGraph;
		createL4FromAvgPia(avgPia, L4Upper, L4Lower, zAxesFromAvgPia);
		
		std::string zAxesStr(outputFilename);
		std::string L4UStr(outputFilename);
		std::string L4LStr(outputFilename);
		zAxesStr += "_axis_field";
		L4UStr += "_L4Upper";
		L4LStr += "_L4Lower";
		
		std::flush(std::cout << "Writing output files" << std::endl);
		
		Reader * axisWriter = new Reader(zAxesStr.c_str(), zAxesStr.c_str());
		axisWriter->setSpatialGraph(zAxesFromAvgPia);
		axisWriter->writeSpatialGraphFile();
		delete axisWriter;
		
		Reader * L4UWriter = new Reader(L4UStr.c_str(), L4UStr.c_str());
		L4UWriter->writeAmiraSurfaceFile(L4Upper);
		delete L4UWriter;
		
		Reader * L4LWriter = new Reader(L4LStr.c_str(), L4LStr.c_str());
		L4LWriter->writeAmiraSurfaceFile(L4Lower);
		delete L4LWriter;
	}
	
	return 0;
}

// calculate avg surface sampled along global z direction
PolyDataPointerType firstSamplingPass(std::vector< Surface * > surfaceVec)
{
	double * bounds = surfaceVec[0]->ptr()->GetBounds();
	for(int ii = 1; ii < surfaceVec.size(); ++ii)
	{
		double * tmpBounds = surfaceVec[ii]->ptr()->GetBounds();
		if(tmpBounds[0] < bounds[0])
			bounds[0] = tmpBounds[0];
		if(tmpBounds[1] > bounds[1])
			bounds[1] = tmpBounds[1];
		if(tmpBounds[2] < bounds[2])
			bounds[2] = tmpBounds[2];
		if(tmpBounds[3] > bounds[3])
			bounds[3] = tmpBounds[3];
		if(tmpBounds[4] < bounds[4])
			bounds[4] = tmpBounds[4];
		if(tmpBounds[5] > bounds[5])
			bounds[5] = tmpBounds[5];
	}
	
	unsigned int xSamples, ySamples;
	double samplingStep = 50, centerZPlane;
	xSamples = (unsigned int)((bounds[1] - bounds[0])/samplingStep) + 1;
	ySamples = (unsigned int)((bounds[3] - bounds[2])/samplingStep) + 1;
	centerZPlane = 0.5*(bounds[4] + bounds[5]);
	
	std::cout << "Computing triangulation of average surface..." << std::endl;
	Delaunay2DFilterPointerType delTriangulation1 = Delaunay2DFilterPointerType::New();
	PointsPointerType piaPoints = PointsPointerType::New();
	PolyDataPointerType piaSampling = PolyDataPointerType::New(); // Delaunay filter needs DataObject as input...
	piaPoints->SetDataTypeToFloat();
	piaSampling->Allocate(1);
	
	for(int ii = 0; ii < xSamples; ++ii)
		for(int jj = 0; jj < ySamples; ++jj)
		{
			double avgPt[3], samplePt[3], sampleAxis[3], tmpZ = 0;
			sampleAxis[0] = 0, sampleAxis[1] = 0, sampleAxis[2] = 1;
			samplePt[0] = bounds[0] + ii*samplingStep;
			samplePt[1] = bounds[2] + jj*samplingStep;
			samplePt[2] = centerZPlane;
			
			int avgCnt = 0;
			for(int n  = 0; n < surfaceVec.size(); ++n)
			{
				surfaceVec[n]->intersectLine(sampleAxis, samplePt);
				if(surfaceVec[n]->getLastIntersectPoint())
				{
					++avgCnt;
					tmpZ += surfaceVec[n]->getLastIntersectPoint()[2];
				}
			}
			if(avgCnt)
			{
				tmpZ /= avgCnt;
				avgPt[0] = samplePt[0], avgPt[1] = samplePt[1], avgPt[2] = tmpZ;
				piaPoints->InsertNextPoint(avgPt);
			}
		}
	piaSampling->SetPoints(piaPoints);
	delTriangulation1->SetInput(piaSampling);
	delTriangulation1->Update();
	return delTriangulation1->GetOutput();
};

// calculate avg surface sampled along local average surface normal
PolyDataPointerType firstSamplingPassWithNormals(std::vector< Surface * > surfaceVec)
{
	double * bounds = surfaceVec[0]->ptr()->GetBounds();
	for(int ii = 1; ii < surfaceVec.size(); ++ii)
	{
		double * tmpBounds = surfaceVec[ii]->ptr()->GetBounds();
		if(tmpBounds[0] < bounds[0])
			bounds[0] = tmpBounds[0];
		if(tmpBounds[1] > bounds[1])
			bounds[1] = tmpBounds[1];
		if(tmpBounds[2] < bounds[2])
			bounds[2] = tmpBounds[2];
		if(tmpBounds[3] > bounds[3])
			bounds[3] = tmpBounds[3];
		if(tmpBounds[4] < bounds[4])
			bounds[4] = tmpBounds[4];
		if(tmpBounds[5] > bounds[5])
			bounds[5] = tmpBounds[5];
	}
	
	unsigned int xSamples, ySamples;
	double samplingStep = 50, centerZPlane;
	xSamples = (unsigned int)((bounds[1] - bounds[0])/samplingStep) + 1;
	ySamples = (unsigned int)((bounds[3] - bounds[2])/samplingStep) + 1;
	centerZPlane = 0.5*(bounds[4] + bounds[5]);
	
	std::cout << "Computing triangulation of average surface..." << std::endl;
	Delaunay2DFilterPointerType delTriangulation1 = Delaunay2DFilterPointerType::New();
	PointsPointerType piaPoints = PointsPointerType::New();
	PolyDataPointerType piaSampling = PolyDataPointerType::New(); // Delaunay filter needs DataObject as input...
	piaPoints->SetDataTypeToFloat();
	piaSampling->Allocate(1);
	
	for(int ii = 0; ii < xSamples; ++ii)
		for(int jj = 0; jj < ySamples; ++jj)
		{
			double avgPt[3], samplePt[3], sampleAxis[3], avgNormal[3], tmpZ = 0;
			avgPt[0] = 0, avgPt[1] = 0, avgPt[2] = 0;
			sampleAxis[0] = 0, sampleAxis[1] = 0, sampleAxis[2] = 1;
			avgNormal[0] = 0, avgNormal[1] = 0, avgNormal[2] = 0;
			samplePt[0] = bounds[0] + ii*samplingStep;
			samplePt[1] = bounds[2] + jj*samplingStep;
			samplePt[2] = centerZPlane;
			
			int avgCnt = 0;
			for(int n  = 0; n < surfaceVec.size(); ++n)
			{
				surfaceVec[n]->intersectLine(sampleAxis, samplePt);
				if(surfaceVec[n]->getLastIntersectPoint())
				{
					++avgCnt;
					tmpZ += surfaceVec[n]->getLastIntersectPoint()[2];
					double normal[3];
					vtkIdType tmpID = surfaceVec[n]->getLastIntersectCellID();
					vtkPolygon::ComputeNormal(surfaceVec[n]->ptr()->GetCell(tmpID)->GetPoints(), normal);
					avgNormal[0] += normal[0], avgNormal[1] += normal[1], avgNormal[2] += normal[2];
				}
			}
			if(avgCnt)
			{
				tmpZ /= avgCnt;
				avgNormal[0] /= avgCnt, avgNormal[1] /= avgCnt, avgNormal[2] /= avgCnt;
				vtkMath::Normalize(avgNormal);
				samplePt[2] = tmpZ;
				int avgCntN = 0;
				for(int n  = 0; n < surfaceVec.size(); ++n)
				{
					surfaceVec[n]->intersectLine(avgNormal, samplePt);
					double * pt = surfaceVec[n]->getLastIntersectPoint();
					if(pt)
					{
						++avgCntN;
// 						vtkMath::Add(avgPt, pt, avgPt);
						avgPt[0] += pt[0], avgPt[1] += pt[1], avgPt[2] += pt[2];
					}
				}
				if(avgCntN)
				{
					vtkMath::MultiplyScalar(avgPt, 1/(double)avgCntN);
					piaPoints->InsertNextPoint(avgPt);
				}
			}
		}
	piaSampling->SetPoints(piaPoints);
	delTriangulation1->SetInput(piaSampling);
	delTriangulation1->Update();
	return delTriangulation1->GetOutput();
};

/******************************************************************************/
/*sample grid with radius 3 mm around D2 column (origin & z axis) determining */
/*local z-axis and layer borders based on avgPia and closest column layer     */
/*z ratios                                                                    */
/******************************************************************************/
void createL4FromAvgPia(PolyDataPointerType avgPia, PolyDataPointerType L4Upper, PolyDataPointerType L4Lower, AmiraSpatialGraph * zAxes)
{
	std::map< int, Column * > avgColumns, avgBarrels;
	std::map< int, double * > avgAxes, avgCenters;
	PolyDataPointerType avgAxesField = readStandardBarrelField(avgColumns, avgBarrels, avgAxes, avgCenters);
	
	Surface * piaIntersector = new Surface(avgPia);
	
	std::flush(std::cout << "Computing convex hull around already established axis field..." << std::endl);
	
	// compute convex hull around barrel field z axes
	// b/c we don't want to create addtl z axes there
	// but only around the already existing ones
	ConvexHullFilterPointerType convexHullFilter = ConvexHullFilterPointerType::New();
	convexHullFilter->SetInput(avgAxesField);
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
	// just to be sure it's doing what it's supposed to:
	zAxes->addPolyDataObject(barrelFieldHull, Barrel);
	
	std::vector< PlanePointerType > hullPlanes;
	for(int ii = 0; ii < barrelFieldHull->GetNumberOfCells(); ++ii)
	{
		PlanePointerType hullFace = PlanePointerType::New();
		double normal[3], origin[3], pCoords[3], * weights = new double[barrelFieldHull->GetCell(ii)->GetNumberOfPoints()];
		int subId;
		barrelFieldHull->GetCell(ii)->GetParametricCenter(pCoords);
		barrelFieldHull->GetCell(ii)->EvaluateLocation(subId, pCoords, origin, weights);
		vtkPolygon::ComputeNormal(barrelFieldHull->GetCell(ii)->GetPoints(), normal);
		hullFace->SetNormal(normal);
		hullFace->SetOrigin(origin);
		hullPlanes.push_back(hullFace);
		delete [] weights;
	}
	
	std::flush(std::cout << "Triangulating L4 borders and z axis field outside of barrel field..." << std::endl);
	
	Delaunay2DFilterPointerType delaunayL4U = Delaunay2DFilterPointerType::New();
	Delaunay2DFilterPointerType delaunayL4L = Delaunay2DFilterPointerType::New();
	PointsPointerType L4UPoints = PointsPointerType::New();
	PointsPointerType L4LPoints = PointsPointerType::New();
	PolyDataPointerType L4USampling = PolyDataPointerType::New();
	PolyDataPointerType L4LSampling = PolyDataPointerType::New();
	L4UPoints->SetDataTypeToFloat(), L4LPoints->SetDataTypeToFloat();
	L4USampling->Allocate(1), L4LSampling->Allocate(1);
	
	int maxRadius = 3000;
	for(int x = -maxRadius; x <= maxRadius; x += 200)
		for(int y = -maxRadius; y <= maxRadius; y += 200)
		{
			if(sqrt(x*x + y*y) > maxRadius)
				continue;
			
			double samplePt[3];
			bool inside = 1;
			samplePt[0] = x, samplePt[1] = y, samplePt[2] = -0.2*sqrt(x*x + y*y);	// make sure it's always below pia...
			for(int ii = 0; ii < hullPlanes.size(); ++ii)
			{
				if(hullPlanes[ii]->EvaluateFunction(samplePt) > 0)
				{
					inside = 0;
					break;
				}
			}
			if(inside)
				continue;
			
			//now we're sure that this is a spot for a new axis
			std::multimap< double, double * > scoreMap = barrelAxisScores(avgPia, samplePt, 0.5);
			// minimum score necessary
			if(scoreMap.rbegin()->first < 0.7)
				continue;
			
			double * tmpAxis = scoreMap.rbegin()->second;
			vtkMath::Normalize(tmpAxis);
			double tmpRAxis[3];
			// make sure axis points where we think it does...
			if(tmpAxis[2] < 0)
				tmpAxis[0] = -tmpAxis[0], tmpAxis[1] = -tmpAxis[1], tmpAxis[2] = -tmpAxis[2];
			tmpRAxis[0] = -tmpAxis[0], tmpRAxis[1] = -tmpAxis[1], tmpRAxis[2] = -tmpAxis[2];
			int closestBarrelID = closestBarrel(samplePt, avgCenters);
			double * layerThickness = computeLayerThicknesses(closestBarrelID, avgColumns, avgBarrels);
			
			piaIntersector->intersectLine(tmpAxis, samplePt);
			double * intersection = piaIntersector->getLastIntersectPoint();
			if(!intersection)
				continue;
			
			double L4UPt[3], L4LPt[3], axisPt1[3], axisPt2[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				L4UPt[ii] = intersection[ii] + layerThickness[0]*tmpRAxis[ii];
				L4LPt[ii] = intersection[ii] + (layerThickness[0] + layerThickness[1])*tmpRAxis[ii];
				axisPt1[ii] = 0.5*(L4UPt[ii] + L4LPt[ii]) + 600*tmpAxis[ii];
				axisPt2[ii] = 0.5*(L4UPt[ii] + L4LPt[ii]) - 600*tmpAxis[ii];
			}
			
			zAxes->addLine(axisPt1, axisPt2, ZAxis);
			L4UPoints->InsertNextPoint(L4UPt);
			L4LPoints->InsertNextPoint(L4LPt);
			
			delete [] layerThickness;
			std::multimap< double, double * >::iterator scoreMapIt;
			for(scoreMapIt = scoreMap.begin(); scoreMapIt != scoreMap.end(); ++scoreMapIt)
				delete [] scoreMapIt->second;
			scoreMap.clear();
		}
	
	L4USampling->SetPoints(L4UPoints);
	L4LSampling->SetPoints(L4LPoints);
	delaunayL4U->SetInput(L4USampling);
	delaunayL4L->SetInput(L4LSampling);
	delaunayL4U->Update();
	delaunayL4L->Update();
	L4Upper->DeepCopy(delaunayL4U->GetOutput());
	L4Lower->DeepCopy(delaunayL4L->GetOutput());
	
	delete piaIntersector;
};

/******************************************************************************/
/*returns pairs of (double score, double * vec) where score is the score of   */
/*the potential local z-axis in direction of vec from centroid to the Pia     */
/******************************************************************************/
std::multimap< double, double * > barrelAxisScores(PolyDataPointerType piaSurface, double * centroid, double alpha)
{
	if(centroid != NULL)
	{
		if(piaSurface->GetNumberOfCells())
		{
			std::vector< double > dotP;
			std::vector< double > dist;
			std::vector< double * > axis;
			std::multimap< double, double * > score;
			for(int ii = 0; ii < piaSurface->GetNumberOfCells(); ++ii)
			{
				CellPointerType currentCell = piaSurface->GetCell(ii);
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
					tmpAxis[jj] = centerPoint[jj] - centroid[jj];
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
				tmpDotP = fabs(tmpDotP);
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
};

int closestBarrel(double samplePt[3], std::map< int, double * > avgCenters)
{
	int closestID = 0;
	double minDist = 1E08;
	std::map< int, double * >::const_iterator avgCenterIt;
	for(avgCenterIt = avgCenters.begin(); avgCenterIt != avgCenters.end(); ++avgCenterIt)
	{
		double tmpDist = vtkMath::Distance2BetweenPoints(samplePt, avgCenterIt->second);
		if(minDist*minDist > tmpDist)
		{
			minDist = sqrt(tmpDist);
			closestID = avgCenterIt->first;
		}
	}
	return closestID;
};

/******************************************************************************/
/*returns vector with vec[0] = Pia-L4Upper distance, vec[1] = L4Upper-L4Lower */
/******************************************************************************/
double * computeLayerThicknesses(int barrel, std::map< int, Column * > avgColumns, std::map< int, Column * > avgBarrels)
{
	double * thicknessVec = new double[2];
	thicknessVec[0] = vtkMath::Distance2BetweenPoints(avgColumns[barrel]->top, avgBarrels[barrel]->top);
	thicknessVec[0] = sqrt(thicknessVec[0]);
	thicknessVec[1] = avgBarrels[barrel]->getHeight();
	return thicknessVec;
};

PolyDataPointerType readStandardBarrelField(std::map< int, Column * >& avgColumns,
			     std::map< int, Column * >& avgBarrels, std::map< int, double * >& avgAxes, std::map< int, double * >& avgCenters)
{
	const char * bfName = "/home/regger/project_src/BarrelField3D/common/average_barrel_field.am";
	std::list< int > barrelLabels;
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
	
	Reader * amReader = new Reader(bfName);
	amReader->readSpatialGraphFile(0);
	std::cout << "Loading Standard Barrel Field..." << std::endl;
	AmiraSpatialGraph * avgBarrelField = amReader->getSpatialGraph();
	PolyDataPointerType avgAxesField = PolyDataPointerType::New();
	if(avgBarrelField)
	{
		if(avgColumns.size())
			avgColumns.clear();
		if(avgBarrels.size())
			avgBarrels.clear();
		if(avgAxes.size())
			avgAxes.clear();
		if(avgCenters.size())
			avgCenters.clear();
		std::list< int >::const_iterator labelIt;
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrelField->isLabelInSpatialGraph(ID))
			{
				PolyDataPointerType currBarrel = PolyDataPointerType::New();
				if(avgBarrelField->extractLandmark(ID, currBarrel))
				{
					std::map< double, int > cellIDMap;	// sorts cell IDs by z value
					if(currBarrel->GetNumberOfCells() != 4)
					{
						std::cout << "Error! Wrong format of barrel " << ID << " in Standard Barrel Field! Check SpatialGraph file." << std::endl;
						return avgAxesField;
					}
					for(int ii = 0; ii < 4; ++ii)
						cellIDMap.insert(std::pair< double, int >(currBarrel->GetCell(ii)->GetBounds()[4], ii));
					
					std::map< double, int >::const_iterator cellIDMapIt = cellIDMap.begin();
					std::map< double, int >::const_reverse_iterator cellIDMapRIt = cellIDMap.rbegin();
					
					PolyDataPointerType columnContours = PolyDataPointerType::New();
					PolyDataPointerType barrelContours = PolyDataPointerType::New();
					IdListPointerType columnCellIds = IdListPointerType::New();
					IdListPointerType barrelCellIds = IdListPointerType::New();
					columnContours->Allocate();
					barrelContours->Allocate();
					
					columnCellIds->InsertId(0, cellIDMapRIt->second);
					columnCellIds->InsertId(1, cellIDMapIt->second);
					++cellIDMapRIt, ++cellIDMapIt;
					barrelCellIds->InsertId(0, cellIDMapRIt->second);
					barrelCellIds->InsertId(1, cellIDMapIt->second);
					columnContours->CopyCells(currBarrel, columnCellIds);
					barrelContours->CopyCells(currBarrel, barrelCellIds);
					
					double barrelTop[3], barrelBottom[3], columnTop[3], columnBottom[3], paramCenter[3];
					//parametric center
					int subID;
					double pCoords[3], * weights1, * weights2, * weights3, *weights4;
					weights1 = new double[barrelContours->GetCell(0)->GetNumberOfPoints()];
					weights2 = new double[barrelContours->GetCell(1)->GetNumberOfPoints()];
					weights3 = new double[columnContours->GetCell(0)->GetNumberOfPoints()];
					weights4 = new double[columnContours->GetCell(1)->GetNumberOfPoints()];
					barrelContours->GetCell(0)->GetParametricCenter(pCoords);
					barrelContours->GetCell(0)->EvaluateLocation(subID, pCoords, barrelTop, weights1);
					barrelContours->GetCell(1)->GetParametricCenter(pCoords);
					barrelContours->GetCell(1)->EvaluateLocation(subID, pCoords, barrelBottom, weights2);
					columnContours->GetCell(0)->GetParametricCenter(pCoords);
					columnContours->GetCell(0)->EvaluateLocation(subID, pCoords, columnTop, weights3);
					columnContours->GetCell(1)->GetParametricCenter(pCoords);
					columnContours->GetCell(1)->EvaluateLocation(subID, pCoords, columnBottom, weights4);
					
					Column * newCol = new Column(columnContours, columnTop, columnBottom);
					Column * newBarrel = new Column(barrelContours, barrelTop, barrelBottom);
					avgColumns.insert(std::pair< int, Column * >(ID, newCol));
					avgBarrels.insert(std::pair< int, Column * >(ID, newBarrel));
					double * avgAxis = new double[3], * avgCenter = new double[3];
					for(int ii = 0; ii < 3; ++ii)
					{
						avgAxis[ii] = barrelTop[ii] - barrelBottom[ii];
						avgCenter[ii] = 0.5*(barrelTop[ii] + barrelBottom[ii]);
					}
					vtkMath::Normalize(avgAxis);
					avgAxes.insert(std::pair< int, double * >(ID, avgAxis));
					avgCenters.insert(std::pair< int, double * >(ID, avgCenter));
// 					std::cout << "barrel " << int2Labels[ID] << " top @ [" << newBarrel->top[0] << "," << newBarrel->top[1] << "," << newBarrel->top[2] << "]" << std::endl;
// 					std::cout << "barrel " << int2Labels[ID] << " bottom @ [" << newBarrel->bottom[0] << "," << newBarrel->bottom[1] << "," << newBarrel->bottom[2] << "]" << std::endl;
// 					std::cout << "column " << int2Labels[ID] << " top @ [" << newCol->top[0] << "," << newCol->top[1] << "," << newCol->top[2] << "]" << std::endl;
// 					std::cout << "column " << int2Labels[ID] << " bottom @ [" << newCol->bottom[0] << "," << newCol->bottom[1] << "," << newCol->bottom[2] << "]" << std::endl;
					
					delete [] weights1, delete [] weights2, delete [] weights3, delete [] weights4;
				}
			}
		}
		if(avgBarrelField->isLabelInSpatialGraph(ZAxis))
		{
			if(!avgBarrelField->extractLandmark(ZAxis, avgAxesField))
			{
				std::cout << "Error! Could not read barrel field z axis field!" << std::endl;
			}
		}
		else
			std::flush(std::cout << "Error! Barrel field z axis field not found in SpatialGraph!" << std::endl);
		
		delete avgBarrelField;
	}
	else
	{
		std::flush(std::cout << "Error! Could not load Standard Barrel Field from file " << bfName << std::endl);
	}
	
	delete amReader;
	return avgAxesField;
};



