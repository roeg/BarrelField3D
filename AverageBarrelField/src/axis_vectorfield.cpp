/****************************************************************************/
/*                                                                          */
/* Program:   AxisVectorField                                               */
/*                                                                          */
/* File:      axis_vectorfield.cpp                                          */
/*                                                                          */
/* Purpose:   Creates axis vector field with 50micron voxels by linear      */
/*            interpolation of axes in and around SBF.                      */
/*            Best used with SBF for entire S1 as input (no edge problems)  */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
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

// # define DEBUG
#define PT_A_FIXED

ImageDataPointerType interpolateAxisField ( PolyDataPointerType coarseAxisField, PolyDataPointerType pia, PolyDataPointerType wm );
void interpolateAxisAtPoint ( PolyDataPointerType coarseAxisField, double pt[3], double axis[3] );
void computeBarycentricCoords ( double pt[3], double * a, double * b, double * c, double bCoords[3]);

int main( int argc , char * argv[])
{
	if(argc == 5)
	{
		const char * barrelFieldFilename = argv[1];
		const char * piaFilename = argv[2];
		const char * wmFilename = argv[3];
		const char * outputFilename = argv[4];
		
		Reader * barrelFieldReader = new Reader(barrelFieldFilename, barrelFieldFilename);
		Reader * piaReader = new Reader(piaFilename, piaFilename);
		Reader * wmReader = new Reader(wmFilename, wmFilename);
		
		barrelFieldReader->readSpatialGraphFile(0);
		AmiraSpatialGraph * barrelField = barrelFieldReader->getSpatialGraph();
		PolyDataPointerType pia = piaReader->readAmiraSurfaceFile();
		PolyDataPointerType wm = wmReader->readAmiraSurfaceFile();
		
		PolyDataPointerType axisFieldCoarse = PolyDataPointerType::New();
		if(!barrelField->extractLandmark(ZAxis, axisFieldCoarse))
		{
			std::cout << "Error! Need coarse axis field for interpolation" << std::endl;
			delete piaReader, delete wmReader;
			delete barrelField;
			return 0;
		}
// 		axisFieldCoarse->Print(std::cout);
		
		ImageDataPointerType axisFieldFine = interpolateAxisField(axisFieldCoarse, pia, wm);
		
		Reader * fieldWriter = new Reader(outputFilename, outputFilename);
		fieldWriter->writeVectorField(axisFieldFine);
		
		delete piaReader, delete wmReader, delete fieldWriter;
		delete barrelField;
	}
	
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		Reader * fieldReader = new Reader(inputFilename, outputFilename);
		
		fieldReader->writeVectorField(fieldReader->readVectorField());
		
		delete fieldReader;
	}
	
	return 0;
}

/****************************************************************************/
/* Creates vtkImageData with 50x50x50 micron^3 voxels and associate local   */
/* z axis vector with each voxel. Vector is interpolated from three closest */
/* z axes at center of voxel.                                               */
/****************************************************************************/
ImageDataPointerType interpolateAxisField ( PolyDataPointerType coarseAxisField, PolyDataPointerType pia, PolyDataPointerType wm )
{
	// set up voxel grid
	double piaBounds[6], wmBounds[6], maxBounds[6], SPACING[3], origin[3];
	int extent[6], offset;
	
	offset = 3;
	SPACING[0] = SPACING[1] = SPACING[2] = 25;
// 	pia->GetBounds(piaBounds);
// 	wm->GetBounds(wmBounds);
	// hard-coded bounds for re-computing with 25micron resolution:
	piaBounds[0] = -1481.05;
	piaBounds[1] = 1718.95;
	piaBounds[2] = -1095.92;
	piaBounds[3] = 1854.08;
	piaBounds[4] = -1613.32;
	piaBounds[5] = 886.678;
	wmBounds[0] = 0;
	wmBounds[1] = 0;
	wmBounds[2] = 0;
	wmBounds[3] = 0;
	wmBounds[4] = 0;
	wmBounds[5] = 0;
	for(int ii = 0; ii < 3; ++ii)
	{
		maxBounds[2*ii] = std::min(piaBounds[2*ii], wmBounds[2*ii]);
		maxBounds[2*ii+1] = std::max(piaBounds[2*ii+1], wmBounds[2*ii+1]);
		extent[2*ii] = 0;
		extent[2*ii+1] = int((maxBounds[2*ii+1] - maxBounds[2*ii])/SPACING[ii] + 1) + 2*offset;
		origin[ii] = maxBounds[2*ii] - offset*SPACING[ii];
	}
	
	ImageDataPointerType fineAxisField = ImageDataPointerType::New();
	fineAxisField->SetExtent(extent);
	fineAxisField->SetOrigin(origin);
	fineAxisField->SetSpacing(SPACING);
	fineAxisField->SetNumberOfScalarComponents(3);
	fineAxisField->SetScalarTypeToDouble();
	fineAxisField->AllocateScalars();
	fineAxisField->Update();
	
	// main loop:
	// interpolate all points on grid
	for(int x = extent[0]; x <= extent[1]; ++x)
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int z = extent[4]; z <= extent[5]; ++z)
			{
				double pt[3], axis[3];
				pt[0] = origin[0] + x*SPACING[0];
				pt[1] = origin[1] + y*SPACING[1];
				pt[2] = origin[2] + z*SPACING[2];
// 				std::cout << "Checking voxel " << x << "," << y << "," << z << std::endl;
				interpolateAxisAtPoint(coarseAxisField, pt, axis);
				double * vec = static_cast< double * >(fineAxisField->GetScalarPointer(x, y, z));
				vec[0] = axis[0];
				vec[1] = axis[1];
				vec[2] = axis[2];
			}
	
// 	double pt[3], axis[3];
// 	pt[0] = origin[0] + 0*SPACING[0];
// 	pt[1] = origin[1] + 20*SPACING[1];
// 	pt[2] = origin[2] + 4*SPACING[2];
// 	interpolateAxisAtPoint(coarseAxisField, pt, axis);
	
	fineAxisField->Update();
	return fineAxisField;
}

/****************************************************************************/
/* performs the actual interpolation and puts result into axis[3].          */
/* Coefficients for linear interpolation are the barycentric coordinates    */
/* of the voxel center in the triangle obtained by projection onto the      */
/* three nearest axes.                                                      */
/****************************************************************************/
void interpolateAxisAtPoint ( PolyDataPointerType coarseAxisField, double pt[3], double axis[3] )
{
	std::map< double, vtkIdType > axisDistanceMap;
	std::map< vtkIdType, double * > projectedPoints;
	
	#ifdef DEBUG
	std::cout << "Interpolating axis @ [" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl;
	#endif
	
	for(int ii = 0; ii < coarseAxisField->GetNumberOfCells(); ++ii)
	{
		double tmpAxis[3], pt1[3], pt2[3], * projectedPt = new double[3], t;
		coarseAxisField->GetCell(ii)->GetPoints()->GetPoint(0, pt1);
		coarseAxisField->GetCell(ii)->GetPoints()->GetPoint(1, pt2);
		for(int jj = 0; jj < 3; ++jj)
			tmpAxis[jj] = pt1[jj] - pt2[jj];
		vtkMath::Normalize(tmpAxis);
		for(int jj = 0; jj < 3; ++jj)
		{
			pt1[jj] += 5000*tmpAxis[jj];
			pt2[jj] -= 5000*tmpAxis[jj];
		}
		double dist = vtkLine::DistanceToLine(pt, pt1, pt2, t, projectedPt);
		axisDistanceMap.insert(std::pair< double, vtkIdType >(dist, ii));
		projectedPoints.insert(std::pair< vtkIdType, double * >(ii, projectedPt));
	}
	
	// because std::map sorts in ascending order,
	// the first three objects are the three
	// nearest axes
	vtkIdType axis1ID, axis2ID, axis3ID, tmpIDb = -1, tmpIDc = -1, tmpIDa = -1;
	double * projPta, * projPtb, * projPtc;
	std::map< double, vtkIdType >::const_iterator closestAxesIt;
	closestAxesIt = axisDistanceMap.begin();
	axis1ID = closestAxesIt->second;
	axis2ID = (++closestAxesIt)->second;
	axis3ID = (++closestAxesIt)->second;
	projPta = projectedPoints[axis1ID];
	projPtb = projectedPoints[axis2ID];
	projPtc = projectedPoints[axis3ID];
	
	#ifdef DEBUG
	std::cout << "Closest axes IDs: " << axis1ID << ", " << axis2ID << ", " << axis3ID << std::endl;
	std::cout << "projPta @ [" << projPta[0] << "," << projPta[1] << "," << projPta[2] << "]" << std::endl;
	std::cout << "projPtb @ [" << projPtb[0] << "," << projPtb[1] << "," << projPtb[2] << "]" << std::endl;
	std::cout << "projPtc @ [" << projPtc[0] << "," << projPtc[1] << "," << projPtc[2] << "]" << std::endl;
	#endif
	
	double bCoords[3];
	computeBarycentricCoords(pt, projPta, projPtb, projPtc, bCoords);
	
	// make sure the point is in fact inside
	// of the triangle formed by these points
	// first, change projPtc
	// if this does not succeed, change projPtb
	// and so on
	std::map< double, vtkIdType >::const_iterator closestAxesItA;
	std::map< double, vtkIdType >::const_iterator closestAxesItB;
	std::map< double, vtkIdType >::const_iterator closestAxesItC;
	
	closestAxesItA = axisDistanceMap.begin();
	#ifndef PT_A_FIXED
	while( (bCoords[0] < 0 || bCoords[1] < 0 || bCoords[2] < 0)
		&& closestAxesItA != axisDistanceMap.end())
	{
	#endif
		closestAxesItB = axisDistanceMap.begin();
		++closestAxesItB;
		while( (bCoords[0] < 0 || bCoords[1] < 0 || bCoords[2] < 0)
			&& closestAxesItB != axisDistanceMap.end())
		{
			closestAxesItC = closestAxesIt;
			while( (bCoords[0] < 0 || bCoords[1] < 0 || bCoords[2] < 0)
				&& closestAxesItC != axisDistanceMap.end())
			{
				tmpIDc = closestAxesItC->second;
				projPtc = projectedPoints[tmpIDc];
// 				#ifdef DEBUG
// 				std::cout << "Oops! Triangle not good... checking next close point @ [";
// 				std::cout << projPtc[0] << "," << projPtc[1] << "," << projPtc[2] << "]" << std::endl;
// 				#endif
				computeBarycentricCoords(pt, projPta, projPtb, projPtc, bCoords);
				++closestAxesItC;
			}
			
			tmpIDb = closestAxesItB->second;
			projPtb = projectedPoints[tmpIDb];
			computeBarycentricCoords(pt, projPta, projPtb, projPtc, bCoords);
			++closestAxesItB;
		}
		
	#ifndef PT_A_FIXED
		tmpIDa = closestAxesItA->second;
		projPta = projectedPoints[tmpIDa];
		computeBarycentricCoords(pt, projPta, projPtb, projPtc, bCoords);
		++closestAxesItA;
	}
	#endif
	
	// if nothing is found for any axis, then pt
	// was probably very close to triangle border
	// to begin with; then nothing needs to be changed
	if(closestAxesItC != axisDistanceMap.end() && tmpIDc >=0)
	{
		axis3ID = tmpIDc;
		#ifdef DEBUG
		std::cout << "Completed triangle with ID " << axis3ID << std::endl;
		#endif
	}
	if(closestAxesItB != axisDistanceMap.end() && tmpIDb >=0)
	{
		axis2ID = tmpIDb;
		#ifdef DEBUG
		std::cout << "Completed triangle with ID " << axis2ID << std::endl;
		#endif
	}
	#ifndef PT_A_FIXED
	if(closestAxesItA != axisDistanceMap.end() && tmpIDa >=0)
	{
		axis1ID = tmpIDa;
		#ifdef DEBUG
		std::cout << "Completed triangle with ID " << axis1ID << std::endl;
		#endif
	}
	#endif
	else if(closestAxesItC == axisDistanceMap.end()
		|| closestAxesItB == axisDistanceMap.end()
		|| closestAxesItA == axisDistanceMap.end())
	{
		projPtc = projectedPoints[axis3ID];
		projPtb = projectedPoints[axis2ID];
		projPta = projectedPoints[axis1ID];
		computeBarycentricCoords(pt, projPta, projPtb, projPtc, bCoords);
		std::cout << "Warning! Point @ [" << pt[0] << "," << pt[1] << "," << pt[2] << "]";
		std::cout << " very close to triangle border!" << std::endl;
		std::cout << "Closest axes IDs: " << axis1ID << ", " << axis2ID << ", " << axis3ID << std::endl;
		std::cout << "Barycentric coordinates = [" << bCoords[0] << "," << bCoords[1] << "," << bCoords[2] << "]" << std::endl;
		
		#ifdef DEBUG
		double normal[3], diff[3], projPt[3], dist;
		vtkTriangle::ComputeNormal(projPta, projPtb, projPtc, normal);
		PlanePointerType plane = PlanePointerType::New();
		plane->SetNormal(normal);
		plane->SetOrigin(projPta);
		dist = plane->EvaluateFunction(pt);
		vtkMath::MultiplyScalar(normal, dist);
		vtkMath::Subtract(pt, normal, projPt);
		dist = plane->EvaluateFunction(projPt);
		std::cout << "Distance proj. pt to triangle = " << dist << std::endl;
		#endif
	}
	
	
	// interpolate three closest axes
	// with weights given by bCoords
	double axis1[3], axis2[3], axis3[3];
	double axis1Tmp1[3], axis1Tmp2[3];
	double axis2Tmp1[3], axis2Tmp2[3];
	double axis3Tmp1[3], axis3Tmp2[3];
	
	coarseAxisField->GetCell(axis1ID)->GetPoints()->GetPoint(0, axis1Tmp1);
	coarseAxisField->GetCell(axis1ID)->GetPoints()->GetPoint(1, axis1Tmp2);
	coarseAxisField->GetCell(axis2ID)->GetPoints()->GetPoint(0, axis2Tmp1);
	coarseAxisField->GetCell(axis2ID)->GetPoints()->GetPoint(1, axis2Tmp2);
	coarseAxisField->GetCell(axis3ID)->GetPoints()->GetPoint(0, axis3Tmp1);
	coarseAxisField->GetCell(axis3ID)->GetPoints()->GetPoint(1, axis3Tmp2);
	for(int ii = 0; ii < 3; ++ii)
	{
		axis1[ii] = axis1Tmp1[ii] - axis1Tmp2[ii];
		axis2[ii] = axis2Tmp1[ii] - axis2Tmp2[ii];
		axis3[ii] = axis3Tmp1[ii] - axis3Tmp2[ii];
	}
	vtkMath::Normalize(axis1);
	vtkMath::Normalize(axis2);
	vtkMath::Normalize(axis3);
	for(int ii = 0; ii < 3; ++ii)
		axis[ii] = bCoords[0]*axis1[ii] + bCoords[1]*axis2[ii] + bCoords[2]*axis3[ii];
	vtkMath::Normalize(axis);
	
	#ifdef DEBUG
	std::cout << "Barycentric coords = [" << bCoords[0] << "," << bCoords[1] << "," << bCoords[2] << "]" << std::endl;
	std::cout << "Interpolated axis = [" << axis[0] << "," << axis[1] << "," << axis[2] << "]" << std::endl;
	#endif
	
	std::map< vtkIdType, double * >::iterator projectedPointsIt;
	for(projectedPointsIt = projectedPoints.begin(); projectedPointsIt != projectedPoints.end(); ++projectedPointsIt)
		delete [] projectedPointsIt->second;
}

// compute barycentric coordinates in 3D
// (from http://facultyfp.salisbury.edu/despickler/personal/C482Hnd.asp)
void computeBarycentricCoords ( double pt[3], double * a, double * b, double * c, double bCoords[3])
{
	// first project point on triangle plane
	double normal[3], diff[3], projPt[3], dist;
	vtkTriangle::ComputeNormal(a, b, c, normal);
	PlanePointerType plane = PlanePointerType::New();
	plane->SetNormal(normal);
	plane->SetOrigin(a);
	dist = plane->EvaluateFunction(pt);
	vtkMath::MultiplyScalar(normal, dist);
	vtkMath::Subtract(pt, normal, projPt);
	
	double n[3], na[3], nb[3], nc[3], n2;
	double vba[3], vca[3], vcb[3], vpb[3], vac[3], vpc[3], vpa[3];
	vtkMath::Subtract(b, a, vba);
	vtkMath::Subtract(c, a, vca);
	vtkMath::Subtract(c, b, vcb);
	vtkMath::Subtract(projPt, b, vpb);
	vtkMath::Subtract(a, c, vac);
	vtkMath::Subtract(projPt, c, vpc);
	vtkMath::Subtract(projPt, a, vpa);
	vtkMath::Cross(vba, vca, n);
	vtkMath::Cross(vcb, vpb, na);
	vtkMath::Cross(vac, vpc, nb);
	vtkMath::Cross(vba, vpa, nc);
	n2 = vtkMath::Dot(n, n);
	bCoords[0] = vtkMath::Dot(n, na)/n2;
	bCoords[1] = vtkMath::Dot(n, nb)/n2;
	bCoords[2] = vtkMath::Dot(n, nc)/n2;
}











