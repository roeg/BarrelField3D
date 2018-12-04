/****************************************************************************/
/*                                                                          */
/* Program:   LandmarkCorrespondence                                        */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   calculates correspondence between two landmark sets within    */
/*            given tolerance radius                                        */
/*            returns non-corresponding landmarks (FN)                      */
/*                                                                          */
/* EMail:     Robert.Egger@mpfi.org                                         */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../common/typedefs.h"
#include "../common/amiraReader.h"
#include "../common/basics.h"

int main( int argc , char * argv[])
{
	if(argc == 4)
	{
		const char * inputFilename1 = argv[1];
		const char * inputFilename2 = argv[2];
		double tolRadius = atof(argv[3]);
		Reader * landmarkReader1 = new Reader(inputFilename1, inputFilename1);
		PointsPointerType refPoints = landmarkReader1->readLandmarkFile();
		std::string outName(inputFilename2);
		outName += "_FP_radius_";
		outName += argv[3];
		Reader * landmarkReader2 = new Reader(inputFilename2, outName.c_str());
		PointsPointerType matchPoints = landmarkReader2->readLandmarkFile();
		
		std::list< int > assignedRefPts;
		std::list< int > assignedMatchPts;
		bool change = 0;
		do
		{
// 			std::cout << "next iteration..." << std::endl;
			int oldNrPts = assignedMatchPts.size();
			std::map< int, int > candidateCorrespondences;
			std::list< int > candidateMatchPts;
			for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
			{
				if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
					continue;
				
				double pt[3];
				matchPoints->GetPoint(ii, pt);
				bool hasCorrespondence = 0;
				double minDist = 1E10;
				int minID = -1;
				for(int jj = 0; jj < refPoints->GetNumberOfPoints(); ++jj)
				{
					double refPt[3];
					refPoints->GetPoint(jj, refPt);
					if(refPt[0] > 1585)
						continue;
					double dist = vtkMath::Distance2BetweenPoints(pt, refPt);
					if(dist <= tolRadius*tolRadius && dist < minDist && std::find(assignedRefPts.begin(), assignedRefPts.end(), jj) == assignedRefPts.end())
					{
						hasCorrespondence = 1;
						minDist = dist;
						minID = jj;
					}
				}
				if(hasCorrespondence)
				{
					candidateCorrespondences.insert(std::pair< int, int >(ii, minID));
					candidateMatchPts.push_back(minID);
				}
			}
			for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
			{
				if(candidateCorrespondences.find(ii) != candidateCorrespondences.end())
				{
					int candID = candidateCorrespondences[ii];
					int multiplicity = 0;
					std::list< int >::const_iterator candidateMatchPtsIt;
					for(candidateMatchPtsIt = candidateMatchPts.begin(); candidateMatchPtsIt != candidateMatchPts.end(); ++candidateMatchPtsIt)
						if(*candidateMatchPtsIt == candID)
							++multiplicity;
					if(multiplicity != 1)
					{
// 						std::cout << "Pt " << candID << " has multiplicity " << multiplicity << std::endl;
						//pick closest of all ambiguous pts
						double pt[3];
						refPoints->GetPoint(candID, pt);
						double minDist = 1E10;
						int minID = -1;
						for(int jj = 0; jj < matchPoints->GetNumberOfPoints(); ++jj)
						{
							if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), jj) != assignedMatchPts.end())
								continue;
							
							double matchPt[3];
							matchPoints->GetPoint(jj, matchPt);
							double dist = vtkMath::Distance2BetweenPoints(pt, matchPt);
							if(dist <= tolRadius*tolRadius)
							{
								minDist = dist;
								minID = jj;
							}
						}
// 						std::cout << "\tminDist = " << sqrt(minDist) << std::endl;
// 						std::cout << "\tminID = " << minID << std::endl;
						assignedMatchPts.push_back(minID);
						assignedRefPts.push_back(candID);
					}
					else
					{
						assignedMatchPts.push_back(ii);
						assignedRefPts.push_back(candID);
					}
				}
			}
			
			candidateCorrespondences.clear();
			candidateMatchPts.clear();
			int newNrPts = assignedMatchPts.size();
			if(newNrPts != oldNrPts) change = 1;
			else change = 0;
		} while(change);
		
		PointsPointerType fnPts = PointsPointerType::New();
		fnPts->SetDataTypeToFloat();
		int fnPtCount = 0;
		for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
		{
			if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
				continue;
			
			double * refPt = new double[3];
			matchPoints->GetPoint(ii, refPt);
			fnPts->InsertPoint(fnPtCount, refPt);
			++fnPtCount;
		}
		
		landmarkReader2->writeLandmarkFile(fnPts);
		delete landmarkReader1, delete landmarkReader2;
	}
	else
		std::cout << "Usage: ./LandmarkCorrespondence [landmark reference set filename] [landmark set 2 filename] [tolerance radius]" << std::endl;
	return 0;
}