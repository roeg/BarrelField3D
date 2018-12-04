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

PointsPointerType assignCorrespondingPoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius);

int main( int argc , char * argv[])
{
	if(argc == 6)
	{
		const char * inputFilename1 = argv[1];
		const char * inputFilename2 = argv[2];
		const char * inputFilename3 = argv[3];
		const char * outputFilename = argv[4];
		double tolRadius = atof(argv[5]);
		std::string outName(outputFilename);
		outName += "_correspondence_radius_";
		outName += argv[5];
		std::string outName1(outName);
		outName1 += "_set1-2";
		Reader * landmarkReader1 = new Reader(inputFilename1, outName1.c_str());
		PointsPointerType matchPoints1 = landmarkReader1->readLandmarkFile();
		std::string outName2(outName);
		outName2 += "_set1-3";
		Reader * landmarkReader2 = new Reader(inputFilename2, outName2.c_str());
		PointsPointerType matchPoints2 = landmarkReader2->readLandmarkFile();
		std::string outName3(outName);
		outName3 += "_set2-3";
		Reader * landmarkReader3 = new Reader(inputFilename3, outName3.c_str());
		PointsPointerType matchPoints3 = landmarkReader3->readLandmarkFile();
		
		std::cout << "Set 1 size" << "\t" << matchPoints1->GetNumberOfPoints() << std::endl;
		std::cout << "Set 2 size" << "\t" << matchPoints2->GetNumberOfPoints() << std::endl;
		std::cout << "Set 3 size" << "\t" << matchPoints3->GetNumberOfPoints() << std::endl;
		
		PointsPointerType correspondenceSet1Set2 = assignCorrespondingPoints(matchPoints1, matchPoints2, tolRadius);
		PointsPointerType correspondenceSet1Set3 = assignCorrespondingPoints(matchPoints1, matchPoints3, tolRadius);
		PointsPointerType correspondenceSet2Set3 = assignCorrespondingPoints(matchPoints2, matchPoints3, tolRadius);
		PointsPointerType correspondenceSet1Set2Set3 = assignCorrespondingPoints(correspondenceSet1Set2, matchPoints3, tolRadius);
		
		std::string summaryName(outName);
		summaryName += "_summary.csv";
		std::ofstream summaryOutStream(summaryName.c_str());
		summaryOutStream << "Correspondence radius (microns)" << "\t" << tolRadius << std::endl;
		summaryOutStream << "Set 1 size" << "\t" << matchPoints1->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Set 2 size" << "\t" << matchPoints2->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Set 3 size" << "\t" << matchPoints3->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Correspondence set 1/2" << "\t" << correspondenceSet1Set2->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Correspondence set 1/3" << "\t" << correspondenceSet1Set3->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Correspondence set 2/3" << "\t" << correspondenceSet2Set3->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Correspondence set 1/2/3" << "\t" << correspondenceSet1Set2Set3->GetNumberOfPoints() << std::endl;
		summaryOutStream.close();
		
		landmarkReader1->writeLandmarkFile(correspondenceSet1Set2);
		landmarkReader2->writeLandmarkFile(correspondenceSet1Set3);
		landmarkReader3->writeLandmarkFile(correspondenceSet2Set3);
		std::string outName123(outName);
		outName123 += "_set1-2-3";
		Reader * landmarkReader123 = new Reader(outName123.c_str(), outName123.c_str());
		landmarkReader123->writeLandmarkFile(correspondenceSet1Set2Set3);
		delete landmarkReader1, delete landmarkReader2, delete landmarkReader3, delete landmarkReader123;
	}
	else
		std::cout << "Usage: ./LandmarkCorrespondence [landmark reference set filename] [landmark set 2 filename] [tolerance radius]" << std::endl;
	return 0;
}

PointsPointerType assignCorrespondingPoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius)
{
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
				if(std::find(assignedRefPts.begin(), assignedRefPts.end(), jj) != assignedRefPts.end())
					continue;
				
				double refPt[3];
				refPoints->GetPoint(jj, refPt);
// 				double dist = vtkMath::Distance2BetweenPoints(pt, refPt);
				double dist = (pt[0] - refPt[0])*(pt[0] - refPt[0]) + (pt[1] - refPt[1])*(pt[1] - refPt[1]);
				if(dist <= tolRadius*tolRadius && dist < minDist)
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
				if(multiplicity > 1)
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
// 						double dist = vtkMath::Distance2BetweenPoints(pt, matchPt);
						double dist = (pt[0] - matchPt[0])*(pt[0] - matchPt[0]) + (pt[1] - matchPt[1])*(pt[1] - matchPt[1]);
						if(dist <= tolRadius*tolRadius && dist < minDist)
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
				else if(multiplicity == 1)
				{
					assignedMatchPts.push_back(ii);
					assignedRefPts.push_back(candID);
				}
				else
				{
					std::cout << "Warning: candidate corresponding point without multiplicity!" << std::endl;
				}
			}
		}
		
		candidateCorrespondences.clear();
		candidateMatchPts.clear();
		int newNrPts = assignedMatchPts.size();
		if(newNrPts != oldNrPts) change = 1;
		else change = 0;
	} while(change);
	
	std::cout << "Found " << assignedMatchPts.size() << " matching points!" << std::endl;
	
	PointsPointerType correspondingPts = PointsPointerType::New();
	correspondingPts->SetDataTypeToFloat();
	int fnPtCount = 0;
	for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
	{
		if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
		{
			double * refPt = new double[3];
			matchPoints->GetPoint(ii, refPt);
			correspondingPts->InsertPoint(fnPtCount, refPt);
			++fnPtCount;
		}
	}
	
	return correspondingPts;
}




