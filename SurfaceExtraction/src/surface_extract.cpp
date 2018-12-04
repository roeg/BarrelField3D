/****************************************************************************/
/*                                                                          */
/* File:      surface_extract.cpp                                           */
/*                                                                          */
/* Purpose:   this program seeks to extracts landmark contours from          */
/*            brightfield image stacks. Landmarks are extracted from        */
/*            different image types:                                        */
/*            Pia, WM & vessels: 4x air plane images                        */
/*            barrels & vessels: 40x oil z-stacks, 8x downsampling in X-Y   */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/*          some methods are imported from the programs NeuroMorph and      */
/*          NeuroCount with kind permission from the authors:               */
/*                                                                          */
/*            Marcel Oberlaender                                            */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/*            Stefan Reissl                                                 */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
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
#include "typedefs.h"
#include "section.h"
/*old Values*/
float XYSAMPLING = 0;
float ZSAMPLING = 0;
float DOWNSAMPLING = 0;

/*float XYSAMPLING = 0.311;
float ZSAMPLING = 0.6;

float XYSAMPLING = 0.395879;
float ZSAMPLING = 1;*/

float averageSomaRadius = 6.0;
float zScale = 1.5;
unsigned long BINSIZE = 0;
unsigned long BOXSIZE = 150;

void closingRadius(std::vector< float > * radiusVector);
void closingRadiusMinima(std::vector< float > * radiusVector);
void closingRadiusSlope(std::vector< float > * radiusVector);
void calculateLocalMaxima(std::vector< float > * inputVector);
int barrelCenterEstimate(std::vector< float > * snrVector);
int bottomCutoff(std::vector< float > * radiusMinima, std::vector< float > * slopeMaxima, std::vector< float > * closingRadius, int barrelCenter);
int topCutoff(std::vector< float > * radiusMinima, std::vector< float > * closingRadius, int barrelCenter);
void modifyOptFlags(std::vector< bool > * optimizeSectionFlags);

void writeBarrelAttributeFile(const char * output_file, std::vector< Contour * > barrelZContours);
int WriteBarrelSpatialGraphFile(const char * output_file, std::vector< Contour * > barrelZContours);

int main( int argc , char * argv[])
{
	if(argc == 6)
	{
		XYSAMPLING = 2.32284;
		ZSAMPLING = 1.0;
		DOWNSAMPLING = 1.0;
		
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		unsigned int start = 1;
		unsigned int stop  = 1;
		int sliceNo = atoi(argv[3]);
		int WMbegin = atoi(argv[4]);
		int SegmentPia = atoi(argv[5]);
		
		float arg1 = 0/*atof(argv[4])*/;
		float arg2 = 0/*atof(argv[5])*/;
		float arg3 = 0/*atof(argv[6])*/;
		
		Segmentation * segmentation = new Segmentation(start, stop, inputFilename, outputFilename, arg1, arg2, arg3, sliceNo);
		
		segmentation->piaWMExtraction(WMbegin, SegmentPia);
		
		delete segmentation;
	}
	if(argc == 3/*10*/)
	{
		XYSAMPLING = 0.23069;
		ZSAMPLING = 1.0;
		DOWNSAMPLING = 8;
		
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		unsigned int start = 1;
		unsigned int stop  = 1;
		int sliceNo = 0;
		float arg1 = 0/*atof(argv[4])*/;
		float arg2 = 0/*atof(argv[5])*/;
		float arg3 = 0/*atof(argv[6])*/;
		
		Segmentation * segmentation = new Segmentation(start, stop, inputFilename, outputFilename, arg1, arg2, arg3, sliceNo);
		
		segmentation->barrelFieldBvp();
		
		delete segmentation;
	}
	if(argc == 1)
	{
		XYSAMPLING = 0.23069;
		ZSAMPLING = 1.0;
		
		const char * inputFilename;
		const char * outputFilename;
		const char * markerFilename;
		const char * medianProjectFilename;
		char * inputPathname = new char[256];
// 		const char * markerFilename = argv[2];
// 		const char * medianProjectFilename = argv[3];
		char * outputPathname = new char[256];
		char * barrelID = new char[256];
		char * transformFilename = new char[256];
		
		int downSampleRate = 8/*atoi(argv[4])*/;
		DOWNSAMPLING = (float)downSampleRate;
// 		int int2 = atoi(argv[6]);
		
		int nrOfSections;
		int firstSection, lastSection;
		int centerSection;
		bool uniformSections;
		bool withOptimization;
		int * planesPerSection, * sectionThickness;
		
		std::cout << std::endl;
		std::cout << "*****************************************************************************" << std::endl;
		std::cout << "*                        Barrel Segmentation v1.0                           *" << std::endl;
		std::cout << "*                                                                           *" << std::endl;
		std::cout << "*          ****         ****         ****        ****         ****          *" << std::endl;
		std::cout << "*         ******       ******       ******      ******       ******         *" << std::endl;
		std::cout << "*         ******       ******       ******      ******       ******         *" << std::endl;
		std::cout << "*          ****         ****         ****        ****         ****          *" << std::endl;
		std::cout << "*                                                                           *" << std::endl;
		std::cout << "*                    Will probably only ever be used by RE...               *" << std::endl;
		std::cout << "*                    HAVE FUN anyways!                                      *" << std::endl;
		std::cout << "*****************************************************************************" << std::endl;
		std::cout << std::endl;
		std::cout << "Nr. of sections: ";
		std::cin >> nrOfSections;
		std::cout << "First section ID: ";
		std::cin >> firstSection;
		lastSection = firstSection + nrOfSections - 1;
		std::cout << "Starting section ID (landmarks are assumed to exist for this section): ";
		std::cin >> centerSection;
		centerSection -= firstSection;
		std::cout << "Path to section folders (Sxx): ";
		std::cin >> inputPathname;
		// 		inputPathname = "/home/NeuroMorph/NeuroMorphData/Barrel_data/240210A1_B/";
		std::cout << "Barrel ID (e.g. C1): ";
		std::cin >> barrelID;
		std::cout << "Uniform sections? yes = 1, no = 0: ";
		std::cin >> uniformSections;
		if(uniformSections)
		{
			planesPerSection = new int[nrOfSections];
			sectionThickness = new int[nrOfSections];
			int tmp;
			std::cout << "Nr. of planes per section: ";
			std::cin >> tmp;
			for(int ii = 0; ii < nrOfSections; ++ii)
				planesPerSection[ii] = tmp;
			std::cout << "Section thickness: ";
			std::cin >> tmp;
			for(int ii = 0; ii < nrOfSections; ++ii)
				sectionThickness[ii] = tmp;
		}
		else
		{
			planesPerSection = new int[nrOfSections];
			sectionThickness = new int[nrOfSections];
			int tmp;
			for(int ii = 0; ii < nrOfSections; ++ii)
			{
				std::cout << "Thickness of section " << firstSection + ii << ": ";
				std::cin >> tmp;
				sectionThickness[ii] = tmp;
				std::cout << "Starting plane of section " << firstSection + ii << ": ";
				std::cin >> tmp;
				planesPerSection[ii] = tmp;
			}
		}
		unsigned int totalNrofPlanes;
		for(int ii = 0; ii < nrOfSections; ++ii)
			totalNrofPlanes += sectionThickness[ii];
		
		std::cout << "Filename for Amira file containing section transformations: ";
		std::cin >> transformFilename;
		
		std::cout << "Output filename: ";
		std::cin >> outputPathname;
		
		std::cout << "Complete pipeline? yes = 1, no = 0: ";
		std::cin >> withOptimization;
		std::cout << std::endl;
		
		Section ** barrelSections = new Section *[nrOfSections];
		
		/******** begin test block for barrel marker propagation ********/
		bool error = 0;
		for(int ii = 0; ii < nrOfSections; ++ii)
		{
			char * inputFilename2 = new char[256];
			char * outputFilename2 = new char[256];
			char * markerFilename2 = new char[256];
			char * medianProjectFilename2 = new char[256];
			
			std::string input(inputPathname);
			input += "S%02d/Plane";
			sprintf(inputFilename2, input.c_str(), ii+firstSection);
			std::string input2(inputFilename2);
			input2 += "%03d_8xdown.bmp";
// 			input += "S%02d/S%02d";
// 			sprintf(inputFilename2, input.c_str(), ii+firstSection, ii+firstSection);
// 			std::string input2(inputFilename2);
// 			input2 += "_zstack_preproc_0%02d.tif";
			inputFilename = input2.c_str();
			int in2Size = input2.size()+1;	//+1: null-termination
			char * inFileName = new char[in2Size];
			for(int jj = 0; jj < in2Size; ++jj)
				inFileName[jj] = inputFilename[jj];
			
			std::string output(outputPathname);
			output += "_S%02d";
			sprintf(outputFilename2, output.c_str(), ii+firstSection/*, ii+firstSection*/);
// 			std::string output1(outputFilename2);
// // 			output1 += output;
// 			outputFilename = output1.c_str();
// 			int outSize = output1.size()+1;	//+1: null-termination
// 			char * outFileName = new char[outSize];
// 			for(int jj = 0; jj < outSize; ++jj)
// 				outFileName[jj] = outputFilename[jj];
			
			std::string marker(inputPathname);
			std::string barrelIDstr(barrelID);
			marker += "S%02d/S%02d_";
			marker += barrelIDstr;
			sprintf(markerFilename2, marker.c_str(), ii+firstSection, ii+firstSection);
			std::string marker2(markerFilename2);
			marker2 += "_landmarks.landmarkAscii";
			markerFilename = marker2.c_str();
			int marker2Size = marker2.size()+1;	//+1: null-termination
			char * markFileName = new char[marker2Size];
			for(int jj = 0; jj < marker2Size; ++jj)
				markFileName[jj] = markerFilename[jj];
			
			std::string median(inputPathname);
			median += "S%02d/S%02d_median_project.tif";
			sprintf(medianProjectFilename2, median.c_str(), ii+firstSection, ii+firstSection);
			std::string median2(medianProjectFilename2);
			medianProjectFilename = median2.c_str();
			int med2Size = median2.size()+1;	//+1: null-termination
			char * medFileName = new char[med2Size];
			for(int jj = 0; jj < med2Size; ++jj)
				medFileName[jj] = medianProjectFilename[jj];
			
			float arg1 = 0;
			float arg2 = 0;
			float arg3 = 1;
			unsigned int start, stop;
			if(uniformSections)
			{
				start = (planesPerSection[ii] - sectionThickness[ii])/2;
				stop  = start + sectionThickness[ii];
			}
			else
			{
				start = planesPerSection[ii];
				stop  = start + sectionThickness[ii];
			}
			int sectionNr = ii + firstSection;
			int amiraSectionNr = sectionNr - firstSection + 1;
			
			std::cout << "Input filename: " << inputFilename << std::endl;
			std::cout << "Median projection filename: " << medianProjectFilename << std::endl;
			std::cout << "Barrel landmarks filename: " << markerFilename << std::endl;
			std::cout << "Output filename: " << outputFilename2 << std::endl;
			std::cout << "Start z plane = " << start << std::endl;
			std::cout << "Stop z plane = " << stop << std::endl;
			
			barrelSections[ii] = new Section(markFileName, inFileName/*inputFilename*/, outputFilename2, transformFilename, medFileName, start, stop, sectionNr, amiraSectionNr, arg1, arg2, arg3);
			if(ii < centerSection)
			{
				barrelSections[ii]->previousSectionID = ii + 1;
// 				std::cout << "previousSectionID = " << barrelSections[ii]->previousSectionID << std::endl;
			}
			if(ii > centerSection)
			{
				barrelSections[ii]->previousSectionID = ii - 1;
// 				std::cout << "previousSectionID = " << barrelSections[ii]->previousSectionID << std::endl;
			}
			std::cout << "Reading transformation parameters from Amira file!" << std::endl;
			int tryCount = 0;
			while(!barrelSections[ii]->readAmiraSectionTransformations() && tryCount < 100)
			{
				++tryCount;
				if(tryCount == 100)
				{
					std::cout << "Error! Could not read transformation parameters for section " << sectionNr << " (Amira Slice" << amiraSectionNr << ")" << std::endl;
					std::cout << "From file " << transformFilename << "... Aborting." << std::endl;
					error = 1;
				}
			}
		}
		
// 		int centerSection = nrOfSections/2;
		bool upOK = 1, downOK = 1;
		int offset = 0;
// 		for(int ii = 0; ii < nrOfSections; ++ii)
		//go up and down from middle section at the same time
		while((upOK || downOK) && !error)
		{
			int startFactor = -1;
			if(!offset)
				startFactor = 1;
			#pragma omp parallel for
			for(int factor = startFactor; factor <= 1; factor+= 2)	// +-1*offset
			{
				int currSectionID = centerSection + factor*offset;
				std::flush(std::cout << "offset = " << offset << "; factor = " << factor << "; " << "currSectionID = " << currSectionID << std::endl);
				if(currSectionID < 0)
				{
					downOK = 0;
					continue;
				}
				if(currSectionID >= nrOfSections)
				{
					upOK = 0;
					continue;
				}
				if(!offset)
				{
					barrelSections[currSectionID]->readAmiraBarrelMarker();
				}
				else if(offset)
				{
					barrelSections[currSectionID]->calculatePropagatedBarrelMarker(barrelSections[barrelSections[currSectionID]->previousSectionID]);
				}
				float arg3 = 0;
				barrelSections[currSectionID]->writeBarrelMarkerLandmarks();
				barrelSections[currSectionID]->startFirstSegmentation(arg3);
				barrelSections[currSectionID]->setBarrelContours(!withOptimization);
			}
			++offset;
		}
		if(withOptimization && !error)
		{
			int nrOfStartSectionBarrels;
			std::vector< std::vector< Contour * > > allContours;
			for(int ii = 0; ii < nrOfSections; ++ii)
			{
	// 			std::flush(std::cout << "getting contours from section " << ii << " (" << ii+firstSection << ")" << std::endl);
				std::vector< std::vector< Contour * > > thisSectionContours = barrelSections[ii]->getBarrelContourVector();
				for(int jj = 0; jj < thisSectionContours.size(); ++jj)
					allContours.push_back(thisSectionContours[jj]);
				if(ii == centerSection)
					nrOfStartSectionBarrels = thisSectionContours[0].size();
			}
			
// 			bool ** optimizeSectionFlags = new bool *[totalNrofPlanes];
// 			for(int ii = 0; ii < totalNrofPlanes; ++ii)
// 				optimizeSectionFlags[ii] = new bool[nrOfStartSectionBarrels];
			std::vector< std::vector< bool > > optimizeSectionFlags;
			float * barrelSNRAvg = new float[nrOfStartSectionBarrels];
			float * barrelSNRAvgCount = new float[nrOfStartSectionBarrels];
			for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
			{
				std::vector< bool > barrelOptFlags;
				optimizeSectionFlags.push_back(barrelOptFlags);
				barrelSNRAvg[ii] = 0;
				barrelSNRAvgCount[ii] = 0;
			}
			for(int z = 0; z < totalNrofPlanes && z < allContours.size(); ++z)
				for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
					if(allContours[z][ii]->getValid())
					{
						barrelSNRAvg[ii] += (*(allContours[z][ii]->attributesPointer()))[7];
						++barrelSNRAvgCount[ii];
					}
			
			for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
				if(barrelSNRAvgCount[ii])
					barrelSNRAvg[ii] /= barrelSNRAvgCount[ii];
			
			for(int z = 0; z < totalNrofPlanes && z < allContours.size(); ++z)
				for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
				{
					if(allContours[z][ii]->getValid())
					{
						float innerSNR1 = (*(allContours[z][ii]->attributesPointer()))[2];
						float innerSNR2 = (*(allContours[z][ii]->attributesPointer()))[3];
						float innerSNR3 = (*(allContours[z][ii]->attributesPointer()))[4];
						float totalSNR = (*(allContours[z][ii]->attributesPointer()))[7];
						float radius = (*(allContours[z][ii]->attributesPointer()))[6];
						float circSNRSlope = (innerSNR3 - innerSNR1)/40;
						float expectedSNR = circSNRSlope*radius + innerSNR1 - circSNRSlope*45;
						
						if(totalSNR < barrelSNRAvg[ii] && totalSNR/expectedSNR < 1.15 && circSNRSlope < 0.005)
							optimizeSectionFlags[ii].push_back(1);
						else
							optimizeSectionFlags[ii].push_back(0);
					}
				}
			
			for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
				modifyOptFlags(&(optimizeSectionFlags[ii]));
			
			for(int z = 0; z < totalNrofPlanes && z < allContours.size(); ++z)
			{
				for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
				{
					if(allContours[z][ii]->getValid())
					{
						allContours[z][ii]->setOptimizeFlag(optimizeSectionFlags[ii][z]);
					}
				}
			}
			
			unsigned int planeCount = 0;
			for(int ii = 0; ii < nrOfSections; ++ii)
			{
// 				std::flush(std::cout << "Setting contours in section " << ii << " (" << ii+firstSection << ")" << std::endl);
				std::vector< std::vector< Contour * > > thisSectionContours;
				for(int jj = 0; jj < sectionThickness[ii]+1; ++jj)
				{
					thisSectionContours.push_back(allContours[planeCount]);
					++planeCount;
				}
				barrelSections[ii]->setBarrelContours(thisSectionContours);
			}
			
			int sectionCount = 0;
			while(sectionCount < nrOfSections)
			{
				#pragma omp parallel for
				for(int ii = 0; ii < 20; ++ii)
				{
					if(sectionCount + ii < nrOfSections)
					{
						std::flush(std::cout << "Starting optimized segmentation in section " << ii << " (" << sectionCount+ii+firstSection << ")" << std::endl);
						barrelSections[sectionCount+ii]->startOptimizedSegmentation();
						barrelSections[sectionCount+ii]->prepareForOutput();
					}
				}
				sectionCount += 20;
			}
			
			/********  start barrel cut-off  ********/
			
			std::vector< std::vector< Contour * > > finalContours;
			for(int ii = 0; ii < nrOfSections; ++ii)
			{
// 				std::flush(std::cout << "Getting contours from section " << ii << " (" << ii+firstSection << ")" << std::endl);
				std::vector< std::vector< Contour * > > thisSectionContours = barrelSections[ii]->getBarrelContourVector();
				for(int jj = 0; jj < thisSectionContours.size(); ++jj)
					finalContours.push_back(thisSectionContours[jj]);
			}
			
			std::string cutFilename(outputPathname);
			cutFilename += "_barrel_z_borders.csv";
			std::ofstream BarrelBorderData( cutFilename.c_str() );
			BarrelBorderData << "Barrel ID\ttop cutoff\tbottom cutoff" << std::endl;

			std::flush(std::cout << "Writing attribute files" << std::endl);
			for(int ii = 0; ii < nrOfStartSectionBarrels; ++ii)
			{
				std::flush(std::cout << "Cutting & writing barrel " << ii << " to " << cutFilename.c_str() << std::endl);
				std::vector< float > finalRadius;
				std::vector< float > regularSNR;
				for(int z = 0; z < finalContours.size(); ++z)
					if(finalContours[z][ii]->getValid())
					{
						if(finalContours[z][ii]->getOptimizeFlag() && finalContours[z][ii]->attributesPointer()->size() > 10)
							finalRadius.push_back((*(finalContours[z][ii]->attributesPointer()))[10]);	//opt radius
						else
							finalRadius.push_back((*(finalContours[z][ii]->attributesPointer()))[6]);		//reg radius
						regularSNR.push_back((*(finalContours[z][ii]->attributesPointer()))[7]);			//reg SNR
					}
				
				int approxBarrelCenter = barrelCenterEstimate(&regularSNR);
// 				std::cout << "barrel center = " << approxBarrelCenter << std::endl;
// 				std::flush(std::cout << "closing radius" << std::endl);
				closingRadius(&finalRadius);
				std::vector< float > finalRadiusSlope(finalRadius);
				std::vector< float > finalRadiusClosing(finalRadius);
// 				std::flush(std::cout << "closing radius slope" << std::endl);
				closingRadiusSlope(&finalRadiusSlope);
// 				std::flush(std::cout << "closing radius minima" << std::endl);
				closingRadiusMinima(&finalRadius);
// 				std::flush(std::cout << "local slope maxima" << std::endl);
				calculateLocalMaxima(&finalRadiusSlope);
// 				std::flush(std::cout << "local minima" << std::endl);
				calculateLocalMaxima(&finalRadius);
// 				std::flush(std::cout << "cutting top" << std::endl);
				int topCut = topCutoff(&finalRadius, &finalRadiusClosing, approxBarrelCenter);
// 				std::flush(std::cout << "cutting bottom" << std::endl);
				int bottomCut = bottomCutoff(&finalRadius, &finalRadiusSlope, &finalRadiusClosing, approxBarrelCenter);
				std::vector< Contour * > thisBarrelOutput;
				int validZ = 0;
				for(int z = 0; z < finalContours.size(); ++z)
					if(finalContours[z][ii]->getValid())
					{
						if(validZ >= topCut && validZ <= bottomCut)
							finalContours[z][ii]->setInsideBarrel(1);
						else
							finalContours[z][ii]->setInsideBarrel(0);
						thisBarrelOutput.push_back(finalContours[z][ii]);
						++validZ;
					}
				
				for(int z = 0; z < thisBarrelOutput.size(); ++z)
					thisBarrelOutput[z]->prepareForOutput();
				
				writeBarrelAttributeFile(outputPathname, thisBarrelOutput);
				
				int zCenterSection = finalContours.size()/2;
				BarrelBorderData << finalContours[zCenterSection][ii]->getBarrelID() << "\t" << topCut << "\t" << bottomCut << std::endl;
				
				WriteBarrelSpatialGraphFile(outputPathname, thisBarrelOutput);
			}
			BarrelBorderData.close();
			
// 			planeCount = 0;
// 			for(int ii = 0; ii < nrOfSections; ++ii)
// 			{
// 				std::flush(std::cout << "setting contours in section " << ii << " (" << ii+firstSection << ")" << std::endl);
// 				std::vector< std::vector< Contour * > > thisSectionContours;
// 				for(int jj = 0; jj < planesPerSection[ii]; ++jj)
// 				{
// 					thisSectionContours.push_back(finalContours[planeCount]);
// 					++planeCount;
// 				}
// 				barrelSections[ii]->setBarrelContours(thisSectionContours);
// 				barrelSections[ii]->writeOptimizedBarrels();
// 			}
			
			/********   end barrel cut-off   ********/
		}
		/********  end test block for barrel marker propagation  ********/
		
		for(int ii = 0; ii < nrOfSections; ++ii)
			delete barrelSections[ii];
		delete [] barrelSections;
	}
	if(argc == 8)
	{
		XYSAMPLING = 0.23069;
		ZSAMPLING = 1.0;
		
		const char * inputFilename = argv[1];
		const char * markerFilename = argv[2];
		const char * medianProjectFilename = argv[3];
		const char * outputFilename = argv[4]/*argv[2]*/;
		const char * barrelID;
		
		unsigned int start = atoi(argv[5]);
		unsigned int stop  = atoi(argv[6]);
// 		unsigned int start = 1;
// 		unsigned int stop  = 1;
		int sliceNo = atoi(argv[7]);
		int downSampleRate = 8/*atoi(argv[4])*/;
		DOWNSAMPLING = (float)downSampleRate;
// 		int int2 = atoi(argv[6]);
// 		
		float arg1 = 0/*atof(argv[7])*/;
		float arg2 = 0/*atof(argv[8])*/;
		float arg3 = 0/*atof(argv[7])*/;
			
		Section * barrelSection = new Section(markerFilename, inputFilename, outputFilename, outputFilename, medianProjectFilename, start, stop, sliceNo, sliceNo, arg1, arg2, arg3);
		barrelSection->readAmiraSectionTransformations();
// 		barrelSection->startFirstSegmentation(arg3);
// 		barrelSection->setBarrelContours();
// // 		barrelSection->writeRegularBarrels();
		delete barrelSection;
		
// 		Segmentation * segmentation = new Segmentation(start, stop, inputFilename, outputFilename, medianProjectFilename, downSampleRate, sliceNo, arg1, arg2, arg3);
// 		
// 		segmentation->barrelExtraction(markerFilename, arg3);
// 		
// 		delete segmentation;
	}
// 	if(argc == 5)
// 	{
// 		XYSAMPLING = 0.23069;
// 		ZSAMPLING = 1.0;
// 		
// 		const char * inputFilename = ""/*argv[1]*/;
// 		const char * dataFilename = argv[1];
// 		const char * medianProjectFilename = ""/*argv[3]*/;
// 		const char * outputFilename = argv[2];
// 		
// 		unsigned int start = 0/*atoi(argv[5])*/;
// 		unsigned int stop  = 0/*atoi(argv[6])*/;
// 		// 		unsigned int start = 1;
// 		// 		unsigned int stop  = 1;
// 		int sliceNo = 1/*atoi(argv[3])*/;
// 		int downSampleRate = 8/*atoi(argv[4])*/;
// 		DOWNSAMPLING = (float)downSampleRate;
// 		// 		int int2 = atoi(argv[6]);
// 		// 		
// 		float arg1 = atof(argv[3]);
// 		float arg2 = atof(argv[4]);
// 		float arg3 = 0/*atof(argv[5])*/;
// 		
// 		
// 		Segmentation * segmentation = new Segmentation(start, stop, inputFilename, outputFilename, medianProjectFilename, downSampleRate, sliceNo, arg1, arg2, arg3);
// 		
// 		segmentation->barrelAttributeProcessing(dataFilename);
// 		
// 		delete segmentation;
// 	}
// 	//getchar();
	return 0;
}

void closingRadius(std::vector< float > * radiusVector)
{
	if(radiusVector->size())
	{
		float * radiusAvg = new float[radiusVector->size()];
		float * radiusDilate = new float[radiusVector->size()];
		float * radiusErode = new float[radiusVector->size()];
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			float tmpAvg = 0;
			float count = 0;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < radiusVector->size())
				{
					++count;
					tmpAvg += (*radiusVector)[ii+jj];
				}
			if(count)
				tmpAvg /= count;
			radiusAvg[ii] = tmpAvg;
		}
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			float tmpMax = 0;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < radiusVector->size())
					if(radiusAvg[ii+jj] > tmpMax)
						tmpMax = radiusAvg[ii+jj];
			radiusDilate[ii] = tmpMax;
		}
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			float tmpMin = 1E06;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < radiusVector->size())
					if(radiusDilate[ii+jj] < tmpMin)
						tmpMin = radiusDilate[ii+jj];
			radiusErode[ii] = tmpMin;
		}
		for(int ii = 0; ii < radiusVector->size(); ++ii)
			(*radiusVector)[ii] = radiusErode[ii];
		
		delete [] radiusAvg, delete [] radiusDilate, delete [] radiusErode;
	}
};

void closingRadiusMinima(std::vector< float > * radiusVector)
{
	if(radiusVector->size())
	{
		bool * minimumFlag = new bool[radiusVector->size()];
		float * radiusMinima = new float[radiusVector->size()];
		float maxRadius = 0;
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			if((*radiusVector)[ii] > maxRadius)
				maxRadius = (*radiusVector)[ii];
			bool isMin = 1;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < radiusVector->size() && jj != 0)
					if((*radiusVector)[ii+jj] < (*radiusVector)[ii])
					{
						isMin = 0;
						break;
					}
			minimumFlag[ii] = isMin;
		}
		float factor =  std::min(maxRadius, 120.0f);
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			if(minimumFlag[ii])
				radiusMinima[ii] = factor*(1 - (*radiusVector)[ii]/factor);
			else
				radiusMinima[ii] = 0;
		}
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			float tmpAvg = 0;
			float count = 0;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < radiusVector->size())
				{
					++count;
					tmpAvg += radiusMinima[ii+jj];
				}
			if(count)
				tmpAvg /= count;
			(*radiusVector)[ii] = tmpAvg;
		}
		
		delete [] minimumFlag, delete [] radiusMinima;
	}
};

void closingRadiusSlope(std::vector< float > * radiusVector)
{
	if(radiusVector->size())
	{
		float * slope = new float[radiusVector->size()];
		for(int ii = 0; ii < radiusVector->size(); ++ii)
		{
			int low = ii - 5;
			int high = ii + 5;
			while(low < 0)
				++low;
			while(high >= radiusVector->size())
				--high;
			int width = high - low + 1;
			slope[ii] = ((*radiusVector)[high] - (*radiusVector)[low])/(float)width;
		}
		for(int ii = 0; ii < radiusVector->size(); ++ii)
			(*radiusVector)[ii] = slope[ii];
		
		delete [] slope;
	}
};

void calculateLocalMaxima(std::vector< float > * inputVector)
{
	if(inputVector->size())
	{
		bool * maxFlag = new bool[inputVector->size()];
		for(int ii = 0; ii < inputVector->size(); ++ii)
		{
			bool isMax = 1;
			for(int jj = -2; jj <= 2; ++jj)
				if(ii + jj >= 0 && ii + jj < inputVector->size() && jj != 0)
					if((*inputVector)[ii+jj] > (*inputVector)[ii])
					{
						isMax = 0;
						break;
					}
			maxFlag[ii] = isMax;
		}
		for(int ii = 0; ii < inputVector->size(); ++ii)
			if(maxFlag[ii])
			{
				int jj = 1;
				int width = 1;
				while(ii + jj < inputVector->size())
				{
					if(maxFlag[ii+jj])
					{
						++width;
						++jj;
					}
					else
						break;
				}
				if(width > 1)
					for(int kk = 0; kk < width - 1; ++kk)
						maxFlag[ii+kk] = 0;
			}
		for(int ii = 0; ii < inputVector->size(); ++ii)
			if(!maxFlag[ii])
				(*inputVector)[ii] = 0;
		
		delete [] maxFlag;
	}
};

int barrelCenterEstimate(std::vector< float > * snrVector)
{
	if(snrVector->size())
	{
		float avgSNR = 0;
		float avgCount = 0;
		for(int ii = 0; ii < snrVector->size(); ++ii)
			if((*snrVector)[ii])
			{
				avgSNR += (*snrVector)[ii];
				++avgCount;
			}
		if(avgCount)
			avgSNR /= avgCount;
		int centerEstimate = 0;
		int centerEstimateCount = 0;
		for(int ii = 0; ii < snrVector->size(); ++ii)
			if((*snrVector)[ii] >= avgSNR)
			{
				centerEstimate += ii;
				++centerEstimateCount;
			}
		if(centerEstimateCount)
		{
			centerEstimate = (float)centerEstimate/(float)centerEstimateCount;
			return centerEstimate;
		}
		return snrVector->size()/2;
	}
	return 0;
};

int bottomCutoff(std::vector< float > * radiusMinima, std::vector< float > * slopeMaxima, std::vector< float > * closingRadius, int barrelCenter)
{
	if(radiusMinima->size() && slopeMaxima->size() && radiusMinima->size() == slopeMaxima->size())
	{
		int bottomEnd = radiusMinima->size() - 1, bottomSlopeEnd = 0;
		float bottomSlopeAvg = 0, bottomSlopeAvgCount = 0;
		float bottomMinimaAvg = 0, bottomMinimaAvgCount = 0;
		for(int ii = barrelCenter; ii < radiusMinima->size(); ++ii)
			if((*slopeMaxima)[ii] > 0)
			{
				bottomSlopeAvg += (*slopeMaxima)[ii];
				++bottomSlopeAvgCount;
			}
		if(bottomSlopeAvgCount)
			bottomSlopeAvg /= bottomSlopeAvgCount;
		// 		std::cout << "bottom slope avg = " << bottomSlopeAvg << std::endl;
		for(int ii = barrelCenter; ii < radiusMinima->size(); ++ii)
			if((*radiusMinima)[ii] > 10)
			{
				bottomMinimaAvg += (*radiusMinima)[ii];
				++bottomMinimaAvgCount;
			}
		if(bottomMinimaAvgCount)
			bottomMinimaAvg /= bottomMinimaAvgCount;
		// 		std::cout << "bottom minima avg = " << bottomMinimaAvg << std::endl;
	// 		std::cout << "bottom minima avg count = " << bottomMinimaAvgCount << std::endl;
		for(int ii = barrelCenter; ii < slopeMaxima->size(); ++ii)
			if((*slopeMaxima)[ii] > bottomSlopeAvg && (*slopeMaxima)[ii] > 1)
			{
				bottomSlopeEnd = ii;
				bool foundBottom = 0;
				for(int jj = bottomSlopeEnd; jj >= barrelCenter; --jj)
					if((*radiusMinima)[jj])
					{
						if((*radiusMinima)[jj] > bottomMinimaAvg)
						{
							foundBottom = 1;
							bottomEnd = jj;
						}
						break;
					}
				if(foundBottom)
					break;
			}
		// 		std::cout << "curr bottom end = " << bottomEnd << std::endl;
		int otherMinimum = 0;
		for(int ii = bottomEnd + 1; ii < bottomEnd + 50 && ii < radiusMinima->size(); ++ii)
			if((*radiusMinima)[ii] > (*radiusMinima)[bottomEnd])
			{
				otherMinimum = ii;
				break;
			}
		if(otherMinimum)
		{
			float betweenMax = 0;
			for(int ii = bottomEnd; ii < otherMinimum; ++ii)
				if((*closingRadius)[ii] > betweenMax)
					betweenMax = (*closingRadius)[ii];
			if(abs(betweenMax - (*closingRadius)[bottomEnd]) < 1.5*abs((*closingRadius)[bottomEnd] - (*closingRadius)[otherMinimum]))
				bottomEnd = otherMinimum;
		}
		return bottomEnd;
	}
	return -1;
};

int topCutoff(std::vector< float > * radiusMinima, std::vector< float > * closingRadius, int barrelCenter)
{
	if(radiusMinima->size() && closingRadius->size() && radiusMinima->size() == closingRadius->size())
	{
		int topEnd = 0;
		float topMinimaAvg = 0, topMinimaAvgCount = 0;
		for(int ii = barrelCenter; ii >= 0; --ii)
			if((*radiusMinima)[ii])
			{
				topMinimaAvg += (*radiusMinima)[ii];
				++topMinimaAvgCount;
			}
		if(topMinimaAvgCount)
			topMinimaAvg /= topMinimaAvgCount;
		// 		std::cout << "top minima avg = " << topMinimaAvg << std::endl;
	// 		std::cout << "top minima avg count = " << topMinimaAvgCount << std::endl;
		std::vector< float > topMinimaList;
		for(int ii = 0; ii <= barrelCenter; ++ii)
			if((*radiusMinima)[ii] > topMinimaAvg)
				topMinimaList.push_back((*radiusMinima)[ii]);
		
		float topRadiusMedian;
		if(topMinimaList.size())
		{
			int n = topMinimaList.size()/2 - 1;	//first element already counted in next line
			if(n < 0)
				n = 0;
			std::vector< float >::iterator medianIt = topMinimaList.begin() + n;
			std::nth_element(topMinimaList.begin(), medianIt, topMinimaList.end());
			topRadiusMedian = *medianIt;
		}
		else
			topRadiusMedian = topMinimaAvg;
// 		std::cout << "top median = " << topRadiusMedian << std::endl;
		bool foundTop = 0;
		for(int ii = barrelCenter; ii >= 0; --ii)
			if((*radiusMinima)[ii] > 1.1*topRadiusMedian)
			{
				foundTop = 1;
				topEnd = ii;
				break;
			}
		if(!foundTop)
			for(int ii = barrelCenter; ii >= 0; --ii)
				if((*radiusMinima)[ii] > topRadiusMedian)
				{
					topEnd = ii;
					break;
				}
		
		int otherMinimum = 0;
		for(int ii = topEnd - 1; ii > topEnd - 50 && ii >= 0; --ii)
			if((*radiusMinima)[ii] > (*radiusMinima)[topEnd])
			{
				otherMinimum = ii;
				break;
			}
		if(otherMinimum)
		{
			float betweenMax = 0;
			for(int ii = topEnd; ii > otherMinimum; --ii)
				if((*closingRadius)[ii] > betweenMax)
					betweenMax = (*closingRadius)[ii];
			if(abs(betweenMax - (*closingRadius)[topEnd]) < 1.5*abs((*closingRadius)[topEnd] - (*closingRadius)[otherMinimum]))
				topEnd = otherMinimum;
		}
		
		return topEnd;
	}
	return -1;
};

void modifyOptFlags(std::vector< bool > * optimizeSectionFlags)
{
	// 			set small blocks of optimize sections to not optimize (remove interspersed planes that otherwise would be flagged as optimize)
	// 			std::flush(std::cout << "setting single non-opt to opt" << std::endl);
	for(int z = 0; z < optimizeSectionFlags->size(); ++z)
	{
		if(!(*optimizeSectionFlags)[z])
		{
			int low = z;
			int high = z;
			while(low > 0)
			{
				--low;
				if((*optimizeSectionFlags)[low])
				{
// 					if(low-1 >= 0)
// 						if((*optimizeSectionFlags)[low-1])
					break;
				}
			}
			while(high < optimizeSectionFlags->size() - 1)
			{
				++high;
				if((*optimizeSectionFlags)[high])
				{
// 					if(high + 1 < optimizeSectionFlags->size())
// 						if((*optimizeSectionFlags)[high+1])
					break;
				}
			}
			int width = high - low;
			if(width <= 4)
				(*optimizeSectionFlags)[z] = 1;
		}
	}
	unsigned int beginMiddle = 0, endMiddle = optimizeSectionFlags->size() - 1;
	std::vector< int > optSectionIDs, nonOptSectionIDs;
	std::vector< int > optSectionSizes, nonOptSectionSizes;
	std::vector< bool > keepOptFlags, keepNonOptFlags;
	unsigned int optCount = 0, nonOptCount = 0;
	bool inOptSection = (*optimizeSectionFlags)[0];
	unsigned int tmpBegin = 0, tmpEnd = 0;
	bool foundBegin = 0, foundEnd = 0;
	// 	std::flush(std::cout << "setting smaller opt blocks to non-opt" << std::endl);
	for(int z = 0; z < optimizeSectionFlags->size(); ++z)
	{
		if(!optCount && (*optimizeSectionFlags)[z])
		{
			++optCount;
			tmpBegin = z;
			optSectionIDs.push_back(z);
		}
		else if(optCount && (*optimizeSectionFlags)[z])
			++optCount;
		if(optCount && !(*optimizeSectionFlags)[z])
		{
			if(z <= optimizeSectionFlags->size()/2)	//look to the left
			{
				int tmpZ = tmpBegin;
				int tmpLeftSize = 0;
				while(tmpZ > 0)
				{
					--tmpZ;
					if((*optimizeSectionFlags)[tmpZ])
						break;
					++tmpLeftSize;
				}
				if(tmpLeftSize > optCount && optCount <= 25)
					keepOptFlags.push_back(0);
				else
					keepOptFlags.push_back(1);
			}
			else	//look to the right
			{
				int tmpZ = z - 1;
				int tmpRightSize = 0;
				while(tmpZ < optimizeSectionFlags->size() - 1)
				{
					++tmpZ;
					if((*optimizeSectionFlags)[tmpZ])
						break;
					++tmpRightSize;
				}
				if(tmpRightSize > optCount && optCount <= 25)
					keepOptFlags.push_back(0);
				else
					keepOptFlags.push_back(1);
			}
			optCount = 0;
		}
	}
	if(optCount)
		keepOptFlags.push_back(1);
	for(int jj = 0; jj < optSectionIDs.size() && jj < keepOptFlags.size(); ++jj)
	{
		if(keepOptFlags[jj])
			continue;
		else
		{
			int tmpCount = optSectionIDs[jj];
			while((*optimizeSectionFlags)[tmpCount] && tmpCount < optimizeSectionFlags->size())
			{
				(*optimizeSectionFlags)[tmpCount] = 0;
				++tmpCount;
			}
		}
	}
};

void writeBarrelAttributeFile(const char * output_file, std::vector< Contour * > barrelZContours)
{
// 	std::flush(std::cout << "writing attribute file" << std::endl);
	std::string output(output_file);
	output += "_ID%02d_attributes.csv";
	char * outChar = new char[256];
	sprintf(outChar, output.c_str(), barrelZContours[0]->getBarrelID());
	// 	std::flush(std::cout << output.c_str() << std::endl);
	std::ofstream AttributeFile(outChar);
	AttributeFile << "\tPlane\tvoronoi ID\tcirc1 SNR\tcirc2 SNR\tcirc3 SNR\treg. threshold\tradius\treg. SNR\toptimize flag\tbarrelArea/VoronoiArea" << std::endl;
	#ifdef DEBUG
	DebugLog << "writing attribute files" << std::endl;
	DebugLog << std::endl;
	#endif
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
		{
			// 	AttributeFile << "Plane\tregion mean" << std::endl;
			std::vector< float >::const_iterator attrIt;
// 			for(int z = 0; z < barrelZContours.size(); ++z)
// 			{
				// 		AttributeFile << z;
				int attrCount = 0;
				for(attrIt = barrelZContours[z]->readAttributes(); attrIt != barrelZContours[z]->attributesEnd(); ++attrIt, ++attrCount)
				{
					AttributeFile << "\t" << *attrIt;
					if(attrCount == 7)
						AttributeFile << "\t" << barrelZContours[z]->getOptimizeFlag();
				}
				AttributeFile << std::endl;
// 			}
		}
	AttributeFile.close();
};

/************************************************************************************/
/*WriteSpatialGraphFile(spatial_graph, filename) writes an "am"-file                */
/*method for only writing one barrel directly from contour vector                   */
/************************************************************************************/
int WriteBarrelSpatialGraphFile(const char * output_file, std::vector< Contour * > barrelZContours)
{
	std::list<std::vector<float> >::iterator contour_it;
	std::list<std::vector<float> >::iterator edge_list_contour_it;
	
	int number_of_edge_points = 0;
	int number_of_contours = 0;
	
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
		{
			if(barrelZContours[z]->edgeListPointer()->size())
			{
				++number_of_contours;
				number_of_edge_points += barrelZContours[z]->edgeListPointer()->size() + 1;
			}
		}
	
	std::string format(output_file);
	format += "_ID%02d_spatialgraph.am";
	char * outChar = new char[256];
	sprintf(outChar, format.c_str(), barrelZContours[0]->getBarrelID());
	
// 	#ifdef DEBUG
	std::cout << "Writing SpatialGraph: " << outChar << std::endl;
	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
// 	#endif
	
	std::ofstream NeuroMorphData( outChar );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatialGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	NeuroMorphData << "define VERTEX " << number_of_contours << std::endl;
	NeuroMorphData << "define EDGE " << number_of_contours  << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                          	<<std::endl;
	NeuroMorphData << "        Neuron { "                                	<<std::endl;
	NeuroMorphData << "            Dendrite {"                           	<<std::endl;
	NeuroMorphData << "                ApicalDendrite {"                 	<<std::endl;
	NeuroMorphData << "                    Color 1 0.5 0.5,"          	<<std::endl;
	NeuroMorphData << "                    Id 4 }"                     	<<std::endl;
	NeuroMorphData << "                BasalDendrite {"         		<<std::endl;
	NeuroMorphData << "                    Color 0.8 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 5 }"				<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 3 }"				<<std::endl;
	NeuroMorphData << "            Axon {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 6 }"				<<std::endl;
	NeuroMorphData << "            Soma {"					<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 7 }"				<<std::endl;
	NeuroMorphData << "            Color 1 0 0,"				<<std::endl;
	NeuroMorphData << "            Id 2 }"					<<std::endl;
	NeuroMorphData << "        Landmark {"					<<std::endl;
	NeuroMorphData << "            Pia {"					<<std::endl;
	NeuroMorphData << "                Color 0 1 0.5,"			<<std::endl;
	NeuroMorphData << "                Id 9 }"				<<std::endl;
	NeuroMorphData << "            Vessel {"				<<std::endl;
	NeuroMorphData << "                Color 1 0.5 0,"			<<std::endl;
	NeuroMorphData << "                Id 10 }"				<<std::endl;
	NeuroMorphData << "            Barrel {"				<<std::endl;
	NeuroMorphData << "                aRow {"				<<std::endl;
	NeuroMorphData << "                    A1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 13 }"			<<std::endl;
	NeuroMorphData << "                    A2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 14 }"			<<std::endl;
	NeuroMorphData << "                    A3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 15 }"			<<std::endl;
	NeuroMorphData << "                    A4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 16 }"			<<std::endl;
	NeuroMorphData << "                Color 1 0.2 0.2,"			<<std::endl;
	NeuroMorphData << "                Id 12 }"				<<std::endl;
	NeuroMorphData << "                bRow {"				<<std::endl;
	NeuroMorphData << "                    B1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 18 }"			<<std::endl;
	NeuroMorphData << "                    B2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 19 }"			<<std::endl;
	NeuroMorphData << "                    B3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 20 }"			<<std::endl;
	NeuroMorphData << "                    B4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 21 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                    Id 17 }"				<<std::endl;
	NeuroMorphData << "                cRow {"				<<std::endl;
	NeuroMorphData << "                    C1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 23 }"			<<std::endl;
	NeuroMorphData << "                    C2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 24 }"			<<std::endl;
	NeuroMorphData << "                    C3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 25 }"			<<std::endl;
	NeuroMorphData << "                    C4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 26 }"			<<std::endl;
	NeuroMorphData << "                    C5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 27 }"			<<std::endl;
	NeuroMorphData << "                    C6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 28 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                    Id 22 }"				<<std::endl;
	NeuroMorphData << "                dRow {"				<<std::endl;
	NeuroMorphData << "                    D1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 30 }"			<<std::endl;
	NeuroMorphData << "                    D2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 31 }"			<<std::endl;
	NeuroMorphData << "                    D3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 32 }"			<<std::endl;
	NeuroMorphData << "                    D4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 33 }"			<<std::endl;
	NeuroMorphData << "                    D5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 34 }"			<<std::endl;
	NeuroMorphData << "                    D6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 35 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                    Id 29 }"				<<std::endl;
	NeuroMorphData << "                eRow {"				<<std::endl;
	NeuroMorphData << "                    E1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 37 }"			<<std::endl;
	NeuroMorphData << "                    E2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 38 }"			<<std::endl;
	NeuroMorphData << "                    E3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 39 }"			<<std::endl;
	NeuroMorphData << "                    E4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 40 }"			<<std::endl;
	NeuroMorphData << "                    E5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 41 }"			<<std::endl;
	NeuroMorphData << "                    E6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 42 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 36 }"				<<std::endl;
	NeuroMorphData << "                greekRow {"				<<std::endl;
	NeuroMorphData << "                    Alpha {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 44 }"			<<std::endl;
	NeuroMorphData << "                    Beta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 45 }"			<<std::endl;
	NeuroMorphData << "                    Gamma {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 46 }"			<<std::endl;
	NeuroMorphData << "                    Delta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 47 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                    Id 43 }"				<<std::endl;
	NeuroMorphData << "                Color 0 1 0,"			<<std::endl;
	NeuroMorphData << "                Id 11 }"				<<std::endl;
	NeuroMorphData << "            WhiteMatter {"				<<std::endl;
	NeuroMorphData << "                Color 0.5 1 0.75,"			<<std::endl;
	NeuroMorphData << "                Id 48 }"				<<std::endl;
	NeuroMorphData << "            OtherBarrels {"				<<std::endl;
	NeuroMorphData << "                Color 1 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 49 }"				<<std::endl;
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"      	        <<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	NeuroMorphData << "\n@1 # Vertices xyz coordinates" 			<< std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
				NeuroMorphData << barrelZContours[z]->edgeListPointer()->front()[X_COORD] << " " << barrelZContours[z]->edgeListPointer()->front()[Y_COORD]  << " " << barrelZContours[z]->edgeListPointer()->front()[Z_COORD]  << std::endl;
		
	NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
			{
				if(barrelZContours[z]->getInsideBarrel())
					NeuroMorphData << Barrel << std::endl;
				else
					NeuroMorphData << 0 << std::endl;
			}
	
	NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
	for(int i=0; i < number_of_contours; i++)
	{
		NeuroMorphData << i << " " << i << std::endl;
	}
	
	NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
				NeuroMorphData << barrelZContours[z]->edgeListPointer()->size() + 1 <<std::endl;
		
	NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
			{
				if(barrelZContours[z]->getInsideBarrel())
					NeuroMorphData << Barrel << std::endl;
				else
					NeuroMorphData << 0 << std::endl;
			}
	
	NeuroMorphData << "@6 # Point xyz coordinates" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
			{
				for(contour_it = barrelZContours[z]->edgeListPointer()->begin(); contour_it != barrelZContours[z]->edgeListPointer()->end(); ++contour_it)
					NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
				NeuroMorphData << barrelZContours[z]->edgeListPointer()->front()[X_COORD] << " " << barrelZContours[z]->edgeListPointer()->front()[Y_COORD]  << " " << barrelZContours[z]->edgeListPointer()->front()[Z_COORD]  << std::endl;
			}
	
	NeuroMorphData << "@7 # Radius at Point" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		if(barrelZContours[z]->getValid())
			if(barrelZContours[z]->edgeListPointer()->size())
			{
				for(int jj = 0; jj < barrelZContours[z]->edgeListPointer()->size() + 1; ++jj)
				{
					// 					NeuroMorphData << 1 << std::endl;
					if(barrelZContours[z]->getOptimizeFlag())
						NeuroMorphData << (*(barrelZContours[z]->attributesPointer()))[10] << std::endl;
					else
						NeuroMorphData << (*(barrelZContours[z]->attributesPointer()))[6] << std::endl;
				}
			}
	
	return 0;
};