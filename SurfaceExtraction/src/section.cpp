/****************************************************************************/
/*                                                                          */
/* File:      section.cpp 						    */
/*                                                                          */
/* Purpose:   implementation of wrapper class for section-wise  	    */
/*	      barrel segmentation of barrel contours			    */
/*                                                                          */
/* Author:    Robert Egger	                                            */
/*            Max-Planck-Florida Institut		                    */
/*                                                                          */
/*                                                                          */
/* EMail: Robert.Egger@maxplanckflorida.org                                 */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include "section.h"

// #define DEBUG

Section::Section(const char * markerFilename, const char * inputFilename, const char * outputFilename, const char * transformFilename, const char * medianProjFilename, int start, int stop, int sectionID, int amiraSectionID, float arg1, float arg2, float arg3)
{
	this->sectionTranslation = new float *[4];
	this->sectionRotation = new float *[4];
	for(int ii = 0; ii < 4; ++ii)
	{
		this->sectionTranslation[ii] = new float[4];
		this->sectionRotation[ii] = new float[4];
		for(int jj = 0; jj < 4; ++jj)
		{
			this->sectionTranslation[ii][jj] = 0;
			this->sectionRotation[ii][jj] = 0;
		}
		this->sectionTranslation[ii][ii] = 1;
		this->sectionRotation[ii][ii] = 1;
	}
	this->thisSectionID = sectionID;
	this->amiraSectionID = amiraSectionID;
	this->nrOfZPlanes = stop - start + 1;
	this->start = start;
	this->stop = stop;
	downSampleRate = 8;
	this->inputFilename = inputFilename;
	this->outputFilename = outputFilename;
	this->transformFilename = transformFilename;
	this->markerFilename = markerFilename;
	this->imageProcessor = new Segmentation(start, stop, inputFilename, outputFilename, medianProjFilename, downSampleRate, sectionID, arg1, arg2, arg3);
	boundingBox = new long[6];
	boundingBox = imageProcessor->getBounds();
	#ifdef DEBUG
	std::string debugFilename(outputFilename);
	debugFilename += "_debug.log";
	DebugLog.open(debugFilename.c_str());
	#endif
};

Section::~Section()
{
	#ifdef DEBUG
	DebugLog.close();
	#endif
	if(zBarrelContoursVector.size())
	{
		for(int ii = 0; ii < zBarrelContoursVector.size(); ++ii)
			if(zBarrelContoursVector[ii].size())
			{
// 				for(int jj = 0; jj < zBarrelContoursVector[ii].size(); ++jj)
// 					delete zBarrelContoursVector[ii][jj];
				zBarrelContoursVector[ii].clear();
			}
		zBarrelContoursVector.clear();
	}
	delete marker;
// 	delete imageProcessor;
};

void Section::startFirstSegmentation(float arg)
{
	this->imageProcessor->setBarrelContours(zBarrelContoursVector);
// 	this->imageProcessor->writeVoronoiDiagram();
	this->imageProcessor->barrelExtraction(arg);
}


void Section::setBarrelContours(bool write)
{
	this->zBarrelContoursVector = this->imageProcessor->getRegularBarrelContours();
	segImageStack = ImageType::New();
	preprocImageStack = ImageType::New();
	segImageStack = this->imageProcessor->getSegmentedBarrelImages();
	preprocImageStack = this->imageProcessor->getPreprocessedBarrelImages();
// 	vessels = new BloodVesselPattern();
	vessels = this->imageProcessor->getVesselPattern();
	if(write)
	{
// 		transformToWorldCoordinates();
		writeRegularBarrels();
	}
// 	std::flush(std::cout << "delete imageProcessor object..." << std::endl);
	delete this->imageProcessor;
}

void Section::setBarrelContours(std::vector< std::vector< Contour * > > contourVec)
{
	this->zBarrelContoursVector = contourVec;
// 	std::flush(std::cout << "deleting barrel contours..." << std::endl);
// 	if(zBarrelContoursVector.size())
// 	{
// 		for(int ii = 0; ii < zBarrelContoursVector.size(); ++ii)
// 			delete zBarrelContoursVector[ii];
// 		zBarrelContoursVector.clear();
// 	}
// 	std::flush(std::cout << "inserting new barrel contours..." << std::endl);
// 	for(int ii = 0; ii < contourVec.size(); ++ii)
// 	{
// 		std::flush(std::cout << "allocating memory for new barrel contour..." << std::endl);
// // 		Contour * newContour = new Contour(contourVec[ii]);
// 		Contour * newContour = new Contour();
// 		newContour->setEdgeList(*(contourVec[ii]->getEdgeList()));
// 		newContour->setAttributeList(*(contourVec[ii]->attributesPointer()));
// 		zBarrelContoursVector.push_back(newContour);
// 	}
};

void Section::writeRegularBarrels()
{
// 	this->imageProcessor = new Segmentation(start, stop, downSampleRate, segImageStack, preprocImageStack, zBarrelContoursVector, outputFilename);
	this->imageProcessor->writeBarrelContours(/*this->zBarrelContoursVector*/);
}

void Section::writeOptimizedBarrels()
{
// 	std::flush(std::cout << "setting barrel contours in segmentation object" << std::endl);
	#ifdef DEBUG
	DebugLog << "transforming barrel contours to world coordinates" << std::endl;
	DebugLog << std::endl;
	#endif
	transformToWorldCoordinates();
	#ifdef DEBUG
	DebugLog << "barrel vector after Section::transformToWorldCoordinates()" << std::endl;
	DebugLog << "\tzBarrelContoursVector.size() = " << zBarrelContoursVector.size() << std::endl;
	for(int z = 0; z < zBarrelContoursVector.size(); ++z)
	{
		DebugLog << "\t\tzBarrelContoursVector[" << z << "].size() = " << zBarrelContoursVector[z].size() << std::endl;
		for(int ii = 0; ii < zBarrelContoursVector[z].size(); ++ii)
		{
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] edge list size = " << zBarrelContoursVector[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] valid = " << zBarrelContoursVector[z][ii]->getValid() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] barrel ID = " << zBarrelContoursVector[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
	this->barrelOptimizer->setBarrelContours(zBarrelContoursVector);
	#ifdef DEBUG
	DebugLog << "writing barrel contours in Segmentation object barrelOptimizer" << std::endl;
	#endif
// 	std::flush(std::cout << "writing barrel contours in segmentation object" << std::endl);
	this->barrelOptimizer->writeBarrelContours(/*this->zBarrelContoursVector*/);
	delete barrelOptimizer;
}

void Section::prepareForOutput()
{
	fillEmptyBarrelContours();
	transformToWorldCoordinates();
};

//expects the plane IDs for each barrel in which it should be segmented
// with the optimizer calculated elsewhere considering all sections
void Section::startOptimizedSegmentation()
{
// 	std::flush(std::cout << "constructing new barrel segmention object..." << std::endl);
	barrelOptimizer = new Segmentation(start, stop, downSampleRate, segImageStack, preprocImageStack, vessels, zBarrelContoursVector, outputFilename);
// 	std::flush(std::cout << "passing images..." << std::endl);
	barrelOptimizer->setVoronoiMap(this->marker->getVoronoiMap());
	barrelOptimizer->setDistanceMap(this->marker->getDistanceMap());
	barrelOptimizer->setBarrelPoints(this->marker->getNumberOfBarrelMarkers());
// 	std::flush(std::cout << "starting optimized region growing..." << std::endl);
	zBarrelContoursVector = barrelOptimizer->getOptimizedBarrelContours();
	#ifdef DEBUG
	DebugLog << "barrel vector after Segmentation object barrelOptimizer->getOptimizedBarrelContours()" << std::endl;
	DebugLog << "\tzBarrelContoursVector.size() = " << zBarrelContoursVector.size() << std::endl;
	for(int z = 0; z < zBarrelContoursVector.size(); ++z)
	{
		DebugLog << "\t\tzBarrelContoursVector[" << z << "].size() = " << zBarrelContoursVector[z].size() << std::endl;
		for(int ii = 0; ii < zBarrelContoursVector[z].size(); ++ii)
		{
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] edge list size = " << zBarrelContoursVector[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] valid = " << zBarrelContoursVector[z][ii]->getValid() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] barrel ID = " << zBarrelContoursVector[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
// 	barrelOptimizer->writeBarrelContours();
// 	delete barrelOptimizer;
};

/****************** deprecated; does not really help **********************/
//replace single small contours by their preceding contour if the SNR during optimization
//has not significantly improved, but has in the preceding contour
void Section::replaceNotCorrectlyOptimizedContours(int direction)
{
// 	std::flush(std::cout << "replacing non-optimized barrel contours!" << std::endl);
	if(direction == 1)	//up
	{
// 		std::flush(std::cout << "direction up!" << std::endl);
		for(int ii = zBarrelContoursVector.size() - 2; ii >= 0; --ii)
		{
			for(int jj = 0; jj < zBarrelContoursVector[ii].size(); ++jj)
				if(zBarrelContoursVector[ii][jj]->getValid())
				{
		// 			std::flush(std::cout << "flag: " << zBarrelContoursVector[ii]->getOptimizeFlag() << std::endl);
					if(zBarrelContoursVector[ii][jj]->getOptimizeFlag())
					{
						float SNRbefore, SNRafter;
						SNRbefore = (*(zBarrelContoursVector[ii][jj]->attributesPointer()))[7];
						SNRafter = (*(zBarrelContoursVector[ii][jj]->attributesPointer()))[11];
		// 				std::flush(std::cout << "SNRbefore = " << SNRbefore << std::endl);
		// 				std::flush(std::cout << "SNRafter = " << SNRafter << std::endl);
						if(SNRafter)
							if(SNRbefore/SNRafter >= 0.95)
							{
								float SNRbefore2, SNRafter2;
								SNRbefore2 = (*(zBarrelContoursVector[ii+1][jj]->attributesPointer()))[7];
								SNRafter2 = (*(zBarrelContoursVector[ii+1][jj]->attributesPointer()))[11];
								// 						std::flush(std::cout << "SNRbefore2 = " << SNRbefore2 << std::endl);
								// 						std::flush(std::cout << "SNRafter2 = " << SNRafter2 << std::endl);
								if(SNRafter2)
								{
									if(SNRbefore2/SNRafter2 < 0.95)
									{
										// 								std::flush(std::cout << "replacing non-optimized barrel contours!" << std::endl);
										zBarrelContoursVector[ii][jj]->setAttributeList(*(zBarrelContoursVector[ii+1][jj]->attributesPointer()));
										zBarrelContoursVector[ii][jj]->replaceEdgeList((zBarrelContoursVector[ii+1][jj]->edgeListPointer()), -1);
										break;
									}
								}
							}
					}
				}
		}
	}
	else if(direction == 2)	//down
	{
// 		std::flush(std::cout << "direction down!" << std::endl);
		for(int ii = 1; ii < zBarrelContoursVector.size(); ++ii)
		{
			for(int jj = 0; jj < zBarrelContoursVector[ii].size(); ++jj)
				if(zBarrelContoursVector[ii][jj]->getValid())
				{
		// 			std::flush(std::cout << "flag: " << zBarrelContoursVector[ii]->getOptimizeFlag() << std::endl);
					if(zBarrelContoursVector[ii][jj]->getOptimizeFlag())
					{
						float SNRbefore, SNRafter;
						SNRbefore = (*(zBarrelContoursVector[ii][jj]->attributesPointer()))[7];
						SNRafter = (*(zBarrelContoursVector[ii][jj]->attributesPointer()))[11];
						if(SNRafter)
							if(SNRbefore/SNRafter >= 0.95)
							{
								float SNRbefore2, SNRafter2;
								SNRbefore2 = (*(zBarrelContoursVector[ii-1][jj]->attributesPointer()))[7];
								SNRafter2 = (*(zBarrelContoursVector[ii-1][jj]->attributesPointer()))[11];
								if(SNRafter2)
								{
									if(SNRbefore2/SNRafter2 < 0.95)
									{
										// 								std::flush(std::cout << "replacing non-optimized barrel contours!" << std::endl);
										zBarrelContoursVector[ii][jj]->setAttributeList(*(zBarrelContoursVector[ii-1][jj]->attributesPointer()));
										zBarrelContoursVector[ii][jj]->replaceEdgeList((zBarrelContoursVector[ii-1][jj]->edgeListPointer()), +1);
										break;
									}
								}
							}
					}
				}
		}
	}
	#ifdef DEBUG
	DebugLog << "barrel vector after Section::replaceNotCorrectlyOptimizedContours()" << std::endl;
	DebugLog << "\tzBarrelContoursVector.size() = " << zBarrelContoursVector.size() << std::endl;
	for(int z = 0; z < zBarrelContoursVector.size(); ++z)
	{
		DebugLog << "\t\tzBarrelContoursVector[" << z << "].size() = " << zBarrelContoursVector[z].size() << std::endl;
		for(int ii = 0; ii < zBarrelContoursVector[z].size(); ++ii)
		{
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] edge list size = " << zBarrelContoursVector[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] valid = " << zBarrelContoursVector[z][ii]->getValid() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] barrel ID = " << zBarrelContoursVector[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
};

/****************************************************************************/
/*in erronous case where the contour at plane z is empty, set it equal to   */
/*a neighboring (in z) contour                                              */
/****************************************************************************/
void Section::fillEmptyBarrelContours()
{
	if(zBarrelContoursVector.size())
	{
		for(int ii = 0; ii < zBarrelContoursVector[0].size(); ++ii)
		{
			std::list< int > emptyContours;
			for(int z = 0; z < zBarrelContoursVector.size(); ++z)
				if(zBarrelContoursVector[z][ii]->edgeListPointer()->size() == 0)
					emptyContours.push_back(z);
			
			std::list< int >::const_iterator emptyContoursIt;
			for(emptyContoursIt = emptyContours.begin(); emptyContoursIt != emptyContours.end(); ++emptyContoursIt)
			{
				int z = *emptyContoursIt;
				int barrelID = ii;
				int offset = 1;
				if(z >= 0 && z < zBarrelContoursVector.size())
				{
					while(!zBarrelContoursVector[z][barrelID]->edgeListPointer()->size())
					{
						if(z - offset >= 0 && z - offset < zBarrelContoursVector.size())
							if(zBarrelContoursVector[z-offset][barrelID]->getValid() && zBarrelContoursVector[z-offset][barrelID]->edgeListPointer()->size())
							{
								zBarrelContoursVector[z][barrelID]->replaceEdgeList(zBarrelContoursVector[z-offset][barrelID]->edgeListPointer(), offset);
								break;
							}
						if(z + offset >= 0 && z + offset < zBarrelContoursVector.size())
							if(zBarrelContoursVector[z+offset][barrelID]->getValid() && zBarrelContoursVector[z+offset][barrelID]->edgeListPointer()->size())
							{
								zBarrelContoursVector[z][barrelID]->replaceEdgeList(zBarrelContoursVector[z+offset][barrelID]->edgeListPointer(), -offset);
								break;
							}
						if((z - offset < 0) && (z + offset >= zBarrelContoursVector.size()))
							break;
						++offset;
					}
				}
			}
		}
	}
	#ifdef DEBUG
	DebugLog << "barrel vector after Section::fillEmptyBarrelContours()" << std::endl;
	DebugLog << "\tzBarrelContoursVector.size() = " << zBarrelContoursVector.size() << std::endl;
	for(int z = 0; z < zBarrelContoursVector.size(); ++z)
	{
		DebugLog << "\t\tzBarrelContoursVector[" << z << "].size() = " << zBarrelContoursVector[z].size() << std::endl;
		for(int ii = 0; ii < zBarrelContoursVector[z].size(); ++ii)
		{
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] edge list size = " << zBarrelContoursVector[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] valid = " << zBarrelContoursVector[z][ii]->getValid() << std::endl;
			DebugLog << "\t\t\tzBarrelContoursVector[" << z << "][" << ii << "] barrel ID = " << zBarrelContoursVector[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
};

//only to be used for middle sections where all initial
//markers are read from an Amira landmarkAscii
//initial barrel IDs are set here and kept in other sections
void Section::readAmiraBarrelMarker()
{
	ImageType::RegionType planeRegion = imageProcessor->getInputRegion();
	ImageType::SizeType planeSize = planeRegion.GetSize();
	planeSize[2] = 1;
	planeRegion.SetSize(planeSize);
	marker = new BarrelMarker();
	marker->readBarrelMarkerFile(markerFilename, planeRegion);
	
	for(int z = 0; z < nrOfZPlanes; ++z)
	{
		std::vector< Contour * > zPlaneContours;
		for(int ii = 0; ii < marker->getNumberOfBarrelMarkers(); ++ii)
		{
			Contour * newContour = new Contour;
			newContour->setBarrelID(ii+1);
			newContour->setValid(1);
			zPlaneContours.push_back(newContour);
		}
		zBarrelContoursVector.push_back(zPlaneContours);
	}
	
	CalcImageType::Pointer voronoiMap = marker->getVoronoiMap();
	CalcImageType::Pointer distanceMap = marker->getDistanceMap();
	CalcIteratorType voronoiIter(voronoiMap, voronoiMap->GetLargestPossibleRegion());
	ConstCalcIteratorType distanceIter(distanceMap, distanceMap->GetLargestPossibleRegion());
	for(voronoiIter.GoToBegin(), distanceIter.GoToBegin(); !voronoiIter.IsAtEnd() && !distanceIter.IsAtEnd(); ++voronoiIter, ++distanceIter)
		if(distanceIter.Get() > 175)
			voronoiIter.Set(0);
	voronoiMap->Update();
	
	imageProcessor->setVoronoiMap(voronoiMap);
	imageProcessor->setDistanceMap(distanceMap);
	imageProcessor->setBarrelPoints(marker->getNumberOfBarrelMarkers());
	imageProcessor->writeVoronoiDiagram();
};

bool Section::readAmiraSectionTransformations()
{
	std::ifstream inputStream(transformFilename);
// 	this->previousSectionID = thisSectionID - 3;
// 	float ** prevSectionTranslation = new float *[4];
// 	float ** prevSectionRotation = new float *[4];
// 	float ** invPrevSectionTransform = new float *[4];
// 	float ** totalTransform = new float *[4];
// 	for(int ii = 0; ii < 4; ++ii)
// 	{
// 		prevSectionTranslation[ii] = new float[4];
// 		prevSectionRotation[ii] = new float[4];
// 		invPrevSectionTransform[ii] = new float[4];
// 		totalTransform[ii] = new float[4];
// 		for(int jj = 0; jj < 4; ++jj)
// 		{
// 			prevSectionTranslation[ii][jj] = 0;
// 			prevSectionRotation[ii][jj] = 0;
// 			invPrevSectionTransform[ii][jj] = 0;
// 			totalTransform[ii][jj] = 0;
// 		}
// 		prevSectionTranslation[ii][ii] = 1;
// 		prevSectionRotation[ii][ii] = 1;
// 		invPrevSectionTransform[ii][ii] = 1;
// 	}
	if(!inputStream.fail())
	{
		const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
		const char * numbers = "0123456789";
		const char * signs = "+-";
		const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
		const char * whitespace = "\t ";
		
		std::string currentLine;
		unsigned int line = 0;
		
		bool parameters = 0;
		bool transform = 0;
		bool correctSection = 0;
		bool correctPrevSection = 0;
		int sectionID = 0;
		unsigned int brackets = 0, transformBrackets = 0;
		unsigned int currentIndex = 0;
		
		while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
		{
			if(currentLine.size())
			{
				if(currentLine.find("@", 0) == 0)
				{
					char * tmp = new char[currentLine.size() - 1];
					currentLine.copy(tmp, currentLine.size() - 1, 1);
					currentIndex = atoi(tmp);
// 					std::cout << "Reading data section " << currentIndex << std::endl;
					delete [] tmp;
					continue;
				}
				
				if(currentIndex == 0)
				{
					std::string::size_type loc = currentLine.find("Parameters", 0);
					if(loc != std::string::npos)
					{
						parameters = 1;
						if(currentLine.find("{", 0) != std::string::npos)
							brackets = 1;
						continue;
					}
// 					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
// 						continue;
					if(parameters && currentLine.find("{", 0) != std::string::npos)
					{
						++brackets;
// 						continue;
					}
					if(parameters && currentLine.find("}", 0) != std::string::npos)
					{
						--brackets;
						if(!brackets)
							parameters = 0;
// 						continue;
					}
					
					loc = currentLine.find("TransformInfo", 0);
					if(loc != std::string::npos && parameters)
					{
						transform = 1;
						if(currentLine.find("{", 0) != std::string::npos)
							transformBrackets = 1;
							continue;
					}
// 					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
// 						continue;
					if(transform && currentLine.find("{", 0) != std::string::npos)
					{
						++transformBrackets;
// 						continue;
					}
					if(transform && currentLine.find("}", 0) != std::string::npos)
					{
						--transformBrackets;
						if(!transformBrackets)
							transform = 0;
// 						continue;
					}
					if(transform && currentLine.find("Slice", 0) != std::string::npos)
					{
						int tmp = currentLine.find_first_of(numbers);
						char tmp2[4];
						currentLine.copy(tmp2, 4, tmp);
						sectionID = atoi(tmp2);
						if(sectionID == amiraSectionID)
							correctSection = 1;
						else
							correctSection = 0;
// 						if(sectionID == thisSectionID)
// 							correctSection = 1;
// 						else
// 							correctSection = 0;
// 						if(sectionID == previousSectionID)
// 							correctPrevSection = 1;
// 						else
// 							correctPrevSection = 0;
					}
					if(transform && correctSection && currentLine.find("Transform ", 0) != std::string::npos)
					{
// 						const char * thisLine = currentLine.c_str();
// 						double * tmpCoords = new double[3];
// 						char ** endptr = new char*;
// 						tmpCoords[0] = strtod(thisLine, endptr);
// 						tmpCoords[1] = strtod(*endptr, endptr);
// 						tmpCoords[2] = strtod(*endptr, endptr);
						
// 						std::cout << "found correct section transform parameters!" << std::endl;
						unsigned int count = 0;
						std::string::size_type loc1, loc2, loc3;
						loc1 = currentLine.find_first_of(numbers, 0);
						loc2 = currentLine.find_first_of(signs, 0);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);	//ignores last value: is always 1 anyways
						while(loc2 != std::string::npos && count < 16)
						{
							char * tmp1 = new char[loc2 - loc1];
							currentLine.copy(tmp1, loc2 - loc1, loc1);
							float ftmp1 = atof(tmp1);
							sectionRotation[count%4][count/4]= ftmp1;	// amira files are columns after each other
							loc3 = loc2;
							loc1 = currentLine.find_first_of(numbers, loc3);
							loc2 = currentLine.find_first_of(signs, loc3);
							if(loc2 != std::string::npos)
								if(loc2 < loc1)
									loc1 = loc2;
							loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
							++count;
						}
						//remove numeric artifacts from z-axis:
						for(int ii = 0; ii < 2; ++ii)
						{
							sectionRotation[2][ii] = 0;
							sectionRotation[ii][2] = 0;
						}
						sectionRotation[2][2] = 1;
						//amira files store complete transformation matrices as T*R -> simple decomposition
						for(int ii = 0; ii < 3; ++ii)
						{
							sectionTranslation[ii][3] = sectionRotation[ii][3];
							sectionRotation[ii][3] = 0;
						}
						std::cout << "translation matrix:" << std::endl;
						for(int ii = 0; ii < 4; ++ii)
						{
							std::cout << "[";
							for(int jj = 0; jj < 4; ++jj)
							{
								if(jj < 3)
									std::cout << sectionTranslation[ii][jj] << ",\t";
								else
									std::cout << sectionTranslation[ii][jj];
							}
							std::cout << "]" << std::endl;
						}
						std::cout << "rotation matrix:" << std::endl;
						for(int ii = 0; ii < 4; ++ii)
						{
							std::cout << "[";
							for(int jj = 0; jj < 4; ++jj)
							{
								if(jj < 3)
									std::cout << sectionRotation[ii][jj] << ",\t";
								else
									std::cout << sectionRotation[ii][jj];
							}
							std::cout << "]" << std::endl;
						}
						inputStream.close();
						return 1;
					}
// 					if(transform && correctPrevSection && currentLine.find("Transform ", 0) != std::string::npos)
// 					{
// // 						std::cout << "found correct section transform parameters!" << std::endl;
// 						unsigned int count = 0;
// 						std::string::size_type loc1, loc2, loc3;
// 						loc1 = currentLine.find_first_of(numbers, 0);
// 						loc2 = currentLine.find_first_of(signs, 0);
// 						if(loc2 != std::string::npos)
// 							if(loc2 < loc1)
// 								loc1 = loc2;
// 						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);	//ignores last value: is always 1 anyways
// 						while(loc2 != std::string::npos && count < 16)
// 						{
// 							char * tmp1 = new char[loc2 - loc1];
// 							currentLine.copy(tmp1, loc2 - loc1, loc1);
// 							float ftmp1 = atof(tmp1);
// 							prevSectionRotation[count%4][count/4]= ftmp1;	// amira files are columns after each other
// 							loc3 = loc2;
// 							loc1 = currentLine.find_first_of(numbers, loc3);
// 							loc2 = currentLine.find_first_of(signs, loc3);
// 							if(loc2 != std::string::npos)
// 								if(loc2 < loc1)
// 									loc1 = loc2;
// 							loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
// 							++count;
// 						}
// 						//remove numeric artifacts from z-axis:
// 						for(int ii = 0; ii < 2; ++ii)
// 						{
// 							prevSectionRotation[2][ii] = 0;
// 							prevSectionRotation[ii][2] = 0;
// 						}
// 						prevSectionRotation[2][2] = 1;
// 						for(int ii = 0; ii < 3; ++ii)
// 						{
// 							prevSectionTranslation[ii][3] = prevSectionRotation[ii][3];
// 							prevSectionRotation[ii][3] = 0;
// 						}
// 						std::cout << "prev section translation matrix:" << std::endl;
// 						for(int ii = 0; ii < 4; ++ii)
// 						{
// 							std::cout << "[";
// 							for(int jj = 0; jj < 4; ++jj)
// 							{
// 								if(jj < 3)
// 									std::cout << prevSectionTranslation[ii][jj] << ",\t";
// 								else
// 									std::cout << prevSectionTranslation[ii][jj];
// 							}
// 							std::cout << "]" << std::endl;
// 						}
// 						std::cout << "prev section rotation matrix:" << std::endl;
// 						for(int ii = 0; ii < 4; ++ii)
// 						{
// 							std::cout << "[";
// 							for(int jj = 0; jj < 4; ++jj)
// 							{
// 								if(jj < 3)
// 									std::cout << prevSectionRotation[ii][jj] << ",\t";
// 								else
// 									std::cout << prevSectionRotation[ii][jj];
// 							}
// 							std::cout << "]" << std::endl;
// 						}
// // 						break;
// 					}
				}
			}
		}
	}
	
	inputStream.close();
	return 0;
// 	//invert previous section transform matrix:
// 	//use known formula for product Rotation*Translation
// 	for(int ii = 0; ii < 2; ++ii)
// 		for(int jj = 0; jj < 2; ++jj)
// 			invPrevSectionTransform[ii][jj] = prevSectionRotation[jj][ii];
// 	invPrevSectionTransform[0][3] = -1*(prevSectionRotation[0][0]*prevSectionTranslation[0][3] + prevSectionRotation[1][0]*prevSectionTranslation[1][3]);
// 	invPrevSectionTransform[1][3] = -1*(prevSectionRotation[0][1]*prevSectionTranslation[0][3] + prevSectionRotation[1][1]*prevSectionTranslation[1][3]);
// 	invPrevSectionTransform[2][3] = -1*prevSectionTranslation[2][3];
// 	
// 	std::cout << "inverse prev section transformation matrix:" << std::endl;
// 	for(int ii = 0; ii < 4; ++ii)
// 	{
// 		std::cout << "[";
// 		for(int jj = 0; jj < 4; ++jj)
// 		{
// 			if(jj < 3)
// 				std::cout << invPrevSectionTransform[ii][jj] << ",\t";
// 			else
// 				std::cout << invPrevSectionTransform[ii][jj];
// 		}
// 		std::cout << "]" << std::endl;
// 	}
// 	
// 	float ** thisSectionTotalTrans = calculateTransAfterRot();
// 	for(int ii = 0; ii < 4; ++ii)
// 		for(int jj = 0; jj < 4; ++jj)
// 			for(int kk = 0; kk < 4; ++kk)
// 				totalTransform[ii][jj] += invPrevSectionTransform[ii][kk]*thisSectionTotalTrans[kk][jj];
// 	
// 	std::cout << "total transformation matrix:" << std::endl;
// 	for(int ii = 0; ii < 4; ++ii)
// 	{
// 		std::cout << "[";
// 		for(int jj = 0; jj < 4; ++jj)
// 		{
// 			if(jj < 3)
// 				std::cout << totalTransform[ii][jj] << ",\t";
// 			else
// 				std::cout << totalTransform[ii][jj];
// 		}
// 		std::cout << "]" << std::endl;
// 	}
// 	
// 	PointSetType::Pointer transformedMarker = PointSetType::New();
// 	for(int ii = 0; ii < marker->getNumberOfBarrelMarkers(); ++ii)
// 	{
// 		PointType oldPoint, newPoint;
// 		if(marker->getPointFromID(ii, &oldPoint))
// 		{
// 			float tmpPoint[4], tmpPoint2[4];
// 			tmpPoint[0] = oldPoint[0], tmpPoint[1] = oldPoint[1];
// 			tmpPoint[2] = 0, tmpPoint[3] = 1;	//homogeneous coordinates
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				tmpPoint2[jj] = 0;
// 				for(int kk = 0; kk < 4; ++kk)
// 				{
// // 					tmpPoint2[jj] += sectionTransform[jj][kk]*tmpPoint[kk];	//inverse rotation
// 					tmpPoint2[jj] += totalTransform[jj][kk]*tmpPoint[kk];	//inverse rotation
// 				}
// // 				tmpPoint2[3] += sectionTransform[jj][3]*tmpPoint[3];	//inverse translation
// // 				tmpPoint2[jj] /= (DOWNSAMPLING*XYSAMPLING);	//transform is in world coordinates, but markers need to be in pixel coordinates
// 			}
// 			newPoint[0] = tmpPoint2[0], newPoint[1] = tmpPoint2[1];
// 			transformedMarker->SetPoint(ii, newPoint);
// 		}
// 	}
// 	marker->setBarrelMarkers(transformedMarker);
// 	
// 	std::string transFilename(amiraFilename);
// 	transFilename += "_S%02d_to_S%02d_total_transform.landmarkAscii";
// 	char transFilename2[256];
// 	sprintf(transFilename2, transFilename.c_str(), thisSectionID+7, previousSectionID+7);
// 	std::string transFilename3(transFilename2);
// 	std::ofstream LandmarkFile(transFilename3.c_str());
// 	
// 	LandmarkFile << "# AmiraMesh 3D ASCII 2.0" << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << "define Markers " << marker->getNumberOfBarrelMarkers() << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << "Parameters {" << std::endl;
// 	LandmarkFile << "\tNumSets 1," << std::endl;
// 	LandmarkFile << "\tContentType \"LandmarkSet\"" << std::endl;
// 	LandmarkFile << "}" << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << "Markers { float[3] Coordinates } @1" << std::endl;
// 	LandmarkFile << std::endl;
// 	LandmarkFile << "# Data section follows" << std::endl;
// 	LandmarkFile << "@1" << std::endl;
// 	
// 	PointType tmpPoint;
// 	for(int ii = 0; ii < marker->getNumberOfBarrelMarkers(); ++ii)
// 	{
// 		marker->getPointFromID(ii, &tmpPoint);
// 		LandmarkFile << tmpPoint[0] << " " << tmpPoint[1] << " " << 0 << std::endl;
// 	}
// 	LandmarkFile.close();
};

//use with extreme caution!!! assumes memory for sectionTranslation and sectionRotation
//is allocated! calculates T*R
float ** Section::calculateTransAfterRot()
{
	float ** mProduct = new float *[4];
	for(int ii = 0; ii < 4; ++ii)
	{
		mProduct[ii] = new float[4];
		for(int jj = 0; jj < 4; ++jj)
			mProduct[ii][jj] = 0;
	}
	
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			for(int kk = 0; kk < 4; ++kk)
				mProduct[ii][jj] += sectionTranslation[ii][kk]*sectionRotation[kk][jj];
	
	return mProduct;
};

//use with extreme caution!!! assumes memory for sectionTranslation and sectionRotation
//is allocated! calculates (T*R)^-1 = R^-1 * T^-1 using known formula for product Rotation*Translation
float ** Section::calculateTransAfterRotInverse()
{
	float ** mInverse = new float *[4];
	for(int ii = 0; ii < 4; ++ii)
	{
		mInverse[ii] = new float[4];
		for(int jj = 0; jj < 4; ++jj)
			mInverse[ii][jj] = 0;
	}
	for(int ii = 0; ii < 2; ++ii)
		for(int jj = 0; jj < 2; ++jj)
			mInverse[ii][jj] = sectionRotation[jj][ii];
	mInverse[0][3] = -1*(sectionRotation[0][0]*sectionTranslation[0][3] + sectionRotation[1][0]*sectionTranslation[1][3]);
	mInverse[1][3] = -1*(sectionRotation[0][1]*sectionTranslation[0][3] + sectionRotation[1][1]*sectionTranslation[1][3]);
	mInverse[2][3] = -1*sectionTranslation[2][3];
	mInverse[3][3] = 1;
	
	return mInverse;
};

//calculates transformation for barrel marker from previous sectionID
//use only while Section::imageProcessor is still valid!!!
//marker don't need to be changed for optimized region growing
void Section::calculatePropagatedBarrelMarker(Section * previousSection)
{
	std::vector< std::vector< Contour * > > previousContours = previousSection->getBarrelContourVector();
	if(!previousContours.size())
	{
		std::cout << "Error! previousContour with ID = " << previousSectionID << " has no barrel z-contours!" << std::endl;
		std::cout << "Cannot calculate barrel marker; aborting." << std::endl;
		return;
	}
	int maxZ = previousContours.size() - 1;
	int nrOfBarrels = previousContours[0].size();
	if(!nrOfBarrels)
	{
		std::cout << "Error! previousContour with ID = " << previousSectionID << " has no barrel contours!" << std::endl;
		std::cout << "Cannot calculate barrel marker; aborting." << std::endl;
		return;
	}
// 	std::flush(std::cout << "Allocating new contours!" << std::endl);
	for(int z = 0; z < nrOfZPlanes; ++z)
	{
		std::vector< Contour * > zPlaneContours;
		for(int ii = 0; ii < nrOfBarrels; ++ii)
		{
			Contour * newContour = new Contour;
			newContour->setBarrelID(previousContours[0][ii]->getBarrelID());	//ID should be ii+1 b/c of voronoiMap
// 			newContour->setValid(previousContours[0][ii]->getValid());
			zPlaneContours.push_back(newContour);
		}
		zBarrelContoursVector.push_back(zPlaneContours);
	}
	float ** avgPoints = new float *[nrOfBarrels];
	PointSetType::Pointer avgMarker = PointSetType::New();
	//case: this on top of previous
	if(previousSectionID > thisSectionID)
	{
// 		std::flush(std::cout << "avg from previous top planes" << std::endl);
		float * totalEdgePtNr = new float[nrOfBarrels];
		for(int ii = 0; ii < nrOfBarrels; ++ii)
		{
			totalEdgePtNr[ii] = 0;
			avgPoints[ii] = new float[2];
			avgPoints[ii][0] = 0, avgPoints[ii][1] = 0;
		}
		for(int ii = 0; ii <= maxZ && ii < 5; ++ii)
		{
			for(int jj = 0; jj < nrOfBarrels; ++jj)
			{
				std::list< std::vector< float > > * prevEdgeList = previousContours[ii][jj]->edgeListPointer();
				std::list< std::vector< float > >::iterator contourIt;
				for(contourIt = prevEdgeList->begin(); contourIt != prevEdgeList->end(); ++contourIt)
				{
					avgPoints[jj][0] += (*contourIt)[0];
					avgPoints[jj][1] += (*contourIt)[1];
				}
				totalEdgePtNr[jj] += (float)prevEdgeList->size();
			}
		}
		for(int ii = 0; ii < nrOfBarrels; ++ii)
		{
			if(totalEdgePtNr[ii])
			{
				avgPoints[ii][0] /= totalEdgePtNr[ii];
				avgPoints[ii][1] /= totalEdgePtNr[ii];
			}
			PointType newAvgPoint;
			newAvgPoint[0] = avgPoints[ii][0], newAvgPoint[1] = avgPoints[ii][1];
			avgMarker->SetPoint(ii, newAvgPoint);
		}
	}
	//case: previous on top of this
	else
	{
// 		std::flush(std::cout << "avg from previous bottom planes" << std::endl);
		float * totalEdgePtNr = new float[nrOfBarrels];
		for(int ii = 0; ii < nrOfBarrels; ++ii)
		{
			totalEdgePtNr[ii] = 0;
			avgPoints[ii] = new float[2];
			avgPoints[ii][0] = 0, avgPoints[ii][1] = 0;
		}
		for(int ii = maxZ; ii >= 0 && ii > maxZ - 5; --ii)
		{
			for(int jj = 0; jj < nrOfBarrels; ++jj)
			{
				std::list< std::vector< float > > * prevEdgeList = previousContours[ii][jj]->edgeListPointer();
				std::list< std::vector< float > >::iterator contourIt;
				for(contourIt = prevEdgeList->begin(); contourIt != prevEdgeList->end(); ++contourIt)
				{
					avgPoints[jj][0] += (*contourIt)[0];
					avgPoints[jj][1] += (*contourIt)[1];
				}
				totalEdgePtNr[jj] += (float)prevEdgeList->size();
			}
		}
		for(int ii = 0; ii < nrOfBarrels; ++ii)
		{
			if(totalEdgePtNr[ii])
			{
				avgPoints[ii][0] /= totalEdgePtNr[ii];
				avgPoints[ii][1] /= totalEdgePtNr[ii];
			}
			PointType newAvgPoint;
			newAvgPoint[0] = avgPoints[ii][0], newAvgPoint[1] = avgPoints[ii][1];
			avgMarker->SetPoint(ii, newAvgPoint);
		}
	}
// 	writeBarrelMarkerLandmarks(avgMarker);
	this->marker = new BarrelMarker();
// 	std::flush(std::cout << "calculating this section translation after rotation inverse matrix..." << std::endl);
	float ** thisSectionTotalTrans = calculateTransAfterRotInverse();
// 	std::flush(std::cout << "calculating previous section translation after rotation matrix..." << std::endl);
	float ** previousSectionTotalTrans = previousSection->calculateTransAfterRot();
	float ** totalTransform = new float *[4];
	for(int ii = 0; ii < 4; ++ii)
	{
		totalTransform[ii] = new float[4];
		for(int jj = 0; jj < 4; ++jj)
			totalTransform[ii][jj] = 0;
	}
// 	std::flush(std::cout << "calculating total transformation matrix..." << std::endl);
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			for(int kk = 0; kk < 4; ++kk)
				totalTransform[ii][jj] += thisSectionTotalTrans[ii][kk]*previousSectionTotalTrans[kk][jj];
	
// 	std::flush(std::cout << "total transformation matrix:" << std::endl);
// 	for(int ii = 0; ii < 4; ++ii)
// 	{
// 		std::cout << "[";
// 		for(int jj = 0; jj < 4; ++jj)
// 		{
// 			if(jj < 3)
// 				std::cout << totalTransform[ii][jj] << ",\t";
// 			else
// 				std::cout << totalTransform[ii][jj];
// 		}
// 		std::cout << "]" << std::endl;
// 	}
	
	PointSetType::Pointer transformedMarker = PointSetType::New();
	//later, calculate the markers from the top/bottom contours; watch sampling below!!!
// 	for(int ii = 0; ii < previousSection->getBarrelMarker()->getNumberOfBarrelMarkers(); ++ii)
	for(int ii = 0, pointID = 0; ii < nrOfBarrels; ++ii)	//use pointID so that marker labels points internally w/o gaps (due to non-valid contours)
	{
		PointType oldPoint, newPoint;
// 		if(previousSection->getBarrelMarker()->getPointFromID(ii, &oldPoint))
		if(avgMarker->GetPoint(ii, &oldPoint))
		{
			float tmpPoint[4], tmpPoint2[4];
			tmpPoint[0] = oldPoint[0], tmpPoint[1] = oldPoint[1];
			tmpPoint[2] = 0, tmpPoint[3] = 1;	//homogeneous coordinates
			for(int jj = 0; jj < 4; ++jj)
			{
				tmpPoint2[jj] = 0;
				for(int kk = 0; kk < 3; ++kk)
					tmpPoint2[jj] += totalTransform[jj][kk]*tmpPoint[kk];
				tmpPoint2[jj] += totalTransform[jj][3]*tmpPoint[3]/(DOWNSAMPLING*XYSAMPLING);	//correction for real world coordinates used in writing barrel contours
			}
			newPoint[0] = tmpPoint2[0], newPoint[1] = tmpPoint2[1];
			//if transformed marker is out of bounds, don't use and set barrel flag to not valid
			if(isInBounds(newPoint))
			{
				transformedMarker->SetPoint(pointID, newPoint);
				zBarrelContoursVector[0][ii]->setValid(1);
				++pointID;
			}
			else
			{
				zBarrelContoursVector[0][ii]->setValid(0);
				std::cout << "Point " << newPoint << " is out of bounds!" << std::endl;
			}
		}
	}
	
	this->marker->setBarrelMarkers(transformedMarker);
	
// 	std::flush(std::cout << "set valid flag in all " << nrOfZPlanes << " z planes" << std::endl);
	//set valid flag in all z planes
	for(int ii = 1; ii < nrOfZPlanes; ++ii)
		for(int jj = 0 ; jj < nrOfBarrels; ++jj)
			zBarrelContoursVector[ii][jj]->setValid(zBarrelContoursVector[0][jj]->getValid());
	
// 	std::flush(std::cout << "remapping voronoimap!" << std::endl);
	adjustVoronoiMapToBarrelIDs();
	
	for(int ii = 0; ii < 4; ++ii)
		delete [] thisSectionTotalTrans[ii], delete [] previousSectionTotalTrans[ii], delete [] totalTransform[ii];
	delete [] thisSectionTotalTrans, delete [] previousSectionTotalTrans, delete [] totalTransform;
};

//if this section is not a middle section, then this method is used
//to set the voronoi IDs equal to the corresponding barrel IDs
//only to be used together with calculatePropagatedBarrelMarker(),
//b/c it assumes the barrel IDs to be in the same order as the
//marker IDs!!!
void Section::adjustVoronoiMapToBarrelIDs()
{
	ImageType::RegionType planeRegion = this->imageProcessor->getInputRegion();
	ImageType::SizeType planeSize = planeRegion.GetSize();
	planeSize[2] = 1;
	planeRegion.SetSize(planeSize);
	this->marker->computeVoronoiDiagram(planeRegion);
	
	int * voronoiRegionMapping = new int[marker->getNumberOfBarrelMarkers()];
	CalcImageType::Pointer voronoiMap = marker->getVoronoiMap();
	CalcImageType::Pointer distanceMap = marker->getDistanceMap();
	CalcIndexIteratorType voronoiIter(voronoiMap, voronoiMap->GetLargestPossibleRegion());
	ConstCalcIteratorType distanceIter(distanceMap, distanceMap->GetLargestPossibleRegion());
	for(int ii = 0, pointID = 0; ii < zBarrelContoursVector[0].size(); ++ii)	//use pointID so that marker labels points internally w/o gaps (due to non-valid contours)
	{
		if(zBarrelContoursVector[0][ii]->getValid())
		{
			PointType tmpPoint;
			if(marker->getPointFromID(pointID, &tmpPoint))
			{
				CalcImageType::IndexType tmpIndex;
				tmpIndex[0] = tmpPoint[0], tmpIndex[1] = tmpPoint[1], tmpIndex[2] = 0;
				voronoiIter.SetIndex(tmpIndex);
				voronoiRegionMapping[int(voronoiIter.Get()+0.5)-1] = ii;
			}
			++pointID;
		}
	}
	
	CalcImageType::Pointer newVoronoiMap = CalcImageType::New();
	newVoronoiMap->SetRegions(voronoiMap->GetLargestPossibleRegion());
	newVoronoiMap->Allocate();
	CalcIteratorType newVoronoiIter(newVoronoiMap, newVoronoiMap->GetLargestPossibleRegion());
	for(voronoiIter.GoToBegin(), newVoronoiIter.GoToBegin(), distanceIter.GoToBegin();
	!voronoiIter.IsAtEnd() && !newVoronoiIter.IsAtEnd() && !distanceIter.IsAtEnd();
	++voronoiIter, ++newVoronoiIter, ++distanceIter)
	{
		if(distanceIter.Get() > 175)
			newVoronoiIter.Set(0);
		else
			newVoronoiIter.Set(voronoiRegionMapping[int(voronoiIter.Get()+0.5)-1] + 1);
	}
	newVoronoiMap->Update();
	this->marker->setVoronoiMap(newVoronoiMap);
	this->imageProcessor->setVoronoiMap(newVoronoiMap);
	this->imageProcessor->setDistanceMap(distanceMap);
	this->imageProcessor->setBarrelPoints(this->marker->getNumberOfBarrelMarkers());
};

//transform all contours to world coordinates
void Section::transformToWorldCoordinates()
{
	float ** worldTransform = calculateTransAfterRot();
	for(int z = 0; z < zBarrelContoursVector.size(); ++z)
		for(int ii = 0; ii < zBarrelContoursVector[z].size(); ++ii)
			if(zBarrelContoursVector[z][ii]->getValid())
			{
				float tmpPoint[4];
				float transPoint[4];
				std::list< std::vector< float > >::iterator contourIt;
				for(contourIt = zBarrelContoursVector[z][ii]->edgeListPointer()->begin();
				contourIt != zBarrelContoursVector[z][ii]->edgeListPointer()->end(); ++contourIt)
				{
					tmpPoint[0] = (*contourIt)[0]*XYSAMPLING*DOWNSAMPLING;
					tmpPoint[1] = (*contourIt)[1]*XYSAMPLING*DOWNSAMPLING;
					tmpPoint[2] = (*contourIt)[2];
					tmpPoint[3] = 1;
					for(int jj = 0; jj < 4; ++jj)
					{
						transPoint[jj] = 0;
						for(int kk = 0; kk < 4; ++kk)
							transPoint[jj] += worldTransform[jj][kk]*tmpPoint[kk];
					}
					(*contourIt)[0] = transPoint[0];
					(*contourIt)[1] = transPoint[1];
					(*contourIt)[2] = transPoint[2];
				}
			}
};

void Section::writeBarrelMarkerLandmarks()
{
	std::string transFilename(outputFilename);
	transFilename += "_barrel_marker.landmarkAscii";
// 	char transFilename2[256];
// 	sprintf(transFilename2, transFilename.c_str(), thisSectionID+7, previousSectionID+7);
// 	std::string transFilename3(transFilename2);
	std::ofstream LandmarkFile(transFilename.c_str());
	
	LandmarkFile << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "define Markers " << marker->getNumberOfBarrelMarkers() << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "Parameters {" << std::endl;
	LandmarkFile << "\tNumSets 1," << std::endl;
	LandmarkFile << "\tContentType \"LandmarkSet\"" << std::endl;
	LandmarkFile << "}" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "Markers { float[3] Coordinates } @1" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "# Data section follows" << std::endl;
	LandmarkFile << "@1" << std::endl;
	
	PointType tmpPoint;
	for(int ii = 0; ii < marker->getNumberOfBarrelMarkers(); ++ii)
	{
		marker->getPointFromID(ii, &tmpPoint);
		LandmarkFile << tmpPoint[0] << " " << tmpPoint[1] << " " << 0 << std::endl;
	}
	LandmarkFile.close();
};

void Section::writeBarrelMarkerLandmarks(PointSetType::Pointer outMarker)
{
	std::string transFilename(outputFilename);
	transFilename += "_avg_marker_no_trans.landmarkAscii";
	// 	char transFilename2[256];
	// 	sprintf(transFilename2, transFilename.c_str(), thisSectionID+7, previousSectionID+7);
	// 	std::string transFilename3(transFilename2);
	std::ofstream LandmarkFile(transFilename.c_str());
	
	LandmarkFile << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "define Markers " << outMarker->GetNumberOfPoints() << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "Parameters {" << std::endl;
	LandmarkFile << "\tNumSets 1," << std::endl;
	LandmarkFile << "\tContentType \"LandmarkSet\"" << std::endl;
	LandmarkFile << "}" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "Markers { float[3] Coordinates } @1" << std::endl;
	LandmarkFile << std::endl;
	LandmarkFile << "# Data section follows" << std::endl;
	LandmarkFile << "@1" << std::endl;
	
	PointType tmpPoint;
	for(int ii = 0; ii < outMarker->GetNumberOfPoints(); ++ii)
	{
		outMarker->GetPoint(ii, &tmpPoint);
		LandmarkFile << tmpPoint[0] << " " << tmpPoint[1] << " " << 0 << std::endl;
	}
	LandmarkFile.close();
};


bool Section::isInBounds(PointType point)
{
	if(point[0] >= boundingBox[0] && point[0] <= boundingBox[1] && point[1] >= boundingBox[2] && point[1] <= boundingBox[3])
		return 1;
	return 0;
};














