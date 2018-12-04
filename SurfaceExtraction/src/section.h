/****************************************************************************/
/*                                                                          */
/* File:      section.h 						    */
/*                                                                          */
/* Purpose:   wrapper class for section-wise barrel segmentation	    */
/*	      of barrel contours					    */
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
#pragma once
#include "typedefs.h"
#include "bloodvesselpattern.h"
#include "segmentation.h"

#ifndef SECTION
#define SECTION

class Section
{
	public:
		Section(const char * markerFilename, const char * inputFilename, const char * outputFilename, const char * transformFilename, const char * medianProjFilename, int start, int stop, int sectionID, int amiraSectionID, float arg1, float arg2, float arg3);
		~Section();
		
		std::vector< std::vector< Contour * > > getBarrelContourVector(){return zBarrelContoursVector;}
		void startFirstSegmentation(float arg);
		void setBarrelContours(bool write);
		void setBarrelContours(std::vector< std::vector< Contour * > > contourVec);
		void writeRegularBarrels();
		
		void startOptimizedSegmentation();
		void replaceNotCorrectlyOptimizedContours(int direction);
		void prepareForOutput();
		void writeOptimizedBarrels();
		
		//for testing of marker transformation
		BarrelMarker * getBarrelMarker(){return this->marker;}
		
		void readAmiraBarrelMarker();
		//calculates transformation for barrel marker from previous sectionID
		//use only while Section::imageProcessor is still valid!!!
		//marker don't need to be changed for optimized region growing
		void calculatePropagatedBarrelMarker(Section * previousSection);
		bool readAmiraSectionTransformations();
		//use with extreme caution!!! assumes memory for sectionTranslation and sectionRotation
		//is allocated! calculates T*R
		float ** calculateTransAfterRot();
		//use with extreme caution!!! assumes memory for sectionTranslation and sectionRotation
		//is allocated! calculates (T*R)^-1 = R^-1 * T^-1 using known formula for product Rotation*Translation
		float ** calculateTransAfterRotInverse();
		
		void writeBarrelMarkerLandmarks();
		
		int previousSectionID, nextSectionID;	//reflects z-orientation of the slices
		
	private:
		BarrelMarker * marker;
		BloodVesselPattern * vessels;
		//contains barrel contours planewise: zBarrelContoursVector[z][barrelID]
		std::vector< std::vector< Contour * > > zBarrelContoursVector;
		Segmentation * imageProcessor;
		Segmentation * barrelOptimizer;
		
		long * boundingBox;
		
		const char * inputFilename;
		const char * outputFilename;
		const char * transformFilename;
		const char * markerFilename;
		std::ofstream DebugLog;
		
		int start, stop;
		int downSampleRate;
		
		ImageType::Pointer segImageStack;
		ImageType::Pointer preprocImageStack;
		
		//transforminfo: simply 4x4 float matrix
		//in amira style: total transformation w.r.t. first slice
		//of the form T*R
		float ** sectionTranslation;
		float ** sectionRotation;
		int nrOfZPlanes;
		int thisSectionID;	//section Sxx
		int amiraSectionID;	//for compatibility with Amira transformation: Slice0001 is always first section w/ transform info
		
		bool isInBounds(PointType point);
		void adjustVoronoiMapToBarrelIDs();
		void transformToWorldCoordinates();
		
		void fillEmptyBarrelContours();
		
		void writeBarrelMarkerLandmarks(PointSetType::Pointer outMarker);
};


#endif