/****************************************************************************/
/*                                                                          */
/* File:      basics.h	                                                    */
/*                                                                          */
/* Purpose:   class for basic functions and common variables                */
/*            especially writing and reading functions                      */
/*            extended functionality (classes Contour & BarrelMarker)       */
/*            added for barrel segmentation pipeline                        */
/*                                                                          */
/* Authors:   Marcel Oberlaender                                            */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/*            Stefan Reissl                                                 */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/*            Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"

#ifndef BASICS
#define BASICS


#define WATERSHED_COLOR_STEP	10
#define KMEANS_COLOR_STEP		13
#define CRITICAL_DENSITY		0.00008	// [N/�m^3]

class Contour
{
	public:
		Contour();
		Contour(std::list<std::vector<float> > _edge_list);
		Contour(Contour * otherContour);
		~Contour();
		void setEdgeList(std::list<std::vector<float> > _edge_list){this->edgeList = _edge_list;}
		void setEdgeList(std::list<std::vector<float> > * _edge_list);
		void replaceEdgeList(std::list<std::vector<float> > * _edge_list, float zOffset);
		std::list<std::vector<float> > * edgeListPointer(){return &edgeList;}
		
		void setOptEdgeList(std::list<std::vector<float> > _edge_list){this->optEdgeList = _edge_list;}
		void setOptEdgeList(std::list<std::vector<float> > * _edge_list);
		std::list<std::vector<float> > * optEdgeListPointer(){return &optEdgeList;}
		
		void setAttributeList(std::vector<float> _attributes){this->attributes = _attributes;}
		void addAttribute(float attr){attributes.push_back(attr);}
		std::vector< float >::const_iterator readAttributes(){return attributes.begin();}
		std::vector< float >::const_iterator attributesEnd(){return attributes.end();}
		std::vector< float > * attributesPointer(){return &attributes;}
		
		void setOptimizeFlag(bool flag){optimize = flag;}
		bool getOptimizeFlag(){return optimize;}
		void setValid(bool flag){validContour = flag;}
		bool getValid(){return validContour;}
		void setBarrelID(int ID){barrelID = ID;}
		int getBarrelID(){return barrelID;}
		void setInsideBarrel(int flag){insideBarrel = flag;}
		int getInsideBarrel(){return insideBarrel;}
		
		void prepareForOutput();
		
	private:
		std::list< std::vector<float> > edgeList;
		std::list< std::vector< float > > optEdgeList;
		std::vector< float > attributes;
		
		bool optimize;
		bool validContour;
		int barrelID;
		bool insideBarrel;	// 1 = keep; 0 = cut away
		
		void contourSmoothing(int radius);
		void contourSampling(int samplingRate);
};

class Basics
{
	public:
		Basics();
		~Basics();
		
	protected:
		const char * inputFilename;						// filename of the input images
		const char * outputFilename;					// filename of CellClusterList.csv, image files (*.sr), landmark file and histo stuff 
		ImageType::RegionType maximum_region;			// contains the maximum region of the input image 
		ImageType::Pointer inputImage;					// main image for all pre processing and cell count processes
		Image2DType::Pointer input2DImage;					// main image for all pre processing and cell count processes
		NeighborhoodOffsetVectorType look_up_table;		// array of for the definition of neighborhood relation for the region growing
		ShapedNeighborhoodOffsetVectorType neighborhood4;	//array with offsets defining a 4-neighborhood in a 2D-plane of a 3D-image
		ShapedNeighborhoodOffsetVectorType neighborhood8;	//array with offsets defining a 8-neighborhood in a 2D-plane of a 3D-image
		ImageType::SizeType radius1;
		long * boundingBox;								// min/max indices of the input image (2D/3D dep. on input)
		unsigned int start;								// start position inside the stack
		unsigned int stop;								// stop position inside the stack
		unsigned int scanType;							// scanType: 1 means confocal, 2 means 2 photon
		bool channel_artifacts;							// if true, artifacts are removed prior to intensity correction; especially important for 2-photon images	
		bool borderStack;								// if false, no splitting of large clusters will be executed
		bool binaryBall;								// array of for the definition of neighborhood relation for the binary erosion
		
		//******** basic input and output functions ********
		
// 		ImageType::Pointer LoadCluster(CellCluster * cluster, const char * path);	// creates an image with color values between 0 and 255 according to the cell cluster list item and the corresponding sr file	
// 		ImageType::Pointer LoadClusterWithBox(CellCluster * cluster);				// creates an black and white image with a surrounding black box according to the cell cluster list item and the corresponding sr file																			// 
// 		void SaveCluster(CellCluster * cluster, ImageType::Pointer image);			// saves the pixel of the image in a sr file with the id number of the cluster		
// 		void DeleteCluster(CellCluster * cluster);									// deletes the sr file with the id number of the cluster
		
		void writeBinaryImage(ImageType::Pointer binaryImage, const char * label, unsigned int objectNumber);	// writes the binary image with a label and a object number in the filename
		void writeBinaryImage(ImageType::Pointer binaryImage, unsigned int objectNumber);						// writes the binary image with a object number in the filename
		void writeImagePlanes();																				// writes the current input image to output filename
		
// 		CellClusterListType readCellClusterList( const char * _filename );			// reads a cell cluster list from the filename and stores it in a STL list
// 		void writeCellClusterList(CellClusterListType clusters, const char * file);	// writes a given sell cluster list to output filename 
		
		void readImage( int start, int stop, const char * inputFilename );	// reads a image stack from start to stop position and stores it in the input image
		void writeInputImage(const char * name);							// saves the current input image to output filename
		
		void read2DImage( const char * inputFilename );
		void writeInput2DImage(const char * name);
		
		std::list<std::list<std::vector<float> > > GetAmiraContours(std::list<std::list< std::vector< float > > > contours);
		std::list<std::list<std::vector<float> > > GetAmiraContours(Contour * contours, int nrOfContours);
		std::list<std::list<std::vector<float> > > GetAmiraContours(std::vector< Contour * > contours, int nrOfContours);
		std::list<std::list<std::vector<float> > > GetAmiraContours(std::vector< std::vector< Contour * > > contours, int nrOfContours);
		std::list<std::vector<float> >  GetAmiraContour(std::list<std::vector<float> >& structure, bool worldCoordinates);
		
		//******** functions for projection images ********
		
		void writeMeanProjectionImage(const char * name);	// writes the mean projection image of the input image
		void writeMaxProjectionImage(const char * name);	// Writes the maximum projection image of the input image
		
		//******** test functions ********
// 		bool isClusterValid(CellCluster *cluster)									// tests the validity of a given cluster (bounds, propability)
// 		{return	isClusterValid(cluster, maximum_region);};
// 		bool isClusterValid(CellCluster *cluster, ImageType::RegionType region);	// tests the validity of a given cluster (bounds, propability)
		
		bool IsInBounds(ImageType::Pointer image, ImageType::IndexType index);		// test whether an index is inside an image region or not
		bool IsNotAtImageBorder(float x, float y, float z);							// tests whether the coordinates are at the bound of the input image
		
		// tests whether an index is inside a box
		bool isInBox(ImageType::RegionType::IndexType index, std::vector<LongVectorType> box) 
		{
			if(index[0]>=(box[0])[0] && index[0]<=(box[0])[1] && index[1]>=(box[1])[0] && index[1]<=(box[1])[1] && index[2]>=(box[2])[0] && index[2]<=(box[2])[1])
				return true;
			
			return false;
		};
		
		// tests whether two indices have a minimum distance from each other
		bool isInBox(ImageType::IndexType pos1, ImageType::IndexType pos2, unsigned int diff)
		{
			
			if( std::abs(pos1[0]-pos2[0]) <= diff && 
				std::abs(pos1[1]-pos2[1]) <= diff &&
				std::abs(pos1[2]-pos2[2]) <= diff)
				return true;
			
			return false;
		};
		
		// tests whether an index is inside a region
		bool isInBox(ImageType::RegionType region, ImageType::IndexType index)
		{
			ImageType::IndexType rIndex = region.GetIndex();
			ImageType::SizeType  rSize  = region.GetSize();
			
			if(rIndex[0]<=index[0] && index[0]< rSize[0] && rIndex[1]<=index[1] && index[1]< rSize[1] && rIndex[2]<=index[2] && index[2]< rSize[2])
				return true;
			
			return false;
		};
		
		
		////******** other functions ********
		
		NeighborhoodOffsetVectorType CreateLookUpTable(); // creates a look up table for the region growing
		void createNeighborhood(); // creates a look up table for 4- and 8-neighborhood
		void computeBounds(ImageType::Pointer image);	//explicit creation of boundingBox in case no image is read in the beginning
		
		void normalizeHistogram(float * histogram, unsigned int noOfBins)	//Normalizes histograms with a certain size
		{
			float maxValue = 0.0;
			
			for(int i = 0; i < noOfBins; ++i)
			{
				if(histogram[i] > maxValue)
					maxValue = histogram[i];
			}
			
			if(maxValue==0) 
				return;
			
			for(int i = 0; i < noOfBins; ++i)
			{
				histogram[i] /= maxValue;
			}
		};
};

class BarrelMarker
{
	public:
		BarrelMarker();
		~BarrelMarker();
		
		void readBarrelMarkerFile(const char * markerFilename, ImageType::RegionType region);
		/*takes marker calculated from previous stack as input*/
		void setBarrelMarkers(PointSetType::Pointer newMarker){this->marker = newMarker;}
		void computeVoronoiDiagram(ImageType::RegionType imageRegion);
		
		CalcImageType::Pointer getVoronoiMap(){return voronoiMap;}
		void setVoronoiMap(CalcImageType::Pointer newVoronoiMap){voronoiMap = newVoronoiMap;}
		CalcImageType::Pointer getDistanceMap(){return distanceMap;}
		unsigned int getNumberOfBarrelMarkers(){return marker->GetNumberOfPoints();}
		bool getPointFromID(int ID, PointType * point){return marker->GetPoint(ID, point);}
		
	private:
		PointSetType::Pointer marker;
		CalcImageType::Pointer voronoiMap;
		CalcImageType::Pointer distanceMap;
};

#endif