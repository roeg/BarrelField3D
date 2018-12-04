/****************************************************************************/
/*                                                                          */
/* File:      segmentation.h 						    */
/*                                                                          */
/* Purpose:   class for segmentation of either Pia/WM images (4x air)	    */
/*            or barrel field images(40x oil)                               */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*                                                                          */
/*                                                                          */
/* EMail: Robert.Egger@maxplanckflorida.org                                 */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/
#include "typedefs.h"
#include "basics.h"
#include "bloodvesselpattern.h"

#ifndef SEGMENTATION
#define SEGMENTATION

struct AmiraContourGraphStruct
{
	AmiraContourGraphStruct(std::list<std::vector<float> > _vertice_list, std::list<std::list<std::vector<float> > > _edge_list);
	AmiraContourGraphStruct();
	std::list<std::vector<float> > vertice_list;
	std::list<std::list<std::vector<float> > > edge_list;
};
typedef struct AmiraContourGraphStruct	AmiraContourGraphType;

class Segmentation : public Basics
{
	public:
		Segmentation( int start, int stop, const char * inputFilename, const char * outputFilename, const char * medProjFilename, int samplingRate, int slice, float arg1, float arg2, float arg3 );
		Segmentation( int start, int stop, const char * inputFilename, const char * outputFilename, float factor, float factor2, float factor3, int slice );
		Segmentation( int start, int stop, int samplingRate, ImageType::Pointer segmentedStack, ImageType::Pointer preprocStack, BloodVesselPattern * vesselPattern, std::vector< std::vector< Contour * > > barrelContours, const char * outputFilename);
		~Segmentation();
		
		void piaWMExtraction(int WMbegin, bool SegmentPia);
		void barrelFieldBvp();
		
		//interface to section class
		void barrelExtraction(bool preProcessed);
		void setVoronoiMap(CalcImageType::Pointer voronoiMap){this->voronoiMap = voronoiMap;}
		void setDistanceMap(CalcImageType::Pointer distanceMap){this->distanceMap = distanceMap;}
		void setBarrelPoints(unsigned int nrOfPoints){this->nrOfBarrelMarkerPoints = nrOfPoints;}
		ImageType::RegionType getInputRegion(){return this->inputImage->GetLargestPossibleRegion();}
		std::vector< std::vector< Contour * > > getRegularBarrelContours();
		std::vector< std::vector< Contour * > > getOptimizedBarrelContours();
		void setBarrelContours(std::vector< std::vector< Contour * > > contours){this->barrelZContours = contours;}
		ImageType::Pointer getSegmentedBarrelImages(){return segmentedImageStack;}
		ImageType::Pointer getPreprocessedBarrelImages(){return preprocImageStack;}
		BloodVesselPattern * getVesselPattern(){return bvPattern;}
		long * getBounds(){return boundingBox;}
		void writeVoronoiDiagram();
		//not used anymore
		void writeBarrelContours(/*std::vector< Contour * > barrelContours*/);
		//not used anymore
		void writeBarrelContours(std::vector< Contour * > barrelContours);
		
	private:
// 		ImageType::Pointer inputImage;
		float factor, factor2, factor3;
		int int1, int2, sliceNumber;
		int XYDownSamplingRate;
		BloodVesselPattern * bvPattern;
		AmiraContourGraphType* amira_contour_graph;
		AmiraContourGraphType* amira_bvp_graph;
		unsigned long nrOfObjects;
		std::vector<unsigned long> sizeOfObjects;
		ImageType::Pointer originalImage;
		ImageType::Pointer medProjImage;
		ImageType::Pointer wmBvpImage;
		ImageType::Pointer sigmoidImage;
		ImageType::Pointer sigmoidImageStack;
		ImageType::Pointer preprocImage;
		ImageType::Pointer preprocImageStack;
		ImageType::Pointer segmentedImage;
		ImageType::Pointer segmentedImageStack;
		ImageType::Pointer fgRemovedImage;
		ObjectImageType::Pointer labeledImage;
		CalcImageType::Pointer calcImage;
		CalcImageType::Pointer voronoiMap;
		CalcImageType::Pointer distanceMap;
		unsigned int nrOfBarrelMarkerPoints;
		std::vector< std::vector< Contour * > > barrelZContours;
		std::ofstream DebugLog;

		/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
		void copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2 );
		/*copy( image1, image2, planeRegion ) performs an iterator copy of planeRegion of image2 (3D) to image1 (2D)     */
		void copy( Image2DType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion );
		/*copy( image1, image2, planeRegion ) performs an iterator copy of planeRegion of image2 (3D) to image1 (3D plane)     */
		void copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion );
		/*copy( image1, image2, planeRegion ) performs an iterator copy of planeRegion of image2 (2D) to image1 (3D)     */
		void copy( ImageType::Pointer workingImage1, Image2DType::Pointer workingImage2, ImageType::RegionType planeRegion );
		/*copy( image1, image2, planeRegion ) performs an iterator copy of image2 (3D plane) to planeRegion of image1 (3D)     */
		void copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion, int second );
		ImageType::Pointer invertImage( ImageType::Pointer workingImage );
		CalcImageType::Pointer invertImage( CalcImageType::Pointer workingImage );
		void medianFilter( unsigned int radius );
		void sigmoidIntensityMapping(float a, float b, float p);
		
		void regionGrowingObjectLabeling( ImageType::Pointer workingImage, int minimum_size, bool neighborhood26 );
		void relabelImage( ObjectImageType::Pointer workingImage, int minimum_size );
		void labeledToGreyscaleImage( ImageType::Pointer& workingImage, unsigned long upperThreshold );
		
		//	calculates the mean and standard deviation of all regions passed in the values vector    
		void calculateRegionalHistograms( float * means, float * sigmas, std::vector< unsigned int > * values, unsigned int divisions );
		void calculateRegionalHistograms( float * means, float * sigmas, std::list< unsigned char > * values, unsigned int divisions );
		//	calculates the histograms of all regions passed in the values vector    
		void calculateRegionalHistograms( int ** histogram, std::list< std::list< unsigned char > > values, unsigned int divisions );
		void calculateRegionalHistograms( int** histogram, std::list< unsigned char > * values, unsigned int divisions );
		void Histogram(float * moments);
		void Histogram(float * moments, ImageType::Pointer workingImage);
		/*completeHistogram() calculates the histogram of the maximum region      */
		int * completeHistogram( ImageType::Pointer tmpImage );
		unsigned int * foregroundHistogram( ImageType::Pointer bgImage, ImageType::Pointer tmpImage );
		unsigned int * objectHistogram( ImageType::Pointer tmpImage );
		unsigned int otsuThreshold(int * histogram);
		
		// accomplishes 3D errosion of the watershed algorithm
		void binaryErosion(ImageType::Pointer workingImage, bool binaryBall, int size);
		void binaryDilation( ImageType::Pointer workingImage, bool binaryBall, int size );
		/*substractImages(a,b,c) gives image c with c=b-a                           */
		void substractImages(ImageType::Pointer image1, ImageType::Pointer image2, ImageType::Pointer image3);
		void binaryOpening( ImageType::Pointer workingImage, bool binaryBall, int size );
		void binaryClosing( ImageType::Pointer workingImage, bool binaryBall, int size );
		
		PointSetType::Pointer readBarrelMarkerFile(const char * markerFilename);
		void computeVoronoiDiagram(PointSetType::Pointer markers);
		void splitWrongBarrelObjects(CalcImageType::Pointer voronoiDiagram, unsigned int noOfCells, bool segProcessing);
		void watershedObjectSplitting(unsigned int * objectSplitFlag);
		void mergeTrueBarrelObjects(CalcImageType::Pointer distanceMap, CalcImageType::Pointer voronoiDiagram, unsigned int noOfCells, std::vector< Contour * > barrelContours, int zCoord, bool segProcessing);
		std::list< std::vector<float> > barrelContourExtraction(ImageType::Pointer barrelImage, long * offset, unsigned int minSize, int zCoord, Contour * contour);
		
		void segmentBarrels(unsigned int noOfPoints, std::vector< Contour * > barrelContours, int z, bool alreadyPreProcessed);
		void morphologicalForegroundRemoval( float rad, float arg2, float arg3 );
		// locates holes resulting from cropping other color channels from the current image and fills them iteratively with the mean of the surrounding pixels  
		void fillHoles2( ImageType::Pointer artifactImage );
		// locates holes in binary inputImage and fills them with foreground value
		void fillBinaryHoles();
		void fillBinaryHoles(ImageType::Pointer workingImage);
		void enlargeHoles( ImageType::Pointer workingImage );
		// calculates the mean of the border
		void calculateBorderPixels( ImageType::Pointer workingImage, std::list< ImageType::IndexType >& objectBorders, std::list< ImageType::IndexType >& objects, std::list< unsigned int >& borderMean );
		void voronoiRegionThreshold(unsigned int noOfRegions, std::vector< Contour * > barrelContours, int z);
		void optimizedVoronoiRegionThreshold(unsigned int noOfRegions, std::vector< Contour * > barrelContours, int z);
		void maxLevelNeighborhoodRegionGrowing(ImageType::Pointer segmentedImage, std::list< ImageType::IndexType > pxIndexList, int startLevel, int maxLevel, bool regular);
		void optimizedNeighborhoodRegionGrowing(ImageType::Pointer segmentedImage, std::list< ImageType::IndexType > pxIndexList, int startLevel, int maxLevel, Contour * contour);
		void fillEmptyBarrelContour(int barrelID, int z);
		
		float averageSphereAroundPxIDs(std::list< ImageType::IndexType > pxIDs);
		float avgSphereSNR(float radius, std::list< ImageType::IndexType > * pxIDs, std::list< ImageType::IndexType > * pxIndexList);
		
		void getWMContour( float alpha_arg, float beta_arg );
		ImageType::Pointer regionGrowingNeighborhoodThreshold( float epsilon, float delta, float gamma, int highThresh, int lowThresh );
		void createBvpImage ( ImageType::Pointer wmImage, ImageType::Pointer fgImage );
		
		void GlobalThreshold(int factor);
		ImageType::Pointer MarkBackground(unsigned int highThreshold);
		std::list< std::vector< float > > contourExtraction( ImageType::Pointer workingImage );
		GraphType createGraphFromContourImage( ImageType::Pointer workingImage );
		std::list< std::vector< int > > minimalDistanceGraph ( GraphType& graph, ImageType::Pointer workingImage, unsigned int samplingRate );
		void contourSmoothing( std::list< std::vector< float > >& contour, int radius );
		void contourSampling ( std::list< std::vector< float > >& contour, int samplingRate );
		
		BloodVesselPattern * ExtractBloodVesselPattern(int imageType, bool WM);
		void bvpPreparation( float alpha_arg, float beta_arg);
		void BloodVesselPatternToAmiraSpatialGraph();
		std::list<std::vector<float> > SetContourVertices(std::list<std::list<std::vector<float> > > contour_edges);
		int WriteSpatialGraphFile(const char * output_file);
		int WriteBVPSpatialGraphFile(const char * output_file);
		//not used anymore
		int WriteBarrelSpatialGraphFile(const char * output_file);
};

#endif

