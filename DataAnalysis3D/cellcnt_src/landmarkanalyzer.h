/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      landmarkanalyzer.h                                            */
/*                                                                          */
/* Purpose:   Class for analysis of 3D neuron soma distributions with       */
/*            respect to barrel columns in cortex/barreloids in VPM         */
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

#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/barrel_field.h"
#include "../../common/amiraReader.h"
#include "../../common/profile.h"
#include "../../common/inputparameters.h"

#ifndef LANDMARKANALYZER_H
#define LANDMARKANALYZER_H

class LandmarkAnalyzer
{
	public:
		LandmarkAnalyzer();
		~LandmarkAnalyzer();
		
		void setPiaSurface(Surface * piaSurface);
		void setWMSurface(Surface * wmSurface);
		void setBarrelField(std::map< int, Column * > barrelField);
		void setBarrelFieldLayer(std::map< int, Column * > barrelField);
		void setSupragranularLayer(std::map< int, Column * > barrelField);
		void setGranularLayer(std::map< int, Column * > barrelField);
		void setInfragranularLayer(std::map< int, Column * > barrelField);
		void setLandmarkSet(PointsPointerType cellLandmarks);
		void setBarreloidField ( std::map< int, ClosedSurface * > barreloidField );
		void standardBFAnalysis();
		
		// main method used for computation of z profiles of all columns
		// taking overlap into account
		void computeColumnProfiles(const char * outputFilename);
		// compute column profiles and apply systematic
		// correction to match IN fraction with HSM 2011 results
		void computeCorrectedColumnProfiles(const char * outputFilename);
		// compute column profiles and apply systematic
		// correction to match IN fraction with HSM 2011 results
		// can be used with layer (supra/gran/infra)
		void computeCorrectedColumnProfilesInLayer(const char * outputFilename);
		// method used for counting neurons in different layers
		// taking overlap into account; layers as in HSM 2010/ MO 2011
		void countLaminarNeuronNumbers(const char* outputFilename, bool IN);
		// main method used for computation of z profiles of all columns
		// not taking overlap into account (i.e., somata may be counted several times)
		void computeSeparateColumnProfiles(const char * outputFilename);
		// main method used for counting somata inside/outside columns
		void countCellsInSeptum(Surface * S1Surface, const char * outputFilename);
		// main method used for counting somata inside/outside columns
		void countInhCellsInSeptumInLayer(Surface * S1Surface, const char * outputFilename);
		// main method used for counting somata inside barreloids
		void countCellsInBarreloids ( const char * outputFilename );
		// method for correcting IN density
		void correctINDensity(ImageDataPointerType density, const char * outputFilename);
		// method for cropping IN density in layers to C row/ arc2
		void cropINDensity(ImageDataPointerType density, const char * outputFilename);
		// method for cropping density to supra/gran/infra layers
		void cropDensityLayers(ImageDataPointerType density, const char * outputFilename);
		// method for cropping density to C row/ arc2
		void cropDensityCRowArc2(ImageDataPointerType density, const char * outputFilename);
		
	private:
		const char * outFilenameGlobal;
		// flags to check for anatomical structures
		bool piaFlag, wmFlag, barrelFlag, barrelLayerFlag, barreloidFlag, cellFlag;
		bool supraFlag, granFlag, infraFlag;
		Surface * pia, * WM;
		std::map< int, Column * > barrelField;
		std::map< int, Column * > barrelFieldLayer;
		std::map< int, Column * > supragranularLayer;
		std::map< int, Column * > granularLayer;
		std::map< int, Column * > infragranularLayer;
		std::map< int, ClosedSurface * > barreloidField;
		PointsPointerType cellSomata;
		BarrelField *  SBF;
		std::map< int, Profile * > columnSomaProfiles;
		// spacing of voxels for volume measurements
		double SPACING;
		
		void assignSomataToColumns(std::map< int, PointsPointerType >& somaColumns,
						std::map< int, Column * >& barrelColumns,
						std::map< int, std::vector< std::vector< double > > >& radialContours,
						std::map< int, double >& minDistances,
						std::map< int, double >& maxDistances);
		// creates z profile of all somata passed as argument
		// z range is given by 'top' and 'bottom'
		Profile * computeColumnSomaProfile(PointsPointerType somata, double top[3], double bottom[3]);
		// creates z profile of all voxels inside column 'ID'
		// taking overlap with other columns into account
		Profile * computeColumnVolumeProfile(std::map< int, Column * > barrelColumns, std::map< int, double > minDistances, 
						     std::map< int, double > maxDistances, 
						     std::map< int, std::vector< std::vector< double > > > radialContours, int ID);
		// creates z profile of all points (somata/voxels)
		// outside of barrel columns
		Profile * computeSeptumProfile(PointsPointerType points, std::map< int, Column * > barrelColumns);
		// creates z profile of all IN somata outside of barrel columns
		// and applies z-correction profile
		Profile * computeCorrectedSeptumProfile(PointsPointerType points, std::map< int, Column * > barrelColumns,
							std::map< int, Column * > layerColumns);
		
		// compute corrected z-profiles
		void correctSomaProfiles(std::map< int, Profile * > columnProfiles);
		
		// helper methods
		int computeLaminarPosition(double x[3], int columnID);
		
		void computeDensityProfile(Profile * rawProfile, Profile * volumeProfile, Profile * densityProfile);
		
		void createCountDataStructures(std::map< int, PointsPointerType >& somaColumns,
					       std::map< int, Column * >& barrelColumns,
					       std::map< int, std::vector< std::vector< double > > >& radialContours,
					       std::map< int, double >& minDistances,
					       std::map< int, double >& maxDistances);
		void createLayerCountDataStructures(std::map< int, PointsPointerType >& somaColumns,
					       std::map< int, Column * >& barrelColumns,
					       std::map< int, std::vector< std::vector< double > > >& radialContours,
					       std::map< int, double >& minDistances,
					       std::map< int, double >& maxDistances);
		
		void writeColumnSomaProfiles(const char * outputFilename);
		void writeZProfile(Profile * zProfile, const char * outputFilename, double globalOffset=0);
		
		Profile * createCorrectionProfile();
		
		ImageDataPointerType createImageVolume(double bounds[6]);
		void calculateExtent(double bounds[6], int extent[6]);
};

#endif // LANDMARKANALYZER_H
