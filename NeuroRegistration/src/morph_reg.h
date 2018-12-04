/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      morph_reg.h                                                   */
/*                                                                          */
/* Purpose:   class providing the methods for registration of 3D neuron     */
/*            morphologies to the standard barrel field                     */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   10.02.2011                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "utilities.h"

// #define PIPELINE_DOC

#ifndef REGISTRATION
#define REGISTRATION

class Registration : public Utilities
{
	public:
		Registration(AmiraSpatialGraph * inputSpatialGraph, bool registerLandmarks, int landmarkMode);
		Registration();
		~Registration();
		
		// main function used for initial registration of 3D neuron morphologies
		int startNeuronRegistration(const char * outputFilename, int axisSelect, int cellType, int goodApical,
									bool piaCorrection, bool onlyAxon, std::string manualHB);
		// DEPRECATED: main function used for computing neuron orientation
		void neuronMorphologyZAxis();
		
		// main function used for re-registration of 3D neuron morphologies
		// to different column (that is, after registration to home column)
		int registerToDifferentColumn(const char * outputFilename, const char * newHomeBarrelLabel, bool onlyAxon);
		
		// returns SpatialGraph object containing neuron morphology
		// and landmarks
		AmiraSpatialGraph * getSpatialGraph() { return spatialGraph; }
		
		// main function used for registration of barrel field contours
		void barrelFieldRegistration(const char * outputFilename, int mode, const char * refBarrel);
		
		// move this somewhere else eventually:
		// 3D reconstruction of manual barrel contours and measurement of
		// barrel parameters (for comparison with automatically reconstructed
		// barrel fields)
		void reconstructionError(const char * filename);
		
	private:
		AmiraSpatialGraph * workingSG;
		
		// global flag: use barrel-pia axis for registration
		// instead of all barrels giving orientation in barrel field
		bool ignoreBarrels;
		// flag: correct for systematic error in Pia distance
		// of artificially created Pia
		// optional for cells that lie below artificial Pia
		bool correctPia;
		double correctPiaDist;
		// correction for systematic error in manual contours
		// holds correction value [micron] in order:
		// supra, gran, infra
		std::map< int, std::vector< double > > contourCorrection;
		void setUpContourCorrection(int mode);
		
		// log file containing all information about transformations
		// and positions computed during registration
		std::ofstream TransformLog;
		
		#ifdef PIPELINE_DOC
		const char * pipelineDocName;
		void writeTransformSG(TransformPointerType transf, const char * label);
		void writeRegStepParameters(const char * label);
		PolyDataPointerType globalPia;
		#endif
		#ifdef REG_ACCURACY
		const char * logFileName;
		double angleUncertainty;
		void regVariability();
		#endif
		
		// methods for 3D reconstruction of barrels
		void barrelReconstruction(PolyDataPointerType pia, double * neuronAxis, std::map< int, double * > * barrelAxes, std::map< int, Column * > * barrels);
		void computeAverageHomeBarrel(PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3], int ID);
		std::vector< double > computeManualBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double * newAxis, double * barrelCenter, std::vector< double * > endPoints, int label, std::map< int, PolyDataPointerType >& avgBarrels);
		
		// methods for registration of landmarks to standard barrel field
		TransformPointerType landmarkRegistration(std::map< int, Column * > * barrels);
		TransformPointerType landmarkRegistration(std::map< int, Column * > * barrels, int mode, int refBarrelID);
		HomogeneousMatrixPointerType alignHomeBarrel(std::map< int, Column * > * barrels);
		HomogeneousMatrixPointerType alignHomeBarrel2(std::map< int, Column * > * barrels);
		HomogeneousMatrixPointerType alignRemainingBarrels(std::map< int, Column * > * barrels);
		void alignHomeBarrelAxis(std::map< int, Column * > * barrels, int HBID);
		
		// methods to automatically determine home barrel after registration
		int registeredHomeBarrel();
		int closestAvgBarrel();
		int closestMorphBarrel();
		
		// transformation of neuron morphology after registration of landmakrks,
		// i.e. rotation of neuron morphology around soma and stepwise z-scaling
		void morphologyRegistration(double * neuronAxis, PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels);
		// transformation of registered neuron morphology to new home barrel
		// includes translation, rotation and stepwise z-scaling
		TransformPointerType transformToNewColumn(double * neuronAxis, int newHomeBarrel);
		
		// various methods to compute position and distances
		// of parts of the neuron morphology during/after registration
		// self-explanatory
		double registeredPiaDistance();
		void registeredSomaPosition(const char * label);
		void registeredSomaPosition(const char * label, PolyDataPointerType pia);
		double somaDistanceToColumnAxis();
		double registeredPiaDendriteDistance();
		double registeredPiaDendriteDistance(PolyDataPointerType pia);
		double unregisteredPiaDendriteDistance(PolyDataPointerType pia);
		double registeredPiaSegmentDistance(int label, int edgeNr);
		double globalPiaSegmentDistance(int label, int edgeNr, PolyDataPointerType pia);
		int getColumnSeptumFlag();
		int getLaminarPosition();
		double morphologyZExtent();
		void maxDendPoint(double maxPoint[3]);
		
		// methods to automatically determine the neuron orientation based on
		// orientation of the apical dendrite where possible, otherwise 
		// using the main axon
		double * detectMorphZAxis(int label, int cellType, AmiraSpatialGraph * morphology);
		// computation of principal axes of dendrite morphology
		// returns axis of smallest moment of inertia
		double * principalAxes(int label, AmiraSpatialGraph * morphology, int biTufted);
		
		// implementation of the various correction factors
		// between manual and automatic landmark contours
		void correctArtificialPia(PolyDataPointerType pia);
		void constantPiaOffset(PolyDataPointerType pia, double offset);
		void constantLandmarkOffset(PolyDataPointerType pia, PolyDataPointerType wm, double offset);
		void correctManualContours(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column* > * barrels);
		void localLandmarkOffset(std::map< int, Column * > * barrels);
		
		// stepwise z-scaling of neuron morphology
		void inhomogeneousZScaling(double zScaling[3], double somaPt[3]);
		// morphology registration from experiment; parallel version
		// OpenMP directives not implemented correctly, but not a major problem
		// since registration usually takes < 1 min
// 		void inhomogeneousZScalingPRL(double zScaling[3], double somaPt[3]);
		
		// implementation of transformation (stepwise z-scaling) of single point coordinates
		void computeNewPosition(double * oldPt, double somaPt[3], double zScaling[3]);
		// DEPRECATED: actual transformation of point coordinates; verbose version for debugging
		void computeNewPosition(double * oldPt, double somaPt[3], double zScaling[3], bool verbose);
		// DEPRECATED: re-registration to different column; soma fixed point
		void inhomogeneousZScaling(double zScaling[3], double relSomaPt[3], double oldSomaPt[3]);
		// DEPRECATED: re-registration to different column; soma fixed point; transformation of point coordinates
		void computeNewPosition(double * oldPt, double relSomaPt[3], double oldSomaPt[3], double zScaling[3]);
		// landmark scaling in one direction only; for visualization
		void inhomogeneousZScaling(double zScaling[3], double refPts[4], double inputRefPts[4], PolyDataPointerType pia, PolyDataPointerType wm);
		void tmpCoordinateSystem(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels);
		void tmpCoordinateSystemInv(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels);
		
		// correction for z-scale in case landmark-based scale factors would lead to 
		// parts of neuron morphology "sticking out of" pia
		std::map< vtkIdType, std::vector< double > > buildLocalZScaleCorrection(double zScaling[3], double somaPt[3]);
		
		// helper methods
		// self-explanatory
		double * axisSurfaceIntersection(PolyDataPointerType surface, double * axis, double * center);
		double * axisSurfaceIntersection(CellLocatorPointerType surfaceLocator, double * axis, double * center);
		double * axisSurfaceIntersection(CellLocatorPointerType surfaceLocator, double * axis, double * center, vtkIdType& intersectCellID);
};



#endif
