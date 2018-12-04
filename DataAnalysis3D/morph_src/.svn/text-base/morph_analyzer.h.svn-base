/****************************************************************************/
/*                                                                          */
/* Program:   MorphAnalyzer                                                 */
/*                                                                          */
/* File:      morph_analyzer.h                                              */
/*                                                                          */
/* Purpose:   Program for analysis of registered neuron morphologies with   */
/*            respect to columns, septa and layers in standardized barrel   */
/*            cortex. E.g., computes distribution of axon in different      */
/*            columns/layers, and computes z-profiles of axons taking local */
/*            orientation into account                                      */
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

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/barrel_field.h"
#include "../../common/amiraReader.h"
#include "../../common/profile.h"
#include "../../common/inputparameters.h"

#ifndef MORPH_ANALYZER
#define MORPH_ANALYZER

class Analyzer
{
	public:
		Analyzer(AmiraSpatialGraph * inputSpatialGraph, InputParameters parameters);
		~Analyzer();
		
		InputParameters parameters;

		// computes 1D profile of structure 'label'
		// along local axis
		std::vector< double >* compute1DProfileLocally(const char * labelStr, double binSize);
		// computes 1D profile of structure 'label' along home barrel axis.
		// not really used (systematic errors!), only for comparison
		// between analysis methods taking/not taking local orientation into account
		std::vector< double >* compute1DProfileGlobally(const char * labelStr, double binSize);
		// Only inside of S1
		std::vector< double >* compute1DProfileInsideS1(const char * labelStr, double binSize);

		
		// this is the standardized way for calculating axon parameters.
		// assumes thay neuron is registered to SBF D2 column
		void computeAxonClusterParameters(const char * outputFilename, double binSize);
		
		AmiraSpatialGraph * getNeuronMorphology() { return neuronMorphology; }
		
		void getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3]);

	private:
		AmiraSpatialGraph * neuronMorphology;
		BarrelField * SBF;
		std::map< std::string, int > neuronLabels2Int;
		PolyDataPointerType avgHBContour;
		
		// determines length assigned to which layer depending on 3D location
		// of the two points defining the line segment
		void assignDistanceToLayers ( double pt1[3], double pt2[3], int pt1Layer, int pt2Layer, std::map< int, Profile * >& pt1Profile, std::map< int, Profile * >& pt2Profile );
		// determines length assigned to which z bin depending on 3D location
		// of the two points defining the line segment
		void assignDistanceToBins(double pt1[3], double pt2[3], double topPt[3], double bottomPt[3], double axis[3], Profile * zProfile);
		// computes edge list indices of branch points
		// branch point is then at end of that edgePtList
		std::list< int > getBranchPointIDs ( int label );
		
		// helper methods
		void initializeConstants();
};

#endif
