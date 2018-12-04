/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      matrixanalyzer.h                                              */
/*                                                                          */
/* Purpose:   Analysis methods for innervation matrices results from the    */
/*            NeuroNet barrel cortex model                                  */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max Planck Institute for Biological Cybernetics               */
/*            Spemannstr. 38-44                                             */
/*            72076 Tuebingen                                               */
/*            Germany                                                       */
/*                                                                          */
/* EMail:     robert.egger@tuebingen.mpg.de                                 */
/*                                                                          */
/* History:   17.06.2014                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/profile.h"

#ifndef MATRIXANALYZER_H
#define MATRIXANALYZER_H

typedef std::pair< unsigned int, unsigned int > ConnectionType;

class TripletMotif
{
	
public:
	TripletMotif(std::list< ConnectionType > connections);
	
	double computeOccurrenceProbability(std::vector< std::vector< double > > innervation);
	unsigned int computeOccurrencesWithStrength(std::vector< std::vector< unsigned int > > synapseRangeMatrix);
	
private:
	bool connections[3][3];
	unsigned int size;
};

class CellTriplet
{
public:
	CellTriplet(unsigned int preCell1, unsigned int preCell2, unsigned int preCell3, unsigned int postCell1, unsigned int postCell2, unsigned int postCell3);
	
	std::vector< unsigned int > preCellIndex;
	std::vector< unsigned int > postCellIndex;
	std::vector< std::vector< double > > innervation;
	
	void setInnervationMatrix(ConnectionMatrix * connectome);
};

class MatrixAnalyzer
{

public:
    MatrixAnalyzer();
    MatrixAnalyzer(ConnectionMatrix * connectome);
    virtual ~MatrixAnalyzer();
	
	void setConnectionMatrix(ConnectionMatrix * connectome){this->connectome = connectome;}
	
	std::vector< unsigned int > parseInputParameters(int argc, char* argv[], int argOffset);
	
	// distribution of innervation I_ij and derived quantities p_ij and n_ij
	// calculated from innervation matrix
	Profile * computeInnervationHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, float binSize);
	Profile * computeProbabilityHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, float binSize);
	Profile * computeSynapseNumberHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType);
	
	// compute distribution and convergence rate (i.e., redundancy)
	// of three-neuron motifs from innervation matrix
	void analyzeTripletMotifs(const char *  outputFilename, unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, unsigned int nrOfTriplets);
	std::vector< Profile * > computeTripletMotifDistribution(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, unsigned int nrOfTriplets);
	
	void writeConnectionMatrix(const char * outputFilename, unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType);
	void writeSynapsesPerCellRows(const char* outputFilename, std::list< unsigned int > postColumns, std::list< unsigned int > postTypes, CellTable* table);
	
private:
	ConnectionMatrix * connectome;
	std::map< unsigned, const char * > int2Labels;
	std::map< std::string, unsigned int > labels2Int;
	std::map< std::string, unsigned int > celltypeLabels2Int;
	std::map< unsigned int, const char * > int2CelltypeLabels;
	
	Profile * computePercentilePoissonHistogram(float innervation, double percentile);
	std::map< unsigned int, std::list< TripletMotif * > > initializeNonRedundantTripletMotifs();
	// randomly select cell triplets
	// constraint: no shared edge
	std::list< CellTriplet * > initializeNonRedundantCellTriplets(SelectionType preSelection, SelectionType postSelection, unsigned int nrOfTriplets);
	double computeHistogramRMSE(Profile * profile1, Profile * profile2);
	void initializeConstants();
};

#endif // MATRIXANALYZER_H
