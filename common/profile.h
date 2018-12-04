/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      profile.h                                                     */
/*                                                                          */
/* Purpose    Class implementing 1D profiles (essentially, histograms) and  */
/*            commonly used methods, e.g. integral, addition of profiles.   */
/*            These profiles are used during analysis of 3D axon/dendrite   */
/*            morphologies or 3D neuron soma distributions.                 */
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
#include "typedefs.h"

#ifndef PROFILE
#define PROFILE

class Profile
{
	public:
		// default constructor: binSize = 50 (micron)
		Profile();
		Profile(double binSize);
		~Profile();
		
		// add double value to bin with index 'bin'
		// automatically creates bins from 0 to 'bin'
		// if they do not exist yet
		void addSegment(double length, unsigned int bin);
		// add 1 to bin with index 'bin'
		// automatically creates bins from 0 to 'bin'
		// if they do not exist yet
		void incrementBin(unsigned int bin);
		// add values of otherProfile to this
		// automatically creates bins from 0 to max. bin
		// of otherProfile if they do not exist yet
		// whether this makes sense is the user's responsibility
		void addProfile(Profile * otherProfile);
		
		double getBinSize(){return binSize;}
		// computes sum of all bins
		double getIntegral(){return integral;}
		// update integral (e.g., after manipulation of bin values)
		void updateIntegral();
		// returns pointer to vector containing values
		// handle responsibly!
		std::vector< double >* getProfile();
		
		void writeProfile(const char * ofName, double binOffset=0);
		
	private:
		double binSize;
		double integral;
		std::vector< double > profile;
};

#endif