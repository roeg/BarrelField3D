/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      inputparameters.h                                             */
/*                                                                          */
/* Purpose:   simple class holding all input flags and parameters           */
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
#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H

class InputParameters
{
	public:
		InputParameters();
		~InputParameters();
		
		// Landmarks
		bool piaFlag, wmFlag;
		double piaSpacing;
		double wmSpacing;
		bool zReversed;
		// Cell structures
		bool somaFlag, dendriteFlag, apicalFlag, basalFlag, axonFlag;
};

#endif // INPUTPARAMETERS_H
