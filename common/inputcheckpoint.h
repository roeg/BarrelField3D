#include "typedefs.h"
#include "basics.h"
#include "barrel_field.h"
#include "amiraReader.h"
#include "inputparameters.h"

#ifndef INPUTCHECKPOINT_H
#define INPUTCHECKPOINT_H

class InputCheckpoint
{
	public:
		InputCheckpoint();
		InputCheckpoint(AmiraSpatialGraph * inputSG);
		~InputCheckpoint();
		
		void run();
		void checkBarrelField();
		void checkNeuronMorphology();
		
		bool getInputOK(){return inputOK;}
		InputParameters getParameters();
	
	private:
		AmiraSpatialGraph * sg;
		std::list< int > barrelLabels;
		
		bool inputOK;
		bool piaFlag, wmFlag;
		double piaSpacing;
		double wmSpacing;
		bool zReversed;
		// Cell structures
		bool somaFlag, dendriteFlag, apicalFlag, basalFlag, axonFlag;
		
		void detectSectionThickness();
		void detectInputZDirection();
// 		void getLandmarkMinMaxIDs(PolyDataPointerType landmark, int& minID, int& maxID);
};

#endif // INPUTCHECKPOINT_H
