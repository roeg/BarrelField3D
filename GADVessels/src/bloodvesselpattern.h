

/*
 *  bloodvesselpattern.h
 *  NeuroMorph
 *
 *  Created by Philip Broser on 4/25/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef BLOODVESSELPATTERN
#define BLOODVESSELPATTERN

#include "typedefs.h"
// #include "compartment.h"
#define VESSEL_X 0
#define VESSEL_Y 1
#define VESSEL_Z 2
#define VESSEL_DIAM 3
#define VESSEL_NR 4

struct AmiraContourGraphStruct
{
	AmiraContourGraphStruct(std::list<std::vector<float> > _vertice_list, std::list<std::list<std::vector<float> > > _edge_list);
	std::list<std::vector<float> > vertice_list;
	std::list<std::list<std::vector<float> > > edge_list;
};

typedef struct AmiraContourGraphStruct AmiraContourGraphType;

class Vessel
{
public:
	Vessel(int x, int y, float diam) {posX=x;posY=y;diameter=diam;};
	float distanceTo(int x, int y) {return sqrt((float) (posX-x)* (posX-x)+ (posY-y)*(posY-y) );};
	int GetX() {return posX;};
	int GetY() {return posY;};
	float GetDiam() {return diameter;};
	
private:
	int posX, posY;
	float diameter;
};

class  BloodVesselPattern
{
	public:
	std::list<Vessel> vesselList;

	BloodVesselPattern(const char * inputFilename, const char * outputFilename);
	BloodVesselPattern(const char * inputFilename, const char * outputFilename, int start, int stop);
	int addVessel(float x, float y, float z, float diam); //Values in physicla coordinates 
	
	int WriteVesselPatternToStream(std::ostream &ostr, unsigned int sizeY);
	int WriteVesselPatternAsHocToStream(std::ostream &ostr);
	std::list<std::list<std::vector<float> > > GetAmiraEdges();

	//vessel extractor Rootines
	void ExtractBloodVesselPattern(int threshold);
	int detectBloodVessel(CalcImage2DType::Pointer, int threshold);
	int plotBloodVesselPatternToImage(CalcImage2DType::Pointer img);
	
	private:
	int start, stop;
	const char * inputFilename;
	const char * outputFilename;
	ImageType::Pointer inputImage;
	void readImage( int start, int stop, const char * inputFilename );
	
	int writeVessel(std::ostream &ostr, std::vector<float> &vessel, unsigned int sizeY);
	int writeVesselAsHoc(std::ostream &ostr, std::vector<float> &vessel);
	std::list<std::vector<float> > GetAmiraEdge(std::vector<float> &vessel);

	std::list<std::vector<float> > blood_vessels;
	
	int idCounter;
	int vesselNr;
};

#endif




