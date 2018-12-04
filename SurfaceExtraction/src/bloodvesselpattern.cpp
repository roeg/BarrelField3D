/*
 *  bloodvesselpattern.cpp
 *  NeuroMorph
 *
 *  Created by Philip Broser on 4/25/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
("vessel"
  (Color RGB (255, 0, 128))
  (Closed)
  (Resolution 0.736015)
  (  405.07     9.40    55.84     0.74)  ;  1, 1
  (  401.85    35.95    55.84     0.74)  ;  1, 2
  (  393.19    45.67    55.84     0.74)  ;  1, 3
  (  382.20    48.81    55.84     0.74)  ;  1, 4
  (  373.29    44.55    55.84     0.74)  ;  1, 5
  (  369.36    29.89    55.84     0.74)  ;  1, 6
  (  371.32    15.87    55.84     0.74)  ;  1, 7
  (  385.10     3.86    55.84     0.74)  ;  1, 8
  (  394.68     4.42    55.84     0.74)  ;  1, 9 Closure point
)  ;  End of contour



("vessel"
  (Color RGB (255, 0, 128))
  (Closed)
  (Resolution 0.736015)
  ( -368.20  -345.86    55.84     0.74)  ;  6, 1
  ( -365.92  -341.48    55.84     0.74)  ;  6, 2
  ( -367.12  -326.00    55.84     0.74)  ;  6, 3
  ( -379.22  -302.23    55.84     0.74)  ;  6, 4
  ( -384.33  -299.94    55.84     0.74)  ;  6, 5
  ( -393.38  -312.29    55.84     0.74)  ;  6, 6
  ( -398.78  -326.18    55.84     0.74)  ;  6, 7
  ( -397.50  -337.25    55.84     0.74)  ;  6, 8
  ( -386.77  -355.10    55.84     0.74)  ;  6, 9
  ( -376.44  -353.81    55.84     0.74)  ;  6, 10
  ( -375.68  -352.35    55.84     0.74)  ;  6, 11 Closure point
)  ;  End of contour
*/



#include "bloodvesselpattern.h"

int verbose=0;

class VesselPoint
{
public:
	VesselPoint(int x, int y, float diam) {posX=x;posY=y;diameter=diam; neigbours=0;};
	float distanceTo(int x, int y) {return sqrt((float) (posX-x)* (posX-x)+ (posY-y)*(posY-y) );};
	int GetX() {return posX;};
	int GetY() {return posY;};
	float GetDiam() {return diameter;};
	int addNeigbour() {neigbours++;return neigbours;}
private:
	int posX, posY;
	float diameter;
	int neigbours;
};

int 
BloodVesselPattern::
writeVessel(std::ostream &ostr, std::vector<float> &vessel, unsigned int sizeY)
{
  float resolution= XYSAMPLING;
  unsigned int Ysize=sizeY;
  ostr<< "(\"vessel\"" << std::endl;
  ostr<< "  (Color RGB (255, 0, 128))" << std::endl;
  ostr<< "  (Closed)" << std::endl;
  ostr<< " (Resolution " << resolution << ")"<< std::endl;
  int Res=12;
  int counter=0;
  for (int i=1;i<=Res;i++)
    {


	float angle=(float) i*(2.0*PI)/Res;
	counter++;
	float posX=(vessel[VESSEL_X]+cos(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING;
	float posY=((sizeY-vessel[VESSEL_Y])+sin(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING;
	float posZ=(vessel[VESSEL_Z])*ZSAMPLING;
	int vesselNr=(int) vessel[VESSEL_NR];
	ostr<< "( " << posX << "\t " << posY << "\t "<< posZ << "\t 0.74) ; "<<  vesselNr << "," << i;
	if (i<Res) {ostr<<std::endl;} else { ostr<< " Closure point" << std::endl;};
    }
    ostr << ")  ;  End of contour" << std::endl << std::endl;

  return 0;
}

std::list<std::list<std::vector<float> > > BloodVesselPattern::GetAmiraEdges()
{
  std::list<std::list<std::vector<float> > > amira_edges;
  std::list<std::vector<float> > tmp_edge;

  std::list<std::vector<float> >::iterator it;
  for(it=blood_vessels.begin();it!=blood_vessels.end();++it)
  {
	tmp_edge = GetAmiraEdge(*it);
	amira_edges.push_back(tmp_edge);
  }

  return amira_edges;
};

std::list<std::vector<float> >  BloodVesselPattern::GetAmiraEdge(std::vector<float> &vessel)
{
  std::list<std::vector<float> > tmp_edge;

  float resolution= XYSAMPLING;
  int Res=12;

  for (int i=1;i<=Res;i++)
  {
  	std::vector<float> contour_point;

	float angle=(float) i*(2.0*PI)/Res;

	contour_point.push_back((vessel[VESSEL_X]+cos(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING*DOWNSAMPLING);
	contour_point.push_back((vessel[VESSEL_Y]+sin(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING*DOWNSAMPLING);
	contour_point.push_back((vessel[VESSEL_Z])*ZSAMPLING);

	tmp_edge.push_back(contour_point);
  }

  std::vector<float> contour_point;
  float angle=(float) (2.0*PI)/Res;

  contour_point.push_back((vessel[VESSEL_X]+cos(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING*DOWNSAMPLING);
  contour_point.push_back((vessel[VESSEL_Y]+sin(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING*DOWNSAMPLING);
  contour_point.push_back((vessel[VESSEL_Z])*ZSAMPLING);

  tmp_edge.push_back(contour_point);

  return tmp_edge;
};

int 
BloodVesselPattern::
writeVesselAsHoc(std::ostream &ostr, std::vector<float> &vessel)
{
int vesselNr=(int) vessel[VESSEL_NR];

ostr<< "{create vessel" << vesselNr << "}"<< std::endl;
ostr<< "{access vessel" << vesselNr << "}"<< std::endl;
ostr<< "{nseg = 1}"<< std::endl;
ostr<< "{strdef color color = \"Blue\"}"<< std::endl;
ostr<< "{pt3dclear()}"<< std::endl;



int Res=12;
int counter=0;
for (int i=1;i<=Res;i++)
{


	float angle=(float) i*(2.0*PI)/Res;
	counter++;
	float posX=(vessel[VESSEL_X]+cos(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING;
	float posY=(vessel[VESSEL_Y]+sin(angle)*(vessel[VESSEL_DIAM]/2.0))*XYSAMPLING;
	float posZ=(vessel[VESSEL_Z])*ZSAMPLING;	
	ostr << "{pt3dadd(" << posX << "," << posY << "," << posZ << "," << XYSAMPLING << ")}" << std::endl;
}
ostr << std::endl;

return 0;
}



int 
BloodVesselPattern::
WriteVesselPatternAsHocToStream(std::ostream &ostr)
{
std::list<std::vector<float> >::iterator it;
for(it=blood_vessels.begin();it!=blood_vessels.end();++it)
	writeVesselAsHoc(ostr, *it);
return 0;
}






int 
BloodVesselPattern::
WriteVesselPatternToStream(std::ostream &ostr, unsigned int sizeY)
{
  unsigned int Ysize = sizeY;
  std::list<std::vector<float> >::iterator it;
  for(it=blood_vessels.begin();it!=blood_vessels.end();++it)
  writeVessel(ostr, *it, Ysize);
  return 0;
}

int 
BloodVesselPattern::
addVessel(float x, float y, float z, float diam)
{
idCounter++;
std::vector<float> vessel;
vessel.push_back(x);
vessel.push_back(y);
vessel.push_back(z);
vessel.push_back(diam);
vessel.push_back(idCounter);
return 0;
}

//VesselExtractor Routines
//// Blood Vessel Detection
int checkPixel (CalcImage2DType::Pointer img, int x, int y, float threshold) 
//true if this is an pixel of a blood vessel, blood vessels are next to zero
{
  
  if (x<0 || x>(img->GetRequestedRegion().GetSize(0)-1) || y<0 || y>(img->GetRequestedRegion().GetSize(1)-1)) return 0;
  
  CalcImage2DType::IndexType index;
  index[0]=x;
  index[1]=y;
  if ( img->GetPixel(index) < threshold ) return 1;
  return 0;
}

int setValue (CalcImage2DType::Pointer img, int x, int y, float value) 
{
  
  if (x<0 || x>img->GetRequestedRegion().GetSize(0) || y<0 || y>img->GetRequestedRegion().GetSize(1)) return 0;
  
  CalcImage2DType::IndexType index;
  index[0]=x;
  index[1]=y;
  img->SetPixel(index,value);
  return 0;
}

std::list<VesselPoint> RegionGrowingFromPosition(CalcImage2DType::Pointer img, int x, int y, float threshold)
{
	std::list<VesselPoint> vesselFill;
	std::list<VesselPoint> fifo;
	int it=0, maxIt=5000000;
	
	fifo.push_back(VesselPoint(x, y, 0.0));
	
	while (fifo.begin()!=fifo.end() && it<maxIt)
	{
	it++;
	VesselPoint vesP=*(fifo.begin());
	fifo.pop_front();
	int presentX=vesP.GetX();
	int presentY=vesP.GetY();
	if ( checkPixel(img, presentX, presentY, threshold))
	{
		setValue(img, presentX, presentY, threshold+1);
		
		for (int dx=-1;dx<2;dx++)
			for (int dy=-1;dy<2;dy++)
			if (dx!=0 || dy!=0)
					if (checkPixel(img, presentX+dx, presentY+dy, threshold))
							{
								vesP.addNeigbour();
								fifo.push_back(VesselPoint(presentX+dx, presentY+dy, 0.0));
							};
		vesselFill.push_back(vesP);
	}


	}

	//Plt Vessel data to File
	std::list<VesselPoint>::iterator vesselIt;
	for (vesselIt=vesselFill.begin();vesselIt!=vesselFill.end();++vesselIt)
				setValue(img, (*vesselIt).GetX(), (*vesselIt).GetY(), 255);
		

	return vesselFill;

}

int CenterOfMassOfList(std::list<VesselPoint> vesselList, float &x, float &y)
{
	int nrOfPoints=0;
	float xPos=0.0;
	float yPos=0.0;
	std::list<VesselPoint>::iterator vesselIt;
	for (vesselIt=vesselList.begin();vesselIt!=vesselList.end();++vesselIt)
	{
		xPos+=(*vesselIt).GetX();
		yPos+=(*vesselIt).GetY();
		nrOfPoints++;
	}
	x=xPos/nrOfPoints;
	y=yPos/nrOfPoints;
	
	return  nrOfPoints;
	
}

int isHereAVesselXLeft(CalcImage2DType::Pointer img, int x, int y, int diameter, float threshold)
{
	for (int dx=0;dx<diameter; dx++)
			if (!checkPixel(img,x+dx,y,threshold)) return 0;
		
	for (int dx=0;dx<10*diameter; dx++)
		if (!checkPixel(img,x+dx,y,threshold)) return dx/2; //return dx Position 
	
	return diameter/2;
}

int isHereAVesselYCenter(CalcImage2DType::Pointer img, int x, int y, int diameter, float threshold)
{
	
	for (int dy=0;dy<diameter/2;dy++)
		{
			
			if (!checkPixel(img,x,y-dy,threshold)) return 0;
			if (!checkPixel(img,x,y+dy,threshold)) return 0;
		}

	int posY=0;
	int updy=0;
	int downdy=0;
	
	for (int dy=0;dy<diameter;dy++)
			if (!checkPixel(img,x,y+dy,threshold)) {updy=dy;break;};
	for (int dy=0;dy<diameter;dy++)
			if (!checkPixel(img,x,y-dy,threshold)) {downdy=dy;break;};
	
	posY=updy-(updy+downdy)/2;
	if (posY==0) posY=1;
	return posY;

}

int isHereAVesselXYCenter(CalcImage2DType::Pointer img, int x, int y, int diameter, float  threshold)
{
	
	for (int dx=0;dx<diameter/2; dx++)
	for (int dy=0;dy<diameter/2;dy++)
		{
			if (dx*dx+dy*dy>diameter/2*diameter/2) continue; 
			if (!checkPixel(img,x+dx,y-dy,threshold)) return 0;
			if (!checkPixel(img,x+dx,y+dy,threshold)) return 0;
			if (!checkPixel(img,x-dy,y-dy,threshold)) return 0;
			if (!checkPixel(img,x-dx,y+dy,threshold)) return 0;
		}

	return 1;
}

float getVesselDiameter(CalcImage2DType::Pointer img, int x, int y, float threshold)
{
	float diamter=0.0;
	int Height=img->GetRequestedRegion().GetSize(0);
	int Width=img->GetRequestedRegion().GetSize(1);
	for (int dx=0;dx+x<Width && x-dx>0  ; dx++)
		for (int dy=0;dy+y<Height && y-dy>0 ;dy++)
		{
			diamter=2.0*sqrt ((float) dx*dx+dy*dy);
			if (!checkPixel(img,x+dx,y-dy,threshold)) return diamter;
			if (!checkPixel(img,x+dx,y+dy,threshold)) return diamter;
			if (!checkPixel(img,x-dy,y-dy,threshold)) return diamter;
			if (!checkPixel(img,x-dx,y+dy,threshold)) return diamter;
		}

	return diamter;


}

#include <list>
int isThereAlreadyAVesselNextToThis(int posX, int posY, int diameter, std::list<BloodVessel> &vesselList)
{
	std::list<BloodVessel>::iterator pos;
    for (pos=vesselList.begin();pos!=vesselList.end();++pos)
      if ( (*pos).distanceTo(posX, posY) < diameter+(*pos).GetDiam() ) return 1;
		
	return 0;

}

int
BloodVesselPattern::
detectBloodVessel(CalcImage2DType::Pointer img, unsigned int contourType, bool WM, int samplingFactor)
{
	
	int Height=img->GetRequestedRegion().GetSize(1);
	int Width=img->GetRequestedRegion().GetSize(0);
	int stepDiameter=0;
	int minDiameter=0;
	int maxDiameter= 0;
	int threshold = 0;
	// 	int threshold=10;//brightfiled decon data, next to zero
	
	//float GradientThreshold=100;
	
	
	//imageDensityData redData;
	//imageDataToDensityDSelectCh(img, &redData, RED);
	int distanceToBoundary=0;

	if(contourType == 1)	//Pia & WM images, 4x air objective
	{
		stepDiameter=2;
		minDiameter=3;		// ~7 mu
		maxDiameter= 50;	// ~115 mu
		if(WM)
			threshold = 100;
		else
			threshold = 50;
	// 	int threshold=10;//brightfiled decon data, next to zero
		
		//float GradientThreshold=100;
		
		
		//imageDensityData redData;
		//imageDataToDensityDSelectCh(img, &redData, RED);
		distanceToBoundary=(int) sqrt((float) Height*Width)/10;
	}
	
	if(contourType == 2)	//Barrel images, 40x oil objective
	{
		stepDiameter=2; //Parameters tbd
		minDiameter=(int)(3*XYSAMPLING*samplingFactor*2 + 0.5);	// ~5.5 mu
		maxDiameter= (int)(50*XYSAMPLING*samplingFactor + 0.5);	// ~92 mu
		threshold = 10;
		// 	int threshold=10;//brightfiled decon data, next to zero
		
		//float GradientThreshold=100;
		
		
		//imageDensityData redData;
		//imageDataToDensityDSelectCh(img, &redData, RED);
		distanceToBoundary=(int) sqrt((float) Height*Width)/10;
	}
	
	for(int y=distanceToBoundary;y<Height-distanceToBoundary;y++)
		for (int x=distanceToBoundary;x<Width-distanceToBoundary;x++)
		{
		  //cout << "*" << flush;
			if (checkPixel(img,x,y,threshold)) 
			  {//cout << "-" << flush;
				
			  int deltaPosX=isHereAVesselXLeft(img, x,y, stepDiameter, threshold);
				if (deltaPosX)
				{
					if (verbose>2) std::cout << "isHereAVesselXLeft gives at:" << x << "," << y << " and deltaPosX:" << deltaPosX << std::endl;
				
					int posX= x+deltaPosX;
					int posY= y;
					//cout << "0" << flush;
					int deltaPosY=isHereAVesselYCenter(img,posX,posY, stepDiameter, threshold);
					//cout << "1" << flush;
					if (deltaPosY) // && ((float) diamX/(2.0*radiusY))<1.5 &&  0.5<((float) diamX/(2.0*radiusY))  )
						{
						posY=y+deltaPosY;
						if (isHereAVesselXYCenter(img,posX,posY, stepDiameter, threshold))
						  if (!isThereAlreadyAVesselNextToThis(posX, posY, minDiameter, vesselList) )
						    //&& VessellGradient (img, posX, posY, GradientThreshold, threshold,diameter))
							{
							  //cout << "2" << flush;
								std::list<VesselPoint> vessList=RegionGrowingFromPosition(img, posX, posY, threshold);
								//cout << "List Size:" << vessList.size() << endl;
								float XCenter;
								float YCenter;
								int nrOfPoints=CenterOfMassOfList(vessList, XCenter, YCenter);
								float vesselDiamter=sqrt( nrOfPoints/3.1415926) *2.0;
								if (verbose) std::cout << "Test Vessel with Diameter="<< vesselDiamter;
								if (verbose) std::cout << " XCenter=" << XCenter << "," <<  "YCenter=" << YCenter << std::endl;
								posX=(int) XCenter;
								posY=(int) YCenter;
											       
							
								if (vesselDiamter>minDiameter && vesselDiamter<maxDiameter)
								{
							
#ifdef DEBUG
								 if (verbose) std::cout << "Blood Vessel at Position x:" << posX << " y:" << posY << " Radius:" << vesselDiamter << std::endl; 
#endif
								BloodVessel newVesselObj(posX, posY, vesselDiamter);
								vesselList.push_back(newVesselObj);
								
								std::vector<float> newVessel;
								newVessel.push_back(posX/**samplingFactor*/);
								newVessel.push_back(posY/**samplingFactor*/);
								float posZ=0.0;
								newVessel.push_back(posZ);
								newVessel.push_back(vesselDiamter/**samplingFactor*/);
								vesselNr++;
								newVessel.push_back(vesselNr);
								blood_vessels.push_back(newVessel);
								};
							}
						};		
				};
					
					//x=x+minDiameter/2;
					
			}
		}
			
			/*
			// Debug
			//Plot Vessel to Image
			for(int y=0;y<img->Height();y++)
				for (int x=0;x<img->Width();x++)
				{
					Pixel myPix=img->GetPixel(x,y);
					myPix.r=0;
					(img->GetPixel(x,y)).Set (myPix.r, myPix.g , myPix.b);
				};
			
			list<BloodVessel>::iterator pos;
			for (pos=vesselList.begin();pos!=vesselList.end();++pos)
			{
				int x=(*pos).GetX();
				int y=(*pos).GetY();
				for (int dx=-4; dx<5; dx++)
				for (int dy=-4; dy<5; dy++)
				{
				Pixel myPix=img->GetPixel(x+dx,y+dy);
				myPix.r=255;
				(img->GetPixel(x+dx,y+dy)).Set (myPix.r, myPix.g , myPix.b);
				}
			}
			
			FILE *fileOut=fopen("vesselFile.bmp","wb");
			BMPWriteImage(img, fileOut);
			fclose(fileOut);
			
			*/
			//AmiraVesselWriter("AmiraTestFile.am", vesselList);
			
			std::list<std::vector<float> > sampleList;
			return 0; //vesselList;
}

int
BloodVesselPattern::
plotBloodVesselPatternToImage(CalcImage2DType::Pointer img)
{

  int Width=img->GetRequestedRegion().GetSize(0);
  int Height=img->GetRequestedRegion().GetSize(1);

  //cout << "Plot Vessel to Image" << endl;
  int i=0;
  std::list<BloodVessel>::iterator pos;
    for (pos=vesselList.begin();pos!=vesselList.end();++pos)
    {
      i++;
      // cout << "Plot Vessel" << i << endl;
		int X=(*pos).GetX(); 
		int Y=(*pos).GetY(); 
		float diam=(*pos).GetDiam();
		//cout << "x=" << X << " y=" << Y << " diam=" << diam << endl;
		int intDiam=(int) diam+10;
		for (int dx=(-1*intDiam);dx<intDiam;dx++)
		  for (int dy=(-1*intDiam);dy<intDiam;dy++)
		    {
		      

		      if ( sqrt((float) dx*dx + (float) dy*dy)>((float) 0.5*diam)  ) continue;
		      
		      int pX=X+dx;
		      int pY=Y+dy;
		      if (pX<0 || pX>Width || pY<0 || pY>Height ) continue; 
		      CalcImage2DType::IndexType index;
		      index[0]=pX;
		      index[1]=pY;
		      img->SetPixel(index,200);
		      //cout << "*";

		    }
		  

		
		




    };



    return 0;



};
