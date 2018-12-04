#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "profile.h"


Profile::Profile()
{
	binSize = 50;
	integral = 0;
};

Profile::Profile(double binSize)
{
	this->binSize = binSize;
	integral = 0;
};

Profile::~Profile()
{
// 	profile.clear();
};

void Profile::addSegment(double length, unsigned int bin)
{
	if(bin < profile.size())
	{
		profile[bin] += length;
		integral += length;
	}
	else
	{
		int diff = bin - profile.size() + 1;
		for(int ii = 0; ii < diff; ++ii)
			profile.push_back(0);
		profile[bin] += length;
		integral += length;
	}
};

void Profile::incrementBin ( unsigned int bin )
{
	if(bin < profile.size())
	{
		profile[bin] += 1;
		integral += 1;
	}
	else
	{
		int diff = bin - profile.size() + 1;
		for(int ii = 0; ii < diff; ++ii)
			profile.push_back(0);
		profile[bin] += 1;
		integral += 1;
	}
}

void Profile::addProfile ( Profile* otherProfile )
{
	std::vector< double > * otherP = otherProfile->getProfile();
	unsigned int otherPLength = otherP->size();
	unsigned int thisPLength = this->profile.size();
	for(int ii = 0; ii < otherPLength; ++ii)
	{
		if(ii < thisPLength)
		{
			this->profile[ii] += (*otherP)[ii];
		}
		else
		{
			this->profile.push_back(0);
			this->profile[ii] += (*otherP)[ii];
		}
	}
	updateIntegral();
}

void Profile::updateIntegral()
{
	integral = 0;
	for(int ii = 0; ii < profile.size(); ++ii)
	{
		integral += profile[ii];
	}
}

std::vector< double >* Profile::getProfile()
{
// 	std::vector< double > outProfile(profile);
// 	return outProfile;
	return &profile;
};

void Profile::writeProfile ( const char* ofName, double binOffset )
{
	std::ofstream ProfileWriter;
	ProfileWriter.open(ofName);
	ProfileWriter << "# 1D Profile" << std::endl;
	ProfileWriter << "Depth [um]\tLength per bin [um]" << std::endl;
	for(int ii = 0; ii < this->profile.size(); ++ii)
		ProfileWriter << ii*binSize+binOffset << "\t" << this->profile[ii] << std::endl;
	ProfileWriter.close();
}

