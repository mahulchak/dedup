#include<iostream>
#include "dedup.h"
#include "seqIO.h"

using namespace std;
using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;

/////////////////////////////////////////////////////////
bool qusort(mI mi1, mI mi2)
{
	return (min(mi1.y1,mi1.y2) < min(mi2.y1,mi2.y2)) ||((min(mi1.y1,mi1.y2) == min(mi2.y1,mi2.y2)) && (max(mi1.y1,mi1.y2)<max(mi2.y1,mi2.y2)));
}
//////////////////////////////////////////////////////
ccov makeChromBucket(int refLen)
{
	ccov v;
	for(int i=0;i<refLen;i++)
	{
		v.push_back(0);
	}
return v;
}	
////////////////////////////////////////////////////////////
bool findInnie(vector<mI> & mums,mI mi)
{
	int i = 0;
	bool found = false;
	//while(!(mi.x2 >mums[i].x2))//until they become just on more than equal
	while((mi.x2 > mums[i].x1) && (i<mums.size()))
	{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
		if((mi.x1 > (mums[i].x1-1)) && (mi.x2 < (mums[i].x2+1)))
		{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
			if(!(mi == mums[i]))
			{
//cout<<"debug\t"<<mi.rn<<'\t'<<mi.x1<<'\t'<<mi.x2<<'\t'<<mums[i].rn<<'\t'<<mums[i].x1<<'\t'<<mums[i].x2<<endl;
				found = true;
				break;
			}
		}
		++i;
	}
	return found;
}
/////////////////////////////////////////////////////////////
bool findInnieQ(vector<mI> & mums,mI mi) // checks if mi query embeds into another query
{
	int i =0;
	bool found = false;
	while((max(mi.y2,mi.y1)) > (min(mums[i].y1,mums[i].y2)) && (i<mums.size()))
	{
		if((min(mi.y1,mi.y2) > (min(mums[i].y1,mums[i].y2)-1)) && (max(mi.y2,mi.y1) <(max(mums[i].y2,mums[i].y1)+1)))
		{
			if(!(mi == mums[i]))
			{
				found = true;
				break;
			}
		}
		++i;
	}
	return found;
}
////////////////////////////////////////////////////////////////////
void storeCords(ccov & masterRef,ccov & masterQ, mI & mi)
{

	int ty1 = 0, ty2 =0;
	
	if(mi.y1 > mi.y2)//if reverse oriented
	{
		ty1 = mi.y2;
		ty2 = mi.y1;
	}
	if(mi.y1 < mi.y2)//forward oriented
	{
		ty1 = mi.y1;
		ty2 = mi.y2;
	}
	for(int i = mi.x1-1; i<mi.x2;i++)
	{
		masterRef[i]++;
	}
	
	for(int j = ty1-1; j<ty2;j++)
	{
		masterQ[j]++;
	}
}
//////////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ, float p)
{
	int d = 0, cov = 0,covCount=0, medCov =0;
	double c;
	vector<double> cc;
	map<int,int> covFreq;//holds coverage frequency for the genomic interval
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t";
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		cov = cov + masterRef[i];	
		covFreq[masterRef[i]]++;
	}
	for(map<int,int>::iterator it = covFreq.begin();it!= covFreq.end();it++)
	{
		if((covCount <int(d*p)) && (d>5))
		{
			covCount = covCount + it->second;
			medCov = it->first;
		}
	}
//cout<<mi.x1<<"\t"<<mi.x2<<"\t"<<mi.y1<<"\t"<<mi.y2<<"\t"<<covCount<<"\t"<<medCov<<"\t";
	c = cov/double(d);
	if(d>5)
	{
		cc.push_back(double(medCov));
	}
	else
	{
		cc.push_back(c);
	}
	
	cov = 0;
	medCov = 0;
	covFreq.erase(covFreq.begin(),covFreq.end());//destroying the previous map
	covCount = 0;
	d = abs(mi.y1-mi.y2);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
		covFreq[masterQ[i]]++;
	}
	for(map<int,int>::iterator it = covFreq.begin();it!= covFreq.end();it++)
	{
		if((covCount <int(d*p)+1) && (d>5))
		{
			covCount = covCount + it->second;
			medCov = it->first;
		}
	}
//cout<<medCov<<endl;
	c = cov/double(d);
	if(d>5)
	{
		cc.push_back(double(medCov));
	}
	else
	{
		cc.push_back(c);	
	}
return cc;
//cout<<c<<endl;
}
		
//////////////////////////////////////////////////////
vector<double> getCoverage(mI & mi, ccov & masterRef, ccov & masterQ)
{
	int d = 0, cov = 0;
	double c;
	vector<double> cc;
	d = mi.x2 - mi.x1;
	for(int i = mi.x1-1;i<mi.x2;i++)
	{
		cov = cov + masterRef[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	cov = 0;
	d = abs(mi.y1-mi.y2);
	for(int i = min(mi.y1,mi.y2)-1;i<max(mi.y1,mi.y2);i++)
	{
		cov = cov + masterQ[i];
	}
	c = cov/double(d);
	cc.push_back(c);
	return cc;
}		
////////////////////////////////////////////////////
int nearestInt(double d)//returns the nearest integer
{	
	int in;
	in = int(d);
	if(in +0.5 >d) //if d was less than in.5
	{
		return in;
	}
	else
	{
		return in+1;
	}
}
////////////////////////////////////////////////////
void readfasta(ifstream & fin,map<string, string> & fastaseq) //reading fasta files
{
	string str,index;
	while(getline(fin,str))
	{
		if(str[0] == '>')
		{
			index = str.substr(1);
		}
		if(str[0] != '>')
		{
			fastaseq[index].append(str);
		}
	}
}

