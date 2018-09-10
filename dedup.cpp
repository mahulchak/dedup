#include<iostream>
#include "dedup.h"
#include "seqIO.h"

using namespace std;

using chroms = map<string,chromPair>;
using ccov = vector<int>;
using vq = vector<qord>;
int main(int argc, char *argv[])
{
	if(argc <2)
	{
		cerr<<"Usage: "<<argv[0]<<" foo.delta asm.fasta"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	
	map<string,ccov> masterRef; //stores sequence coverage but it can also be used to find reference chromosome lengths
	map<string,ccov>masterQ; //stores sequence coverage but it can also be used to find query chromosome lengths
	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
	map<string,string> seq;

	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	map<string,bool> qStrand; //stores whether query strand is forward strand or reverse strand
	mI tempmi;

	string line, chromName,refName,qName,indexAln;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,indelPos =0;
	
	vector<double> vd;
	vector<mI> qvm;//qvm is query sorted vm
	vector<string> altSeq;//store the seqname of alt haplotype
	size_t pos1,pos2,namePos;
	double hcr =0,hcq =0;
	ifstream fin, refFasta;
	ofstream fout,pctg,actg;
	fin.open(argv[1]);
	while(getline(fin,line))
	{
		
		if(line.find('>') != string::npos)//start of an aligning chromosome description
		{
						
			refName = line.substr(1,line.find(' ')-1);
			pos1 = line.find(' '); //position of the first space
			pos2 = line.find(' ',pos1+1);//position of the second space
			qName = line.substr(pos1+1, pos2-pos1-1); //up to the second space
//cout<<qName<<endl;
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;
			seqLen[indexAln].push_back(refLen);
			seqLen[indexAln].push_back(qLen);
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
			if(masterRef[refName].size() == 0)//if they have not been created
			{
				masterRef[refName] = makeChromBucket(refLen);
			}
			if(masterQ[qName].size() == 0)//if they have not been created
			{
				masterQ[qName] = makeChromBucket(qLen);
			}
		}
		if((line.size() <10) && (refName != "") && (count > -1))
		{
			indelPos = stoi(line);		
			if(indelPos ==0) //reached the end of the indel description
			{
				storeCords(masterRef[refName],masterQ[qName],tempmi);
				allChrom[indexAln].mums.push_back(tempmi);
			}
				
			count++;
			
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
//cout<<line<<endl;
				tempmi.rn = refName;
				tempmi.qn = qName;		
				refStart = stoi(line,&pos1);
				refEnd = stoi(line.substr(pos1),&pos2);
				qStart = stoi(line.substr(pos1+pos2), &namePos);
				qEnd = stoi(line.substr(pos1+pos2+namePos));
				tempmi.x1 = refStart;
				tempmi.x2 = refEnd;
				tempmi.y1 = qStart;
				tempmi.y2 = qEnd;
				count = 0;

		}
	}
//cout<<"storage finished\t"<<allChrom.size()<<endl;
	fin.close();
	for(chroms::iterator it = allChrom.begin();it!= allChrom.end();it++)
	{
		qvm.clear();
		indexAln = it->first;
//cout<<indexAln<<endl;
		sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
		for(unsigned int i = 0; i<allChrom[indexAln].mums.size();i++)
		{
			tempmi = allChrom[indexAln].mums[i];
			vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn],0.3);
cout<<"ncm\t"<<indexAln<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<'\t'<<vd[0]<<'\t'<<vd[1]<<endl;
			if((vd[0] == 2) && (vd[1]==2))
			{
				qvm.push_back(tempmi);
				count = count + (tempmi.x2 - tempmi.x1);
cout<<"cm\t"<<indexAln<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<'\t'<<vd[0]<<'\t'<<vd[1]<<endl;
			}
		}
		if(qvm.size()>0)
		{
			sort(qvm.begin(),qvm.end());
			hcr = double(count)/double(qvm[qvm.size()-1].x2-qvm[0].x1);
			//hcq = double(count)/double(abs(qvm[qvm.size()-1].y1-qvm[0].y2));
			hcq = double(count)/double(seqLen[indexAln][1]);
			//if((hcr > 0.5) && (hcq > 0.5))//covers at least 50% of alignment and query length
			if(hcr > 0.5)
			{
				//cout<<indexAln<<"\tLength is\t"<<count<<"\t"<<qvm[qvm.size()-1].x2-qvm[0].x1<<"\t"<<hcr<<"\t"<<qvm[qvm.size()-1].y2<<'\t'<<qvm[0].y1<<'\t'<<hcq<<'\t'<<seqLen[indexAln][1]<<endl;
				if((find(altSeq.begin(),altSeq.end(),qvm[0].qn) == altSeq.end()) && \
				(qvm[0].qn != qvm[0].rn))
				{
					altSeq.push_back(qvm[0].qn);
					//cout<<indexAln<<"\tLength is\t"<<count<<"\t"<<qvm[qvm.size()-1].x2-qvm[0].x1<<"\t"<<hcr<<"\t"<<qvm[qvm.size()-1].y2<<'\t'<<qvm[0].y1<<'\t'<<hcq<<'\t'<<seqLen[indexAln][1]<<endl;
				}
			}
		}		
		qvm.clear();
		count = 0; //reset count for the next alignment
		allChrom[indexAln].mums.clear();
	}
	refFasta.open(argv[2]);
	readfasta(refFasta,seq);
	pctg.open("p_ctg.fasta");
	actg.open("a_ctg.fasta");
	for(map<string,string>::iterator it = seq.begin();it!=seq.end();it++)
	{
		chromName = it->first;
		if(find(altSeq.begin(),altSeq.end(),chromName) != altSeq.end())
		{
		//	actg<<">"<<chromName<<endl;
		//	actg<<it->second<<endl;
		}
		if(find(altSeq.begin(),altSeq.end(),chromName) == altSeq.end())
		{
		//	pctg<<">"<<chromName<<endl;
		//	pctg<<it->second<<endl;
		}
	}				
	pctg.close();
	actg.close();
	refFasta.close();
return 0;
}
			


