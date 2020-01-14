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
		cerr<<"Usage: "<<argv[0]<<" foo.delta ref.fasta query.fasta cutoff mode(h/l) last_out.txt prefix"<<endl;
		exit(EXIT_FAILURE);
	}

	chroms allChrom;
	
	map<string,ccov> masterRef; //stores sequence coverage but it can also be used to find reference chromosome lengths
	map<string,ccov>masterQ; //stores sequence coverage but it can also be used to find query chromosome lengths
	map<string,ccov>masterHQ; //same as masterQ but records coverage only for homologous pairs

	map<string,int> lookUpRef;//encodes reference names in integers to save space
	map<string,int> lookUpQ;//same as above for query

	map<string,ccov> chromDensityRef;//stores whether a position in ref maps to >2 queries
	map<string,ccov> chromDensityQ;//stores whether a position in query maps to >2 refs

	map<string,vector<string> > cp; //cp is an alias for Chromosome partner. Each reference name index has a vector of unqiue alignments which are part of these
	map<string,vector<string> > hcp;//hcp stands for homologous cp
	
	map<string,map<int,vq> > umRef;//stores the coordinates of unique reference to query map; requires re-reading the file
	map<string,string> refseq;
	map<string,string> qseq;
	map<string,vector<int> > seqLen;//length of sequences.first element is ref and second is query
	
	mI tempmi,prevmi,nextmi,temprmi,tempmi2;

	string foo = string(argv[1]);
	string line, chromName,refName,qName,indexAln,fileName;
	int refStart = 0, refEnd = 0, qStart = 0, qEnd = 0, refLen =0, qLen =0, count = -1,qGap =0, indelPos =0, refChromCount=0, qChromCount = 0;
	unsigned int index = 0;
	
	vector<double> vd(2),vn(2);
	vector<int> vi;
	vector<mI> vmi,tempVmi,vm,qvm,gapmi,mastermi,transmi,invmi;//mastermi holds all mums in a delta file
	size_t pos1,pos2,namePos;
	
	ifstream fin, refFasta, qFasta,flast;//flast for opening lastz output
	ofstream fout,fcords,fcm;
	fin.open(argv[1]);
	fileName = "cords."+string(argv[6])+".txt";
	fcords.open(fileName);
	while(getline(fin,line))
	{
		
		if(line.find('>') != string::npos)//start of an aligning chromosome description
		{
						
			refName = line.substr(1,line.find(' ')-1);
			pos1 = line.find(' '); //position of the first space
			pos2 = line.find(' ',pos1+1);//position of the second space
			qName = line.substr(pos1+1, pos2-pos1-1); //up to the second space
			pos1 = line.find(' ',pos2+1); //recycling pos1 to find pos3
			refLen = stoi(line.substr(pos2+1,pos1-pos2));//reference length
			qLen = stoi(line.substr(pos1));//from last space till end 
			indexAln = refName + qName;
			count = -1;
			seqLen[indexAln].push_back(refLen);
			seqLen[indexAln].push_back(qLen);
			cp[refName].push_back(indexAln); //adding the alignment to the list of refName alignments
			if(lookUpRef[refName] == 0)
			{
				lookUpRef[refName] = ++refChromCount;
			}
			if(lookUpQ[qName] == 0)
			{
				lookUpQ[qName] = ++qChromCount;
			}
			if(masterRef[refName].size() == 0)//if they have not been created
			{
				masterRef[refName] = makeChromBucket(refLen);
				chromDensityRef[refName] = makeChromBucket(refLen);	
			}
			if(masterQ[qName].size() == 0)//if they have not been created
			{
				masterQ[qName] = makeChromBucket(qLen);
				chromDensityQ[qName] = makeChromBucket(qLen);
			}
		}
		if((line.size() <10) && (refName != "") && (count > -1))
		{
			indelPos = stoi(line);
			if(indelPos ==0) //reached the end of the indel description
			{
				mastermi.push_back(tempmi);//add to the master list
				storeCords(masterRef[refName],masterQ[qName],tempmi);
				if(tempmi.rn != tempmi.qn)//only record the non-self alignments
				allChrom[indexAln].mums.push_back(tempmi);
				storeNameCount(chromDensityRef[refName],chromDensityQ[qName],lookUpRef,lookUpQ,tempmi);
			}
			count++;
		}
		if((line.find('>') == string::npos) && (line.size() >10) && (refName != "")) //when describing alignment segments
		{
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
				tempmi.l = refEnd - refStart;
				count = 0;
		}
	}//complete reading the delta file
	fin.close();
	fileName = "cm."+string(argv[6])+".txt";
	fcm.open(fileName);
	//it takes 1min and 17s for dros genome to finish this 
/////////////////////LASTZ input//////////////////////
	for(chroms::iterator it = allChrom.begin();it!= allChrom.end();it++)
	{
		vm.clear();
		qvm.clear();
		tempVmi.clear();
		invmi.clear();
		indexAln = it->first;
		count = 0;
		sort(allChrom[indexAln].mums.begin(),allChrom[indexAln].mums.end());
		if(allChrom[indexAln].mums.size() >2)
		{
			for(unsigned int i = 0; i<allChrom[indexAln].mums.size();i++)
			{
				tempmi = allChrom[indexAln].mums[i];
				findInnie(allChrom[indexAln].mums,tempmi);
				vn = getChromCount(tempmi,chromDensityRef[tempmi.rn],chromDensityQ[tempmi.qn]);
				if(tempmi.c != 'r' && tempmi.c != 'q' && tempmi.c != 'd')
				{
					if(vn[0] <0.01 && vn[1] <0.01)
					{
						transmi.push_back(tempmi);
					}
					else
					{
						allChrom[indexAln].ncm.push_back(tempmi);
					}
				}
				else
				{
					allChrom[indexAln].ncm.push_back(tempmi);
//cout<<"ncm\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<qvm.size()<<endl;					
				}
				vd = getCoverage(tempmi,masterRef[tempmi.rn],masterQ[tempmi.qn],0.3);
				fcords<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<tempmi.c<<'\t'<<vn[0]<<'\t'<<vn[1]<<'\t'<<vd[0]<<'\t'<<vd[1]<<endl;
			}
//######################filtered all shadowed mums#########################//
			sort(allChrom[indexAln].ncm.begin(),allChrom[indexAln].ncm.end());
///////////////////////////////testing this part/////////////////////////		
			tempVmi = transmi;
			sort(tempVmi.begin(),tempVmi.end(),qusort);
			qvm.clear();//experimental
			sort(transmi.begin(),transmi.end());
			qvm.push_back(transmi[0]);//adding the first 
			indelPos =0;
			for(unsigned int i = 0; i<transmi.size();i++)
			{
//cout<<transmi[i].rn<<'\t'<<transmi[i].x1<<'\t'<<transmi[i].x2<<'\t'<<transmi[i].qn<<'\t'<<transmi[i].y1<<'\t'<<transmi[i].y2<<'\t'<<transmi.size()<<'\t'<<i<<endl;
//cout<<tempVmi[i].rn<<'\t'<<tempVmi[i].x1<<'\t'<<tempVmi[i].x2<<'\t'<<tempVmi[i].qn<<'\t'<<tempVmi[i].y1<<'\t'<<tempVmi[i].y2<<'\t'<<tempVmi.size()<<'\t'<<i<<'\t'<<double(tempVmi[i].x2)/double(tempVmi[i].y2)<<endl;
				//indelPos = findTrans(tempVmi,transmi[i]);
//cout<<transmi[i].rn<<'\t'<<transmi[i].x1<<'\t'<<transmi[i].x2<<'\t'<<transmi[i].qn<<'\t'<<transmi[i].y1<<'\t'<<transmi[i].y2<<'\t'<<transmi.size()<<'\t'<<i<<'\t'<<indelPos<<endl;				
				tempmi = transmi[i];
				nextmi = returnMumByQ1(tempmi.y1,tempVmi);//or nextmi
				prevmi = returnMumByQ2(tempmi.y1,tempVmi);//or prevmi
				if((i>0) && (tempmi.y1 < tempmi.y2 && nextmi.y1 < nextmi.y2) && (nextmi == transmi[i+1]) && (i!= transmi.size()-1) && (prevmi == transmi[i-1]))
				{
					tempmi.c = 'c';//syntenic mum
					qvm.push_back(tempmi);
//cout<<"this-cm\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<endl;
				}
				if((i>0) && (i!= transmi.size()-1) && (tempmi.y1 > tempmi.y2 && prevmi.y1 >prevmi.y2) && (prevmi == transmi[i+1]) && (nextmi == transmi[i-1]))
				{
					tempmi.c = 'c';
					qvm.push_back(tempmi);
				}
//cout<<"prev\t"<<prevmi.rn<<'\t'<<prevmi.x1<<'\t'<<prevmi.x2<<'\t'<<prevmi.qn<<'\t'<<prevmi.y1<<'\t'<<prevmi.y2<<'\t'<<endl;
//cout<<"this\t"<<tempmi.rn<<'\t'<<tempmi.x1<<'\t'<<tempmi.x2<<'\t'<<tempmi.qn<<'\t'<<tempmi.y1<<'\t'<<tempmi.y2<<'\t'<<qvm.size()<<endl;
//cout<<"next\t"<<nextmi.rn<<'\t'<<nextmi.x1<<'\t'<<nextmi.x2<<'\t'<<nextmi.qn<<'\t'<<nextmi.y1<<'\t'<<nextmi.y2<<'\t'<<endl;
			}
			transmi.clear();
///////////////////////////////end of test//////////////////////////////
			tempVmi = qvm;
			if(qvm.size()>1)
			{
				sort(allChrom[indexAln].ncm.begin(),allChrom[indexAln].ncm.end());
				sort(tempVmi.begin(),tempVmi.end(),qusort);
				sort(qvm.begin(),qvm.end());
				for(unsigned int i = 0;i<qvm.size()-1;i++)
				{
					prevmi = qvm[i];
					temprmi = qvm[i+1];
//cout<<tempVmi[i].rn<<'\t'<<tempVmi[i].x1<<'\t'<<tempVmi[i].x2<<'\t'<<tempVmi[i].qn<<'\t'<<tempVmi[i].y1<<'\t'<<tempVmi[i].y2<<endl;
					//tempmi = returnMumByQ2(prevmi.y1,tempVmi);
					allChrom[indexAln].cm.push_back(qvm[i]);
					if(temprmi.x1 - prevmi.x2 > 100)
					{
						if((prevmi.y1 > prevmi.y2) && (temprmi.y1 > temprmi.y2) && \
						 (prevmi.y2 - temprmi.y1 >100)) // both inverted
						{
							tempmi2.rn = prevmi.rn;
							tempmi2.x1 = prevmi.x2;
							tempmi2.x2 = temprmi.x1;
							tempmi2.qn = prevmi.qn;
							tempmi2.y1 = temprmi.y1;//maintain the forward strand style
							tempmi2.y2 = prevmi.y2;
							tempmi2.c = 'i';
						}
						if((prevmi.y1 < prevmi.y2) && (temprmi.y1 < temprmi.y2) && \
						(temprmi.y1 - prevmi.y2 > 100)) //both forward
						{
							tempmi2.rn = prevmi.rn;
							tempmi2.x1 = prevmi.x2;
							tempmi2.x2 = temprmi.x1;
							tempmi2.qn = prevmi.qn;
							tempmi2.y1 = prevmi.y2;
							tempmi2.y2 = temprmi.y1;
						}
						if((prevmi.y1 < prevmi.y2) && (temprmi.y1 > temprmi.y2))//first is forward but next is inverted
						{
							tempmi = returnMumByQ1(prevmi.y1,tempVmi);//get the next query
							if(min(tempmi.y1,tempmi.y2) - prevmi.y2 >100)
							{
								tempmi2.rn = prevmi.rn;
								tempmi2.x1 = prevmi.x2;
								tempmi2.x2 = temprmi.x1;
								tempmi2.qn = prevmi.qn;
								tempmi2.y1 = prevmi.y2;
								tempmi2.y2 = min(tempmi.y1,tempmi.y2);
							}
						}
					}//this is for x>100
					if((i>0) && (prevmi.y1 < prevmi.y2) && (qvm[i-1].y1 > qvm[i-1].y2) && \
						(prevmi.x1 - qvm[i-1].x2>100))
					{
						tempmi = returnMumByQ2(prevmi.y1,tempVmi);//get the previous query
//cout<<prevmi.rn<<'\t'<<prevmi.x1<<'\t'<<prevmi.x2<<'\t'<<prevmi.qn<<'\t'<<prevmi.y1<<'\t'<<prevmi.y2<<'\t'<<tempmi.y1<<endl;
						if(prevmi.y1 - max(tempmi.y1,tempmi.y2) > 100)
						{
							tempmi2.rn = prevmi.rn;
							tempmi2.x1 = qvm[i-1].x2;
							tempmi2.x2 = prevmi.x1;
							tempmi2.qn = prevmi.qn;
							tempmi2.y1 = max(tempmi.y1,tempmi.y2);
							tempmi2.y2 = prevmi.y1;
						}
					}	
					if(tempmi2.x1 != 0)
					{
//cout<<"GAP\t"<<tempmi2.rn<<"\t"<<tempmi2.x1<<"\t"<<tempmi2.x2<<"\t"<<tempmi2.qn<<"\t"<<tempmi2.y1<<"\t"<<tempmi2.y2<<endl;
						gapmi.push_back(tempmi2);
						tempmi2.x1 = 0; // reset it
					}				
				}//this is for the for loop
				allChrom[indexAln].cm.push_back(qvm[qvm.size()-1]);//add the last element
				if(allChrom[indexAln].cm.size() > 0)
				{
					for(unsigned int i =0;i<gapmi.size();i++)
					{	
						tempmi = gapmi[i];
						gapCloser(tempmi,allChrom[indexAln].ncm,allChrom[indexAln].cm);
						//gapCloser(tempmi,transmi,allChrom[indexAln].cm);
					}
					//hcp[allChrom[indexAln].cm[0].rn].push_back(indexAln);//homologous alignment
					sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());//sorting by length
					}
					sort(allChrom[indexAln].cm.begin(),allChrom[indexAln].cm.end());//sorting by length
					//add the lastz alignments to the ncm vector
//					allChrom[indexAln].ncm.insert(allChrom[indexAln].ncm.end(),allChrom[indexAln].last.begin(),allChrom[indexAln].last.end());
					transmi.clear();
					vm = allChrom[indexAln].cm;
					sort(vm.begin(),vm.end(),qusort);
					for(unsigned int i = 0;i<allChrom[indexAln].cm.size();i++)
					{
						tempmi = allChrom[indexAln].cm[i];
						count = count + (tempmi.x2-tempmi.x1);
						if(i>0)
						{
							qGap = qGap + abs(max(vm[i-1].y1,vm[i-1].y2) - min(vm[i].y1,vm[i].y2));
					//		avgDensity = avgDensity + ((tempmi.x2-tempmi.x1) * (double(tempmi.y2)/double(tempmi.x2))); 
						}
									}
					//avgDensity = avgDensity/count;
					cout<<indexAln<<'\t'<<count<<'\t'<<qGap<<'\t'<<double(count)/double(qGap)<<endl;
					if((double(count)/double(qGap) > 1) && (allChrom[indexAln].cm.size()>1))//if more than 1 element is present then use it for variant calling
					{
						for(unsigned int i = 0;i<allChrom[indexAln].cm.size();i++)
						{
							tempmi = allChrom[indexAln].cm[i];
							fcm<<tempmi.rn<<"\t"<<tempmi.x1<<"\t"<<tempmi.x2<<"\t"<<tempmi.qn<<"\t"<<tempmi.y1<<"\t"<<tempmi.y2<<'\t'<<endl;
						}
						hcp[allChrom[indexAln].cm[0].rn].push_back(indexAln);//homologous alignment
					}
				}					
				
			}//qvm.size()>0
			qvm.clear();
			tempVmi.clear();
		}//vm.size()>2
		count = 0; //reset count for the next alignment
		qGap = 0;
		allChrom[indexAln].mums.clear();
	}
	fcm.close();
	fcords.close();
	
return 0;
}
			


