// KMC program
// A3-B2 code

#include<iostream>
#include<cstdio>
#include<vector>
#include<utility>
#include<map>
#include<queue>
#include<algorithm>
#include<cmath>
#include<ctime>
#include<omp.h>
#include"tree.hh"
#define MATRIX_A 0x9908b0dfUL //for RNG
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

#define DIST_TYPE short

using namespace std;

typedef pair<size_t,size_t> TUPLE;

// logging options
char logfn[50];
bool LOG = true;
size_t logFreq = 5000;
FILE *flog;

// Global variables relating to KMC
//const size_t MAX=255;
const size_t MAX=(2<<sizeof(DIST_TYPE)*8-1);
const double Nbsq=46.50993;
double ca0=0.09221; // can be overwritten w/ argument input into main
double PAB,la; // initialized in initialize()
size_t MWfreq=10; // output/record frequency for DPw
double conversion=0.835;
const bool LOOPFREE=false; // loop free or not, affects probability and final conversion

 /* A4B4 settings (comment this line to use this setting)
const size_t 	fA=4,fB=4;
const char 		dA=1,dB=1;
const size_t 	mA=1,mB=1;
size_t 			NRA=1000; // can be overwritten w/ argument input into main
size_t			NRB=NRA*fA/fB;
// */

// /* A4B2 settings (comment this line to use this setting)
const size_t 	fA=4,fB=2;
const char 		dA=0,dB=1;
const size_t 	mA=0,mB=1;
size_t 			NRA=1000; // can be overwritten w/ argument input into main
size_t			NRB=NRA*fA/fB;
// */

char PATH[]="";
size_t suffix;

// KMC functions to call
void initialize();
void KMCstep();
void KMCconv(double conversion);
void output();

// RNG related
unsigned long mt[624];
int mti=625;
double fn[128], wn[128];
int kn[128];
void RNG_initialize(unsigned long);
unsigned long rand32();
size_t seed = 1000;
bool RndSeed = true;

// variables relating to KMC
size_t tmpJA,tmpJB; // temporary containers
vector<double> p;
//vector<vector<unsigned DIST_TYPE> > JunctionDist;
map< pair<size_t,size_t>, DIST_TYPE> newJunctionDist;
map< pair<size_t,size_t>, DIST_TYPE>::iterator mapIt;
vector<vector<size_t> > neighA,neighB;
vector<size_t> urA,urB; // unreacted A and B junctions, size changing
vector<double> sumA; // sum of probability relating to each A junction. size=NRA
double sum;

vector<size_t> molSize; // size of each molecule
vector<size_t> mol; // molecule idx for each junction; size=number of Junctions
size_t largestMol; // index of largest molecule
vector<double> DPw,DPwr; // holds DPw and reduced DPw for each step
vector<double> DPz,DPwX;
vector<double> DP4,DP5,DP6,DP7,DP8,DP9,DP10;
vector<double> Conv;
vector<double> SolFrac;
vector<double> nBranchs; // number of branching rxns

size_t loopVecSize = 255;
size_t writeLoopSize = 16;
vector<size_t> loop;
vector<double> loopfrac;
vector<size_t> loop1;
vector<size_t> loop2;
vector<size_t> loop3;
vector<size_t> loop4;

double DPs[] = {0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.735,0.74,0.745,0.75,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.8,0.805,0.806,0.807,0.808,0.809,0.81,0.811,0.812,0.813,0.814,0.815,0.816,0.817,0.818,0.819,0.82,0.821,0.822,0.823,0.824,0.825,0.826,0.827,0.828,0.829,0.83,0.831,0.832,0.833,0.834,0.835};
double NWPs[] = {0.73,0.74,0.75,0.76,0.77,0.78,0.785,0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83,0.835};
size_t nDPdist = sizeof(DPs)/sizeof(double);
vector<double> DPdistP(DPs,DPs+sizeof(DPs)/sizeof(double));
size_t nNW = sizeof(NWPs)/sizeof(double);
vector<double> NWP(NWPs,NWPs+sizeof(NWPs)/sizeof(double));
vector<vector<size_t> > DPcumdist;
vector<vector<size_t> > DDdist;
vector<vector<double> > DPcumdistNorm;

//vector<size_t> DDdist; // Distances between connected pairs
//vector<size_t> DPdist; // Distances between connected pairs given a primary loop is located at origin
//vector<size_t> DSdist; // Distances between connected pairs given a secondary loop is located at origin
//vector<size_t> PPdist; // Distances between pairs of primary loop containing junctions
//vector<double> fPPdist; // Distances between pairs of primary loop containing junctions
//vector<size_t> PSdist; // Distances between primary loop - secondary loop containing junctions
//vector<size_t> SPdist; // Distances between secondary loop - primary loop containing junctions
//vector<size_t> SSdist; // Distances between secondary loop containing junctions

/*
vector<size_t> DistSnapShot;
vector<vector<size_t> > Distances;
vector<size_t> DegreeSnapShot;
vector<vector<size_t> > Degrees;
*/

// helper functions
size_t dist(unsigned char);
int getJunctionDistAA(size_t,size_t);
int getJunctionDistBB(size_t,size_t);
void SelectJunct(size_t &,size_t &,size_t &,size_t &);
void UpdateSum(const size_t,const size_t,const size_t,const size_t);
void UpdateLoop(const size_t,const size_t);
//void CollectConnected(const size_t,const size_t,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,
//																vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&);
void CollectConnectedA(const size_t,const size_t,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&);
void CollectConnectedB(const size_t,tree<TUPLE>&);
//void UpdateJuncDist(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&);
//void UpdateJuncDist(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);
//void UpdateJuncDist(const size_t,const size_t,const tree<TUPLE>&,const tree<size_t>&);
void UpdateJuncDist(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const tree<TUPLE>&);
//void UpdateMol(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);
void UpdateMol(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&);
void updateWriteData(double);
void updateDPcum();
DIST_TYPE getNewJunctionDist(size_t,size_t);

// MAIN FUNCTION
int main(int argc,char* argv[])
{
    // handle argument input
	if(argc<2) { }
	else { ca0=atof(argv[1]); }
	if(argc==3) {
        NRA = atoi(argv[2]);
        NRB = atoi(argv[2])*fA/fB;
	} else if(argc==4) {
        NRA = atoi(argv[2]);
        NRB = atoi(argv[2])*fA/fB;
		suffix = atoi(argv[3]);
	}
sprintf(logfn,"%d.log",suffix); 
clock_t c;
c=clock();
	// RUN KMC simulation
    initialize();
//flog=fopen(logfn,"a");
//fprintf(stderr,"initialization done!\n");
//fprintf(stderr,"initialization time = %f\n",(double)(clock()-c)/CLOCKS_PER_SEC);
//    KMCconv(conversion);
//fprintf(stderr,"SEED = %d\n",seed);
cerr << seed << "\t";
    KMCconv(conversion);
//fprintf(stderr,"KMC time = %f\n",(double)(clock()-c)/CLOCKS_PER_SEC);
    output();
//fprintf(stderr,"output time = %f\n",(double)(clock()-c)/CLOCKS_PER_SEC);
cout << seed << "\t";
cout<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;
//	fclose(flog);
    return 0;
}

void output()
{
    FILE *fp;
    double c_star = ca0 * (6.02214e23/1.0e24)*pow(Nbsq,1.5);

    // set file name
    char fn[50]; // no more than 50 characters fore file name
//    char PATH[]="output/";
    char PATH[]="";
    char prefixLoopFrac[]="lpfrac_", prefixMW[]="MW_",prefixMWdist[]="cumDist_",prefixDegree[]="deg_",prefixDist[]="Dist_";

	// write loop fraction
    sprintf(fn,"%s%scs=%1.4fA%dB%d.txt",PATH,prefixLoopFrac,c_star,NRA,NRB);
	fp=fopen(fn,"a");
    fprintf(fp,"%.8f\t",c_star);
	for(size_t l=0;l<writeLoopSize;++l) 
		fprintf(fp,"%.8f\t",loopfrac[l]);
	// write loop fraction to screen
    printf("%.8f\t",c_star);
	for(size_t l=0;l<9;++l) printf("%.8f\t",loopfrac[l]);
//	printf("\n");
	fprintf(fp,"\n");
    fclose(fp);

    // write MW
    sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.txt",PATH,prefixMW,c_star,NRA,NRB,suffix);
    fp=fopen(fn,"a");
    fprintf(fp,"Conv\tDPw\tDPwr\tDPz\tX\tSolFrac\tBranchFrac\tLoop1\tLoop2\tLoop3\tLoop4\t");
	fprintf(fp,"DP4\tDP5\tDP6\tDP7\tDP8\tDP9\tDP10\n");
    for(size_t l=0;l<DPw.size();++l) {
		fprintf(fp,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%d\t%d\t%d\t%d\t",Conv[l],DPw[l],DPwr[l],DPz[l],DPwX[l],SolFrac[l],nBranchs[l],loop1[l],loop2[l],loop3[l],loop4[l]);
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",DP4[l],DP5[l],DP6[l],DP7[l],DP8[l],DP9[l],DP10[l]);
	}
    fprintf(fp,"\n");
    fclose(fp);


	// write MW distribution
	sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.txt",PATH,prefixMWdist,c_star,NRA,NRB,suffix);
	fp=fopen(fn,"a");
	fprintf(fp,"DP");
//	for(size_t i=(DPdistP.size());i>0;--i) fprintf(fp,"\t%1.3f",DPdistP[i-1]);
	for(size_t i=0;i<DPdistP.size();++i) fprintf(fp,"\t%1.3f",DPdistP[i]);
//	for(size_t i=(DPdistP.size());i>0;--i) fprintf(fp,"\t%1.3f",DPdistP[i-1]);
//	fprintf(fp,"\n");
//	fprintf(fp,"0");
//	for(size_t i=0;i<2*DPdistP.size();++i) fprintf(fp,"\t0");
	fprintf(fp,"\n");
	for(size_t i=0;i<DPcumdist[0].size()-1;++i) {
		fprintf(fp,"%d",i+1);
		for(size_t j=0;j<DPcumdist.size();++j) fprintf(fp,"\t%d",DPcumdist[j][i+1]);
//		for(size_t j=0;j<DPcumdistNorm.size();++j) fprintf(fp,"\t%f",DPcumdistNorm[j][i]);
		fprintf(fp,"\n");
	}
    fclose(fp);

	FILE *fp1;
	sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.txt",PATH,"DD_",c_star,NRA,NRB,suffix);
	fp1=fopen(fn,"a");
	fprintf(fp1,"D");
	for(size_t i=0;i<NWP.size();++i) fprintf(fp1,"\t%1.3f",NWP[i]);
	fprintf(fp1,"\n");
	for(size_t i=0;i<DDdist[0].size()-1;++i) {
		fprintf(fp1,"%d",i);
		for(size_t j=0;j<DDdist.size();++j) fprintf(fp1,"\t%d",DDdist[j][i]);
		fprintf(fp1,"\n");
	}
	fclose(fp1);

}

void writeDD(double currConv)
{ 
	size_t maxSize = 100;
	vector<size_t> dist(maxSize+1,0);//,cumdist(molSize.size(),0);
	for(size_t i=0;i<dist.size();++i) dist[i] = 0;
	size_t LL;
	for(size_t l=0;l<NRA;++l) {
		if(mol[l]!=largestMol) continue;
		for(size_t m=0;m<NRA;++m) {
			if(mol[m]!=largestMol) continue;
			if(m==l) continue;
				LL=getJunctionDistAA(l,m)/2;
				if(LL<maxSize+1) dist[LL]++;
		}
	}
	DDdist.push_back(dist);
}

void writeNW(double currConv)
{
	char fn[50];
    double c_star = ca0 * (6.02214e23/1.0e24)*pow(Nbsq,1.5);
	sprintf(fn,"%s%scs=%1.4fx=%1.2fA%dB%d_%02d.csv",PATH,"NW_",c_star,currConv,NRA,NRB,suffix);
	FILE *fp;
	fp=fopen(fn,"a");
	
	fprintf(fp,"Target,Source,Type\n");	

	for(size_t i=0;i<neighA.size();++i) {
//		fprintf(fp,"%d",i);
		for(size_t j=0;j<neighA[i].size();++j) {
			size_t AA = neighA[i][j];
			for(size_t k=0;k<neighB[AA].size();++k) {
				if(neighB[AA][k]>i)
				fprintf(fp,"A%d,A%d,Undirected\n",i,neighB[AA][k]);
//				fprintf(fp,",%d",neighB[AA][k]);
			}
		}
//		fprintf(fp,"\n");
	}
/*
	for(size_t i=0;i<neighA.size();++i) {
		fprintf(fp,"A%d",i);
		for(size_t j=0;j<neighA[i].size();++j) {
			fprintf(fp,",B%d",neighA[i][j]);
		}
		fprintf(fp,"\n");
	}
*/
}

void KMCstep()
{
    // JA JB holds the selected junction for this step
    // idxA idxB holds the index of JA JB in urA urB
    size_t JA,JB,idxA,idxB;

    // select the junction A and junction B to react
    SelectJunct(JA,JB,idxA,idxB);
//cout<<endl<<"JA,JB="<<JA<<","<<JB<<endl;
    // Updates the sum of relative probabilities of unreacted A-B pairs, by deleting the probabilities for the A and B reacting this step
    UpdateSum(JA,JB,idxA,idxB);

    // Update Loop information
    UpdateLoop(JA,JB);

    // Find all junctions connected to JA and JB

    // AcA holds all A junctions that is connected to JA
    // BcA holds all B junctions that is connected to JA
    // AcB holds all A junctions that is connected to JB
    // BcB holds all B junctions that is connected to JB
//    vector<size_t> AcA,BcA,AcB,BcB;
//	 vector<size_t> dAcA,dBcB,dAcB,dBcA;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
//    CollectConnected(JA,JB,AcA,AcB,BcA,BcB,dAcA,dAcB,dBcA,dBcB);
	vector<size_t> AcA,BcA,dAcA,dBcA;
	tree<TUPLE> trB;
	CollectConnectedA(JA,JB,AcA,BcA,dAcA,dBcA);
	CollectConnectedB(JB,trB);

    // Updates connectivity by updating neighbor list of JA and JB
    neighA[JA].push_back(JB);
    neighB[JB].push_back(JA);

    // update JunctionDist
    UpdateJuncDist(JA,JB,AcA,BcA,dAcA,dBcA,trB);
//    UpdateJuncDist(JA,JB,AcB,BcA,dAcB,dBcA);
//    UpdateJuncDist(JA,JB,AcA,BcB,dAcA,dBcB);
/*
cout<<"JA,JB = "<<JA<<","<<JB<<endl;
for(size_t i=0;i<NRA;++i) {
	for(size_t j=0;j<NRB;++j) cout<<JunctionDist[i][j]<<"\t";
	cout<<endl;
}
*/

    // Update the molecule grouping info and molecule sizes
//    UpdateMol(JA,JB,AcA,AcB,BcA,BcB);
    UpdateMol(JA,JB,AcA,BcA);
	
	trB.clear();
}

void KMCconv(double conversion)
{
    // do KMCstep until designated conversion
    size_t NA=NRA*fA,NB=NRB*fB;
    size_t start = NA>NB ? NB : NA;
    //double end = NA-((NRB-1)*fA+(NRA-NRB+1));a
    vector<double> tmpDPdistP,tmpNWP;
	for(size_t i=DPdistP.size();i>0;--i) tmpDPdistP.push_back(DPdistP[i-1]);
	for(size_t i=NWP.size();i>0;--i) tmpNWP.push_back(NWP[i-1]);
    double End = (1-conversion)*start;
	if(LOOPFREE) End = 1 + start - ( start/fA+start/fB );
//fprintf(stderr,"SEED = %d\n",seed);
    for(size_t i = start; i > End; i--) {
		double currConv = (double)(start-i+1)/(double)start;
//if((start-i) % logFreq == 0) fprintf(stderr,"reached i=%d \n",i);
//cout<<"i= "<<i;
        KMCstep();
//cout<<endl;
        if(!((start-i+1)%MWfreq)) updateWriteData(currConv);
		if(!tmpDPdistP.empty()) 
			if(currConv>=tmpDPdistP.back()) {
				updateDPcum();
				tmpDPdistP.pop_back();
			}
		if(!tmpNWP.empty())
			if(currConv>=tmpNWP.back()) {
				writeNW(currConv);
				writeDD(currConv);
				tmpNWP.pop_back();
			}

		if(sum < 1) break;
    }
//fprintf(stderr,"end KMC, into output\n");
    // calculate loop fraction from loop[]
    size_t totalbond = start - ceil(End);
	for(size_t l=0;l<loop.size();++l)
		loopfrac[l]=((double)loop[l])/totalbond;
}

void updateDPcum()
{
	size_t maxSize = 5000;
//	vector<size_t> dist(molSize.size()+1,0);//,cumdist(molSize.size(),0);
	vector<size_t> dist(maxSize+1,0);//,cumdist(molSize.size(),0);
//	vector<double> cumdistnorm(molSize.size(),0);
	size_t sum = 0;
	for(size_t i=0;i<molSize.size();++i) {
		if(molSize[i] <= maxSize) dist[molSize[i]]++;
//		sum += molSize[i];
	}
	DPcumdist.push_back(dist);
//	for(size_t i=0;i<molSize.size();++i) 
//		cumdistnorm[i] = (double)dist[i]/(double)sum;
//	DPcumdistNorm.push_back(cumdistnorm);
/*	
	size_t sum=0;
	for(size_t i=1;i<dist.size();++i) {
		sum += dist[i];
		cumdist[i-1] = sum;
	}
	for(size_t i=0;i<cumdist.size();++i) 
		cumdistnorm[i] = (double)cumdist[i]/(double)sum;
	DPcumdist.push_back(cumdist);
	DPcumdistNorm.push_back(cumdistnorm);
*/
}

void updateWriteData(double currConv)
{
    size_t NA=NRA*fA,NB=NRB*fB;
	double totalWWW=0,totalWW=0,totalW=0,W4=0,W5=0,W6=0,W7=0,W8=0,W9=0,W10=0;
	size_t maxSize = molSize[largestMol];
	for(size_t l=0;l<molSize.size();++l) {
		W4 += pow(molSize[l],4);
		W5 += pow(molSize[l],5);
		W6 += pow(molSize[l],6);
		W7 += pow(molSize[l],7);
		W8 += pow(molSize[l],8);
		W9 += pow(molSize[l],9);
		W10 += pow(molSize[l],10);
		totalWWW += molSize[l]*molSize[l]*molSize[l];
		totalWW += molSize[l]*molSize[l];
		totalW += molSize[l];
	}
	loop1.push_back(loop[1]);
	loop2.push_back(loop[2]);
	loop3.push_back(loop[3]);
	loop4.push_back(loop[4]);
	Conv.push_back(currConv);
	DPw.push_back(totalWW/totalW);
	DPz.push_back(totalWWW/totalWW);
	DP4.push_back(W4/totalWWW);
	DP5.push_back(W5/W4);
	DP6.push_back(W6/W5);
	DP7.push_back(W7/W6);
	DP8.push_back(W8/W7);
	DP9.push_back(W9/W8);
	DP10.push_back(W10/W9);
	DPwX.push_back((totalWW-(double)maxSize*maxSize)/((double)totalW)); // Critical factor
	DPwr.push_back((totalWW-(double)maxSize*maxSize)/(totalW-(double)maxSize)); // Critical factor
	SolFrac.push_back((totalW-(double)maxSize)/totalW);
//	DPwr.push_back((totalWW-maxSize*maxSize)/(totalW-maxSize)); // strictly DPwr
	nBranchs.push_back((double)loop[0]/(double)(NA));
//	Distances.push_back(DistSnapShot);
//	Degrees.push_back(DegreeSnapShot);
 
}

void initialize()
{
    // initialize PAB, la
    PAB=(1.0e24/6.02214e23)*pow((3.0/(2.0*3.14159*Nbsq)),1.5);
    la=PAB/ca0;
    // initialize matrix p for determining probability
    p = vector<double>(MAX+1,0);
	if(!LOOPFREE) {
    	for(size_t i=1;i<MAX;++i)
       		p[i] = 1.0 + la*NRA*fA*pow(i,-1.5);
    	p[MAX] = 1.0;
	}
    p[0] = 1.0;

    // initialize JunctionDist array
//    for(size_t i=0;i<NRA;++i)
//        JunctionDist.push_back(vector<unsigned DIST_TYPE>(NRB,0));
		

    // initialize neighA neighB
    for(size_t i=0;i<NRA;++i)
        neighA.push_back(vector<size_t>());
    for(size_t i=0;i<NRB;++i)
        neighB.push_back(vector<size_t>());

    // initialize unreacted A and unreacted B vectors
    for(size_t i=0;i<NRA;++i) urA.push_back(i);
    for(size_t i=0;i<NRB;++i) urB.push_back(i);

    // initialize cumulative probability
    sum = NRA*fA * NRB*fB;
    sumA = vector<double>(NRA,fA*NRB*fB);

    // initialize cluster size vectors
    molSize = vector<size_t>(NRA,mA);
	molSize.insert(molSize.end(),NRB,mB);
    // for A4B2 system, A junction has size 0
//    for(size_t i=0;i<NRA;++i) molSize.push_back(0);
//    for(size_t i=0;i<NRB;++i) molSize.push_back(1);
    // for A4B4 system, A,B has size
//    molSize = vector<size_t>(NRA+NRB,1);

    mol = vector<size_t>(NRA+NRB,0);
    largestMol=0;
	for(size_t i=0;i<mol.size();++i) mol[i] = i;
	size_t NNN = NRA*fA>NRB*fB ? NRB*fB : NRA*fA;
    DPw.reserve(NNN);
    DPz.reserve(NNN);
    DP4.reserve(NNN);
    DP5.reserve(NNN);
    DP6.reserve(NNN);
    DP7.reserve(NNN);
    DP8.reserve(NNN);
    DP9.reserve(NNN);
    DP10.reserve(NNN);
    DPwr.reserve(NNN);
    DPwX.reserve(NNN);
    SolFrac.reserve(NNN);
    nBranchs.reserve(NNN);
    loop1.reserve(NNN);
    loop2.reserve(NNN);
    loop3.reserve(NNN);
    loop4.reserve(NNN);

	DPcumdist.reserve(nDPdist);
	DPcumdistNorm.reserve(nDPdist);
	DDdist.reserve(nNW);
//	for(size_t i=nDPdist;i>0;--i) DPdistP.push_back(1.0/nDPdist*i);
//	for(size_t i=nNW;i>0;--i) NWP.push_back(1.0/nNW*i);
/*
	DistSnapShot = vector<size_t>(MAX+1,0);
	DistSnapShot[0] = (NRA+NRB)*(NRA+NRB-1)/2;
	Distances.reserve(NNN);
	DegreeSnapShot = vector<size_t>((fA>fB?fA:fB)+1,0);
	DegreeSnapShot[0] = NRA+NRB;
	Degrees.reserve(NNN);
*/
    // initialize loop vector
    loop = vector<size_t>(loopVecSize,0);
    loopfrac = vector<double>(loop.size(),0);

	// correlation relations vectors
//	DDdist = vector<size_t>(MAX+1,0);
//	PPdist = vector<size_t>(MAX+1,0);
//	fPPdist = vector<double>(MAX+1,0);
//	PSdist = vector<size_t>(MAX+1,0);
//	SPdist = vector<size_t>(MAX+1,0);
//	SSdist = vector<size_t>(MAX+1,0);
//	DPdist = vector<size_t>(MAX+1,0);
//	DSdist = vector<size_t>(MAX+1,0);

    // initialize RNG
    if(RndSeed == true)
	    seed = time(NULL)+clock()+suffix;
//seed = 100;
    RNG_initialize(seed);
}

void SelectJunct(size_t &JA,size_t &JB,size_t &idxA,size_t &idxB)
{
    double stop = rand32()/4294967296.0*sum; //rand32()/4294967296.0 is a uniformly distributed random number between 0 and 1
    double cump = 0;
    for(size_t i=0;i<urA.size();++i) {
        if(cump+sumA[urA[i]] >= stop) {
            JA = urA[i];
            idxA = i;
            break;
        } else {
            cump += sumA[urA[i]];
        }
    }
    for(size_t i=0;i<urB.size();++i) {
        JB = urB[i];
//        double tmpP = p[dist(JunctionDist[JA][JB])] * (double)(fB - neighB[JB].size()) * (double)(fA - neighA[JA].size());
        double tmpP = p[dist( getNewJunctionDist(JA,JB) )] * (double)(fB - neighB[JB].size()) * (double)(fA - neighA[JA].size());
        if(cump+tmpP >= stop) {
            idxB = i;
            break;
        } else {
            cump += tmpP;
        }
    }
}

void UpdateSum(const size_t JA,const size_t JB,const size_t idxA,const size_t idxB)
{
    double probAB=0;
    for(size_t i=0;i<urA.size();++i) {
        tmpJA = urA[i];
//        probAB = p[dist(JunctionDist[tmpJA][JB])] * (double)(fA - neighA[tmpJA].size());
        probAB = p[dist( getNewJunctionDist(tmpJA,JB) )] * (double)(fA - neighA[tmpJA].size());
        sum -= probAB;
        sumA[tmpJA] -= probAB;
    }
    for(size_t i=0;i<urB.size();++i) {
        tmpJB = urB[i];
//        probAB = p[dist(JunctionDist[JA][tmpJB])] * (double)(fB - neighB[tmpJB].size());
        probAB = p[dist( getNewJunctionDist(JA,tmpJB) )] * (double)(fB - neighB[tmpJB].size());
        sum -= probAB;
        sumA[JA] -= probAB;
    }
//    probAB = p[dist(JunctionDist[JA][JB])];
    probAB = p[dist(getNewJunctionDist(JA,JB))];
    sum += probAB;
    sumA[JA] += probAB;

    // erase element for urA/urB
    // if size == fA-1 or fB-1, then all nodes have been reacted for the junction
    // neighbor list not updated here because it affects the calculation of original distance between junctions
    
//	DegreeSnapShot[neighA[JA].size()]--;DegreeSnapShot[neighA[JA].size()+1]++;
//	DegreeSnapShot[neighB[JB].size()]--;DegreeSnapShot[neighB[JB].size()+1]++;

    if(neighA[JA].size() == fA-1) urA.erase(urA.begin()+idxA);
    if(neighB[JB].size() == fB-1) urB.erase(urB.begin()+idxB);
}

void UpdateLoop(const size_t JA,const size_t JB)
{
    for(size_t i=0;i<loop.size();++i) {
//        if(dist(JunctionDist[JA][JB]) == i) {
        if(dist(getNewJunctionDist(JA,JB)) == i) {
        // loop[0] is for new branch, not loop
        // for system with fA or fB =2, loop[n] is for n/2 order loop
        // for other system, loop[n] is for nth order loop
        // loop[n] should be zero for odd n
            loop[i]++;
            break;
        }
    }
}

void CollectConnectedA(const size_t JA,const size_t JB,vector<size_t> &AcA,vector<size_t> &BcA,vector<size_t> &dAcA,vector<size_t> &dBcA)
{
	for(size_t i=0;i<NRA;++i) 
		if(getJunctionDistAA(i,JA)>=0) {
			AcA.push_back(i);
			dAcA.push_back(getJunctionDistAA(i,JA));
		} 
	for(size_t i=0;i<NRB;++i) 
//		if(JunctionDist[JA][i]>0) {
		if(newJunctionDist.find(make_pair(JA,i)) != newJunctionDist.end() ) {
			BcA.push_back(i);
//			dBcA.push_back(JunctionDist[JA][i]);
			dBcA.push_back(newJunctionDist[make_pair(JA,i)]);
		}
/*
cout<<endl;
cout<<"AcA: ";
for(size_t i=0;i<AcA.size();++i) cout<<AcA[i]<<" ";
cout<<endl;
cout<<"BcA: ";
for(size_t i=0;i<BcA.size();++i) cout<<BcA[i]<<" ";
cout<<endl;
*/
}

void CollectConnectedB(const size_t J,tree<TUPLE> &tr)
{
	vector<bool> visitedA(NRA,false);
	vector<bool> visitedB(NRB,false);
	
	TUPLE t(0,J);
	tr=tree<TUPLE>(t);
	
	visitedB[J]=true;

	queue<tree<TUPLE>::iterator> Q;
	Q.push(tr.begin());
	tree<TUPLE>::iterator tmpIt;

	while(!Q.empty()) {
		tmpIt = Q.front();
		TUPLE tmpT = *tmpIt;
		Q.pop();
		size_t Jtype = (1 + tmpT.first) %2;
		size_t tmpJ = tmpT.second;
		if(Jtype == 0) {
			for(size_t i=0;i<neighA[tmpJ].size();++i) {
				if(!visitedB[neighA[tmpJ][i]]) Q.push(tr.append_child(tmpIt,TUPLE(tmpT.first+1,neighA[tmpJ][i])));
				visitedB[neighA[tmpJ][i]]=true;
			}
		}
		if(Jtype == 1) {
			for(size_t i=0;i<neighB[tmpJ].size();++i) {
				if(!visitedA[neighB[tmpJ][i]]) Q.push(tr.append_child(tmpIt,TUPLE(tmpT.first+1,neighB[tmpJ][i])));
				visitedA[neighB[tmpJ][i]]=true;
			}
		}
	}
/*
tmpIt=tr.begin();
while(tmpIt!=tr.end()) {
	for(size_t i=0;i<tr.depth(tmpIt);++i) cout<<" ";
	if(tr.depth(tmpIt)%2==0) cout<<"B";
	else cout<<"A";
	cout<<*tmpIt<<" ";
	cout<<endl;
	++tmpIt;
}
*/
}

void UpdateJuncDist(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &BcA,const vector<size_t> &dAcA,const vector<size_t> &dBcA,const tree<TUPLE> &trB)
{
size_t counter=0,counter2=0;
	size_t olddist,newdist,tmpJA,tmpJB;
	double probAB;
	tree<TUPLE>::iterator ItA,ItB;
	ItA = trB.begin(); 
	++ItA;
	if(ItA!=trB.end()) {
	for(size_t i=0;i<BcA.size();++i) {
		tmpJB = BcA[i];
		ItA = trB.begin();
		while(ItA != trB.end()) {
counter2++;
			if( (*ItA).first % 2 == 0 ) {++ItA;continue;}
			tmpJA = (*ItA).second;
//			olddist = JunctionDist[tmpJA][tmpJB];
			olddist = getNewJunctionDist(tmpJA,tmpJB);
			newdist = dBcA[i] + (*ItA).first + 1;
			if(olddist == 0 || newdist < olddist) {
counter++;
				if(newdist > MAX) newdist = MAX;
//				JunctionDist[tmpJA][tmpJB] = newdist;
				newJunctionDist[make_pair(tmpJA,tmpJB)] = newdist;
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                sum += probAB;
                sumA[tmpJA] += probAB;
			}
			else { ItA.skip_children(); }
			++ItA;	
		}
	}
	}
//cout<<"\t"<<counter<<"/"<<counter2;

counter=0,counter2=0;
	ItB = trB.begin();
	for(size_t i=0;i<AcA.size();++i) {
		tmpJA = AcA[i];
		ItB = trB.begin();
		while(ItB != trB.end()) {
counter2++;
			if( (*ItB).first % 2 == 1 ) {++ItB;continue;}
			tmpJB = (*ItB).second;
//			olddist = JunctionDist[tmpJA][tmpJB];
			olddist = getNewJunctionDist(tmpJA,tmpJB);
			newdist = dAcA[i] + (*ItB).first + 1;
			if(olddist == 0 || newdist < olddist) {
counter++;
				if(newdist > MAX) newdist = MAX;
//				JunctionDist[tmpJA][tmpJB] = newdist;
				newJunctionDist[make_pair(tmpJA,tmpJB)] = newdist;
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                sum += probAB;
                sumA[tmpJA] += probAB;
			}
			else { ItB.skip_children(); }
			++ItB;	
		}
	}
//cout<<"\t"<<counter<<"\t"<<counter2<<endl;
/*
	while(ItA != trA.end()) {
		size_t tmpJ1 = *ItA;
		ItB = trB.begin();
		while(ItB != trB.end()) {
counter2++;
			// trA.depth(ItA) is even if J1 is A type
			// trB.depth(ItB) is odd if J2 is A type
			// if of same type (trA.depth+trB.depth is odd), continue
			if( (trA.depth(ItA) + trB.depth(ItB) ) %2 == 1) {++ItB;continue;}
			size_t tmpJ2 = *ItB;
			if( trA.depth(ItA)%2 == 0) {tmpJA=tmpJ1;tmpJB=tmpJ2;}
			else {tmpJA=tmpJ2;tmpJB=tmpJ1;}
			olddist = JunctionDist[tmpJA][tmpJB];
			newdist = trA.depth(ItA) + trB.depth(ItB) + 1;
			if(olddist == 0 || newdist < olddist) {
counter++;
				if(newdist > MAX) newdist = MAX;
				JunctionDist[tmpJA][tmpJB] = newdist;
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                sum += probAB;
                sumA[tmpJA] += probAB;
			}
			else { ItB.skip_children(); }
			++ItB;	
		}
		++ItA;
	}
*/
}

/*
void UpdateJuncDist(const size_t JA,const size_t JB,const vector<size_t> &Ac,const vector<size_t> &Bc,
																	const vector<size_t> &dAc,const vector<size_t> &dBc)
{
size_t counter=0,counter2=0;
    size_t olddist,newdist;
    double probAB;
//    #pragma omp parallel for
    for(size_t i=0;i<Ac.size();++i) {
        tmpJA = Ac[i];
        for(size_t j=0;j<Bc.size();++j) {
            tmpJB = Bc[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = dAc[i] + dBc[j] + 1;
++counter2;
            if(olddist == 0 || newdist < olddist) {
++counter;
                if(newdist > MAX) newdist = MAX;
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
//                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
//                }
                continue;
            }
        }
    }
//cout<<"\t"<<counter<<"/"<<counter2;
}
*/
void UpdateMol(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &BcA)
{
    size_t mA,mB,newMol,delMol;
	mA = mol[JA];
	mB = mol[NRA+JB];
	// new molecule number should be the smaller of the two molecule numbers
	if(mA != mB) { // if mA==mB then this is intramolecular rxn, no change is needed
//        newMol = (mA<mB)?mA:mB;
		newMol = mB;
//        delMol = (mA<mB)?mB:mA;
		delMol = mA;
        for(size_t j=0;j<AcA.size();++j) mol[AcA[j]] = newMol;
        for(size_t j=0;j<BcA.size();++j) mol[NRA+BcA[j]] = newMol;
//        for(size_t j=0;j<AcB.size();++j) mol[AcB[j]] = newMol;
//        for(size_t j=0;j<BcB.size();++j) mol[NRA+BcB[j]] = newMol;
        molSize[newMol] = molSize[mA] + molSize[mB];
        molSize[delMol] = 0;
        if(molSize[newMol]> molSize[largestMol]) largestMol = newMol;
    }
}

int getJunctionDistAA(size_t J1,size_t J2)
{ // return -1 if J1 J2 are not connected, otherwise return minimal distance between them
    if(J1 == J2) return 0;
    if(neighA[J2].size()<1) return -1;
//    if(JunctionDist[J1][neighA[J2][0]] == 0) return -1;
    if( newJunctionDist.find(make_pair(J1,neighA[J2][0])) == newJunctionDist.end() ) return -1;
    size_t minD = MAX,JB;
    for(size_t i=0;i<neighA[J2].size();++i) {
        JB = neighA[J2][i];
//        if(JunctionDist[J1][JB] < minD)
        if(newJunctionDist[make_pair(J1,JB)] < minD)
            minD = newJunctionDist[make_pair(J1,JB)];
//            minD = JunctionDist[J1][JB];
    }
    return minD+1;
}

int getJunctionDistBB(size_t J1,size_t J2)
{
    if(J1 == J2) return 0;
    if(neighB[J2].size()<1) return -1;
//    if(JunctionDist[neighB[J2][0]][J1] == 0) return -1;
    if(newJunctionDist.find(make_pair(neighB[J2][0],J1)) == newJunctionDist.end() ) return -1;
    size_t minD = MAX,JA;
    for(size_t i=0;i<neighB[J2].size();++i) {
        JA = neighB[J2][i];
//        if(JunctionDist[JA][J1] < minD)
        if(newJunctionDist[make_pair(JA,J1)] < minD)
            minD = newJunctionDist[make_pair(JA,J1)];
//            minD = JunctionDist[JA][J1];
    }
    return minD+1;
}

size_t dist(unsigned char JunctionDist)
{
    // JunctionDist should always be odd, since only different types of junctions react
    if(JunctionDist == 0) return 0;
    else return (JunctionDist+1)/2*(dA+dB);
}

DIST_TYPE getNewJunctionDist(size_t JA,size_t JB)
{
	mapIt = newJunctionDist.find(make_pair(JA,JB));
	if (mapIt == newJunctionDist.end()) return 0;
	else return mapIt->second;
}

//Initializes random number generator with seed
//RNG is Mersenne Twister MT19937 algorithm
void RNG_initialize(unsigned long seed)
{
	mt[0]= seed & 0xffffffffUL;
    for(mti=1; mti<624; mti++)
	{
        mt[mti] = (1812433253UL*(mt[mti-1]^(mt[mti-1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
	double dn = 3.442619855899;
	int i;
	const double m1 = 2147483648.0;
	double q;
	double tn = 3.442619855899;
	const double vn = 9.91256303526217E-03;

	q = vn/exp(-0.5*dn*dn);

	kn[0] = (dn/q)*m1;
	kn[1] = 0;

	wn[0] = q/m1;
	wn[127] = dn/m1;

	fn[0] = 1.0;
	fn[127] = exp(-0.5*dn*dn);

	for(i = 126; i >= 1; i--)
	{
		dn = sqrt(-2*log(vn/dn + exp(-0.5*dn*dn)));
		kn[i+1] = (dn/tn)*m1;
		tn = dn;
		fn[i] = exp(-0.5*dn*dn);
		wn[i] = dn/m1;
	}
}

//Returns a random long between 0 and 4294967295
unsigned long rand32()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if(mti >= 624)
	{
        int kk;

        for(kk=0;kk<227;kk++)
		{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+397] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for(;kk<623;kk++)
		{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk-227] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[623]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}




