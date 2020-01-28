#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

int Instance::writeILP(unsigned int fctObj, unsigned int T){
	vector< vector<int> > invCorresp(m); // Gives indices of jobs such that machines are qualified
	for(int i=0;i<n;i++){
		for(unsigned int j=0;j<corresp[i].size();j++)
			invCorresp[corresp[i][j]-1].push_back(i+1);	
	}

	/*******Open file to .lp format*******/
	ostringstream pl;
	pl<<"./LPFiles/"<<instanceName<<".lp";
	ofstream fpl;
	cout<<"Writing files "<<pl.str()<<endl;
	fpl.open(pl.str().c_str());		

	// Writing objective function
	switch(fctObj){
		case 1 :fpl<<"Maximize"<<endl;
			for(int i=1;i<=n;i++){
				for(unsigned int j=1;j<=corresp[i-1].size();j++){
					int bound = (H<(int)(T-exec[i-1][j-1]))?H:T-exec[i-1][j-1];
					for(int t=0;t<=bound;t++){
						fpl<<(((t==0)&&(i*j==1))?" ":" +")<<((t<=H-exec[i-1][j-1])?(w[i-1]):(w[i-1]*(H-t)/exec[i-1][j-1]))<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;					
					}
				}		
			}
			break;
		case 2 :fpl<<"Minimize"<<endl;
			for(int i=1;i<=n;i++){
				for(unsigned int j=1;j<=corresp[i-1].size();j++){
					for(unsigned int t=0;t<=T-exec[i-1][j-1];t++){
						fpl<<(((i*j==1)&&(t==0))?" ":" +")<<c[i-1]*(t+exec[i-1][j-1])<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;	
					}
				}
			}
			break;
		default :fpl<<"Minimize"<<endl;
			for(int i=1;i<=l;i++)
				for(int j=0;j<=m;j++)
					fpl<<" +"<<"u"<<i<<"_"<<j;
			fpl<<" -"<<l;
			break;
	}
	fpl<<endl<<"Subject To"<<endl;

	// Writing constraints
	if(fctObj != 3){
		// 1: A job runs only once
		cout<<"Writing assignment constraints..."<<endl;
		for(int i=1;i<=n;i++){
			for(unsigned int j=1;j<=corresp[i-1].size();j++)
				for(unsigned int t=0;t<=T-exec[i-1][j-1];t++)
					fpl<<((j+t==1)?" ":" +")<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;
			if(fctObj==1)
				fpl<<" = 1"<<endl;
			else
				fpl<<" = 1"<<endl;
		}
		//2: Machine capacity, sequence-dependent setup times and mask transport constraints
		cout<<"Writing machine capacity constraints..."<<endl;
		int nbTerm = 0;
		for(int i=1;i<=n;i++){
			for(unsigned int j=1;j<=corresp[i-1].size();j++){
				int machineIndex = corresp[i-1][j-1];
				for(unsigned int t=0;t<=T-exec[i-1][j-1];t++){// For each (i,j,t)
					nbTerm = 0;
					for(unsigned int i0=1;i0<=invCorresp[machineIndex-1].size();i0++){
						int currentJobIndex = invCorresp[machineIndex-1][i0-1];
						int indexJForI0 = indexMachineQualif(currentJobIndex,machineIndex);
						if(currentJobIndex != i && indexJForI0!=-1){ // For each job != i going on the machine
							int bound = ((int)(T-exec[currentJobIndex-1][indexJForI0])<=(int)(t+exec[i-1][j-1]-1+setup[families[machineIndex-1]-1][batch[i-1][j-1]-1][batch[currentJobIndex-1][indexJForI0]-1]))?T-exec[currentJobIndex-1][indexJForI0]:t+exec[i-1][j-1]-1+setup[families[machineIndex-1]-1][batch[i-1][j-1]-1][batch[currentJobIndex-1][indexJForI0]-1];
							for(int t0=t;t0<=bound;t0++){
								fpl<<" + u"<<currentJobIndex<<"_"<<machineIndex<<"_"<<t0;
								nbTerm++;
							}
						}
					}
					for(int i0=1;i0<=n;i0++){
						if((phi[i0-1] == phi[i-1])&&(i0!=i)){
							for(unsigned int j0=1;j0<=corresp[i0-1].size();j0++){
								if(corresp[i-1][j-1]!=corresp[i0-1][j0-1]){
									int bound = ((int)(T-exec[i0-1][j0-1])<=(int)(t+exec[i-1][j-1]))?T-exec[i0-1][j0-1]:t+exec[i-1][j-1];
									for(int t0=t;t0<=bound;t0++){
							  			fpl<<" + u"<<i0<<"_"<<corresp[i0-1][j0-1]<<"_"<<t0;
										nbTerm++;
									}
								}
							}
						}
					}
					if(nbTerm != 0)
						fpl<<" +"<<nbTerm<<" u"<<i<<"_"<<machineIndex<<"_"<<t<<" <= "<<nbTerm<<endl;
				}
			}
		}
		// Impossible to start at t=0 when the mask is not initially on the machine
		cout<<"Writing initial mask transport constraint..."<<endl;
		for(int i=1;i<=n;i++)
			for(unsigned int j=1;j<=corresp[i-1].size();j++)
				if(initMask[phi[i-1]-1]!=corresp[i-1][j-1])
					fpl<<" + u"<<i<<"_"<<corresp[i-1][j-1]<<"_0";
		fpl<<" = 0"<<endl;
	}
	else{
		cout<<endl<<endl<<"Formulation SC_POP"<<endl<<endl;	
		cout<<"Writing Set Cover constraints..."<<endl;
		vector< vector<int> > jobM(l);	// Jobs lists per mask
		for(int i=1;i<=n;i++)
			jobM[phi[i-1]-1].push_back(i);
		// For each mask A_i, adding the jobM[i-1].size() constraints
		for(int i=1;i<=l;i++){
			// Determining the M_j(i) sets
			vector< vector<int> > machineI(m+1);
			for(unsigned int k=1;k<=jobM[i-1].size();k++)
				for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
					machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);
			machineI[initMask[i-1]].push_back(0);
			jobM[i-1].push_back(0);
			for(unsigned int k1=0;k1<jobM[i-1].size();k1++){
				for(int j=0;j<=m;j++)
					if(member(machineI[j],jobM[i-1][k1])!=-1)
						fpl<<" +"<<"u"<<i<<"_"<<j;
				fpl<<" >= 1"<<endl;
			}
		}
	}

	// Specification of binaries
	fpl<<"Binaries"<<endl;
	if(fctObj!=3){
		for(unsigned int t=0;t<=T;t++)
			for(int i=1;i<=n;i++)		
				for(unsigned int j=1;j<=corresp[i-1].size();j++)
					if(t <= T - exec[i-1][j-1])
						fpl<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t<<endl;
	}
	else{
		for(int i=1;i<=l;i++)		
			for(int j=0;j<=m;j++)
					fpl<<"u"<<i<<"_"<<j<<endl;
	}
	fpl<<"End";
	fpl.close();
	return(EXIT_SUCCESS);
}

vector<vector<vector<int> > > Instance::makeScheduleFromILP(int fctObj){
	vector<vector<vector<int> > > schedule(m);
	ostringstream solutionPath;
	if(fctObj == 4)
		solutionPath<<"./SolFiles/"<<instanceName<<"_ILPMO.sol";
	else
		solutionPath<<"./SolFiles/"<<instanceName<<".sol";
	ifstream fst(solutionPath.str().c_str(),ios::in);
	cout<<"Reading file "<<solutionPath.str()<<endl;
	char buff[512];	
	if(fctObj != 3){
		fst.getline(buff,512);	
		fst.getline(buff,512);
		while(fst.good()){
			char car = fst.get();
			if(car == 'u'){
				int indM, indJ, tem;
				fst>>indJ;
				fst.get();
				fst>>indM;
				fst.get();
				fst>>tem;
				vector<int> elt(2);
				elt[0] = indJ;
				elt[1] = tem;
				schedule[indM-1].push_back(elt);
				fst.getline(buff,512);		
			}
		}
	}
	else{
		fst.getline(buff,512);	
		fst.getline(buff,512);
		// Form partial solution
		vector<vector<vector<int> > > partialSol(l);
		vector<vector<int> > u(l);
		int indl,indM;
		while(fst.good()){
			fst.get();
			fst>>indl;
			fst.get();
			fst>>indM;
			fst.getline(buff,512);
			u[indl-1].push_back(indM);
		}
		// To fix a f.good()-related bug
		u[indl-1].pop_back();
		vector< vector<int> > jobM(l); // Each item contains the list of jobs for mask i	
		for(int i=1;i<=n;i++)
			jobM[phi[i-1]-1].push_back(i);
		for(int i=1;i<=l;i++){
			// Start by determining sets M_j(i)
			vector< vector<int> > machineI(m+1);
			for(unsigned int k=1;k<=jobM[i-1].size();k++)
				for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
					machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);		
			partialSol[i-1].resize(2);
			if(initMask[i-1] != 0){
				for(unsigned int dd=0;dd<machineI[initMask[i-1]].size();dd++){
					if(machineI[initMask[i-1]][dd] != 0){
						partialSol[i-1][0].push_back(initMask[i-1]);
						partialSol[i-1][1].push_back(machineI[initMask[i-1]][dd]);
					}
			 	}
			}
			for(unsigned int d=0;d<u[i-1].size();d++){
				if(u[i-1][d] != initMask[i-1]){
					for(unsigned int dd=0;dd<machineI[u[i-1][d]].size();dd++){
						if(member(partialSol[i-1][1],machineI[u[i-1][d]][dd])==-1){
							partialSol[i-1][0].push_back(u[i-1][d]);
							partialSol[i-1][1].push_back(machineI[u[i-1][d]][dd]);
						}
					}
				}
			}
		}
		schedule = completeSchedule(partialSol);
	}
	return schedule;
}

int Instance::writeILPmasks(int infBound){
	vector<vector<vector<int> > > schedule(m);
	// File opening
	ostringstream pl;
	pl<<"./LPFiles/"<<instanceName<<"_bound_"<<infBound<<".lp";
	ofstream fpl;
	cout<<"Writing file "<<pl.str()<<endl;
	fpl.open(pl.str().c_str());
	// Objective function
	fpl<<"Minimize"<<endl;
	for(int i=1;i<=l;i++)
		for(int j=0;j<=m;j++)
			fpl<<" +"<<"u"<<i<<"_"<<j;
	fpl<<" -"<<l;

	fpl<<endl<<"Subject To"<<endl;
	
	// Constraints
	cout<<endl<<endl<<"Formulation SC_POP"<<endl<<endl;	
	cout<<"Writing Set Cover constraints..."<<endl;
	vector< vector<int> > jobM(l);	// Each item contains the jobs for mask i	
	for(int i=1;i<=n;i++)
		jobM[phi[i-1]-1].push_back(i);

	// For each mask A_i introduce jobM[i-1].size() constraints
	for(int i=1;i<=l;i++){
		// Start by determining sets M_j(i)
		vector< vector<int> > machineI(m+1);
		for(unsigned int k=1;k<=jobM[i-1].size();k++){
			for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
				machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);		
		}
		machineI[initMask[i-1]].push_back(0);
		jobM[i-1].push_back(0);
		for(unsigned int k1=0;k1<jobM[i-1].size();k1++){
			for(int j=0;j<=m;j++)
				if(member(machineI[j],jobM[i-1][k1])!=-1)
					fpl<<" +"<<"u"<<i<<"_"<<j;
			fpl<<" >= 1"<<endl;
		}
	}
	// Additional constraint
	for(int i=1;i<=l;i++)
		for(int j=0;j<=m;j++)
			fpl<<" +"<<"u"<<i<<"_"<<j;
	fpl<<" >= "<<l+infBound<<endl;

	// Specification of binaries
	fpl<<"Binaries"<<endl;

	for(int i=1;i<=l;i++)
		for(int j=0;j<=m;j++)
			fpl<<"u"<<i<<"_"<<j<<endl;

	fpl<<"End";
	fpl.close();
	return(EXIT_SUCCESS);
}


