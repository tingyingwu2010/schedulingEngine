#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

vector<vector<int> > Instance::setCover(int i, vector<int> JJ, vector<vector<int> > MM){
	vector<vector <int> > cover(2);
	int s = JJ.size();
	vector<int> indices;
	vector< vector<int> > sol;
	unsigned int index = 0;
	while(s > 0){
		int max = initMask[i-1];
		for(int k=1;k<=m;k++){
			if(MM[max].size()<MM[k].size())
				max = k;
		}
		indices.push_back(max);		
		sol.resize(sol.size()+1);
		sol[sol.size()-1].resize(MM[max].size());
		sol[sol.size()-1] = MM[max];
		if(max == initMask[i-1])
			index = sol.size()-1;		
		s-=MM[max].size();
		for(unsigned int k=0;k<MM[max].size();k++){
			for(int j=0;j<=m;j++){
				if(j != max){
					unsigned int jj = 0;
					while(jj<MM[j].size()){
						if(MM[j][jj] == MM[max][k]){
							MM[j].erase(MM[j].begin()+jj);
							jj--;
						}
						jj++;
					}
				}
			}
		}
		MM[max].clear();
	}
	// We build the solution
	for(unsigned int dd=0;dd<sol[index].size();dd++){
		if(sol[index][dd] != 0){
			cover[0].push_back(initMask[i-1]);
			cover[1].push_back(sol[index][dd]);
		}
	}
	for(unsigned int d=0;d<sol.size();d++){
		if(d != index){
			for(unsigned int dd=0;dd<sol[d].size();dd++){
				if(sol[d][dd] != 0){
					cover[0].push_back(indices[d]);
					cover[1].push_back(sol[d][dd]);
				}
			}
		}
	}
	return cover;
}

vector<vector<vector<int> > > Instance::completeSchedule(vector<vector<vector<int> > > partialSolution){
	vector<vector<vector<int> > > finalSchedule(m);
	vector<int> S(l); // Current machine per mask
	vector<int> D(m,0); // Last job per machind
	for(int i=1;i<=l;i++)
		S[i-1] = initMask[i-1];
	for(int i=1;i<=l;i++){
		for(unsigned int j=1;j<=partialSolution[i-1][0].size();j++){
			int currentJob = partialSolution[i-1][1][j-1];
			int currentMachine = partialSolution[i-1][0][j-1];
			int tau = 0;
			if(S[i-1] != currentMachine)
				tau = 1;
			S[i-1] = currentMachine;
			int tyijm1 = (j == 1)?0:finalSchedule[partialSolution[i-1][0][j-2]-1][finalSchedule[partialSolution[i-1][0][j-2]-1].size()-1][1];
			int rhom1 = (j == 1)?0:exec[partialSolution[i-1][1][j-2]-1][indexMachineQualif(partialSolution[i-1][1][j-2],partialSolution[i-1][0][j-2])];
			int tDx = (D[currentMachine-1] == 0)?0:finalSchedule[currentMachine-1][finalSchedule[currentMachine-1].size()-1][1];
			int imQ;
			if(D[currentMachine-1]!=0)
				imQ = indexMachineQualif(D[currentMachine-1],currentMachine);
			int rhoDx = (D[currentMachine-1] == 0)?0:exec[D[currentMachine-1]-1][imQ];
			int Beta = (D[currentMachine-1] == 0)?0:setup[families[currentMachine-1]-1][batch[D[currentMachine-1]-1][imQ]-1][batch[currentJob-1][indexMachineQualif(currentJob,currentMachine)]-1];
			// Formula
			vector< int > elt(2);
			elt[0] = currentJob;
			elt[1] = (tyijm1+rhom1+tau > tDx+rhoDx+Beta)?tyijm1 + rhom1 + tau:tDx + rhoDx + Beta;
			finalSchedule[currentMachine-1].push_back(elt);
			D[currentMachine-1] = currentJob;
		}
	}
	return finalSchedule;
}

double Instance::harmonicSerie(int a){
	double alpha = 0.0;
	for(int i=1;i<=a;i++)
		alpha+=(double)(1.0/i);
	return alpha;
}

vector<vector<vector<int> > > Instance::numberMaskMoves(){
	vector<vector<vector<int> > > finalSchedule;
	vector<vector<vector<int> > > partialSolution(l);
	// Insert the jobs requiring each mask in an array
	vector< vector<int> > jobM(l);
	for(int i=1;i<=n;i++)
		jobM[phi[i-1]-1].push_back(i);
	unsigned int Delta = 0;
	double lowerBound = 0.0;
	double upperBound = 0.0;
	// For each mask A_i
	for(int i=1;i<=l;i++){
		// Start with the sets M_j(i)
		vector< vector<int> > machineI(m+1);
		for(unsigned int k=1;k<=jobM[i-1].size();k++)
			for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
				machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);		
		// Check if M_R_i(i) has greatest cardinality
		bool MriGrand = true;
		int cp = 0;		
		while(cp <= m && MriGrand){
			if(cp != initMask[i-1] && machineI[cp].size() > machineI[initMask[i-1]].size())
				MriGrand = false;
			cp++;
		}
		if(!MriGrand){
			jobM[i-1].push_back(0);
			machineI[initMask[i-1]].push_back(0);	
		}		
		unsigned int Deltai = 0;
		double Bi = 1.0;
		for(unsigned int ii=0;ii<machineI.size();ii++)
			if(Deltai < machineI[ii].size())
				Deltai = machineI[ii].size();
		if(Delta<Deltai)
			Delta = Deltai;
		// Min cover
		partialSolution[i-1] = setCover(i,jobM[i-1],machineI);
		for(unsigned int ii=0;ii<partialSolution[i-1][0].size()-1;ii++)
			if(partialSolution[i-1][0][ii] != partialSolution[i-1][0][ii+1])
				Bi++;
		if(initMask[i-1] != partialSolution[i-1][0][0])
			Bi++;
		upperBound+=Bi-1.0;
		lowerBound+=Bi/harmonicSerie(Deltai) - 1.0;
	}
	// We have the subproblems solutions. Now complete final schedule
	finalSchedule = completeSchedule(partialSolution);
	cout<<endl<<endl<<"Approximation ratio is "<<1/harmonicSerie(Delta)<<endl<<endl;
	cout<<"We have the inequalities: "<<ceil(lowerBound)<<" <= OPT <= "<<upperBound<<endl<<endl;
	return finalSchedule;
}


