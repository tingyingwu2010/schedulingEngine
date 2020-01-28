#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

double Instance::minorantCriterion2(){
	double nb = 0.0;
	for(int i=0;i<n;i++){
		int min = 0;
		for(unsigned int j=0;j<corresp[i].size();j++)
			if(exec[i][j]<exec[i][min])
				min = j;
		nb+=(c[i]*exec[i][min]);
	}
	return nb;
}

double Instance::computeAntiIdeal3D(int obj){
	if(obj == 1){
		double optIn2ForCriterion1(computeAntiIdeal(1,2));
		double optIn3ForCriterion1(computeAntiIdeal(1,3));
		return min(optIn2ForCriterion1,optIn3ForCriterion1);
	}
	if(obj == 2){
		double optIn1ForCriterion2(computeAntiIdeal(2,1));
		double optIn3ForCriterion2(computeAntiIdeal(2,3));
		return max(optIn1ForCriterion2,optIn3ForCriterion2);
	}
	else{
		double optIn1ForCriterion3(computeAntiIdeal(3,1));
		double optIn2ForCriterion3(computeAntiIdeal(3,2));
		return max(optIn1ForCriterion3,optIn2ForCriterion3);
	}
}

double Instance::computeAntiIdeal(int obj, int toOpen){
	ostringstream path;
	path.str("");
	double nadir(0);
	path<<"../Instances/"<<instanceName<<"/Best_Known_Sol_"<<toOpen<<".txt";
	ifstream init(path.str().c_str(),ios::in);
	char buff[512];
	vector<vector<vector<int> > > schedule(m);
	if(init){
		init.getline(buff,512);
		init.getline(buff,512);
		for(int j=1;j<=m;j++){
			int lot(0),temps(0);
			init.getline(buff,512,':');
			init.get();
			init.getline(buff,512,' ');
			while(buff[0] == 'T'){ // There is a job for the current machine
				init>>lot;
				init.getline(buff,512,'s');
				init.get();
				init>>temps;
				vector<int> couple(2);
				couple[0] = lot;
				couple[1] = temps;
				schedule[j-1].push_back(couple);
				init.getline(buff,512,'|');
				init.get();
				init.getline(buff,512,' ');
			}
		}
		nadir = evaluateSchedule(schedule,obj,0,0,0,0,0,0,0);
		init.close();
	}
	else{
		if(obj==2)
			nadir = majorantCriterion2();
	}
	return nadir;
}

double Instance::majorantCriterion2(){ 
	double nb = 0.0;
	double ci = 0.0;
	double setupi = 0.0;
	for(int i=0;i<n;i++){
		int max = 0;
		for(unsigned int j=0;j<corresp[i].size();j++){
			if(exec[i][j]>exec[i][max])
				max = j;		
		}
		setupi = setup[families[max]-1][batch[i][max]-1][0];
		for(unsigned int k=1;k<setup[families[max]-1][batch[i][max]-1].size();k++)
			if(setup[families[max]-1][batch[i][max]-1][k]>setupi)
				setupi = setup[families[max]-1][batch[i][max]-1][k];
		if(i == n-1)
			setupi = 0;
		if(i == 0)
			ci++;
		ci+=exec[i][max];
		nb+=(c[i]*ci);
		ci+=setupi+1;
	}
	return nb;
}

int Instance::writeMIP_BCT(double l1, double l2, double r1, double r2, unsigned int T){	
	vector< vector<int> > invCorresp(m); // Gives indices of lots for which the machine is qualified
	for(int i=0;i<n;i++)
		for(unsigned int j=0;j<corresp[i].size();j++)
			invCorresp[corresp[i][j]-1].push_back(i+1);

	/*******Opening file to .lp format*******/
	ostringstream pl;
	pl<<"./LPFiles/"<<instanceName<<"_ILPMO.lp";
	ofstream fpl;
	cout<<"Writing file "<<pl.str()<<endl;
	fpl.open(pl.str().c_str());		

	// Writing objective function

	fpl<<"Minimize z"<<endl;
		
	fpl<<endl<<"Subject To"<<endl;

	// Writing constraints
	// Max linearizing
	fpl<<"z";
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			int bound = (H<(int)(T-exec[i-1][j-1]))?H:T-exec[i-1][j-1];
			for(int t=0;t<=bound;t++){
				fpl<<" + "<<((t<=H-exec[i-1][j-1])?(l1*w[i-1]):(l1*w[i-1]*(H-t)/exec[i-1][j-1]))<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;					
			}				
		}		
	}
	fpl<<" >= "<<l1*r1<<endl;

	fpl<<"z";
	for(int i=1;i<=n;i++)
		for(unsigned int j=1;j<=corresp[i-1].size();j++)
			for(int t=0;t<=(int)(T-exec[i-1][j-1]);t++)
				fpl<<" - "<<l2*c[i-1]*(t+exec[i-1][j-1])<<" u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;
	fpl<<" >= "<<-l2*r2<<endl;

	// Assignment constraints
	cout<<"Writing constraints of type 1..."<<endl;

	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++)
			for(int t=0;t<=(int)(T-exec[i-1][j-1]);t++)
				fpl<<((j+t==1)?" ":" +")<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;
		fpl<<" = 1"<<endl;	
	}

	// Capacity constraints, sequence-dependent setup times and mask transports
	cout<<"Writing constraints of type 2..."<<endl;
	int nbTerm = 0;
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			int machineIndex = corresp[i-1][j-1];
			for(unsigned int t=0;t<=T-exec[i-1][j-1];t++){ // For each (i,j,t)
				nbTerm = 0;
				for(unsigned int i0=1;i0<=invCorresp[machineIndex-1].size();i0++){
					int currentJobIndex = invCorresp[machineIndex-1][i0-1];
					int indexJForI0 = indexMachineQualif(currentJobIndex,machineIndex);
					if(currentJobIndex != i && indexJForI0!=-1){ // For each job != i going on machine
						int bound = (T-(exec[currentJobIndex-1][indexJForI0])<=(t+exec[i-1][j-1]-1+setup[families[machineIndex-1]-1][batch[i-1][j-1]-1][batch[currentJobIndex-1][indexJForI0]-1]))?T-(exec[currentJobIndex-1][indexJForI0]):(t+exec[i-1][j-1]-1+setup[families[machineIndex-1]-1][batch[i-1][j-1]-1][batch[currentJobIndex-1][indexJForI0]-1]);
						for(int t0=t;t0<=bound;t0++){
							fpl<<" + u"<<currentJobIndex<<"_"<<machineIndex<<"_"<<t0;
							nbTerm++;
						}
					}			
				}
				for(int i0=1;i0<=n;i0++){
					if((phi[i0-1] == phi[i-1])&&(i0!=i)){
						for(unsigned int j0=1;j0<=corresp[i0-1].size();j0++){
							if(corresp[i-1][j-1] != corresp[i0-1][j0-1]){
								int bound = (T-(int)(exec[i0-1][j0-1])<=(t+(int)(exec[i-1][j-1])))?T-(int)(exec[i0-1][j0-1]):(t+(int)(exec[i-1][j-1]));
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
	// Initial state constraint
	cout<<"Writing constraint of type 3..."<<endl;
	for(int i=1;i<=n;i++)
		for(unsigned int j=1;j<=corresp[i-1].size();j++)
			if(initMask[phi[i-1]-1] != corresp[i-1][j-1])
				fpl<<" + u"<<i<<"_"<<corresp[i-1][j-1]<<"_0";			
	fpl<<" = 0"<<endl;
	fpl<<"Binaries"<<endl;
	for(unsigned int t=0;t<=T;t++)
		for(int i=1;i<=n;i++)		
			for(unsigned int j=1;j<=corresp[i-1].size();j++)
				if(t <= T - exec[i-1][j-1])
					fpl<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t<<endl;
	fpl<<"End";
	fpl.close();
	return(EXIT_SUCCESS);
}


