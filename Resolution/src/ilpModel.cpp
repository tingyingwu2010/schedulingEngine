#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

int Instance::ecrirePLNE(unsigned int fctObj, unsigned int T){
	vector< vector<int> > invcorresp(m); //Donne les indices des lots pour lesquels chaque machine est qualifiee
	for(int i=0;i<n;i++){
		for(unsigned int j=0;j<corresp[i].size();j++)
			invcorresp[corresp[i][j]-1].push_back(i+1);	
	}

	/*******Ouverture du fichier au format .lp*******/
	ostringstream pl;
	pl<<"./LPFiles/"<<instanceName<<".lp";
	ofstream fpl;
	cout<<"Ecriture du fichier "<<pl.str()<<endl;
	fpl.open(pl.str().c_str());		

	// Writing objective function
	switch (fctObj){
		case 1 :fpl<<"Maximize"<<endl;
			for(int i=1;i<=n;i++){
				for(unsigned int j=1;j<=corresp[i-1].size();j++){
					int borne = (H<(int)(T-exec[i-1][j-1]))?H:T-exec[i-1][j-1];
					for(int t=0;t<=borne;t++){
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
			for(int i=1;i<=l;i++){
				for(int j=0;j<=m;j++){
					fpl<<" +"<<"u"<<i<<"_"<<j;
				}
			}
			fpl<<" -"<<l;
			break;
	}		

	fpl<<endl<<"Subject To"<<endl;

	//Ecriture des contraintes
	
	if(fctObj!=3){
	
	//1 Un lot s'execute une unique fois
	cout<<"Ecriture des contraintes de type 1..."<<endl;
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			for(unsigned int t=0;t<=T-exec[i-1][j-1];t++){
				//cout<<T<<" "<<exec[i-1][j-1]<<endl;cin.get();
				fpl<<((j+t==1)?" ":" +")<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;			
			}
		}
		if(fctObj==1)
			fpl<<" = 1"<<endl;
		else
			fpl<<" = 1"<<endl;	
	}
/**/
	//2 Un seul lot en cours d'execution a la fois sur une machine et temps de setup et non simultaneite d'execution de lots 
	// ayant le meme masque requis
	cout<<"Ecriture des contraintes de type 2..."<<endl;
	int nbterm = 0;
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			int indicemachine = corresp[i-1][j-1];
			for(unsigned int t=0;t<=T-exec[i-1][j-1];t++){//Pour tout triplet (i,j,t)
				nbterm = 0;
				for(unsigned int i0=1;i0<=invcorresp[indicemachine-1].size();i0++){
					int indiceLotCourant = invcorresp[indicemachine-1][i0-1];
					int indicejPouri0 = indexMachineQualif(indiceLotCourant,indicemachine);
					if(indiceLotCourant != i && indicejPouri0!=-1){//pour tout lot != i pouvant aller sur la machine
						int borne = ((int)(T-exec[indiceLotCourant-1][indicejPouri0])<=(int)(t+exec[i-1][j-1]-1+setup[families[indicemachine-1]-1][batch[i-1][j-1]-1][batch[indiceLotCourant-1][indicejPouri0]-1]))?T-exec[indiceLotCourant-1][indicejPouri0]:t+exec[i-1][j-1]-1+setup[families[indicemachine-1]-1][batch[i-1][j-1]-1][batch[indiceLotCourant-1][indicejPouri0]-1];
						for(int t0=t;t0<=borne;t0++){
							fpl<<" + u"<<indiceLotCourant<<"_"<<indicemachine<<"_"<<t0;
							nbterm++;
						}
					}			
				}
				for(int i0=1;i0<=n;i0++){
						if((phi[i0-1]==phi[i-1])&&(i0!=i)){
							for(unsigned int j0=1;j0<=corresp[i0-1].size();j0++){
							 if(corresp[i-1][j-1]!=corresp[i0-1][j0-1]){
							  int borne= ((int)(T-exec[i0-1][j0-1])<=(int)(t+exec[i-1][j-1]))?T-exec[i0-1][j0-1]:t+exec[i-1][j-1];
							   for(int t0=t;t0<=borne;t0++){
							  	fpl<<" + u"<<i0<<"_"<<corresp[i0-1][j0-1]<<"_"<<t0;
								nbterm++;
							  }
							 }
							}					
						}				
				}
				if(nbterm!=0)
					fpl<<" +"<<nbterm<<" u"<<i<<"_"<<indicemachine<<"_"<<t<<" <= "<<nbterm<<endl;
			}
		}
	}

	/**/
		//On ne peut pas commencer a t=0 avec un job dont le masque requis n'est pas initialement sur la machine
		cout<<"Ecriture de la contrainte de type 3..."<<endl;
		for(int i=1;i<=n;i++){
			for(unsigned int j=1;j<=corresp[i-1].size();j++){
				if(initMask[phi[i-1]-1]!=corresp[i-1][j-1])
					fpl<<" + u"<<i<<"_"<<corresp[i-1][j-1]<<"_0";			
			}
		}
		fpl<<" = 0"<<endl;
	/**/
	}

	else{
		cout<<endl<<endl<<"Formulation SC_POP"<<endl<<endl;	
		cout<<"Ecriture des contraintes de Set Cover..."<<endl;
	
		vector< vector<int> > jobM(l);	//Chaque case contient la liste des jobs pour le masque Ai	
		for(int i=1;i<=n;i++){
			jobM[phi[i-1]-1].push_back(i);
		}
		
	/*	for(int i=1;i<=l;i++){cout<<"A"<<i<<" : ";
			for(int j=1;j<=jobM[i-1].size();j++)
				cout<<jobM[i-1][j-1]<<" ";
			cout<<endl;
		}
		cin.get();
		cin.get();
	*/

		//Pour tout masque A_i on va introduire jobM[i-1].size() contraintes
		for(int i=1;i<=l;i++){
			//On commence par determiner les ensembles M_j(i)
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
	}

	// Specification of binaries
	fpl<<"Binaries"<<endl;
	if(fctObj!=3){
		for(unsigned int t=0;t<=T;t++)
			for(int i=1;i<=n;i++)		
				for(unsigned int j=1;j<=corresp[i-1].size();j++)
					if(t<= T - exec[i-1][j-1])
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
	vector<vector<vector<int> > > ordo(m);
	ostringstream cheminSol;
	if(fctObj==4)
		cheminSol<<"./SolFiles/"<<instanceName<<"_ILPMO.sol";
	else
		cheminSol<<"./SolFiles/"<<instanceName<<".sol";
	ifstream fst(cheminSol.str().c_str(),ios::in);
	cout<<"Lecture du fichier "<<cheminSol.str()<<endl;
	char buff[512];	
	if(fctObj!=3){
		fst.getline(buff,512);	
		fst.getline(buff,512);
		while(fst.good()){
			char car = fst.get();
			if(car=='u'){
				int indM, indJ, tem;
				fst>>indJ;
				fst.get();
				fst>>indM;
				fst.get();
				fst>>tem;
				vector<int> elt(2);
				elt[0] = indJ;
				elt[1] = tem;
				ordo[indM-1].push_back(elt);
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
			for(unsigned int k=1;k<=jobM[i-1].size();k++){
				for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
					machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);		
			}
			
			partialSol[i-1].resize(2);
			if(initMask[i-1]!=0){
				for(unsigned int dd=0;dd<machineI[initMask[i-1]].size();dd++){
					if(machineI[initMask[i-1]][dd]!=0){
						partialSol[i-1][0].push_back(initMask[i-1]);
						partialSol[i-1][1].push_back(machineI[initMask[i-1]][dd]);
					}
			 	}
			}
			for(unsigned int d=0;d<u[i-1].size();d++){
				if(u[i-1][d]!=initMask[i-1]){
					for(unsigned int dd=0;dd<machineI[u[i-1][d]].size();dd++){
						if(member(partialSol[i-1][1],machineI[u[i-1][d]][dd])==-1){
							partialSol[i-1][0].push_back(u[i-1][d]);
							partialSol[i-1][1].push_back(machineI[u[i-1][d]][dd]);
						}
					}
				}
			}
		}
		ordo = completerOrdonnancement(partialSol);
	}
	return ordo;
}

int Instance::ecrirePLNEmasques(int infBound){
	vector<vector<vector<int> > > ordo(m);
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

	//Pour tout masque A_i on va introduire jobM[i-1].size() contraintes
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
			for(int j=0;j<=m;j++){
				if(member(machineI[j],jobM[i-1][k1])!=-1)
					fpl<<" +"<<"u"<<i<<"_"<<j;
			}
			fpl<<" >= 1"<<endl;
		}
	}
	// Additional constraint
	for(int i=1;i<=l;i++)
		for(int j=0;j<=m;j++)
			fpl<<" +"<<"u"<<i<<"_"<<j;
	fpl<<" >= "<<l+infBound<<endl;
	
	//Specification des variables binaires
	fpl<<"Binaries"<<endl;

	for(int i=1;i<=l;i++)		
		for(int j=0;j<=m;j++)
			fpl<<"u"<<i<<"_"<<j<<endl;

	fpl<<"End";
	fpl.close();
	return(EXIT_SUCCESS);
}


