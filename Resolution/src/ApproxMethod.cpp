#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

vector<vector<int> > Instance::couvertureMin(int i, vector<int> JJ, vector<vector<int> > MM){
	vector<vector <int> > couverture(2);
	int s = JJ.size();
	vector<int> indices;
	vector< vector<int> > sol;
	unsigned int indice = 0;
	while(s > 0){
		//cout<<"Iteration courante : nombre d'elements non couverts egal a "<<s<<endl;
		
		int max=initMask[i-1];
		for(int k=1;k<=m;k++){
			if(MM[max].size()<MM[k].size())
				max=k;
		}
		indices.push_back(max);		
		/**
		cout<<"Ensemble courant trouve pour couvrir : "<<endl<<"{ ";
		for(int kk=0;kk<MM[max].size()-1;kk++)
			cout<<"J"<<MM[max][kk]<<" , ";
		cout<<"J"<<MM[max][MM[max].size()-1]<<" }"<<endl<<endl;
		cout<<"Indice de la machine en question : "<<max<<endl;
		**/
		sol.resize(sol.size()+1);
		sol[sol.size()-1].resize(MM[max].size());
		sol[sol.size()-1] = MM[max];
		if(max==initMask[i-1])
			indice = sol.size()-1;		
		s-=MM[max].size();
		for(unsigned int k=0;k<MM[max].size();k++){
			for(int j=0;j<=m;j++){
				if(j!=max){
					unsigned int jj=0;
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
	/**
	cout<<"Affichage de la sol : "<<endl;
	for(int d=0;d<sol.size();d++){
	cout<<"{ ";
	for(int kk=0;kk<sol[d].size()-1;kk++)
			cout<<"J"<<sol[d][kk]<<" , ";
	cout<<"J"<<sol[d][sol[d].size()-1]<<" }"<<endl<<endl;		
	}	
	cout<<"L'indice de la machine de depart est "<<indice+1<<endl;
	cin.get();
	**/
	
	//On construit couverture a partir de indices et sol
	for(unsigned int dd=0;dd<sol[indice].size();dd++){
		if(sol[indice][dd]!=0){
			couverture[0].push_back(initMask[i-1]);
			couverture[1].push_back(sol[indice][dd]);
		}
	}
	for(unsigned int d=0;d<sol.size();d++){
		if(d!=indice){
			for(unsigned int dd=0;dd<sol[d].size();dd++){
				if(sol[d][dd]!=0){
					couverture[0].push_back(indices[d]);
					couverture[1].push_back(sol[d][dd]);
				}
			}
		}
	}
	//if(i==2){	
	/**
		cout<<endl;
		for(int ii=0;ii<couverture[0].size();ii++)
			cout<<couverture[0][ii]<<" ";
		cout<<endl;
		for(int ii=0;ii<couverture[0].size();ii++)
			cout<<couverture[1][ii]<<" ";
		cout<<endl;
		cin.get();
	**/
	//}
	return couverture;
}

vector<vector<vector<int> > > Instance::completerOrdonnancement(vector<vector<vector<int> > > solpartielle){
	vector<vector<vector<int> > > ordoFinal(m);
	vector<int> S(l); //Machine courante de chaque masque
	vector<int> D(m,0); //Dernier job sur chaque machine
	/**
	for(int i=0;i<l;i++){
	cout<<"Masque "<<i+1<<" : ";
		for(int j=0;j<solpartielle[i].size();j++)
			cout<<"(M"<<solpartielle[i][0][j]<<" , J"<<solpartielle[i][1][j]<<" )";
	cout<<endl;
	}
	cin.get();cin.get();
	}**/
	for(int i=1;i<=l;i++)
		S[i-1]=initMask[i-1];
	for(int i=1;i<=l;i++){
		for(unsigned int j=1;j<=solpartielle[i-1][0].size();j++){
			int jobact = solpartielle[i-1][1][j-1];
			int machact = solpartielle[i-1][0][j-1];
			//cout<<"Masque A"<<i<<", Job a la position "<<j<<" : job J"<<jobact<<" pour la machine M"<<machact<<" : "<<endl;
			int tau = 0;
			if(S[i-1]!=machact)
				tau=1;
			S[i-1] = machact;
			//cout<<"Le masque A"<<i<<" est maintenant sur M"<<S[i-1]<<endl;
			int tyijm1 = (j==1)?0:ordoFinal[solpartielle[i-1][0][j-2]-1][ordoFinal[solpartielle[i-1][0][j-2]-1].size()-1][1];
			//cout<<"La date de debut d'execution du job precedent pour le masque A"<<i<<" est "<<tyijm1<<endl;
			int rhom1 = (j==1)?0:exec[solpartielle[i-1][1][j-2]-1][indexMachineQualif(solpartielle[i-1][1][j-2],solpartielle[i-1][0][j-2])];
			//cout<<"Le temps d'execution de ce job sur la machine concernee est de"<<rhom1<<endl;
			int tDx = (D[machact-1]==0)?0:ordoFinal[machact-1][ordoFinal[machact-1].size()-1][1];
			//if(tDx!=0)
			//cout<<"La date de debut d'execution du job precedent, jobJ"<<ordoFinal[machact-1][ordoFinal[machact-1].size()-1][0]<<", sur la machine M"<<machact<<" est "<<tDx<<endl;
			int imQ;
			if(D[machact-1]!=0)
				imQ = indexMachineQualif(D[machact-1],machact);
			int rhoDx = (D[machact-1]==0)?0:exec[D[machact-1]-1][imQ];
			//cout<<"Le temps d'execution de ce job sur cette machine est de"<<rhoDx<<endl;
			int Beta = (D[machact-1]==0)?0:setup[families[machact-1]-1][batch[D[machact-1]-1][imQ]-1][batch[jobact-1][indexMachineQualif(jobact,machact)]-1];
			//cout<<"Le temps de setup pour passer du job precedent sur la machine M"<<machact<<" a J"<<jobact<<" est de "<<Beta<<endl;
			//formule
			vector< int > elt(2);
			elt[0] = jobact;
			elt[1] = (tyijm1+rhom1+tau > tDx+rhoDx+Beta)?tyijm1+rhom1+tau:tDx+rhoDx+Beta;
			ordoFinal[machact-1].push_back(elt);				
			D[machact-1] = jobact;
			//cout<<"La date de debut d'execution du job J"<<jobact<<" sur la machine M"<<machact<<" est donc "<<elt[1]<<endl;
		}
	}
	return ordoFinal;
}

double Instance::nombreHarmonique(int a){
	double alpha=0.0;
	for(int i=1;i<=a;i++)
		alpha+=(double)(1.0/i);
	return alpha;
}

vector<vector<vector<int> > > Instance::nombreDeplacementsMasques(){
	vector<vector<vector<int> > > ordoFinal;
	vector<vector<vector<int> > > solpartielle(l);
	//On range dans un tableau les jobs necessitant chaque masque
	vector< vector<int> > jobM(l);
	
	for(int i=1;i<=n;i++){
		jobM[phi[i-1]-1].push_back(i);
	}

	/**
	for(int i=1;i<=l;i++){
		cout<<"Masque "<<i<<" : ";
		for(int j=1;j<=jobM[i-1].size();j++){
			cout<<jobM[i-1][j-1]<<" ";
		}
		cout<<endl;
	}
	vector< vector<int> > machineJ(m);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=corresp[i-1].size();j++){
			machineJ[corresp[i-1][j-1]-1].push_back(i);
		}
	}
	
	for(int i=1;i<=m;i++){
		cout<<"Machine "<<i<<" : ";
		for(int j=1;j<=machineJ[i-1].size();j++){
			cout<<machineJ[i-1][j-1]<<" ";
		}
		cout<<endl;
	cin.get();cin.get();
	}**/
	unsigned int Delta = 0;
	double borneinf = 0.0;
	double bornesup = 0.0;
	//Pour tout masque A_i
	for(int i=1;i<=l;i++){
		//On commence par determiner les ensembles M_j(i)
		vector< vector<int> > machineI(m+1);
		for(unsigned int k=1;k<=jobM[i-1].size();k++){
			for(unsigned int d=1;d<=corresp[jobM[i-1][k-1]-1].size();d++)
				machineI[corresp[jobM[i-1][k-1]-1][d-1]].push_back(jobM[i-1][k-1]);		
		}
		
		/**
		
		if(i==5){
		 for(int ii=0;ii<=m;ii++){cout<<"Machine "<<ii<<" : ";
		  for(int j=1;j<=machineI[ii].size();j++){
			cout<<machineI[ii][j-1]<<" ";
		}cout<<endl;}
		
		for(int ii=1;ii<=m;ii++){
			for(int j=1;j<=machineI[ii].size();j++){
				cout<<machineI[ii][j-1]<<" ";
			}
			cout<<endl;
		}
		}
		cin.get();
		**/

		//On verifie si M_R_i(i) a la plus grande cardinalite
		int MriGrand=true;
		int cp=0;		
		while(cp<=m && MriGrand==true){
			if(cp!=initMask[i-1] && machineI[cp].size() > machineI[initMask[i-1]].size())
				MriGrand=false;
			cp++;
		}
		if(MriGrand==false){
			jobM[i-1].push_back(0);
			machineI[initMask[i-1]].push_back(0);	
		}
/**
		if(i==5){
		for(int ii=1;ii<=jobM[i-1].size();ii++){
			cout<<jobM[i-1][ii-1]<<" ";
		}
		cout<<endl;
		for(int ii=0;ii<=m;ii++){
			for(int j=1;j<=machineI[ii].size();j++){
				cout<<machineI[ii][j-1]<<" ";
			}
			cout<<endl;
		}
		cin.get();
		}
**/
		
		unsigned int Deltai = 0;
		double Bi = 1.0;
	
		for(unsigned int ii=0;ii<machineI.size();ii++)
			if(Deltai<machineI[ii].size())
				Deltai = machineI[ii].size();
		if(Delta<Deltai)
			Delta = Deltai;
		//Couverture Min
		solpartielle[i-1]=couvertureMin(i,jobM[i-1],machineI);
		for(unsigned int ii=0;ii<solpartielle[i-1][0].size()-1;ii++)
			if(solpartielle[i-1][0][ii] != solpartielle[i-1][0][ii+1])
				Bi++;
		if(initMask[i-1] != solpartielle[i-1][0][0])
			Bi++;
		//cout<<"Nombre de deplacements du masque A"<<i<<" : "<<Bi-1<<endl;
		bornesup+=Bi-1.0;
		//cout<<"TOTAL : "<<bornesup<<endl;
		borneinf+=Bi/nombreHarmonique(Deltai) - 1.0;
	}
	//On a les solutions de tous les sous-problemes; il faut maintenant former l'ordo final.
	//afficheTemps();

	ordoFinal = completerOrdonnancement(solpartielle);
	cout<<endl<<endl<<"Le rapport d'approximation de l'algorithme sur ce probleme est d'environ "<<1/nombreHarmonique(Delta)<<endl<<endl;
	cout<<"On a l'encadrement de l'optimum suivant : "<<ceil(borneinf)<<" <= OPT <= "<<bornesup<<endl<<endl;
	return ordoFinal;
}



