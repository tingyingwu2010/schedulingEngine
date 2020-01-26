#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>

using namespace std;

double Instance::minorantCritere2(){
	double nb = 0.0;
	for(int i=0;i<n;i++){
		int min=0;
		for(unsigned int j=0;j<corresp[i].size();j++){
			if(exec[i][j]<exec[i][min])
				min=j;		
		}
		nb=nb+(c[i]*exec[i][min]);
	}
	return nb;
}

double Instance::calculAntiIdeal3D(int obj){
	if(obj==1){
		double optEn2SelonCritere1(calculAntiIdeal(1,2));
		double optEn3SelonCritere1(calculAntiIdeal(1,3));
		return min(optEn2SelonCritere1,optEn3SelonCritere1);
	}
	if(obj==2){
		double optEn1SelonCritere2(calculAntiIdeal(2,1));
		double optEn3SelonCritere2(calculAntiIdeal(2,3));
		return max(optEn1SelonCritere2,optEn3SelonCritere2);
	}
	else{
		double optEn1SelonCritere3(calculAntiIdeal(3,1));
		double optEn2SelonCritere3(calculAntiIdeal(3,2));
		return max(optEn1SelonCritere3,optEn2SelonCritere3);
	}
}

double Instance::calculAntiIdeal(int obj, int toOpen){
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
			nadir=majorantCritere2();
	}
	return nadir;
}

//On calcule la solution comme si on mettait tout sur la meme machine avec des temps d'execution maximaux et le temps de setup max entre
//chaque couple de taches
double Instance::majorantCritere2(){ 
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
		if(i==n-1)
			setupi=0;
		if(i==0)
			ci++;
		ci += exec[i][max];
		nb = nb + c[i]*ci;
		ci+=setupi+1;
	}
	
	return nb;
}

int Instance::ecrireMIP_BCT(double l1, double l2, double r1, double r2, unsigned int T){	
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

	//Ecriture des contraintes
	//Contraintes de max
	fpl<<"z";
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			int borne = (H<(int)(T-exec[i-1][j-1]))?H:T-exec[i-1][j-1];
			for(int t=0;t<=borne;t++){
				fpl<<" + "<<((t<=H-exec[i-1][j-1])?(l1*w[i-1]):(l1*w[i-1]*(H-t)/exec[i-1][j-1]))<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;					
			}				
		}		
	}
	fpl<<" >= "<<l1*r1<<endl;
	
	fpl<<"z";
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			for(int t=0;t<=(int)(T-exec[i-1][j-1]);t++){
				fpl<<" - "<<l2*c[i-1]*(t+exec[i-1][j-1])<<" u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;	
			}			
		}		
	}
	fpl<<" >= "<<-l2*r2<<endl;

	//1 Un lot s'execute une unique fois
	cout<<"Ecriture des contraintes de type 1..."<<endl;
	
	for(int i=1;i<=n;i++){
		for(unsigned int j=1;j<=corresp[i-1].size();j++){
			for(int t=0;t<=(int)(T-exec[i-1][j-1]);t++){
				//cout<<T<<" "<<exec[i-1][j-1]<<endl;cin.get();
				fpl<<((j+t==1)?" ":" +")<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t;			
			}
		}
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
				for(unsigned int i0=1;i0<=invCorresp[indicemachine-1].size();i0++){
					int indiceLotCourant = invCorresp[indicemachine-1][i0-1];
					int indicejPouri0 = indexMachineQualif(indiceLotCourant,indicemachine);
					if(indiceLotCourant != i && indicejPouri0!=-1){//pour tout lot != i pouvant aller sur la machine
						int borne = (T-(exec[indiceLotCourant-1][indicejPouri0])<=(t+exec[i-1][j-1]-1+setup[families[indicemachine-1]-1][batch[i-1][j-1]-1][batch[indiceLotCourant-1][indicejPouri0]-1]))?T-(exec[indiceLotCourant-1][indicejPouri0]):(t+exec[i-1][j-1]-1+setup[families[indicemachine-1]-1][batch[i-1][j-1]-1][batch[indiceLotCourant-1][indicejPouri0]-1]);
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
							  int borne= (T-(int)(exec[i0-1][j0-1])<=(t+(int)(exec[i-1][j-1])))?T-(int)(exec[i0-1][j0-1]):(t+(int)(exec[i-1][j-1]));
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

	fpl<<"Binaries"<<endl;
	for(unsigned int t=0;t<=T;t++)
		for(int i=1;i<=n;i++)		
			for(unsigned int j=1;j<=corresp[i-1].size();j++)
				if(t<= T - exec[i-1][j-1])
					fpl<<"u"<<i<<"_"<<corresp[i-1][j-1]<<"_"<<t<<endl;

	fpl<<"End";
	fpl.close();
	return(EXIT_SUCCESS);
}


