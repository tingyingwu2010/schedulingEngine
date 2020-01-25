#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<queue>

using namespace std;

void Instance::printVector(vector<unsigned int> v){
	for(unsigned int i=0;i<v.size();i++)
		cout<<v[i]<<"\t";
	cout<<endl;
}

//Procédure de croisement
Individu Instance::CroisementMemetique(Individu parent1, Individu parent2){
	vector<int> motif(n);
	vector<int> buf;
	vector<int> indicesDe1;
	Individu enfant(parent1.fctObj,n,corresp);
	for(int i=0;i<n;i++){
		motif[i]=rand()%2;
		if(motif[i]==1)
			indicesDe1.push_back(i);
	}	
	int ii=0;
	int j=0;
	while(ii<n){
		if(motif[ii]==0){
			enfant.machines[ii]=parent1.machines[ii];
			enfant.jobs[ii]=parent1.jobs[ii];
		}
		else{
			enfant.machines[ii]=-1;
			buf.push_back(parent1.jobs[ii]);
		}
		ii++;
	}
	j=0;
	for(int k=0;k<n;k++){
		int mem = member(buf,parent2.jobs[k]);
		if(mem!=-1){
			enfant.machines[indicesDe1[j]] = parent2.machines[k];
			enfant.jobs[indicesDe1[j]] = parent2.jobs[k];
			buf.erase(buf.begin() + mem);
			j++;
		}
	}
	return enfant;
}

//Procédure spéciale de calcul de la nouvelle fitness pour la somme des Ci apres le voisinage V1
double Instance::nouvelleFitness(vector<vector<unsigned int> > &enfant, vector<unsigned int> &completionTimes, int indice, double fitness){	
	vector<int> tDispo(m,0); //Date de disponibilité au plus tôt de chaque machine (t=0 au début)
	vector<int> D(m,0); //Dernier job sur chaque machine (nécessaire pour calculer le setup) : 0 si aucun
	for(int j=0;j<m;j++){
		int k(indice-1);
		while(k>=0){
			if(enfant[0][k]-1==(unsigned int)(j)){
				tDispo[j]=completionTimes[k];	
				D[j]=enfant[1][k];		
				k=0;			
			}
			k--;
		}
	}
	vector<unsigned int> S(l); //Machine courante de chaque masque
	vector<int> tMasque(l,0); //Date de dispo au plus tôt de chaque masque
	for(int i=0;i<l;i++){
		S[i]=initMask[i];
		int k(indice-1);
		while(k>=0){
			if(phi[enfant[1][k]-1]-1==i){
				S[i]=enfant[0][k];
				tMasque[i]=completionTimes[k];				
				k=0;			
			}
			k--;
		}
	}
	///////////////////////////////////////////////////////////////
	//Les données sont désormais préparées
	//On commence à l'indice donné
	for(int j=indice;j<n;j++){//On parcourt l'individu
		//cout<<"Masque A"<<phi[jobact-1]<<", Job à la position "<<j+1<<" : job J"<<jobact<<" pour la machine M"<<machact<<" : "<<endl;
		if(S[phi[enfant[1][j]-1]-1]!=enfant[0][j]){
			tMasque[phi[enfant[1][j]-1]-1]++; //Si il faut déplacer le masque
			S[phi[enfant[1][j]-1]-1]=enfant[0][j];
		}
		int imQ(0);
		int imQ1=indexMachineQualif(enfant[1][j],enfant[0][j]);
		if(D[enfant[0][j]-1]!=0)
			imQ=indexMachineQualif(D[enfant[0][j]-1],enfant[0][j]);
		int Beta=(D[enfant[0][j]-1]==0)?0:setup[families[enfant[0][j]-1]-1][batch[D[enfant[0][j]-1]-1][imQ]-1][batch[enfant[1][j]-1][imQ1]-1];
		//formule
		int Fin=max(tMasque[phi[enfant[1][j]-1]-1],tDispo[enfant[0][j]-1]+Beta)+exec[enfant[1][j]-1][imQ1];
		//On actualise toutes les infos, sauf si on est à la fin, auquel cas ça ne sert à rien
		if(j-n+1!=0){
			D[enfant[0][j]-1]=enfant[1][j];
			tDispo[enfant[0][j]-1]=Fin;
			tMasque[phi[enfant[1][j]-1]-1]=Fin;
		}
		//Actualisation de la fitness!!
		fitness=fitness+c[enfant[1][j]-1]*(Fin-completionTimes[j]);
		//La date de fin est conservée
		completionTimes[j]=Fin;
	}
	return fitness;
}

//Méthode de décodage
vector<vector<vector<int> > > Instance::HeuristiqueMemetique(Individu individu/**, vector<unsigned int> &completionTimes**/){
	vector<vector<vector<int> > > scheduleFinal(m);
	vector<int> tDispo(m,0); //Date de disponibilité au plus tôt de chaque machine (t=0 au début)
	vector<int> S(l); //Machine courante de chaque masque
	vector<int> tMasque(l,0); //Date de dispo au plus tôt de chaque masque
	vector<int> D(m,0); //Dernier job sur chaque machine (nécessaire pour calculer le setup) : 0 si aucun
	for(int i=1;i<=l;i++)
		S[i-1]=initMask[i-1];
	for(int j=0;j<n;j++){//On parcourt l'individu
		int jobact=individu.jobs[j];
		int machact=individu.machines[j];

		vector<int> elt(2);
		elt[0]=jobact;

		//cout<<"Masque A"<<phi[jobact-1]<<", Job à la position "<<j+1<<" : job J"<<jobact<<" pour la machine M"<<machact<<" : "<<endl;
		if(S[phi[jobact-1]-1]!=machact){
			tMasque[phi[jobact-1]-1]++;
			S[phi[jobact-1]-1] = machact;
		}
		int imQ;
		if(D[machact-1]!=0)
			imQ=indexMachineQualif(D[machact-1],machact);
		int Beta = (D[machact-1]==0)?0:setup[families[machact-1]-1][batch[D[machact-1]-1][imQ]-1][batch[jobact-1][indexMachineQualif(jobact,machact)]-1];
		//formule
		int Debut =(tMasque[phi[jobact-1]-1] > tDispo[machact-1]+Beta)?tMasque[phi[jobact-1]-1]:tDispo[machact-1]+Beta;
		int Fin = Debut + exec[jobact-1][indexMachineQualif(jobact,machact)];
		elt[1]=Debut;
		scheduleFinal[machact-1].push_back(elt);	
		//On actualise toutes les infos
		D[machact-1] = jobact;
		tDispo[machact-1] = Fin;
		tMasque[phi[jobact-1]-1] = Fin;
		//La date de fin est conservée
		/**completionTimes[j]=Fin;**/
	}
	return scheduleFinal;
}


double Instance::evaluateSchedule(vector< vector< vector<int> > > schedule, int fct, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){ 
	//fct : 1->wafer, 2->wiCi, 3->mask moves, 4->Chebychev tricriteria
	if(fct == 1){ // Number of wafers
		double sum = 0;
		for(int j=0;j<m;j++){
			for(unsigned int i=0;i<schedule[j].size();i++){
				int job = schedule[j][i][0];
				int ti = schedule[j][i][1];
				int h = exec[job-1][indexMachineQualif(job,j+1)];
				double theta =(H<ti+h)?H:ti+h;
				theta = (theta-ti)/h;
				sum = sum+((H<ti)?0:w[job-1]*theta);		
			}
		}
		obj1 = sum;
		return sum;
	}
	if(fct == 2){ // Weighted sum of completion times
		double sum=0;
		for(int j=0;j<m;j++){
			for(unsigned int i=0;i<schedule[j].size();i++){
				int job = schedule[j][i][0];
				int Ci = schedule[j][i][1] + exec[job-1][indexMachineQualif(job,j+1)];
				sum = sum+c[job-1]*Ci;		
			}
		}
		obj2 = sum;
		return sum;
	}
	else if(fct == 3){ // Number of mask moves
		vector< vector<int> > maskM(l);
		vector< vector<int> > lots(n);
		double deplmask = 0;
		int nbl = 0;
		for(int i=0;i<l;i++)
			maskM[i].push_back(initMask[i]);
		int k = 0;
		for(int j=0;j<m;j++){
			nbl+=schedule[j].size();
			for(unsigned int i=0;i<schedule[j].size();i++){
				lots[k].push_back(schedule[j][i][0]);
				lots[k].push_back(schedule[j][i][1]);
				lots[k].push_back(j+1);
				k++;
			}
		}
		// Sorting in increasing order of starting times
		for(int i=0;i<n;i++){
			int min = i;
			for(int j=i+1;j<n;j++)
				if(lots[j][1] < lots[min][1])
					min = j;		
			lots[i].swap(lots[min]);
		}
		// It is sorted. We can now find the number of mask moves
		for(int i=0;i<n;i++){
			if(maskM[phi[lots[i][0]-1]-1][maskM[phi[lots[i][0]-1]-1].size()-1]!=lots[i][2]){
				maskM[phi[lots[i][0]-1]-1].push_back(lots[i][2]);
				deplmask=deplmask+1;
			}
		}
		obj3 = deplmask;
		return deplmask;
	}
	else{ // Chebychev tricriteria
		double c1(evaluateSchedule(schedule,1,r1,r2,r3,l1,l2,l3,epsilon));
		double c2(evaluateSchedule(schedule,2,r1,r2,r3,l1,l2,l3,epsilon));
		double c3(evaluateSchedule(schedule,3,r1,r2,r3,l1,l2,l3,epsilon));
		c1=l1*(r1-c1);
		c2=l2*(c2-r2);
		c3=l3*(c3-r3);
		return (max(c1,max(c2,c3)) + epsilon*(c1+c2+c3));
	}
}


//Calcul de la qualité de la solution (fitness)
double Instance::fitnessMemetique(Individu individu, int fctObj, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	if(fctObj==4){
		double c1(fitnessMemetique(individu,1,r1,r2,r3,l1,l2,l3,epsilon));
		double c2(fitnessMemetique(individu,2,r1,r2,r3,l1,l2,l3,epsilon));
		double c3(fitnessMemetique(individu,3,r1,r2,r3,l1,l2,l3,epsilon));
		c1=l1*(r1-c1);
		c2=l2*(c2-r2);
		c3=l3*(c3-r3);
		double comp(epsilon*(c1+c2+c3));
		if(c1>c2){
			if(c1>c3){
				return c1+comp;
			}
			else{
				return c3+comp;
			}
		}
		if(c2>c3){
			return c2+comp;
		}
		else{
			return c3+comp;
		}
	}
	double fitn(0.0);
	vector<int> tDispo(m,0); //Date de disponibilité au plus tôt de chaque machine (t=0 au début)
	vector<int> S(l); //Machine courante de chaque masque
	vector<int> tMasque(l,0); //Date de dispo au plus tôt de chaque masque
	vector<int> D(m,0); //Dernier job sur chaque machine (nécessaire pour calculer le setup) : O si aucun
	for(int i=1;i<=l;i++)
		S[i-1]=initMask[i-1];
	for(int j=0;j<n;j++){//On parcourt l'individu
		int jobact=individu.jobs[j];
		int machact=individu.machines[j];
		//cout<<"Masque A"<<phi[jobact-1]<<", Job a la position "<<j+1<<" : job J"<<jobact<<" pour la machine M"<<machact<<" : "<<endl;
		if(S[phi[jobact-1]-1]!=machact){
			tMasque[phi[jobact-1]-1]++;
			S[phi[jobact-1]-1]=machact;
			if(fctObj==3)
				fitn++;
		}
		//cout<<"Le masque A"<<phi[jobact-1]<<" est maintenant sur M"<<S[phi[jobact-1]-1]<<endl;
		//cout<<"La date de disponibilite du masque A"<<phi[jobact-1]<<" est "<<tMasque[phi[jobact-1]-1]<<endl;
		//cout<<"La date de de disponibilite de la machine M"<<machact<<" est "<<tDispo[machact-1]<<endl;
		int imQ(0);
		int imQ1(indexMachineQualif(jobact,machact));
		if(D[machact-1]!=0)
			imQ=indexMachineQualif(D[machact-1],machact);
		int Beta=(D[machact-1]==0)?0:setup[families[machact-1]-1][batch[D[machact-1]-1][imQ]-1][batch[jobact-1][imQ1]-1];
		//cout<<"Le temps de setup pour passer du job precedent sur la machine M"<<machact<<" a J"<<jobact<<" est de "<<Beta<<endl;
		//formule
		int Debut=(tMasque[phi[jobact-1]-1]>tDispo[machact-1]+Beta)?tMasque[phi[jobact-1]-1]:tDispo[machact-1]+Beta;
		int Fin=Debut + exec[jobact-1][imQ1];		
		//On actualise toutes les infos
		D[machact-1] = jobact;
		tDispo[machact-1] = Fin;
		tMasque[phi[jobact-1]-1] = Fin;
		if(fctObj==1)
			if(H>Debut){
				if(H>=Fin)
				fitn+=w[jobact-1];
			else
				fitn+=(float)(w[jobact-1]*(H-Debut)/exec[jobact-1][imQ1]);
		}		
		if(fctObj==2)
			fitn+=c[jobact-1]*Fin;
	}
	return fitn;
}

//Procédure de recherche locale.
double Instance::rechercheLocaleMemetique(Individu &enfant, int nbiter, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	//Variables pour quantifier l'utilité du test de voisinage significatif
	int fctObj(enfant.fctObj);
	double nbTentatives(0.0);
	double delta(0.0);
	double fit(0.0);
	enfant.fitness=fitnessMemetique(enfant,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
	fit=enfant.fitness;
	//Meilleure solution trouvée dans la recherche locale
	Individu best(enfant);
	//cout<<"Debut : "<<enfant<<endl;
	for(int i=0;i<2*nbiter;i++){
		//cout<<"     "<<fit<<endl;
		//On crée un vecteur temporaire de Completion Times, initialisé au vecteur principal par constructeur de copie
		/**vector<unsigned int> tmpCompletion(completionTimes);**/
		int p=rand()%100;
		if(p<p1){//Voisinage V1 : changement de machine

			/**
			//Checkpoint
			cout<<"INDIVIDU DE DEPART : "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Fitness = "<<fit<<endl;
			///////////////////////////////////////////////////////////////
			**/
			int indiceChgt=rand()%n;
			int indiceAncMach=enfant.machines[indiceChgt];
			int indNouvMach=rand()%(corresp[enfant.jobs[indiceChgt]-1].size()); //La nouvelle machine doit etre qualifiée
			enfant.machines[indiceChgt]=corresp[enfant.jobs[indiceChgt]-1][indNouvMach];

			double fitact(0.0);
			
			fitact=fitnessMemetique(enfant,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			/////////////////////////////////////////////////////////////////

			/**
			//Checkpoint
			cout<<endl<<"........................"<<endl<<"INDIVIDU MODIFIE DONT LA MACHINE CHANGE EN POSITION "<<indiceChgt<<": "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Nouvelle fitness = "<<fitact<<endl;
			cout<<"Delta : "<<fitact-fit<<endl;
			cin.get();
			///////////////////////////////////////////////////////////////
			**/

			delta=fitact-fit;
			double pr = (double)(rand()/((double)(RAND_MAX)));
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0)) /**/ && ( (fctObj==1 && pr>exp(delta/T)) || (fctObj!=1 && pr>exp(-delta/T)) ) /**/){//On n'améliore pas
				enfant.machines[indiceChgt]=indiceAncMach;
			}
			else{ //On améliore, on actualise aussi les C_i
				fit=fitact;
				enfant.fitness=fitact;
			}
		}
		else if(p<p1+p2){ //Voisinage V2
			/**
			//Checkpoint
			cout<<"INDIVIDU DE DEPART : "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Fitness = "<<fit<<endl;
			///////////////////////////////////////////////////////////////
			**/
			nbTentatives++;
			int indiceBoug=rand()%n;
			int indDep=rand()%(n+1);//On peut aller à n+1 pour le mettre en dernière position
			/**int aa(indiceBoug);**/
			/**int bb(indDep);**/
			
			/**
			//Utilisation d'une propriété
			int indice1=indiceBoug-1;
			int indice2=indiceBoug+1;	
			while((indice1>=0)&&(enfant[0][indice1]!=enfant[0][indiceBoug] && phi[enfant[1][indice1]-1]!=phi[enfant[1][indiceBoug]-1]))
				indice1--;
			while((indice2<n)&&(enfant[0][indice2]!=enfant[0][indiceBoug] && phi[enfant[1][indice2]-1]!=phi[enfant[1][indiceBoug]-1]))
				indice2++;
			ratioIntervalleInterdit+=100.0*(indice2-indice1+1)/n;
			//Choix de l'indice de déplacement
			if(indice1>=0 && indice2<n){
			        indDep = rand()%(indice1+n-indice2+2);
				if(indDep>indice1)
					indDep = indice2 + indDep - indice1 - 1;	
			}
			else if(indice1>=0){
				indDep = rand()%(indice1+1);
			}
			else if(indice2<n){
				indDep = indice2+rand()%(n-indice2+1);
			}
			**/
			//On a l'indice, on effectue maintenant le déplacement
			enfant.machines.insert(enfant.machines.begin()+indDep,enfant.machines[indiceBoug]);
			enfant.jobs.insert(enfant.jobs.begin()+indDep,enfant.jobs[indiceBoug]);
			if(indDep<=indiceBoug)
				indiceBoug++;
			else
				indDep--;
			enfant.machines.erase(enfant.machines.begin()+indiceBoug);
			enfant.jobs.erase(enfant.jobs.begin()+indiceBoug);

			//On teste si le changement n'est pas améliorant
			double fitact(0.0);
			/**if(fctObj==2){
				fitact=nouvelleFitness(enfant,tmpCompletion,min(aa,bb),fit);
			}
			else**/
			fitact=fitnessMemetique(enfant,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			///////////////////////////////////////////////////////////////////////


			/**
			//Checkpoint
			cout<<endl<<"........................"<<endl<<"INDIVIDU DONT ON A DEPLACE LA POSITION "<<aa<<" A LA POSITION "<<bb<<": "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Nouvelle fitness = "<<fitact<<endl;
			cout<<"Delta : "<<fitact-fit<<endl;
			cin.get();
			///////////////////////////////////////////////////////////////
			**/


			delta=fitact-fit;
/**...**/		double pr = (double)(rand()/((double)(RAND_MAX)));
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0)) /**/ && ((fctObj==1 && pr>exp(delta/T)) || (fctObj!=1 && pr>exp(-delta/T))) /**/){//On n'améliore pas. On revient à la solution d'origine
				enfant.machines.insert(enfant.machines.begin()+indiceBoug,enfant.machines[indDep]);
				enfant.jobs.insert(enfant.jobs.begin()+indiceBoug,enfant.jobs[indDep]);
				if(indDep>indiceBoug)
					indDep++;
				enfant.machines.erase(enfant.machines.begin()+indDep);
				enfant.jobs.erase(enfant.jobs.begin()+indDep);
				enfant.fitness=fit;
			}
			else{ //On améliore, on actualise aussi les C_i
				fit=fitact;
				enfant.fitness=fitact;
			}
		}
		else{ //Voisinage V3

			/**
			//Checkpoint
			cout<<"INDIVIDU DE DEPART : "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Fitness = "<<fit<<endl;
			///////////////////////////////////////////////////////////////
			**/

			int indiceBoug=rand()%n;
			int indDep=rand()%n;
			//On a l'indice, on effectue maintenant l'échange
			swap(enfant.machines[indiceBoug],enfant.machines[indDep]);
			swap(enfant.jobs[indiceBoug],enfant.jobs[indDep]);
			//On teste si le changement n'est pas améliorant
			double fitact(0.0);
			/**if(fctObj==2){
				fitact=nouvelleFitness(enfant,tmpCompletion,min(indiceBoug,indDep),fit);
			}
			else**/
			fitact=fitnessMemetique(enfant,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			/**
			//Checkpoint
			cout<<endl<<"........................"<<endl<<"INDIVIDU DONT ON A ECHANGE LES POSITIONS "<<indiceBoug<<" ET "<<indDep<<": "<<endl;
			printVector(enfant[1]);
			cout<<endl;
			printVector(enfant[0]);
			cout<<endl;
			cout<<"Liste des dates de fin : "<<endl;
			printVector(completionTimes);
			cout<<endl<<"Nouvelle fitness = "<<fitact<<endl;
			cout<<"Delta : "<<fitact-fit<<endl;
			cin.get();
			///////////////////////////////////////////////////////////////
			**/

			delta=fitact-fit;
			double pr=(double)(rand()/((double)(RAND_MAX)));
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0)) /**/ && ((fctObj==1 && pr>exp(delta/T)) || (fctObj!=1 && pr>exp(-delta/T))) /**/){//On n'améliore pas. On revient à la solution d'origine
				swap(enfant.machines[indiceBoug],enfant.machines[indDep]);
				swap(enfant.jobs[indiceBoug],enfant.jobs[indDep]);
				enfant.fitness=fit;
			}
			else{ //On améliore, on actualise aussi les C_i
				fit=fitact;
				enfant.fitness=fitact;
			}
		}
		if(i%palier==0) T=T*a; //actualisation de la température
		//On actualise la meilleure solution trouvée dans la recherche locale			
		if(enfant>best){
			best=enfant;
		}
	}
	enfant=best;
	return best.fitness;
}

//Recuit simulé
double Instance::recuitSimule(Individu &i1, int nbIter, int conv, int stagnation, double tauxRedemarrage, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	int fctObj(i1.fctObj);
	double Tredemarrage(T/10.0);
	i1.fitness=fitnessMemetique(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
	double fit(i1.fitness);
	double delta(0.0);
	int nbSansAm(0);
	int meilleureSansAm(0);
	//Meilleure solution trouvée dans la recherche locale
	Individu best(i1);
	//cout<<"Solution de départ : "<<fit<<endl;
	for(int i=0;i<nbIter;i++){
		int p=rand()%100;
		if(p<p1){//Voisinage V1 : changement de machine
			int indiceChgt=rand()%n;
			int indiceAncMach=i1.machines[indiceChgt];
			int indNouvMach=rand()%(corresp[i1.jobs[indiceChgt]-1].size()); //La nouvelle machine doit etre qualifiée
			i1.machines[indiceChgt]=corresp[i1.jobs[indiceChgt]-1][indNouvMach];
			double fitact(0.0);
			fitact=fitnessMemetique(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta=fitact-fit;
			double pr=(double)(rand()/((double)(RAND_MAX)));
			if((fctObj==1 && delta<=0) || (fctObj!=1 && delta>=0))
				nbSansAm++;
			else
				nbSansAm=0;
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0)) && ( (fctObj==1 && pr>exp(delta/T)) || (fctObj!=1 && pr>exp(-delta/T)) )){ //On n'améliore pas
				i1.machines[indiceChgt]=indiceAncMach;
				i1.fitness=fit;
			}
			else{ //On améliore, on actualise aussi les C_i
				fit=fitact;
				i1.fitness=fitact;
			}
		}
		else if(p<p1+p2){ //Voisinage V2
			int indiceBoug=rand()%n;
			int indDep=rand()%(n+1);//On peut aller à n+1 pour le mettre en dernière position
			//On a l'indice, on effectue maintenant le déplacement
			i1.machines.insert(i1.machines.begin()+indDep,i1.machines[indiceBoug]);
			i1.jobs.insert(i1.jobs.begin()+indDep,i1.jobs[indiceBoug]);
			if(indDep<=indiceBoug)
				indiceBoug++;
			else
				indDep--;
			i1.machines.erase(i1.machines.begin()+indiceBoug);
			i1.jobs.erase(i1.jobs.begin()+indiceBoug);
			//On teste si le changement n'est pas améliorant
			double fitact(0.0);
			fitact=fitnessMemetique(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta=fitact-fit;
			double pr=(double)(rand()/((double)(RAND_MAX)));
			if((fctObj==1 && delta<=0) || (fctObj!=1 && delta>=0))
				nbSansAm++;
			else
				nbSansAm=0;
			//cout<<"delta="<<delta<<" | pr="<<pr<<" | "<<exp(-delta/T)<<endl;cin.get();
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0))&&((fctObj==1 && pr>exp(delta/T)) || (fctObj!=1 && pr>exp(-delta/T)))){ //On n'améliore pas. On revient à la solution d'origine
				i1.machines.insert(i1.machines.begin()+indiceBoug,i1.machines[indDep]);
				i1.jobs.insert(i1.jobs.begin()+indiceBoug,i1.jobs[indDep]);
				if(indDep>indiceBoug)
					indDep++;
				i1.machines.erase(i1.machines.begin()+indDep);
				i1.jobs.erase(i1.jobs.begin()+indDep);
				i1.fitness=fit;
			}
			else{ //On améliore,
				fit=fitact;
				i1.fitness=fitact;
			}
		}
		else{ //Voisinage V3
			int indiceBoug=rand()%n;
			int indDep=rand()%n;
			//On a l'indice, on effectue maintenant l'échange
			swap(i1.machines[indiceBoug],i1.machines[indDep]);
			swap(i1.jobs[indiceBoug],i1.jobs[indDep]);
			//On teste si le changement n'est pas améliorant
			double fitact(0.0);
			fitact=fitnessMemetique(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta=fitact-fit;
			double pr=(double)(rand()/((double)(RAND_MAX)));
			if((fctObj==1 && delta<=0) || (fctObj!=1 && delta>=0))
				nbSansAm++;
			else
				nbSansAm=0;
			if(((fctObj==1 && delta<0) || (fctObj!=1 && delta>0)) && ((fctObj==1 && pr>exp(delta/T))||(fctObj!=1 && pr>exp(-delta/T)))){//On n'améliore pas. On revient à la solution d'origine
				swap(i1.machines[indiceBoug],i1.machines[indDep]);
				swap(i1.jobs[indiceBoug],i1.jobs[indDep]);
				i1.fitness=fit;
			}
			else{ //On améliore,
				fit=fitact;
				i1.fitness=fitact;
			}
		}
		if(i%palier==0) T=T*a; //actualisation de la température
		if(nbSansAm>=conv){
			T=Tredemarrage;
			nbSansAm=0;
			//cout<<"Convergence!!!!!!"<<endl;
			//cin.get();
		}
		//On actualise la meilleure solution trouvée dans la recherche locale
		if(i1>best){
			best=i1;
			meilleureSansAm=0;
		}
		else
			meilleureSansAm++;
		//Lorsqu'on a beaucoup stagné au niveau du best, ça veut souvent dire qu'on dégrade trop au redémarrage de T
		if(meilleureSansAm>=stagnation){ //Soit on diminue Tredemarrage, soit on le fait descendre plus vite en diminuant a
			Tredemarrage=Tredemarrage*tauxRedemarrage;
			//a=a*tauxRedemarrage;
			meilleureSansAm=0;
		}
		//if(i%1000==0)
		//	cout<<i1.fitness<<"\t"<<best.fitness<<"\t"<<exp(-1/T)<<endl;
		//cin.get();
	}
	i1=best;
	return best.fitness;
}

//Algorithme principal
vector<vector<vector<int> > > Instance::POPMemetique(int Npop, int nbIter, int conv, int stagnation, double tauxRedemarrage, int fctObj, int Ktour, double Rep, double Ts, int ReinIter, int nbIterSA, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, unsigned int duree, double epsilon){
	int Desc=Rep*Npop;
	int Surv=Desc*Ts;	
	int nbIterSansAm=0;	
	cout<<"Nombre d'iterations sans amelioration avant reinitialisation : "<<ReinIter<<endl<<endl;
	int sansAmelioration=0;	

	int iterLoc=nbIter;

	priority_queue<Individu> population;
	vector<Individu> pop;
	for(int i=0;i<Npop;i++){
		Individu i1(fctObj,n,corresp);
		//rechercheLocaleMemetique(i1,iterLoc,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3);
		recuitSimule(i1,iterLoc,conv,stagnation,tauxRedemarrage,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3,epsilon);
		pop.push_back(i1);
	}
	for(int i=0;i<Npop;i++)
		population.push(pop[i]);

	Individu bestSol(population.top());

	int nIt(0);
	//Il faut avoir une gestion de la bestSol depuis la dernière réinitialisation pour ne pas réinitialiser et re-converger bêtement juste parce qu'on n'a pas encore modifié la bestSol; on appelle la variable bestCourant;
	double bestCourant(bestSol.fitness);
	double bestCout(bestSol.fitness);
	clock_t dureeArret=clock()+(duree*CLOCKS_PER_SEC); //Date limite

	while(nbIterSansAm<nbIterSA && clock()<dureeArret){
		//Petite procédure d'affichage		
		if(nIt%20 == 0)
			cout<<"Cout de la meilleure solution depuis : "<<endl<<"le debut\t(la derniere reinitialisation)"<<endl;
		priority_queue<Individu> enfants;
		priority_queue<Individu> temp;
		vector<Individu> enf;

		//Sélection puis Croisement pour former les Npop*Rep individus de la génération suivante
		for(int j=0;j<Desc;j++){
			//Sélection par K-Tournoi
			int selec1=rand()%Npop;
			int ind(0);
			for(int i=0;i<Ktour-1;i++){
				ind=rand()%Npop;
				if(pop[ind]>pop[selec1]){
					selec1=ind;
				}
			}
			//selec1 est l'indice du premier parent sélectionné par K-Tournoi
			int selec2=rand()%Npop;
			ind=0;
			for(int i=0;i<Ktour;i++){
				ind=rand()%Npop;
				if(pop[ind]>pop[selec2]){
					selec2=ind;
				}
			}
			//selec2 est l'indice du second parent sélectionné par K-Tournoi
			//On effectue un croisement (crossover) pour former un enfant
			enf.push_back(CroisementMemetique(pop[selec1],pop[selec2]));
			//rechercheLocaleMemetique(enf[j],iterLoc,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3);
			recuitSimule(enf[j],iterLoc,conv,stagnation,tauxRedemarrage,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3,epsilon);
			enfants.push(enf[j]);
		}
		if(Desc!=0){
			for(int j=0;j<Surv;j++){ //On chope les Surv meilleurs enfants
				temp.push(enfants.top());
				enfants.pop();
			}
			for(int j=0;j<Npop-Surv;j++){ //On chope les Npop-Surv meilleurs parents
				temp.push(population.top());
				population.pop();
			}
			while(!population.empty()){ //On vide la population
				population.pop();		
			}
			for(int j=0;j<Npop;j++){ //On réunit le tout dans la population
				population.push(temp.top());
				pop[j]=temp.top();
				temp.pop();		
			}
		}

/**		for(int j=0;j<Npop;j++){
			rechercheLocaleMemetique(pop[j],iterLoc,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3);
			if(pop[j]>bestSol){
				population.push(pop[j]);
				bestSol=pop[j];
			}
		}
**/
		//Gestion de bestCout = meilleure solution globale
		if(fctObj==1 && population.top().fitness>bestCout){
			bestCout=population.top().fitness;		
			bestSol=population.top();
			nbIterSansAm=0;
		}
		else if(fctObj!=1 && population.top().fitness<bestCout){
			bestCout=population.top().fitness;
			bestSol=population.top();
			nbIterSansAm=0;
		}
		else{
			nbIterSansAm++;
		}
		//Gestion de bestCourant = meilleure solution depuis la dernière réinitialisation
		if(fctObj==1 && population.top().fitness>bestCourant){
			bestCourant=population.top().fitness;
			sansAmelioration=0;
		}
		else if(fctObj!=1 && population.top().fitness<bestCourant){
			bestCourant=population.top().fitness;
			sansAmelioration=0;
		}
		else{
			sansAmelioration++;
		}
		//En cas de convergence précoce, on réinitialise : 
		if(sansAmelioration==ReinIter){
			cout<<"On reinitialise la population : "<<endl<<endl<<endl;
			while(!population.empty()){ //On vide la population
				population.pop();		
			}
			pop[0]=bestSol; //--> On garde ou pas la bestSol dans la population réinitialisée?
			population.push(pop[0]);
			for(int ii=1;ii<Npop;ii++){
				Individu i1(fctObj,n,corresp);
				//rechercheLocaleMemetique(i1,iterLoc,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3);
				recuitSimule(i1,iterLoc,conv,stagnation,tauxRedemarrage,p1,p2,T,palier,a,r1,r2,r3,l1,l2,l3,epsilon);
				pop[ii]=i1;
			}
			for(int ii=1;ii<Npop;ii++)
				population.push(pop[ii]);
			//Actualisation de bestCourant
			bestCourant=population.top().fitness;
			sansAmelioration=0;
		}
		nIt++;
		cout<<bestCout<<"\t\t("<<bestCourant<<")"<<endl;
	}
	
	//On construit une solution à partir de la meilleure trouvée, qu'on conserve
	return HeuristiqueMemetique(bestSol);
}


