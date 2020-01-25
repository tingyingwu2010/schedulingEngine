#include"ProblemePhoto.h"
#include<sstream>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<cmath>

using namespace std;

vector<vector<vector<vector<int> > > > ProblemePhoto::tricritere(double r1, double r2, double l1, double l2, int ecart, int nbiter, int p1, double T, int palier, double a){
	//La premiere etape consiste a recuperer une solution qui optimise le critere 3
	ecrirePLNE(3,0);
	ostringstream commande; 
 	commande.str("");	
	commande<<"scip -c \"read ./FichiersPL/"<<nomInstance<<".lp" << " opt write solution "<<"./FichiersSol/"<<nomInstance<<".sol"<<" quit \""<<endl;
	cout<<commande.str().c_str();
	if (system(commande.str().c_str()))
		cout<<"Fin de l'appel au solveur"<<endl;
	//Fonction pour lire le fichier .sol et construire un ordo
	vector<vector<vector<int > > > ordo=formerOrdoPLNE(3);
	afficheOrdoConsole(ordo);
	EvaluationOrdonnancement(ordo,3,0,0,0,0,0,0,0);
	double c3=getcout(3);
	cout<<endl<<endl<<"La solution trouvee est de cout "<<c3<<endl<<endl<<endl;
	if((getMSol3() == -1) || (c3<getMSol3()))
		ecrireMeilleureConnue(3,c3,ordo);
	//Ici, on a le cout de la meilleure solution dans c3 et la solution elle-meme dans ordo
	//L'etape suivante consiste a construire un individu, sous la forme des individus du memetique, a partir de ordo
	Individu chromosome = creerIndividuInitial(ordo);
	rechercheLocaleTricritere(chromosome,nbiter,p1,T,palier,a,r1,r2,l1,l2);
	
	vector< vector< vector< vector<int> > > > solutions(ecart+1);
	solutions[0]=ordo; //On a range une solution faiblement Pareto-optimale
	//On reproduit le meme processus en creant des solutions initiales moins bonnes sur f3
	for(int k=1;k<=ecart;k++){
		//Fonction qui construit une solution a partir du PLNE auquel on ajoute la contrainte qui impose la borne inf
		ecrirePLNEmasques(c3+k);
		ostringstream commande1; 
 		commande1.str("");
		commande1<<"scip -c \"read ./FichiersPL/"<<nomInstance<<"_borne_"<<c3+k<<".lp" << " opt write solution "<<"./FichiersSol/"<<nomInstance<<".sol"<<" quit \""<<endl;
		cout<<commande1.str().c_str();
		if (system(commande1.str().c_str()))
			cout<<"Fin de l'appel au solveur"<<endl;
		vector<vector<vector<int > > > ordo1=formerOrdoPLNE(3);
		Individu chromosome1=creerIndividuInitial(ordo1);
		//On lui applique de la recherche locale pour ameliorer les 2 autres criteres
		rechercheLocaleTricritere(chromosome1,nbiter,p1,T,palier,a,r1,r2,l1,l2);
		//On le stocke dans un tableau de solutions
		solutions[k]=ordo1;
	}
	//A la fin de la boucle, on a ecart+1 solutions qu'on peut comparer et évaluer sur les 3 critères. L'une d'elles est faiblement
	//Pareto-optimale : la solutions[0].
	return solutions;
}

Individu ProblemePhoto::creerIndividuInitial(vector<vector<vector<int > > > ordo){
	Individu individu(n);
	for(int i=1;i<=n;i++){
		int j=0;
		while(ordo[j].size()==0 && j<=m-1)
			j++;
		int indmin=j;
		
		for(int k=j+1;k<m;k++)
			if(ordo[k].size()>=1 && ordo[k][0][1] < ordo[indmin][0][1])
				indmin=k;
		individu.jobs[i-1]=indmin+1;
		individu.machines[i-1]=ordo[indmin][0][0];
		ordo[indmin][0].erase(ordo[indmin][0].begin());
		ordo[indmin][0].erase(ordo[indmin][0].begin());
		ordo[indmin].erase(ordo[indmin].begin());
	}
	return individu;
}

double ProblemePhoto::rechercheLocaleTricritere(Individu individu, int nbiter, int p1, double T, int palier, double a, double r1, double r2, double l1, double l2){
	//Pour savoir quel voisinage utiliser, on met en place une distribution de probas sur (V1,V2)
	double pV1=p1/100.0;//Probabilité d'opérer le voisinage V1
	double delta=0.0;

	double fit = fitnessMemetique(individu,4,r1,r2,0,l1,l2,0,0);

	//Meilleure solution trouvee dans la recherche locale	
	double cbest = fit;
	Individu best(individu);
	for(int i=0;i<2*nbiter;i++){
		//cout<<"     "<<fit<<endl;
		double p = ((float)rand())/((float)RAND_MAX);
		if(p<pV1){//Voisinage V1 : changement de machine
			int indiceChgt = rand()%n;
			int indiceAncMach = individu.jobs[indiceChgt];
			bool voisinage1=false;
			bool restreint=false;
			int indRestreint=0;
			/**
			//Le bout de code qui suit calcule la machine qui accueille le masque juste avant
			int maskA = phi[(*individu)[1][indiceChgt]-1];
			int numJob = (*individu)[1][indiceChgt];
			int machineAvant=0;
			int curseur=indiceChgt-1;
			while(curseur>=0 && phi[(*individu)[1][curseur]-1]!=maskA)
				curseur--;
			if(curseur==-1)
				machineAvant=initMask[maskA-1];
			else
				machineAvant=(*individu)[0][curseur];
			//Le bout de code qui suit calcule la machine qui accueille le masque juste apres			
			int machineApres=0;
			curseur=indiceChgt+1;
			while(curseur<=n-1 && phi[(*individu)[1][curseur]-1]!=maskA)
				curseur++;
			if(curseur<=n-1)
				machineApres=(*individu)[0][curseur];
			//Ici on va determiner les cas ou on peut eventuellement operer le changement
			if(machineApres==0){
				if(machineAvant!=indiceAncMach)			
					voisinage1=true;
			}
			else{
				if(machineApres!=indiceAncMach && machineAvant!=indiceAncMach)
					voisinage1=true;
				else{
					if(machineApres==indiceAncMach && estQualifiee(numJob,machineAvant)){
						voisinage1=true;
						restreint=true;
						indRestreint=machineAvant;
					}
					if(machineAvant==indiceAncMach && estQualifiee(numJob,machineApres)){
						voisinage1=true;
						restreint=true;
						indRestreint=machineApres;
					}							
				}	
			}
			**/
			voisinage1=true;
			if(voisinage1){
				int indNouvMach=rand()%(corresp[individu.machines[indiceChgt]-1].size()); //La nouvelle machine doit être qualifiée
				if(restreint)
					individu.jobs[indiceChgt]=indRestreint;
				else
					individu.jobs[indiceChgt]=corresp[individu.machines[indiceChgt]-1][indNouvMach];
				double fitact=fitnessMemetique(individu,4,r1,r2,0,l1,l2,0,0);
				delta = fitact-fit;
				double pr = (double)(rand()/((double)(RAND_MAX)));
				if((delta>0) && (pr>exp(-delta/T)) ){//On n'améliore pas
					individu.jobs[indiceChgt]=indiceAncMach;			
				}
				else{ //On améliore
					fit=fitact;				
				}
			}
		}
		else{//Voisinage V2
			int indiceBoug=rand()%n; //L'indice du couple (job,machine) a deplacer
			int indDep=rand()%(n+1);
			/**
			int indice1=indiceBoug-1;
			int indice2=indiceBoug+1;
			int mask1=phi[(*individu)[1][indiceBoug]-1]; //On garde le masque concerne
			int machineAct=(*individu)[0][indiceBoug]; //On garde la machine
			//On ne va autoriser le deplacement que dans un rayon bien defini : pas au dela des utilisations par d'autres machines 				//du meme masque			
			while(indice1>=0 && (phi[(*individu)[1][indice1]-1]!=mask1 || (*individu)[0][indice1]==machineAct))
				indice1--;
			while(indice2<=n-1 && (phi[(*individu)[1][indice2]-1]!=mask1 || (*individu)[0][indice2]==machineAct))
				indice2++;
			indice1++;
			indDep=indice1 + rand()%(indice2-indice1);//On peut aller a n+1 pour le mettre en derniere position
			**/
			//On a l'indice, on effectue maintenant le deplacement
			individu.jobs.insert(individu.jobs.begin()+indDep,individu.jobs[indiceBoug]);
			individu.machines.insert(individu.machines.begin()+indDep,individu.machines[indiceBoug]);
			if(indDep<=indiceBoug)
				indiceBoug++;
			else
				indDep--;
			individu.jobs.erase(individu.jobs.begin()+indiceBoug);
			individu.machines.erase(individu.machines.begin()+indiceBoug);

			//On teste si le changement n'est pas ameliorant
			double fitact=fitnessMemetique(individu,4,r1,r2,0,l1,l2,0,0);
			delta = fitact-fit;
			double pr = (double)(rand()/((double)(RAND_MAX)));
			if((delta>0) && (pr>exp(-delta/T)) ){//On n'ameliore pas. On revient a la solution d'origine
				individu.jobs.insert(individu.jobs.begin()+indiceBoug,individu.jobs[indDep]);
				individu.machines.insert(individu.machines.begin()+indiceBoug,individu.machines[indDep]);
				if(indDep>indiceBoug)
					indDep++;
				individu.jobs.erase(individu.jobs.begin()+indDep);
				individu.machines.erase(individu.machines.begin()+indDep);
			}
			else //Si on ameliore
				fit=fitact;
		}
		if(i%palier==0) 
			T=T*a; //actualisation de la température
		//double total=fitnessMemetique(individu,3,0,0,0,0,0,0);
		//cout<<"Cout de la solution actuelle : "<<total<<endl;
		cout<<"Cout de la solution actuelle : "<<fit<<endl;
		//On actualise la meilleure solution trouvée dans la recherche locale			
		if(fit<cbest){
			cbest=fit;
			for(unsigned int k=0;k<individu.jobs.size();k++){
				best.jobs[k]=individu.jobs[k];
				best.machines[k]=individu.machines[k];
			}
		}
	}
	for(unsigned int k=0;k<individu.jobs.size();k++){
		individu.jobs[k]=best.jobs[k];
		individu.machines[k]=best.machines[k];
	}

	return cbest;
}


