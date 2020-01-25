#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<ctime>
#include<fstream>
#include<sstream>
#include"Instance.h"
#include"POP.h"
#include"scip_exception.hpp"


using namespace std;

int main(int argc, char **argv){
	if(argc!=2){
		cout<<"usage : ./exec instancePath"<<endl;
		return(EXIT_FAILURE); 
	}
 
	string path_file(argv[1]);
 	string instanceName;
	instanceName.replace(instanceName.begin(), instanceName.end(), path_file.substr(path_file.find_last_of("/")+1));
	srand(time(NULL));
	
	cout<<endl<<endl<<"Scheduling Optimizer: "<<endl<<endl; 
  	cout<<instanceName<<endl;

 	cout<<"Loading files... "<<endl;
 	Instance instance(path_file);

	string askObjectiveFunction("Objective function (1 for the number of wafers before H, 2 for the weighted sum of completion times,\n3 for the number of mask moves):");

 	bool terminated = false;
 	while(!terminated){
		cout<<endl<<endl<<"Which method do you want to use to solve the problem?"<<endl<<endl;
		cout<<"b: Bicriteria PLMO Chebychev"<<endl;
		cout<<"a: Non-exact method"<<endl;
		cout<<"e: Exact method"<<endl;
		cout<<"g: Approximation algorithm for mask move minimization"<<endl;
		cout<<"t: Multicriteria non-exact method"<<endl;
		cout<<"p: Bicriteria non-exact method"<<endl;
		cout<<"Other: Quit."<<endl<<endl<<endl;

		char c;
		cin>>c;
	
		if(c == 'a'){
			cout<<"Please enter: "<<endl<<endl;
			cout<<askObjectiveFunction;
			int fctObj;
			cin>>fctObj;
			if(fctObj == 1){
				cout<<"For which value of horizon H?"<<endl;
				int hor;
				cin>>hor;
				instance.modifierHorizon(hor);
			}
			cout<<endl<<"The size of the memetic algorithm population"<<endl;
			int Npop;
			cin>>Npop;

			int n(instance.jobsCount());
			//Settings
			//Paramètres du génétique
			int K(Npop/3);
			double Rep((n>70)?0:0.4);
			double Ts(1.0);
			int ReinIter(300);
			// Simulated annealing parameters
			int nbIterSA(1000);
			int nbIter((n>400)?10000000:((n>70)?25000000+(200-n)*100000:20000+(70-n)*1000));
			int p1(40); //p1/100 = Probability of neighbourhood V1
			int p2(50); //p2/100 = Probability of neighbourhood V2
			double T((n>70)?15000:2000), a((n>70)?0.99:(max(0.7,0.2+n/10))); // Temperature and decreasing coefficient
			int step(1); // Step to update T
			double resetRatio((n>70)?0.9:0.8); // Ration of temperature decrease when it is reset (compared to its initial value)
			int conv((n>70)?5000:1000); // Number of iterations without improvement before temperature reset
			int stagnation((n>70)?600000:4000); // Stagnation of best solution before reset

			//////////////////////////////////////////////////////////////////////////////////////////////////
	
			srand(time(NULL));
			cout<<"Thanks!"<<endl<<endl;
			cout<<"Start solving of instance "<<instanceName<<"."<<endl;
			cout<<"Objective function is ";
			if(fctObj == 1) cout<<"the number of processed wafers before horizon H, to maximize";
			else if(fctObj == 2) cout<<"the weighted sum of completion times, to minimize";
			else cout<<"the number of mask moves, to minimize";
			cout<<"."<<endl;
			cout<<"Population size is "<<Npop<<"."<<endl;

			clock_t ti,tf;
			ti = clock();

			vector< vector<vector<int> > > ordo=instance.POPMemetique(Npop,
									    nbIter,
									    conv,
									    stagnation,
									    resetRatio,
									    fctObj,
									    K,
									    Rep,
									    Ts,
									    ReinIter,
									    nbIterSA,
									    p1,
									    p2,
									    T,
									    step,
									    a,0,0,0,0,0,0,120,0);
			
			tf=clock();
			instance.evaluateSchedule(ordo,fctObj,0,0,0,0,0,0,0);
			double cc = instance.getCost(fctObj);
	
			instance.printScheduleToConsole(ordo);
			instance.visualization(ordo,fctObj,2);

			cout<<endl<<endl;

			cout<<"Solution cost: "<<cc<<endl<<endl;

			cout<<"Solution found in "<<(tf-ti)*1e-6<<" seconds."<<endl;
	
			if(fctObj == 1 && (instance.getMSol1() == -1 || cc > instance.getMSol1()))
				instance.writeBestKnownSolution(fctObj, cc, ordo);
			if(fctObj==2 && (instance.getMSol2() == -1 || cc < instance.getMSol2()))
					instance.writeBestKnownSolution(fctObj, cc, ordo);
			if(fctObj==3 && (instance.getMSol3() == -1 || cc < instance.getMSol3()))
				instance.writeBestKnownSolution(fctObj, cc, ordo);
		}
		else if(c=='e'){
			cout<<"Please enter: "<<endl<<endl;
			cout<<askObjectiveFunction;
			unsigned int fctObj, T;
			char reut;
			cin>>fctObj;
			if(fctObj == 2 || fctObj == 1){
				cout<<"Value of T such that there is an optimal solution with makespan lower than T (for binary variables indices)"<<endl;
				cin>>T;
			}
			if(fctObj == 2){
				cout<<"Do you want to use the best known sol? [o/N]"<<endl;
				cin>>reut;
			}
			if(fctObj == 1){
				cout<<"For which H?"<<endl;
				int hor;
				cin>>hor;
				instance.modifierHorizon(hor);
			}
			if(fctObj == 3){
				T=0;
				instance.ecrirePLNE(fctObj,T);
				ostringstream command; 
 				command.str("");
				command<<"scip -c \"read ./LPFiles/"<<instanceName<<".lp" << " opt write solution "<<"./SolFiles/"<<instanceName<<".sol"<<" quit \""<<endl;
			  	cout<<command.str().c_str();
	  			if(system(command.str().c_str()))
	      				cout<<"End of solver call"<<endl;
				// Reading the .sol file and creating a schedule
				vector<vector<vector<int > > > schedule = instance.makeScheduleFromILP(fctObj);
				instance.printScheduleToConsole(schedule);
				instance.evaluateSchedule(schedule,fctObj,0,0,0,0,0,0,0);
				double cc = instance.getCost(fctObj);
				cout<<endl<<endl<<"Solution cost: "<<cc<<endl<<endl<<endl;
				instance.visualization(schedule, fctObj, 1);
				if(instance.getMSol3() == -1 || cc < instance.getMSol3())
					instance.writeBestKnownSolution(fctObj, cc, schedule);		
			}
			else{ //fctObj = 1 or 2
				POP Instance(path_file, T, fctObj, &instance, reut);
				Instance.initialize_scip();
				Instance.create_ILP();
				cout<<"printed"<<endl;
				Instance.solve();
				if(Instance.write_best_solution()){
					cout<<"Solution OK and saved."<<endl;
					cout<<"End of solver call."<<endl;
					// Reading the .sol file and creating a schedule
					vector<vector<vector<int > > > schedule = instance.makeScheduleFromILP(fctObj);
					instance.printScheduleToConsole(schedule);
					instance.evaluateSchedule(schedule,fctObj,0,0,0,0,0,0,0);
					double cc = instance.getCost(fctObj);
					cout<<endl<<endl<<"Solution cost: "<<cc<<endl<<endl<<endl;
					instance.visualization(schedule, fctObj, 1);
					if(fctObj == 1 && (instance.getMSol1() == -1 || cc > instance.getMSol1()))
						instance.writeBestKnownSolution(fctObj, cc, schedule);
					if(fctObj == 2 && (instance.getMSol2() == -1 || cc < instance.getMSol2()))
						instance.writeBestKnownSolution(fctObj, cc, schedule);
				}
				else
    					cout<<"OUCH, solution not valid. Increase the value of T..."<<endl;
			}
		}
		else if(c == 't'){
			cout<<"Please enter : "<<endl<<endl;
			cout<<"Horizon H: "<<endl;
			int hor;
			cin>>hor;
			instance.modifierHorizon(hor);

	//Point de reference
	double r1,r2,r3;
	r1=(instance.getBound1()==-1)?(instance.nombrePlaquettes()):(instance.getBound1());
	r2=(instance.getBound2()==-1)?(instance.minorantCritere2()):(instance.getBound2());
	r3=instance.getMSol3(); 
	
	//On calcule le point ideal
	double nadir1(instance.calculAntiIdeal3D(1));
	double nadir2(instance.calculAntiIdeal3D(2));
	double nadir3(instance.calculAntiIdeal3D(3));

	//On evite certains cas
	if(nadir1==r1)
		nadir1=nadir1-1;
	if(nadir2==r2)
		nadir2=nadir2+1;
	if(nadir3==r3)
		nadir3=nadir3+1;
	

	//Ponderations
	double l1(1.0),l2((r1-nadir1)/(nadir2-r2)),l3((r1-nadir1)/(nadir3-r3));
	double coeff1(0.0), coeff2(0.1);

	cout<<endl<<"La taille de la population de l'algorithme memetique : "<<endl;
	int Npop;
	cin>>Npop;

	cout<<endl<<"Duree limite de resolution (en secondes) : "<<endl;
	unsigned int duree;
	cin>>duree;

	int n(instance.jobsCount());
	//Paramétrages /////////////////////////////////////////////////////////////////////////////////////
	
	//Paramètres du génétique
	int Ktour(Npop/3);
	double Rep((n>70)?0:0.4);
	double Ts(1.0);
	int ReinIter(300);
	int nbIterSA(1000);
	//Paramètres du recuit simulé
	int nbIter((n>400)?10000000:((n>70)?25000000+(200-n)*100000:20000+(70-n)*1000));
	int p1(40); //p1/100 = Probabilité du voisinage V1
	int p2(50); //p2/100 = Probabilité du voisinage V2
	double T((n>70)?15000:2000), a((n>70)?0.99:(max(0.7,0.2+n/10))); //température et coeff de décroissance
	int palier(1); //palier pour actualiser T
	double tauxRedemarrage((n>70)?0.9:0.8); //Taux de diminution de la température quand on la redémarre (par rapport à la valeur initiale)
	int conv((n>70)?5000:1000); //Largeur de la convergence, en nombre d'itérations, avant de redémarrer la température
	int stagnation((n>70)?600000:4000); //Stagnation du best avant de diminuer le coefficient de redémarrage

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////

	ostringstream outputParetoTricritere;
	outputParetoTricritere<<path_file<<"/outputParetoTricritere.dat";
	ofstream srt(outputParetoTricritere.str());
	double epsilon(0.001);

	while(coeff2<6.0){
		coeff1=0.1;
		while(coeff1<6.0){
			double ll2(coeff1*l2), ll3(coeff2*l3);
			//Affichage des parametres de la fonction d'agregation
			cout<<endl<<endl<<"Point de reference choisi : ("<<r1<<" , "<<r2<<" , "<<r3<<")."<<endl;
			cout<<"Approximation du Point anti-ideal : ("<<nadir1<<" , "<<nadir2<<" , "<<nadir3<<")."<<endl;
			cout<<"Ponderations : ("<<l1<<" , "<<ll2<<" , "<<ll3<<")."<<endl<<endl;

			clock_t ti,tf;
			ti = clock();

			vector< vector<vector<int> > > ordo = instance.POPMemetique(Npop, nbIter, conv, stagnation, tauxRedemarrage, 4, Ktour, Rep, Ts, ReinIter, nbIterSA, p1, p2, T, palier, a, r1, r2, r3, l1, ll2, ll3, duree, epsilon);
	
			tf=clock();
	
			double cc(instance.evaluateSchedule(ordo,4, r1, r2, r3, l1, ll2, ll3, epsilon));
			double c1(instance.evaluateSchedule(ordo,1, 0, 0, 0, 0, 0, 0, epsilon));
			double c2(instance.evaluateSchedule(ordo,2, 0, 0, 0, 0, 0, 0, epsilon));
			double c3(instance.evaluateSchedule(ordo,3, 0, 0, 0, 0, 0, 0, epsilon));

			cout<<endl<<endl;
	
			cout<<"Le cout de la solution est de "<<cc<<" ( "<<c1<<" , "<<c2<<" , "<<c3<<" ) "<<endl<<endl;

			cout<<"Solution trouvee en "<<(tf-ti)*1e-6<<" secondes."<<endl;
	
			if((instance.getMSol1() == -1) || (c1>instance.getMSol1()))
				instance.writeBestKnownSolution(1,c1,ordo);
			if((instance.getMSol2() == -1) || (c2<instance.getMSol2()))
				instance.writeBestKnownSolution(2,c2,ordo);
			if((instance.getMSol3() == -1) || (c3<instance.getMSol3()))
				instance.writeBestKnownSolution(3,c3,ordo);
			if(coeff1<1.0)
				coeff1+=0.2;
			else
				coeff1+=1.0;
			srt<<c1<<"\t"<<c2<<"\t"<<c3<<endl;
		}
		srt<<endl;
		if(coeff2<1.0)
			coeff2+=0.2;
		else
			coeff2+=1.0;
	}
	srt.close();
}
else if(c=='g'){
	vector<vector<vector<int> > > ordo = instance.nombreDeplacementsMasques();
	instance.printScheduleToConsole(ordo);
	instance.evaluateSchedule(ordo,3,0,0,0,0,0,0,0);
	double cc = instance.getCost(3);
	cout<<endl<<endl<<"Le nombre de deplacements de masque est de "<<cc<<endl<<endl<<endl;
	instance.visualization(ordo,3,3);
	if((instance.getMSol3() == -1) || (cc < instance.getMSol3()))
		instance.writeBestKnownSolution(3, cc, ordo);
}
	else if(c=='b'){
		cout<<"Please enter: "<<endl<<endl;
		cout<<"The value T for time indices of the model variables:"<<endl;
		unsigned int T;
		cin>>T;
		cout<<"For which value of horizon H?"<<endl;
		int hor;
		cin>>hor;
		instance.modifierHorizon(hor);
	
		// Reference points computed with ideal values on criteria 1 and 2
		int C((instance.getBound1()==-1)?instance.nombrePlaquettes():instance.getBound1());
		double r1(C);
		double r2((instance.getBound2()==-1)?instance.minorantCritere2():instance.getBound2());
	
		double nadir1(instance.calculAntiIdeal(1,2));	
		double nadir2(instance.calculAntiIdeal(2,1));
	
		// Avoid some cases
		if(nadir1==r1)
			nadir1=nadir1-1;
		if(nadir2==r2)
			nadir2=nadir2+1;

		// Ref point can either be manually computed or approached
		double l1(1.0);
		double l2((r1-nadir1)/(nadir2-r2));
	
		ostringstream sortie;
		sortie<<path_file<<"/outputPareto.par";
		ofstream outputPareto(sortie.str());
		// Compensation factor epsilon
		double epsilon(0.00);

		double coefficientL2(0.0);
		while(coefficientL2<100.0){
			double ll2(coefficientL2*l2);
			cout<<endl<<endl<<"Reference point: ("<<r1<<" , "<<r2<<")."<<endl;
			cout<<"Approximation of anti-ideal point: ("<<nadir1<<" , "<<nadir2<<")."<<endl;
			cout<<"Weights: ("<<l1<<" , "<<ll2<<")."<<endl<<endl;

			POP Instance(path_file,T,4,&instance,l1,ll2,r1,r2,epsilon);
			Instance.initialize_scip();
			Instance.create_ILPMO();
			cout<<"printed"<<endl;
			Instance.solve();
			if(Instance.write_best_solution())
				cout<<"Solution OK and saved"<<endl;
  			else
    				cout<<"OUCH, solution not valid"<<endl;
			cout<<"End of call to solver"<<endl;
			cout<<endl<<endl<<"Reference point: ("<<r1<<" , "<<r2<<")"<<endl<<endl;
			// Function for reading .sol file and build a schedule
			vector<vector<vector<int > > > ordo = instance.makeScheduleFromILP(4);
			instance.printScheduleToConsole(ordo);
			instance.evaluateSchedule(ordo,1,0,0,0,0,0,0,epsilon);
			instance.evaluateSchedule(ordo,2,0,0,0,0,0,0,epsilon);
			double cc1 = instance.getCost(1);
			double cc2 = instance.getCost(2);
			cout<<endl<<endl<<"Found solution has value, for each criterion: ("<<cc1<<" , "<<cc2<<")"<<endl<<endl<<endl;
			outputPareto<<cc1<<" "<<cc2<<endl;
			if(coefficientL2<1.0)
				coefficientL2+=0.02;
			else
				coefficientL2+=2.0;
			if((instance.getMSol1() == -1) || (cc1>instance.getMSol1()))
				instance.writeBestKnownSolution(1,cc1,ordo);
			if((instance.getMSol2() == -1) || (cc2<instance.getMSol2()))
				instance.writeBestKnownSolution(2,cc2,ordo);
		}
		outputPareto.close();
	}
	else if(c=='p'){
		cout<<"Pour quelle valeur de l'horizon H ?"<<endl;
		int hor;
		cin>>hor;
		instance.modifierHorizon(hor);
	
		cout<<endl<<"Taille de la population de l'algorithme memetique : "<<endl;
		int Npop;
		cin>>Npop;

		cout<<endl<<"Duree limite de resolution (en secondes) : "<<endl;
		unsigned int duree;
		cin>>duree;

		int n(instance.jobsCount());
		//Settings
		int Ktour(Npop/3);
		double Rep((n>70)?0:0.4);
		double Ts(1.0);
		int ReinIter(300);
		int nbIterSA(1000);
		int nbIter((n>400)?10000000:((n>70)?25000000+(200-n)*100000:1000));
		int p1(40); //p1/100 = Probabilité du voisinage V1
		int p2(50); //p2/100 = Probabilité du voisinage V2
		double T((n>70)?15000:100), a((n>70)?0.99:0.9); //température et coeff de décroissance
		int palier(1); //palier pour actualiser T
		double tauxRedemarrage((n>70)?0.9:0.8); //Taux de diminution de la température quand on la redémarre (par rapport à la valeur initiale)
		int conv((n>70)?5000:1000); //Largeur de la convergence, en nombre d'itérations, avant de redémarrer la température
		int stagnation((n>70)?600000:4000); //Stagnation du best avant de diminuer le coefficient de redémarrage

		//////////////////////////////////////////////////////////////////////////////////////////////////
	
		//on determine le point de reference en recuperant les meilleures bornes trouvees sur les criteres 1 et 2
		int C((instance.getBound1()==-1)?instance.nombrePlaquettes():instance.getBound1());
		double r1(C);
		double r2((instance.getBound2()==-1)?instance.minorantCritere2():instance.getBound2());
		double nadir1(instance.calculAntiIdeal(1,2));	
		double nadir2(instance.calculAntiIdeal(2,1));
	
		// Avoid some cases
		if(nadir1==r1)
			nadir1=nadir1-1;
		if(nadir2==r2)
			nadir2=nadir2+1;

		// Reference point can be computed thanks to an approximation	
		double l1(1.0);
		double l2((r1-nadir1)/(nadir2-r2));

		ostringstream sortie;
		sortie<<path_file<<"/outputParetoBicritere.par";
		ofstream outputPareto(sortie.str());
		// Compensation factor
		double epsilon(0.001);
		double coefficientL2(0.0);
		while(coefficientL2<100.0){
			double ll2(coefficientL2*l2);
			cout<<endl<<endl<<"Point de reference choisi : ("<<r1<<" , "<<r2<<")."<<endl;
			cout<<"Approximation du Point anti-ideal : ("<<nadir1<<" , "<<nadir2<<")."<<endl;
			cout<<"Ponderations : ("<<l1<<" , "<<ll2<<")."<<endl<<endl;
			vector< vector<vector<int> > > ordo=instance.POPMemetique(Npop, nbIter, conv, stagnation, tauxRedemarrage, 4, Ktour, Rep, Ts, ReinIter, nbIterSA, p1, p2, T, palier, a, r1, r2, 0, l1, ll2, 0, duree, epsilon);	
			cout<<endl<<endl<<"Point de reference : ("<<r1<<" , "<<r2<<")"<<endl<<endl;
			instance.printScheduleToConsole(ordo);
			instance.evaluateSchedule(ordo,1,0,0,0,0,0,0,epsilon);
			instance.evaluateSchedule(ordo,2,0,0,0,0,0,0,epsilon);
			double cc1=instance.getCost(1);
			double cc2=instance.getCost(2);
			cout<<endl<<endl<<"La solution trouvee a pour cout sur chacun des criteres : ("<<cc1<<" , "<<cc2<<")"<<endl<<endl<<endl;
			outputPareto<<cc1<<" "<<cc2<<endl;
			if(coefficientL2<1.0)
				coefficientL2+=0.02;
			else
				coefficientL2+=2.0;
			if((instance.getMSol1() == -1) || (cc1>instance.getMSol1()))
				instance.writeBestKnownSolution(1,cc1,ordo);
			if((instance.getMSol2() == -1) || (cc2<instance.getMSol2()))
				instance.writeBestKnownSolution(2,cc2,ordo);
		}
		outputPareto.close();
	}
	else
		terminated = true;
	}
	return(EXIT_SUCCESS);
}



