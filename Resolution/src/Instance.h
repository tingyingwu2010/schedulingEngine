#ifndef ______PBMPHOTO______
#define ______PBMPHOTO______

#include<cstdio>
#include<vector>
#include<string>
#include"Individu.h"

#define MAXPLAQUETTE 25 

class Instance{
	//Attributes of the class
 private :
	int n; // Jobs count
	int m; // Machines count
	int l; // Masks count
	int nB; // Setup family types count
	int f; // Setup families count
	int H; // Horizon
	std::vector<int> phi; // Required mask for each job
	std::vector<int> initMask; // Initial mask position
	std::vector<double> w; // Number of wafers for each job
	std::vector<double> c; // Priority of each job
	std::vector<int> families; // Setup family type for each machine
	std::vector< std::vector<int> > corresp; // Qualified machines for each job
	std::vector< std::vector<int> > exec; // Processing time for each couple job/machine for qualified machines
	std::vector< std::vector<int> > batch; // Setup family (batch-config) for each couple job/machine for qualified machines
	std::vector< std::vector< std::vector<int> > > setup; // Setup time for each machine and couple of setup family
	std::string instanceName; // Name of instance file name
	std::vector<std::vector< float > > colors; // Masks colors for visualization
	double obj1; // Solution cost for the number of processed wafers before H
	double obj2; // Solution cost for the weighted sum of completion times
	double obj3; // Solution cost for the number of mask moves
	double bestKnownSol1; // Best known value for objective 1
	double bestKnownSol2; // Best known value for objective 2
	double bestKnownSol3; // Best known value for objective 3
	double bestBound1; // Best upper bound for objective 1
	double bestBound2; // Best lower bound for objective 2

	// Methods of the class
 public :
	
	Instance(std::string name); // Constructor that reads file name and creates the instance
	
	// Utils
	bool estQualifiee(int i, int j); //Indique si la machine j est qualifiee pour le lot i
	bool estTrainQualifiee(std::vector<int> train, int j); //Indique si la machine est qualifiee pour le train
	unsigned int indexMachineQualif(int job, int machine); // Returns index of qualification of machine j for job i
	bool estBatchCompatible(std::vector<int> train , int batchMachine, int j); //Indique si une machine est batch-config-convenable pour 											   //tous les lots d'un train
	double nombreHarmonique(int a); //Calcule 1 + 1/2 + 1/3 + ... + 1/a
	int member(std::vector<int> tab, int j); // Return index of item j in tab, -1 if it does not belong to it	
	int nombrePlaquettes(); //Retourne le nombre total de plaquettes du probleme
	void writeBestKnownSolution(int fctObj, double cc, std::vector< std::vector<std::vector<int> > > schedule); // Updates best known solution

	// Getters
	int nombreMachines();
	int jobsCount(); // Gets the number of jobs
	int nombreMasques();
	double getCost(int obj); // Gets the cost of a schedule given an objective function code
	std::vector< std::vector<int> > getCorresp();
	std::vector< std::vector<int> > getExec(); 
	std::vector< std::vector<int> > getBatch();
	std::vector< std::vector< std::vector<int> > > getSetup();
	std::vector<int> getPhi();
	std::vector<int> getInitMask();
	std::vector<double> getW();
	std::vector<double> getC();
	std::vector<int> getFamilies();
	int getH();
	std::string getInstanceName(); // Returns the instance name, as it appears on the folder containing the input
	double getMSol1();
	double getMSol2();
	double getMSol3();
	double getBound1(); // Gets best upper bound for the first criterion
	double getBound2(); // Gets best lower bound for the second criterion

	// Setters
	void modifierHorizon(int nH); //Modifie l'horizon de temps H en lui donnant la valeur nH
	void setBound(int obj, double val); // Updates the best bound on given objective 

	// Solver methods
	
	// Metaheuristics	
	std::vector< std::vector<std::vector<int> > > POPMemetique(int Npop, int nbIter, int conv, int stagnation, double tauxRedemarrage, int fctObj, int Ktour, double Rep, double Ts, int ReinIter, int nbIterSA, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, unsigned int duree, double epsilon); 
						//Algorithme memetique avec operations sur la sequence de lots affectes aux machines; 
						//parametres : Npop, taille de la population; fctObj, objectif de l'optimisation(1->moves, 2->ciCI, 
						//3->masques), Ktour est le nombre de parents selectionnes pour le tournoi de la selection, Rep est 
						//le taux de reproduction de la generation, soit le pourcentage de la population que represente la 
						//descendance, Ts est le taux de survivants de la nouvelle generation, ReinIter le nombre d'iterations sans ameliorations avant de reinitialiser la population et nbIterSA le nombre d'iterations sans ameliorations avant d'arreter, p1 et p2 les probabilites d'utilisation respectives des voisinages V1 et V2 dans la recherche locale, T est la temperature dans le recuit simule, palier est le nombre d'iterations sans diminuer la temperature T, a est le coefficient par lequel on multiplie T, les coordonnees r1, r2 et r3 sont celles du point du reference, et l1, l2 et l3 sont les coefficients. Enfin, on borne la durée de résolution en secondes et la valeur epsilon correspond au coefficient du terme de compensation dans la fonction d'agregation utilisée
	double fitnessMemetique(Individu individu, int fctObj, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); //Calcule l'évaluation (ou fitness) d'un individu dans l'algorithme mémétique
	Individu CroisementMemetique(Individu parent1, Individu parent2); //Croisement des parents parent1 et parent2
	double rechercheLocaleMemetique(Individu &enfant, int nbiter, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon);//Lance un algorithme de recherche locale sur un individu/enfant
	std::vector<std::vector<std::vector<int> > > HeuristiqueMemetique(Individu individu); //Complète un individu du mémétique en une solution
	double nouvelleFitness(std::vector<std::vector<unsigned int> > &enfant, std::vector<unsigned int> &completionTimes, int indice, double fitness); //Calcule la nouvelle fitness dans la recherche locale, en économisant quelques calculs
	double evaluateSchedule(std::vector< std::vector< std::vector<int> > > schedule, int fct, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); // Given a criterion, evaluate the schedule
	double recuitSimule(Individu &i1, int nbIter, int conv, int stagnation, double tauxRedemarrage, int p1, int p2, double T, int palier, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); //Recuit simulé

	// Linear Programming
	int ecrirePLNE(unsigned int fctObj, unsigned int T); //Ecrit le PLNE dans le dossier FichiersPL
	std::vector<std::vector<std::vector<int > > > makeScheduleFromILP(int fctObj); // Transfer obtained solution form LP to the internal structure

	// Approximation algorithm
	std::vector< std::vector< std::vector<int> > > nombreDeplacementsMasques(); // Approximation algorithm to minimize the number of mask moves
	std::vector< std::vector<int> > couvertureMin(int i, std::vector<int> JJ,std::vector< std::vector<int> > MM);//Resoud un probleme de 														      //couverture minimale par des 														      //ensembles par l'algorithme 														      //glouton de Lovasz
	std::vector<std::vector<std::vector<int> > > completerOrdonnancement(std::vector<std::vector<std::vector<int> > > solpartielle); //Complete une solution (x,y) en une solution (x,y,t) de facon gloutonne 
	
	// Display methods
	void printVector(std::vector<unsigned int> v); // Prints the array
	void afficheTemps(); //Affiche les temps d'execution sur la console (en lignes les jobs et en colonnes les machines)
	void printScheduleToConsole(std::vector<std::vector<std::vector<int> > > schedule); // Prints the solution to the console
	int visualization(std::vector<std::vector<std::vector<int> > > schedule, int fctObj, int meth); // Builds a .pdf Gantt chart, displaying the solution returned by the optimization method. Here, fctObj equals 1 for the number of wafers, 2 for the weighted sum of completion times, 3 for the number of mask moves and 4 for the Chebychev aggregation function of the 3 previous objective functions. The meth parameter equals 1 for the ILP, 2 for the metaheuristic, 3 for the approximation algorithm, and 4 for the multi-objective ILP

	//PLMO pour la norme de Tchebycheff (version bicritere)
	int ecrireMIP_BCT(double l1, double l2, double r1, double r2, unsigned int T); //(l1,l2) est le jeu de poids, (r1,r2) le point de reference et T une valeur entière telle qu'il existe une solution optimale dont le makespan lui est inferieur
	double minorantCritere2();//Donne un minorant de la somme ponderee des dates de fin d'execution des jobs pour le point de reference
	double majorantCritere2();//Donne un majorant de la somme ponderee des dates de fin d'execution des jobs pour le poids du critere 2
	double calculAntiIdeal(int obj, int aouvrir);//Retourne une approximation de la coordonnee obj du point Nadir en bicritere
	double calculAntiIdeal3D(int obj);//Retourne une approximation de la coordonnee obj du point Nadir en tricritere

	//Methode tricritere pour la garantie de faible Pareto-optimalite
	std::vector<std::vector<std::vector<std::vector<int> > > > tricritere(double r1, double r2, double l1, double l2, int ecart, int nbiter, int p1, double T, int palier, double a); //(l1,l2) est le jeu de poids, (r1,r2) le point de reference, ecart est l'ecart maximal par rapport a l'optimum pour le troisieme critere, les autres parametres sont ceux de la recherche locale recuit simule
	Individu creerIndividuInitial(std::vector<std::vector<std::vector<int > > > ordo); //Crée un individu à partir de ordo. Correspond à la fonction réciproque du codage mémétique qui à un individu associe une solution
	double rechercheLocaleTricritere(Individu individu, int nbiter, int p1, double T, int palier, double a, double r1, double r2, double l1, double l2); //Algorithme de recherche locale qui ne dégrade jamais le nombre de déplacements de masques

	int ecrirePLNEmasques(int borninf); //Crée un PLNE pour le critère 3 avec borne inf imposée
};


#endif
