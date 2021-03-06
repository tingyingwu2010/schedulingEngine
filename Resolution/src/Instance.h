#ifndef ______PBMPHOTO______
#define ______PBMPHOTO______

#include<cstdio>
#include<vector>
#include<string>
#include"Individual.h"

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
	bool isQualified(int job, int machine); // Tells whether given machine is qualified for given job
	unsigned int indexMachineQualif(int job, int machine); // Returns index of qualification of machine j for job i
	bool estBatchCompatible(std::vector<int> train , int batchMachine, int j); //Indique si une machine est batch-config-convenable pour 											   //tous les lots d'un train
	double harmonicSerie(int a); // Computes 1 + 1/2 + 1/3 + ... + 1/a
	int member(std::vector<int> tab, int j); // Return index of item j in tab, -1 if it does not belong to it	
	int numberWafers(); // Returns the total number of wafers
	void writeBestKnownSolution(int fctObj, double cc, std::vector< std::vector<std::vector<int> > > schedule); // Updates best known solution

	// Getters
	int machinesCount(); // Gets the number of machines
	int jobsCount(); // Gets the number of jobs
	int masksCount(); // Gets the number of masks
	double getCost(int obj); // Gets the cost of a schedule given an objective function code
	std::vector< std::vector<int> > getCorresp(); // Gets qualification of machines
	std::vector< std::vector<int> > getExec();  // Gets processing times
	std::vector< std::vector<int> > getBatch(); // Gets setup families
	std::vector< std::vector< std::vector<int> > > getSetup(); // Gets setup times
	std::vector<int> getPhi(); // Gets required masks per job
	std::vector<int> getInitMask(); // Gets initial position of each mask
	std::vector<double> getW(); // Gets number of wafers per job
	std::vector<double> getC(); // Gets priorities per job
	std::vector<int> getFamilies(); // Gets family per machine
	int getH(); // Returns horizon H
	std::string getInstanceName(); // Returns the instance name, as it appears on the folder containing the input
	double getBestSol1(); // Returns best solution on the number of wafers before horizon
	double getBestSol2(); // Returns best solution on the weighted sum of completion times
	double getBestSol3(); // Returns best solution on the number of mask moves
	double getBound1(); // Gets best upper bound for the first criterion
	double getBound2(); // Gets best lower bound for the second criterion

	// Setters
	void modifyHorizon(int nH); // Sets horizon H to a new value
	void setBound(int obj, double val); // Updates the best bound on given objective 

	// Solver methods
	
	// Metaheuristics	
	std::vector< std::vector<std::vector<int> > > POPMemetic(int Npop, int nbIter, int conv, int stagnation, double resetRatio, int fctObj, int Ktour, double Rep, double Ts, int ReinIter, int nbIterSA, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, unsigned int duration, double epsilon); //Memetic algorithm. Parametres : Npop is the population size; fctObj is the objective function, Ktour  is the number of parents for the selection tournament, Rep is the reproduction ratio, Ts is the offspring ratio, ReinIter is the number of iterations without improvements before resetting the population, nbIterSA is the number of iterations without improvement before stopping, p1 and p2 are the neighbourhood probabilities for V1 and V2, T is the simulated annealing temperature, step is the number of iterations without decreasing temperature T, a is the multiplication coefficient for T, coordinates r1, r2 and r3 are those of reference point, and l1, l2, l3 are the coefficients. Then, duration is the solving time limit and epsilon is the compensation coefficient in the aggregation function
	double fitnessMemetic(Individual individual, int fctObj, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); //Compute an individual fitness
	Individual crossoverMemetic(Individual parent1, Individual parent2); //Crossover parents parent1 and parent2
	double localSearchMemetic(Individual &child, int nbiter, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); // Local search algorithm
	std::vector<std::vector<std::vector<int> > > HeuristicMemetic(Individual individu); // Complete an individual to a solution
	double newFitness(std::vector<std::vector<unsigned int> > &child, std::vector<unsigned int> &completionTimes, int indice, double fitness); // Computes the new fitness in local search, saving some computations
	double evaluateSchedule(std::vector< std::vector< std::vector<int> > > schedule, int fct, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); // Given a criterion, evaluate the schedule
	double simulatedAnnealing(Individual &i1, int nbIter, int conv, int stagnation, double resetRatio, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon); // Simulated annealing

	// Linear Programming
	int writeILP(unsigned int fctObj, unsigned int T); // Writes ILP in LPFiles folder
	std::vector<std::vector<std::vector<int > > > makeScheduleFromILP(int fctObj); // Transfer obtained solution form LP to the internal structure

	// Approximation algorithm
	std::vector< std::vector< std::vector<int> > > numberMaskMoves(); // Approximation algorithm to minimize the number of mask moves
	std::vector< std::vector<int> > setCover(int i, std::vector<int> JJ,std::vector< std::vector<int> > MM); // Solves a Set Cover problem using Lovasz approximation algorithm
	std::vector<std::vector<std::vector<int> > > completeSchedule(std::vector<std::vector<std::vector<int> > > partialSolution); // Completes an (x,y) solution to an (x,y,t) solution 

	// Display methods
	void printVector(std::vector<unsigned int> v); // Prints the array
	void displayTime(); // Displays processing times on the console (in lines, the jobs, in columns, the machines)
	void printScheduleToConsole(std::vector<std::vector<std::vector<int> > > schedule); // Prints the solution to the console
	int visualization(std::vector<std::vector<std::vector<int> > > schedule, int fctObj, int meth); // Builds a .pdf Gantt chart, displaying the solution returned by the optimization method. Here, fctObj equals 1 for the number of wafers, 2 for the weighted sum of completion times, 3 for the number of mask moves and 4 for the Chebychev aggregation function of the 3 previous objective functions. The meth parameter equals 1 for the ILP, 2 for the metaheuristic, 3 for the approximation algorithm, and 4 for the multi-objective ILP

	// ILPMO for Chebychev norm
	int writeMIP_BCT(double l1, double l2, double r1, double r2, unsigned int T); // (l1,l2) is a weight couple, (r1,r2) the reference point and T an integer such that there exists an optimal solution with makespan lower than T
	double minorantCriterion2(); // Returns a lower bound of the weighted sum of completion times for the reference point
	double majorantCriterion2(); // Returns un upper bound of the weighted sum of completion times
	double computeAntiIdeal(int obj, int aouvrir); // Returns an approximation of coordinate obj of point Nadir in bicriteria
	double computeAntiIdeal3D(int obj); // Returns an approximation of coordinate obj of point Nadir in tricriteria

	// Tricriteria method with Pareto-optimality guarantee
	std::vector<std::vector<std::vector<std::vector<int> > > > tricritere(double r1, double r2, double l1, double l2, int ecart, int nbiter, int p1, double T, int step, double a); //(l1,l2) est le jeu de poids, (r1,r2) le point de reference, ecart est l'ecart maximal par rapport a l'optimum pour le troisieme critere, les autres parametres sont ceux de la recherche locale recuit simule
	Individual creerIndividuInitial(std::vector<std::vector<std::vector<int > > > ordo); //Crée un Individual à partir de ordo. Correspond à la fonction réciproque du codage mémétique qui à un Individual associe une solution
	double rechercheLocaleTricritere(Individual individual, int nbiter, int p1, double T, int step, double a, double r1, double r2, double l1, double l2); //Algorithme de recherche locale qui ne dégrade jamais le nombre de déplacements de masques
	int writeILPmasks(int lowerBound); // Creates an ILP for criterion 3 with lower bound
};


#endif
