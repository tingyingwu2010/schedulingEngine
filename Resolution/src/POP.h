#ifndef ___POP_H_________
#define ___POP_H_________

#include"Instance.h"
#include<scip/scip.h>
#include<scip/scipdefplugins.h>
#include"objscip/objscip.h"
#include"scip/scip.h" 
#include"scip/cons_linear.h"
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<sys/stat.h>
#include<vector>
#include<list>

using namespace std;

class POP{

	public:
		
		SCIP *scip; // The SCIP handler instance
  		Instance *K; // The instance
		int T; //The T in the ILP
		int primalHeuristic; // Tells whether we must use primal heuristic (1) or not (0)
		char reut; // 'o' if we re-use the best known solution
		int obj; // Objective function code: 1 for the wafers count, 2 for the weighted sum of completion times, and 4 for the Chebychev bicriteria version
		double l1; // Reference point for criterion 1 in the multi-objective problem
		double l2; // Reference point for criterion 2 in the multi-objective problem
		double r1; // Nadir point for criterion 1 in the multi-objective problem
		double r2; // Nadir point for criterion 2 in the multi-objective problem
		double epsilon; // Epsilon for the augmented Chebychev multi-objective function

		// Variables
  		vector< vector< vector<SCIP_VAR *> > > v_varu;  // u_ijt variables with time index 
		vector<SCIP_VAR *> v_varz; // For multi-objective version

	  	// Constructors
		POP(string path_file, int T, int fctObj, Instance *pbm, char reutt); // Mono-criterion version
		POP(string path_file, int T, int fctObj, Instance *pbm, double l1, double l2, double r1, double r2, double eps); // Bi-criteria version

		// Other methods
  		void initialize_scip(); // Initializes the model 
		void create_ILP(); // Creates variables and constraints
		void create_ILPMO(); // Creates variables and constraints for the multi-objective problem
		void initialMaskConstraint(SCIP_CONS *Cte); // Constraint for initial location of auxiliary resources
		void addBestBoundConstraint(SCIP_CONS * Cte); // Additional constraints if there is an already known solution
		void movingAuxiliaryResourcesConstraints(SCIP_CONS *Cte); // Additional constraints on auxiliary resources
		void bigMConstraints(SCIP_CONS *Cte, std::vector< std::vector<int> > invCorresp); // Big-M Constraints for non-simultaneous processes on machines constraints
		void generalizedValidInequalities(SCIP_CONS *Cte, std::vector< std::vector<int> > invCorresp); // Constraints C(u,l,j,k) generalized valid inequalities
		void auxiliaryResourcesValidInequalities(SCIP_CONS *Cte); // Valid inequalities related to auxiliary resources transport times
		void solve(); // Solving method
		int write_best_solution(); // Writing the optimal solution by testing validity
		SCIP_RETCODE write_solution(SCIP_SOL *sol); // Writing the optimal solution
};

#endif
