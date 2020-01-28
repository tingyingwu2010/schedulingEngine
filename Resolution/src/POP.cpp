#include"POP.h"
#include"POPCte.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<ostream>
#include<sstream>
#include<math.h>
#include"scip_exception.hpp"

using namespace std;

POP::POP(string path_file, int TT, int fctObj, Instance *pbm, char reutt){
	T = TT;
	obj = fctObj;
	reut = reutt;
	primalHeuristic = 0;
	K=pbm;
	int i;
 	// Creating variables u_ijt
	v_varu.resize(K->jobsCount());
	for(i=0;i<K->jobsCount();i++){
		v_varu[i].resize(K->getCorresp()[i].size());
		for(unsigned int j=0;j<K->getCorresp()[i].size();j++)
			v_varu[i][j].resize(1+T-K->getExec()[i][j]);
	}
	v_varz.resize(1);
}

POP::POP(string path_file, int TT, int fctObj, Instance *pbm, double ll1, double ll2, double rr1, double rr2, double eps){
	l1 = ll1;
	l2 = ll2;
	r1 = rr1;
	r2 = rr2;
	T = TT;
	epsilon = eps;
	obj = fctObj;
	K=pbm;
 	int i;
 	// Creating variables u_ijt
	v_varu.resize(K->jobsCount());
	for(i=0;i<K->jobsCount();i++){
		v_varu[i].resize(K->getCorresp()[i].size());
		for(unsigned int j=0;j<K->getCorresp()[i].size();j++)
			v_varu[i][j].resize(1+T-K->getExec()[i][j]);
 	} 
	v_varz.resize(1);
}

void POP::initialize_scip(){
	/* initialize SCIP */
	SCIP_CALL_EXC( SCIPcreate(&scip) );

	/* include default SCIP plugins */
	SCIP_CALL_EXC( SCIPincludeDefaultPlugins(scip) );

        SCIP_CALL_EXC( SCIPsetIntParam(scip, "display/verblevel", 5) );

	 /*
	SCIP_CALL_EXC( SCIPincludeConshdlrLinear(scip) );
        SCIP_CALL_EXC( SCIPincludeConshdlrIntegral(scip) );
	SCIP_CALL_EXC( SCIPincludeNodeselBfs(scip) );
	SCIP_CALL_EXC( SCIPincludeDispDefault(scip) );


	SCIP_CALL_EXC( SCIPincludeHeurTrivial(scip) );
	SCIP_CALL_EXC( SCIPincludeHeurTrySol(scip) );

	SCIP_CALL_EXC( SCIPincludeBranchruleAllfullstrong(scip) );
	SCIP_CALL_EXC( SCIPincludeNodeselBfs(scip) );
	SCIP_CALL_EXC( SCIPincludeNodeselEstimate(scip) );
 	SCIP_CALL_EXC( SCIPincludeBranchruleFullstrong(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchruleInference(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchruleMostinf(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchruleLeastinf(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchrulePscost(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchruleRandom(scip) );
	SCIP_CALL_EXC( SCIPincludeBranchruleRelpscost(scip) );

	SCIP_CALL_EXC( SCIPincludePresolBoundshift(scip) );
	SCIP_CALL_EXC( SCIPincludePresolDualfix(scip) );
	SCIP_CALL_EXC( SCIPincludePresolImplics(scip) );
	SCIP_CALL_EXC( SCIPincludePresolInttobinary(scip) );
	SCIP_CALL_EXC( SCIPincludePresolProbing(scip) );
	SCIP_CALL_EXC( SCIPincludePresolTrivial(scip) );
	 */

	/* include CBCol constraint handlers */
       	SCIP_CALL_EXC( SCIPincludeObjConshdlr(scip, new POPCte(scip, this), TRUE) );

	/* include Primal Heuristic *****************************************************************************/
	//SCIP_CALL_EXC( SCIPincludeObjHeur(scip, new POPHeur(scip, this), TRUE) );

	/* include a Branching Rule *****************************************************************************/
	//SCIP_CALL_EXC( SCIPincludeObjBranchrule(scip, new AffBranchRule(scip,this), FALSE));

	// Create an empty problem
	SCIP_CALL_EXC( SCIPcreateProb(scip, "POP", NULL, NULL, NULL, NULL, NULL, NULL , NULL));
}

void POP::addBestBoundConstraint(SCIP_CONS * Cte){
	ostringstream namebuff;	
	//Here, add the best known bound if it exists
	if(obj == 1 && K->getBound1() != -1){
		double coeff;
		namebuff.str(""); namebuff <<"BestKnownBound";
       		SCIP_CALL_EXC(SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), K->getBound1(), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
		for(int ii=0;ii<K->jobsCount();ii++){
			for(unsigned int jj=0;jj<K->getCorresp()[ii].size();jj++){
				for(int tt=0;tt<=T-K->getExec()[ii][jj];tt++){
					coeff = (tt>K->getH())?0:( (tt<=K->getH()-K->getExec()[ii][jj])?K->getW()[ii]:K->getW()[ii]*(K->getH()-tt)/K->getExec()[ii][jj]);					
					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[ii][jj][tt], coeff));
		              	}
        		}
   		}
		SCIP_CALL_EXC(SCIPaddCons(scip, Cte));
		SCIP_CALL_EXC(SCIPreleaseCons(scip, &Cte));
	}

	if(obj == 2 && K->getBound2() != -1){
		double coeff;
		namebuff.str(""); namebuff <<"BestKnownBound";
       		SCIP_CALL_EXC(SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, K->getBound2(), SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
		for(int ii=0;ii<K->jobsCount();ii++){
			for(unsigned int jj=0;jj<K->getCorresp()[ii].size();jj++){
				for(int tt=0;tt<=T-K->getExec()[ii][jj];tt++){
					coeff = K->getC()[ii]*(tt+K->getExec()[ii][jj]);	
					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[ii][jj][tt], coeff));
		              	}
        		}
   		}
		SCIP_CALL_EXC(SCIPaddCons(scip, Cte));
		SCIP_CALL_EXC(SCIPreleaseCons(scip, &Cte));
	}
}

void POP::bigMConstraints(SCIP_CONS *Cte, vector< vector<int> > invCorresp){
	ostringstream namebuff;
	int i,t;
	unsigned int j;	
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			int indexMachine = K->getCorresp()[i-1][j-1];
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){ //For each (i,j,t)
				// Compute the corresponding big-M
				unsigned int nbterm=0;
				// Add the number of elements of Qj - 1
				nbterm += invCorresp[indexMachine-1].size() - 1;
				// Then add number of elements of E_phi_i minus those of M_j(phi_i)
				for(int nb=0;nb<K->jobsCount();nb++)
					if(K->getPhi()[nb] == K->getPhi()[i-1] && !K->isQualified(nb+1,indexMachine))
						nbterm++;
				// M uijt + sum of uiojoto <= M
				namebuff.str(""); namebuff <<"Big_" <<i<<"_"<<j<<"_"<<t;
				SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), nbterm, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
				for(unsigned int i0=1;i0<=invCorresp[indexMachine-1].size();i0++){
					int indexCurrentJob = invCorresp[indexMachine-1][i0-1];
					int indexjFori0 = K->indexMachineQualif(indexCurrentJob,indexMachine);
					if(indexCurrentJob != i && indexjFori0!=-1){ // For each job != i for which machine is qualified
						int bound = (T-K->getExec()[indexCurrentJob-1][indexjFori0]<=t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1])?T-K->getExec()[indexCurrentJob-1][indexjFori0]:t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1];
						int lowerBound = (t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1]<0)?0:t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1];
						for(int t0=lowerBound;t0<=bound;t0++)
							SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[indexCurrentJob-1][indexjFori0][t0], 1) );
					}
				}
				for(int i0=1;i0<=K->jobsCount();i0++){
					if((K->getPhi()[i0-1]==K->getPhi()[i-1])&&(i0!=i)){
						for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++){
							if(K->getCorresp()[i-1][j-1]!=K->getCorresp()[i0-1][j0-1]){
								int bound = (T-K->getExec()[i0-1][j0-1]<=t+K->getExec()[i-1][j-1])?T-K->getExec()[i0-1][j0-1]:t+K->getExec()[i-1][j-1];
								int lowerBound = (t-K->getExec()[i0-1][j0-1]<0)?0:t-K->getExec()[i0-1][j0-1];
								for(int t0=lowerBound;t0<=bound;t0++)
									SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t0], 1) );
							}
						}
					}
				}
				SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t], nbterm) );
				SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
				SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
			}
		}
	}
}

void POP::initialMaskConstraint(SCIP_CONS *Cte){
	ostringstream namebuff;
	int i;
	unsigned int j;
	namebuff.str(""); namebuff <<"Init";
	SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	for(i=1;i<=K->jobsCount();i++)
			for(j=1;j<=K->getCorresp()[i-1].size();j++)
				if(K->getInitMask()[K->getPhi()[i-1]-1]!=K->getCorresp()[i-1][j-1])
					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][0], 1) );			
	SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
        SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
}

void POP::movingAuxiliaryResourcesConstraints(SCIP_CONS *Cte){
	ostringstream namebuff;
	int i,t;	
	unsigned int j;
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){ // For each (i,j,t)
				namebuff.str(""); namebuff <<"Cut3_"<<i<<"_"<<j<<"_"<<t;
		        	SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)); // Create a constraint
				int bound = (t+K->getExec()[i-1][j-1]-1>T-K->getExec()[i-1][j-1])?T-K->getExec()[i-1][j-1]:t+K->getExec()[i-1][j-1]-1;
				for(int t0=t;t0<=bound;t0++)
					SCIP_CALL_EXC(SCIPaddCoefLinear(scip,Cte,v_varu[i-1][j-1][t0],1));					
				for(int i0=1;i0<=K->jobsCount();i0++)
					if((K->getPhi()[i0-1]==K->getPhi()[i-1])&&(i0!=i))
						for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++)
							if(K->getCorresp()[i-1][j-1]!=K->getCorresp()[i0-1][j0-1] && t+K->getExec()[i-1][j-1]<=T-K->getExec()[i0-1][j0-1])
								SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t+K->getExec()[i-1][j-1]],1));
				SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
				SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
			}
		}
	}
}

void POP::generalizedValidInequalities(SCIP_CONS *Cte, vector< vector<int> > invCorresp){
	ostringstream namebuff;
	// Constraints C_(u,l,j,k) generalizing the previous ones, cuts given by van den Akker and adapted by Bitar
	for(int gap(1);gap<4;gap++){
		for(int lower(1);lower<T/2;lower++){
			for(int upper(lower+1);upper<=lower+gap;upper++){
				for(unsigned int j=1;j<=(unsigned int)(K->machinesCount());j++){
					for(unsigned int k=1;k<=invCorresp[j-1].size();k++){ 
						// For each job that machine j can process
						namebuff.str("");
						namebuff <<"Cut_C_"<<lower<<"_"<<upper<<"_"<<j<<"_"<<invCorresp[j-1][k-1];
						// Add coeffs
						SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
						// There are two distinct sums in this constraint: One for k and one for the others. We make out both with a test
						for(unsigned int i0=1;i0<=invCorresp[j-1].size();i0++){
							int lowerBound = 0;
							int bound = 0;
							int mm = K->indexMachineQualif(invCorresp[j-1][i0-1],j);
							if(i0 == k){ // First sum, for k
								lowerBound = (lower-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1<0)?0:lower-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1;
								bound = (upper>T-K->getExec()[invCorresp[j-1][i0-1]-1][mm])?T-K->getExec()[invCorresp[j-1][i0-1]-1][mm]:upper;
							}
							else{ // Second sum, for the others
								lowerBound = (upper-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1<0)?0:upper-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1;
								bound = (lower>T-K->getExec()[invCorresp[j-1][i0-1]-1][mm])?T-K->getExec()[invCorresp[j-1][i0-1]-1][mm]:lower;
							}
							for(int t0=lowerBound;t0<=bound;t0++)
								SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[invCorresp[j-1][i0-1]-1][mm][t0], 1));
						}
						// End of constraint writing, validation
						SCIP_CALL_EXC( SCIPaddCons(scip, Cte));
						SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte));
					}
				}
	// For each couple (j,k) of machine * job we have a constraint of type (l,u).
			}
		}
	}
}

void POP::auxiliaryResourcesValidInequalities(SCIP_CONS *Cte){
	ostringstream namebuff;	
	for(int k(1);k<=K->masksCount();k++){ // For each mask
		for(int time(0);time<=T;time++){ // For each time index
			namebuff.str("");
			namebuff <<"Cut2_"<<k<<"_"<<time;
		        SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)); // Create a constraint
			for(int i0(1);i0<=K->jobsCount();i0++){
				if(K->getPhi()[i0-1] == k){ // For each job requiring the mask
					for(unsigned int j0(1);j0<=K->getCorresp()[i0-1].size();j0++){ // For each qualified machine
						int lowerBound(max(time-K->getExec()[i0-1][j0-1]+1,0));
						int bound(min(time,T-K->getExec()[i0-1][j0-1]));
						for(int t0(lowerBound);t0<=bound;t0++){
							SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t0], 1));
						}
					}
				}
			}
			SCIP_CALL_EXC( SCIPaddCons(scip, Cte));
			SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte));
		}
	}
}

void POP::create_ILP(){
	cout<<"before presolve..."<<endl;
	
	SCIP_CALL_EXC( SCIPsetPresolving(scip, SCIP_PARAMSETTING_DEFAULT, false) );
	SCIP_CALL_EXC( SCIPsetObjsense(scip, (obj == 1)?SCIP_OBJSENSE_MAXIMIZE:SCIP_OBJSENSE_MINIMIZE));
	cout<<"Creating LP"<<endl;
	ostringstream namebuff;


	/******************************************************/
	/*		Create and add variables	      */
	/******************************************************/

	// Variables Uijt
	for(int i=0;i<K->jobsCount();i++){
		for(unsigned int j=0;j<K->getCorresp()[i].size();j++){
			for(int t=0;t<=T-K->getExec()[i][j];t++){
				namebuff.str(""); namebuff << "u"<<"_"<<i+1<<"_"<<K->getCorresp()[i][j]<<"_"<<t;
	 			double coeff = (t>K->getH())?0:( (t<=K->getH()-K->getExec()[i][j])?K->getW()[i]:K->getW()[i]*(K->getH()-t)/K->getExec()[i][j]);
				SCIP_CALL_EXC( SCIPcreateVar(scip, &(v_varu[i][j][t]), namebuff.str().c_str(), 0.0, 1.0,((obj==1)?coeff:K->getC()[i]*(t+K->getExec()[i][j])),SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) ) ;
				SCIP_CALL_EXC( SCIPaddVar(scip, v_varu[i][j][t]));
			}
        	}
	}

	/***********************************************************/
	/*			Creating constraints		   */
	/***********************************************************/

	SCIP_CONS * Cte;
	int i, t; unsigned int j;
     
	SCIP_CALL_EXC(SCIPcreatePOPCte(scip, &Cte, "POPCte", TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE));
	SCIP_CALL_EXC(SCIPaddCons(scip, Cte));
	SCIP_CALL_EXC(SCIPreleaseCons(scip, &Cte));

 	// Assignment constraints	
	for(i=0;i<K->jobsCount();i++){
		namebuff.str(""); namebuff << "Aff_" << i;
		SCIP_CALL_EXC(SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, 1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
		for(j=0;j<K->getCorresp()[i].size();j++)
			for(t=0;t<1+T-K->getExec()[i][j];t++)
				SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i][j][t], 1));
		SCIP_CALL_EXC(SCIPaddCons(scip, Cte));
		SCIP_CALL_EXC(SCIPreleaseCons(scip, &Cte));
	}

	vector< vector<int> > invCorresp(K->machinesCount()); // Giving indices of qualified machines for each job
	for(i=0;i<K->jobsCount();i++)
		for(j=0;j<K->getCorresp()[i].size();j++)
			invCorresp[K->getCorresp()[i][j]-1].push_back(i+1);


	//bigMConstraints(Cte, invCorresp);

	cout<<"Assignment constraints"<<endl;
	// Constraints replacing the big-M
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			int indexMachine = K->getCorresp()[i-1][j-1];
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){ // For each (i,j,t)
				// Masks handling
				for(int i0(1);i0<=K->jobsCount();i0++){
						if((K->getPhi()[i0-1]==K->getPhi()[i-1])&&(i0!=i)){
							namebuff.str(""); namebuff <<"Mask_" <<i<<"_"<<indexMachine<<"_"<<t<<"_"<<i0;
							SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
							// Adding the u_ijt
							SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t], 1) );
							// Remaining u_i'j't'
							for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++){
								int dep = (K->getCorresp()[i0-1][j0-1]==K->getCorresp()[i-1][j-1])?1:0;
								int bound = (T-K->getExec()[i0-1][j0-1]<=t+K->getExec()[i-1][j-1]-dep)?T-K->getExec()[i0-1][j0-1]:t+K->getExec()[i-1][j-1]-dep;
								int lowerBound = (t-K->getExec()[i0-1][j0-1]+dep<0)?0:t-K->getExec()[i0-1][j0-1]+dep;
								for(int t0=lowerBound;t0<=bound;t0++)
								  	SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t0],1) );
							}
							SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
							SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
						}
				}
				// Sequence-dependent setup times
				for(unsigned int i0(1);i0<=invCorresp[indexMachine-1].size();i0++){
					int indexCurrentJob = invCorresp[indexMachine-1][i0-1]; // It is the i'
					int indexjFori0 = K->indexMachineQualif(indexCurrentJob,indexMachine);
					if(indexCurrentJob != i && indexjFori0!=-1){// For each job != i processable by the machine
						namebuff.str(""); namebuff <<"Setup_" <<i<<"_"<<indexMachine<<"_"<<t<<"_"<<indexCurrentJob;
						SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
						// Add u_ijt
						SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t], 1));
						int bound = (T-K->getExec()[indexCurrentJob-1][indexjFori0]<=t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1])?T-K->getExec()[indexCurrentJob-1][indexjFori0]:t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1];
						int lowerBound = (t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1]<0)?0:t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1];
						for(int t0=lowerBound;t0<=bound;t0++)
							SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[indexCurrentJob-1][indexjFori0][t0], 1));
						SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
						SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );			
					}
	  			}
			}
      		}
     	}

	// Last constraint: Cannot start at t=0 if the mask is not initially on the machine
	initialMaskConstraint(Cte);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Valid inequalities, to test their impact on the relaxation value

	//generalizedValidInequalities(Cte);
	//auxiliaryResourcesValidInequalities(Cte);
	//movingAuxiliaryResourcesConstraints(Cte);
	//addBestBoundConstraint(Cte);
	
	cout<<"LP done"<<endl;
}

void POP::create_ILPMO(){

	cout<<"before presolve..."<<endl;

	SCIP_CALL_EXC( SCIPsetPresolving(scip, SCIP_PARAMSETTING_AGGRESSIVE, false) );
  
  SCIP_CALL_EXC( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

  cout<<"Creating LP"<<endl;

  ostringstream namebuff;

  /******************************************************/
  /*		Création et ajout des variables		*/
  /******************************************************/

  ////////////////
  // Variables Uijt
   double coeffc1, coeffc2;
   for(int i=0;i<K->jobsCount();i++){
	for(unsigned int j=0;j<K->getCorresp()[i].size();j++){
		for(int t=0;t<=T-K->getExec()[i][j];t++){
		        namebuff.str(""); namebuff << "u"<<"_"<<i+1<<"_"<<K->getCorresp()[i][j]<<"_"<<t;
			coeffc1 = (t>K->getH())?0:((t<=K->getH()-K->getExec()[i][j])?K->getW()[i]:K->getW()[i]*(K->getH()-t)/K->getExec()[i][j]);
			coeffc2 = K->getC()[i]*(t+K->getExec()[i][j]);
	 SCIP_CALL_EXC( SCIPcreateVar(scip, &(v_varu[i][j][t]), namebuff.str().c_str(), 0.0, 1.0, epsilon*(l2*coeffc2-l1*coeffc1),SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) ) ;

         SCIP_CALL_EXC( SCIPaddVar(scip, v_varu[i][j][t]));
		}
        }
   }

   
	

   //////////Variable Z
	SCIP_CALL_EXC( SCIPcreateVar(scip, &(v_varz[0]), namebuff.str().c_str(), 0.0, SCIPinfinity(scip), 1.0,SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) ) ;

         SCIP_CALL_EXC( SCIPaddVar(scip, v_varz[0]));
		
	

//   /***********************************************************/
//   /*			Création des contraintes		*/
//   /***********************************************************/


    SCIP_CONS * Cte;
    int i,t; unsigned int j;
     // On cree les contraintes

     SCIP_CALL_EXC( SCIPcreatePOPCte(scip, &Cte, "POPCte", TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE ) );

     SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
     SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );

     // On cree les ctes d'affectation

     for (i=0;i<K->jobsCount();i++){

       namebuff.str(""); namebuff << "Aff_" << i;
       SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, 1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

       for(j=0;j<K->getCorresp()[i].size();j++){
	for(t=0;t<1+T-K->getExec()[i][j];t++){
	SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i][j][t], 1) );
        }
       }

     SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
     SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
     }

     vector< vector<int> > invCorresp(K->machinesCount()); //Donne les indices des taches pour lesquels chaque machine est qualifiee
     for(i=0;i<K->jobsCount();i++){
		for(j=0;j<K->getCorresp()[i].size();j++)
			invCorresp[K->getCorresp()[i][j]-1].push_back(i+1);	
     }

/**********************************************************************************************************************************/
     // On cree les contraintes de setup/ressources auxiliaires/capacite de machines
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			int indexMachine = K->getCorresp()[i-1][j-1];
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){//Pour tout triplet (i,j,t)
				//On commence par calculer le bigM correspondant
				unsigned int nbterm=0;
				//On rajoute le nombre d'elements de Qj - 1
				nbterm += invCorresp[indexMachine-1].size() - 1;
				//Puis on rajoute le nombre d'elements de E_phi_i moins ceux de M_j(phi_i)
				for(int nomb = 0;nomb<K->jobsCount();nomb++)
					if(K->getPhi()[nomb] == K->getPhi()[i-1] && !K->isQualified(nomb+1,indexMachine))
						nbterm++;
	
       // M uijt + somme des uiojoto <= M
       namebuff.str(""); namebuff <<"Big_" <<i<<"_"<<j<<"_"<<t;
       SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), nbterm, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

        for(unsigned int i0=1;i0<=invCorresp[indexMachine-1].size();i0++){
					int indexCurrentJob = invCorresp[indexMachine-1][i0-1];
					int indexjFori0 = K->indexMachineQualif(indexCurrentJob,indexMachine);
					if(indexCurrentJob != i && indexjFori0!=-1){//pour tout lot != i pouvant aller sur la machine
						int bound = (T-K->getExec()[indexCurrentJob-1][indexjFori0]<=t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1])?T-K->getExec()[indexCurrentJob-1][indexjFori0]:t+K->getExec()[i-1][j-1]-1+K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[i-1][j-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1];
						int lowerBound=(t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1]<0)?0:t-K->getExec()[indexCurrentJob-1][indexjFori0]+1-K->getSetup()[K->getFamilies()[indexMachine-1]-1][K->getBatch()[indexCurrentJob-1][indexjFori0]-1][K->getBatch()[i-1][j-1]-1];
						for(int t0=lowerBound;t0<=bound;t0++){
							//cout<<indexCurrentJob-1<<","<<j-1<<","<<t0<<endl;
							SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[indexCurrentJob-1][indexjFori0][t0], 1) );
						}
					}			
				}
				for(int i0=1;i0<=K->jobsCount();i0++){
						if((K->getPhi()[i0-1]==K->getPhi()[i-1])&&(i0!=i)){
							for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++){
							 if(K->getCorresp()[i-1][j-1]!=K->getCorresp()[i0-1][j0-1]){
							  int bound= (T-K->getExec()[i0-1][j0-1]<=t+K->getExec()[i-1][j-1])?T-K->getExec()[i0-1][j0-1]:t+K->getExec()[i-1][j-1];
							  int lowerBound=(t-K->getExec()[i0-1][j0-1]<0)?0:t-K->getExec()[i0-1][j0-1];
							   for(int t0=lowerBound;t0<=bound;t0++){
							  	SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t0], 1) );
							  }
							 }
							}					
						}				
				}
     SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t], nbterm) );

     SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
     SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
       }
      }
     }

	//Derniere contrainte : on ne peut commencer a t=0 si le masque n'est pas initialement sur la machine
	namebuff.str(""); namebuff <<"Init";
       SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	for(i=1;i<=K->jobsCount();i++){
			for(j=1;j<=K->getCorresp()[i-1].size();j++){
				if(K->getInitMask()[K->getPhi()[i-1]-1]!=K->getCorresp()[i-1][j-1])
				 SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][0], 1) );			
			}
	}
	SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
        SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );

	//Les deux contraintes de definition du max

	// l1*r1 <= z + l1*f1  <= +scipinfinity
	namebuff.str(""); namebuff <<"Z1";
       SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, l1*r1, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varz[0], 1));
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){//Pour tout triplet (i,j,t)
				if(t<=K->getH())
					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t],(t<=K->getH()-K->getExec()[i-1][j-1])?l1*K->getW()[i-1]:l1*K->getW()[i-1]*(K->getH()-t)/K->getExec()[i-1][j-1]) );
			}
		}
	}

	SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
        SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );


	//-z + l2*f2 <= l2*r2
	namebuff.str(""); namebuff <<"Z2";
        SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), l2*r2, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varz[0], -1));
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){//Pour tout triplet (i,j,t)
					SCIP_CALL_EXC( SCIPaddCoefLinear(scip, Cte, v_varu[i-1][j-1][t],l2*K->getC()[i-1]*(t+K->getExec()[i-1][j-1])));
			}
		}
	}

	SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
        SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Inegalites valides supplementaires pour tester l'impact sur la relaxation

	//Capacite des machines

	/**/
	for(j=1;j<=(unsigned int)(K->machinesCount());j++){
		int mm = K->indexMachineQualif(invCorresp[j-1][0],j);
		int minprocess = K->getExec()[invCorresp[j-1][0]-1][mm];
		for(unsigned int iter=1;iter<invCorresp[j-1].size();iter++){
			mm=K->indexMachineQualif(invCorresp[j-1][iter],j);
			if(K->getExec()[invCorresp[j-1][iter]-1][mm]<minprocess)
				minprocess = K->getExec()[invCorresp[j-1][iter]-1][mm];
		}
		for(int temps=0;temps<=T-minprocess;temps++){ //Pour tout couple (j,t)
			namebuff.str(""); namebuff <<"Cut_"<<j<<"_"<<temps;
		        SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
			
			for(unsigned int i0=1;i0<=invCorresp[j-1].size();i0++){
				mm = K->indexMachineQualif(invCorresp[j-1][i0-1],j);
				int lowerBound=(temps-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1<0)?0:temps-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1;
				int bound=(temps>T-K->getExec()[invCorresp[j-1][i0-1]-1][mm])?K->getExec()[invCorresp[j-1][i0-1]-1][mm]:temps;
				for(int t0=lowerBound;t0<=bound;t0++){
					SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[invCorresp[j-1][i0-1]-1][mm][t0], 1));
				}
			}
			SCIP_CALL_EXC( SCIPaddCons(scip, Cte));
			SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte));
		}
	}

	/**/

	/**/

	//Les contraintes C_(u,l,j,k) qui  generalisent les precedentes, coupes donnees par van den Akker et adaptees par Bitar
	int gap=1;
	//On va faire le tout pour des valeurs de (l,u) variant de (1,2) a (4,5)
	for(int lower=1;lower<5;lower++){
	for(int upper=lower+1;upper<=lower+gap;upper++){
	for(j=1;j<=(unsigned int)(K->machinesCount());j++){//On veut le faire pour tout j et tout k
		for(unsigned int k=1;k<=invCorresp[j-1].size();k++){//Pour tout job qui peut etre execute par la machine j choisie
			namebuff.str("");
			namebuff <<"Coupe_C_"<<lower<<"_"<<upper<<"_"<<j<<"_"<<invCorresp[j-1][k-1];
			//OK on a nomme la contrainte, et maintenant on ajoute ses coeffs
			SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
			//Il y a deux sommes distinctes dans cette contrainte : celle pour k et celle pour les autres. On distingue
			//les 2 cas par un test
			for(unsigned int i0=1;i0<=invCorresp[j-1].size();i0++){
				int lowerBound=0;
				int bound=0;
				int mm = K->indexMachineQualif(invCorresp[j-1][i0-1],j);
				if(i0==k){//On a la premiere somme, celle qui concerne k
					lowerBound=(lower-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1<0)?0:lower-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1;
					bound=(upper>T-K->getExec()[invCorresp[j-1][i0-1]-1][mm])?T-K->getExec()[invCorresp[j-1][i0-1]-1][mm]:upper;			
				}
				else{//La deuxieme somme, celle qui concerne les autres
					lowerBound=(upper-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1<0)?0:upper-K->getExec()[invCorresp[j-1][i0-1]-1][mm]+1;
					bound=(lower>T-K->getExec()[invCorresp[j-1][i0-1]-1][mm])?T-K->getExec()[invCorresp[j-1][i0-1]-1][mm]:lower;
				}
				for(int t0=lowerBound;t0<=bound;t0++){
					//if(invCorresp[j-1][k-1]==3 && j==1)
					//cout<<"ajout de variable u_"<<invCorresp[j-1][i0-1]<<"_"<<K->getCorresp()[invCorresp[j-1][i0-1]-1][mm]<<"_"<<t0<<endl;
					SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[invCorresp[j-1][i0-1]-1][mm][t0], 1));
				}
			}//Fin de l'ecriture de la contrainte, on la valide
			SCIP_CALL_EXC( SCIPaddCons(scip, Cte));
			SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte));
		}
	}
	//Pour tout couple (j,k) de machine * job qualifie pour cette machine, on a ecrit une contrainte de type (l,u).
	}}


	// Auxiliary resources
	for(int k=1;k<=K->masksCount();k++){//Pour tout masque
		for(int temps=0;temps<=T;temps++){ //Pour tout indice de temps
			namebuff.str(""); namebuff <<"Cut2_"<<k<<"_"<<temps;
		        SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));//On cree une contrainte
			for(int i0=1;i0<=K->jobsCount();i0++){
				if(K->getPhi()[i0-1]==k){ //Pour toute tache qui necessite le masque fixe
					for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++){//Pour toute machine qualifiee pour cette tache
						int lowerBound=(temps-K->getExec()[i0-1][j0-1]+1<0)?0:temps-K->getExec()[i0-1][j0-1]+1;
						int bound=(temps>T-K->getExec()[i0-1][j0-1])?T-K->getExec()[i0-1][j0-1]:temps;
						for(int t0=lowerBound;t0<=bound;t0++){
							SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t0], 1));
						}
					}
				}
			}
			SCIP_CALL_EXC( SCIPaddCons(scip, Cte));
			SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte));
		}
	}

	// Moving auxiliary resources
	for(i=1;i<=K->jobsCount();i++){
		for(j=1;j<=K->getCorresp()[i-1].size();j++){
			for(t=0;t<=T-K->getExec()[i-1][j-1];t++){//Pour tout triplet (i,j,t)
				namebuff.str(""); namebuff <<"Cut3_"<<i<<"_"<<j<<"_"<<t;
		        	SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &Cte, namebuff.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));//On cree une contrainte
				int bound = (t+K->getExec()[i-1][j-1]-1>T-K->getExec()[i-1][j-1])?T-K->getExec()[i-1][j-1]:t+K->getExec()[i-1][j-1]-1;
				for(int t0=t;t0<=bound;t0++){
					SCIP_CALL_EXC(SCIPaddCoefLinear(scip,Cte,v_varu[i-1][j-1][t0],1));					
				}
				for(int i0=1;i0<=K->jobsCount();i0++){
					if((K->getPhi()[i0-1]==K->getPhi()[i-1])&&(i0!=i)){
						for(unsigned int j0=1;j0<=K->getCorresp()[i0-1].size();j0++){
							if(K->getCorresp()[i-1][j-1]!=K->getCorresp()[i0-1][j0-1] && t+K->getExec()[i-1][j-1]<=T-K->getExec()[i0-1][j0-1]){
								SCIP_CALL_EXC(SCIPaddCoefLinear(scip, Cte, v_varu[i0-1][j0-1][t+K->getExec()[i-1][j-1]],1));
							}
						}					
					}				
				}
		     	SCIP_CALL_EXC( SCIPaddCons(scip, Cte) );
     			SCIP_CALL_EXC( SCIPreleaseCons(scip, &Cte) );
			}
		}
	}

	/**/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////


     cout<<"PL fait"<<endl;
 }


/////////////////////
void POP::solve()
{
   // this tells scip to start the solution process
   SCIP_CALL_EXC( SCIPsolve(scip) );
}

SCIP_RETCODE POP::write_solution(SCIP_SOL *sol){
	if(sol == NULL)
		cerr<<"No solution"<<endl;
	else{
		ostringstream fileName;
		fileName.str("");
		if(obj == 4)
			fileName <<"./SolFiles/"<<K->getInstanceName()<<"_PLMO.sol";
		else
			fileName <<"./SolFiles/"<<K->getInstanceName()<<".sol";
		ofstream f(fileName.str().c_str());
		int i,t; unsigned int j;
		f<<"optimal solution of cost "<<SCIPgetSolOrigObj(scip,sol)<<endl<<endl;
		for(i=1;i<=K->jobsCount();i++)
			for(j=1;j<=K->getCorresp()[i-1].size();j++)
				for(t=0;t<=T-K->getExec()[i-1][j-1];t++)
					if(SCIPgetSolVal(scip, sol, v_varu[i-1][j-1][t])>=0.99)
						f<<"u"<<i<<"_"<<K->getCorresp()[i-1][j-1]<<"_"<<t<<" = 1"<<endl;
		f.close();
	}
	// Here write the best bound
	double bound = SCIPgetDualbound(scip);
	if(obj==1 && (bound<K->getBound1() || K->getBound1()==-1)){ // Written in the txt file
		K->setBound(1,bound);
		ostringstream fileNamebound;
		fileNamebound.str("");
		fileNamebound<<"../Instances/"<<K->getInstanceName()<<"/Best_upper_bound_1.txt";
		ofstream g(fileNamebound.str().c_str());
		g<<bound<<endl;
		g.close();
	}
	if(obj==2 && bound>K->getBound2()){//On ecrit dans le fichier meilleure_bound_inf_2
		ostringstream fileNamebound1;
		K->setBound(2,bound);
		fileNamebound1.str("");
		fileNamebound1<<"../Instances/"<<K->getInstanceName()<<"/Best_lower_bound_2.txt";
		ofstream g1(fileNamebound1.str().c_str());
		g1<<bound<<endl;
		g1.close();
	}
	return SCIP_OKAY;
}

int POP::write_best_solution(){
	SCIP_SOL *sol=SCIPgetBestSol(scip);
	if(sol == NULL){
		cerr<<"There was no optimization, not known solution. Launch optim"<<endl;
		return false;
	}
	else{
		write_solution(sol);
		return true;
	}
}



