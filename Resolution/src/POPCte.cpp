#include "POPCte.h"
#include "scip_exception.hpp"
#include<iostream>
#include "scip/cons_setppc.h"

using namespace std;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Sort of Constructor


/** creates and captures a Cte constraint */
SCIP_RETCODE SCIPcreatePOPCte(
				   SCIP*                 scip,               /**< SCIP data structure */
				   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
				   const char*           name,               /**< name of constraint */
				   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
				   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
				   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
				   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
				   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
				   SCIP_Bool             local,              /**< is constraint only valid locally? */
				   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
				   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
				   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
				   )
{ 
  /* find the Cte constraint handler */
    SCIP_CONSHDLR* conshdlr=NULL;
    conshdlr = SCIPfindConshdlr(scip, "POPCte");
    if( conshdlr == NULL )
     {
         SCIPerrorMessage("J'ai pas trouve le handler de POPCte\n");
         return SCIP_PLUGINNOTFOUND;
     } 
	
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, NULL, initial, separate, enforce, check, propagate,
 			    local, modifiable, dynamic, removable, FALSE) );

  return SCIP_OKAY;
}



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/// Check if an integer solution is a solution of the problem

SCIP_RETCODE POPCte::scip_check(
				    SCIP*              scip,               /**< SCIP data structure */
				    SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				    SCIP_CONS**        conss,              /**< array of constraints to process */
				    int                nconss,             /**< number of constraints to process */
				    SCIP_SOL*          sol,                /**< the solution to check feasibility for */
				    SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
				    SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
				    SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
				    SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
				    ){


  *result = SCIP_FEASIBLE;
	
  return SCIP_OKAY;
}



////////////////////////////////////////////////
/////////////////////////////////////////////////
//////////// SEPARATION POPCte



int  separ_POPCte(
		      SCIP*              scip,               /**< SCIP data structure */
		      SCIP_SOL*          sol,                /**< primal solution that should be separated */
		      SCIP_RESULT*       result,              /**< pointer to store the result of the separation call */
		      POP *Instance
		      )
{	

  *result = SCIP_DIDNOTFIND;
 
/**
  if (Instance->ChercheOwaCte(sol, Lmin, Lmax)>0){

    *result = SCIP_CONSADDED;
}
**/
  //SCIP_CALL(SCIPwriteLP(scip,"afftest.lp"));
  return 0;
	
}





/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// Four functions for separation


/** separation method of constraint handler for LP solution */
SCIP_RETCODE POPCte::scip_sepalp(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				     )
{
  //cout << "SEPALP"<<endl;
	
  separ_POPCte(scip, NULL, result, Instance);
	
  //cout<<"end SEPALP"<<endl;
  return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution 
    The method is called outside the LP solution loop (e.g., by a relaxator or a primal heuristic), 
    which means that there is no valid LP solution. Instead, the method should produce cuts that 
    separate the given solution. **/
SCIP_RETCODE POPCte::scip_sepasol(
				      SCIP*              scip,               /**< SCIP data structure */
				      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				      SCIP_CONS**        conss,              /**< array of constraints to process */
				      int                nconss,             /**< number of constraints to process */
				      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				      SCIP_SOL*          sol,                /**< primal solution that should be separated */
				      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				      )
{
  //cout << "SEPASOL\n";
	
  separ_POPCte(scip, sol, result, Instance);
	
  return SCIP_OKAY;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// ENFOLP is launched when there is no constraint left to add
/// The method is called at the end of the node processing loop for a node where the LP was solved. 
/// The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved 
/// by branching, reducing a variable's domain to exclude the solution or separating the solution with 
/// a valid cutting plane.

SCIP_RETCODE POPCte::scip_enfolp(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				     SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				     ){
 
 // Il faudrait tester si on est bien dans le cas suivant:

  // On est ici car le sous-probleme ne possede plus qu'une solution integere
  // On nous demande ce qu'on en fait: on ne peut ni la couper, ni brancher
  // On va calculer le y qui en decoule et maj la borne sup si besoin

  //  cout << "ENFOLP "<<endl;
  
   
	SCIP_SOL* newsol;

	double val;
	int i,t; unsigned int j;
	int integer;
  
	SCIP_CALL( SCIPcreateSol (scip, &newsol, NULL) );
	integer = true;
	i = 0;
	while((integer)&&(i<Instance->K->jobsCount())){
		j=0;
		while((integer)&&(j<Instance->K->getCorresp()[i].size())){
			t=0;
			while((integer)&&(t<=Instance->T - Instance->K->getExec()[i][j])){
				val = SCIPgetSolVal(scip, NULL, Instance->v_varu[i][j][t]);
				SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varu[i][j][t], val) );
   				if((val<=0.99999)&&(val>=0.000001))
					integer=false;
				t++;
			}
			j++;
		}
		i++;
	}
	if(integer){
		if(separ_POPCte(scip, NULL, result, Instance)==0)
			*result = SCIP_FEASIBLE;
		else
			*result=SCIP_CONSADDED;
	}
	else
		*result=SCIP_BRANCHED;

	//SCIP_CALL(SCIPwriteLP(scip, "essai.lp"));
	return SCIP_OKAY;	
}

/*** The method is called at the end of the node processing loop for a node where the LP was not solved.
     The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be 
     resolved by branching, reducing a variable's domain to exclude the solution or adding an additional 
     constraint. Separation is not possible, since the LP is not processed at the current node.  ***/

SCIP_RETCODE POPCte::scip_enfops(
				     SCIP*              scip,               /**< SCIP data structure */
				     SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				     SCIP_CONS**        conss,              /**< array of constraints to process */
				     int                nconss,             /**< number of constraints to process */
				     int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				     SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				     SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
				     SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				     ){
	
  *result = SCIP_FEASIBLE; // Indicates that all the constraints of the handler are feasible

 cout << "ENFOPS\n";

  return SCIP_OKAY;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// This method is called, after a constraint is added or removed from the transformed problem.
// If the constraint may get violated by changing the variable in any direction, it should call SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg)

SCIP_RETCODE POPCte::scip_lock(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
									   *   constraint handler does not need constraints */
				   int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
				   int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
				   ){

cout << "LOCK\n";

// int k,j;

// for (k=0;k<Instance->K->nbcrit;k++){
//  for(j=0;j<Instance->K->nbcrit;j++){
 
//  SCIP_CALL( SCIPaddVarLocks(scip, Instance->v_varx[k][j] , nlockspos, nlocksneg));  //If the constraint may become violated by decreasing the value of a variable

			      //nlocksneg, nlockspos)); //If the constraint may become violated by increasing the value of a variable,


   //, nlocksneg, nlockspos)); //If the constraint may become violated by increasing the value of a variable,

   //nlockspos, nlocksneg));  If the constraint may become violated by decreasing the value of a variable
   //nlockspos + nlocksneg, nlockspos + nlocksneg)); If the constraint may become violated by changing the variable in any direction, 
//  }
// }

 
  return SCIP_OKAY;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// ??

/** transforms constraint data into data belonging to the transformed problem */
SCIP_RETCODE POPCte::scip_trans(
				    SCIP*              scip,               //**< SCIP data structure *
				    SCIP_CONSHDLR*     conshdlr,           //**< the constraint handler itself *
				    SCIP_CONS*         sourcecons,         //**< source constraint to transform *
				    SCIP_CONS**        targetcons          //**< pointer to store created target constraint *
				    )
{
	
  SCIP_CALL_EXC( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, NULL,
			    SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
			    SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
			    SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
			    SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

 cout << "TRANS??\n";
	
  return SCIP_OKAY;
}




/** constraint display method of constraint handler
 *
 *  The constraint handler should store a representation of the constraint into the given text file.
 */
SCIP_RETCODE POPCte::scip_print(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         cons,               /**< the constraint that should be displayed */
   FILE*              file                /**< the text file to store the information into */
   )
{
 



   SCIPinfoMessage(scip, file, "POP Cte...\n");
   
   return SCIP_OKAY;
}


