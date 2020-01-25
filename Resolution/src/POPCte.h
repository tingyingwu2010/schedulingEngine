#ifndef POPCte_H_
#define POPCte_H_

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>

#include"POP.h"

using namespace std;


class POPCte : public scip::ObjConshdlr
{
 public:

  POP *Instance;

  POPCte(SCIP* scip, POP *IInstance)
   : scip::ObjConshdlr(scip,"POPCte", "POPCte Constraint separator",
		       2000000, -2000000, -2000000, //	int sepapriority, int enfopriority, int checkpriority, 
		       1, -1, 1, 0,//	int sepafreq, int propfreq, int eagerfreq, int maxprerounds, 
		       FALSE, FALSE, FALSE, TRUE, SCIP_PROPTIMING_BEFORELP ) //SCIP_Bool  	delaysepa, SCIP_Bool  	delayprop, SCIP_Bool  	delaypresol, SCIP_Bool  	needscons, unsigned int  	timingmask 

    {
      Instance=IInstance;

    };
	
  /** destructor */
  virtual ~POPCte()
    {
    };



	
  virtual SCIP_RETCODE scip_check(
				  SCIP*              scip,               /**< SCIP data structure */
				  SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				  SCIP_CONS**        conss,              /**< array of constraints to process */
				  int                nconss,             /**< number of constraints to process */
				  SCIP_SOL*          sol,                /**< the solution to check feasibility for */
				  SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
				  SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
				  SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
				  SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
				  );
	
  virtual SCIP_RETCODE scip_enfolp(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				   );
	
  virtual SCIP_RETCODE scip_enfops(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
				   SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
				   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
				   );
  virtual SCIP_RETCODE scip_lock(
				 SCIP*              scip,               /**< SCIP data structure */
				 SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */

				 SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
									 *   constraint handler does not need constraints */
				 int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
				 int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
				 );
	
  /** transforms constraint data into data belonging to the transformed problem */
  virtual SCIP_RETCODE scip_trans(
				  SCIP*              scip,               //**< SCIP data structure *
				  SCIP_CONSHDLR*     conshdlr,           //**< the constraint handler itself *
				  SCIP_CONS*         sourcecons,         //**< source constraint to transform *
				  SCIP_CONS**        targetcons          //**< pointer to store created target constraint *
				  );
	
  /** separation method of constraint handler for LP solution
   *  possible return values for *result (if more than one applies, the first in the list should be used):
   *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
   *  - SCIP_CONSADDED  : an additional constraint was generated
   *  - SCIP_REDUCEDDOM : a variable's domain was reduced
   *  - SCIP_SEPARATED  : a cutting plane was generated
   *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
   *  - SCIP_DIDNOTRUN  : the separator was skipped
   *  - SCIP_DELAYED    : the separator was skipped, but should be called again
   */
	
  virtual SCIP_RETCODE scip_sepalp(
				   SCIP*              scip,               /**< SCIP data structure */
				   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				   SCIP_CONS**        conss,              /**< array of constraints to process */
				   int                nconss,             /**< number of constraints to process */
				   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				   );
	
  /** separation method of constraint handler for arbitrary primal solution
   *  possible return values for *result (if more than one applies, the first in the list should be used):
   *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
   *  - SCIP_CONSADDED  : an additional constraint was generated
   *  - SCIP_REDUCEDDOM : a variable's domain was reduced
   *  - SCIP_SEPARATED  : a cutting plane was generated
   *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
   *  - SCIP_DIDNOTRUN  : the separator was skipped
   *  - SCIP_DELAYED    : the separator was skipped, but should be called again
   */
  virtual SCIP_RETCODE scip_sepasol(
				    SCIP*              scip,               /**< SCIP data structure */
				    SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
				    SCIP_CONS**        conss,              /**< array of constraints to process */
				    int                nconss,             /**< number of constraints to process */
				    int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
				    SCIP_SOL*          sol,                /**< primal solution that should be separated */
				    SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
				    );
	
	
SCIP_RETCODE scip_print(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         cons,               /**< the constraint that should be displayed */
   FILE*              file                /**< the text file to store the information into */
				);
};

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
				  );


#endif
