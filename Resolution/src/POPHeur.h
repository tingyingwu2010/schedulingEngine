#ifndef __POPHEUR_H__
#define __POPHEUR_H__

#include "objscip/objscip.h"
#include "POP.h"


class POPHeur : public scip::ObjHeur
   {
      
   public:

      POP * Instance;
	
      /** default constructor http://scip.zib.de/doc/html/HEUR.html */
       POPHeur(SCIP* scip, POP *);

      /** destructor */
      virtual ~POPHeur()
      {
      }


      /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
      virtual SCIP_RETCODE scip_free(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** initialization method of primal heuristic (called after problem was transformed) */
      virtual SCIP_RETCODE scip_init(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** deinitialization method of primal heuristic (called before transformed problem is freed) */
      virtual SCIP_RETCODE scip_exit(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
       *
       *  This method is called when the presolving was finished and the branch and bound process is about to begin.
       *  The primal heuristic may use this call to initialize its branch and bound specific data.
       *
       */
      virtual SCIP_RETCODE scip_initsol(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
       *
       *  This method is called before the branch and bound process is freed.
       *  The primal heuristic should use this call to clean up its branch and bound data.
       */
      virtual SCIP_RETCODE scip_exitsol(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** execution method of primal heuristic
       *
       *  Searches for feasible primal solutions. The method is called in the node processing loop.
       *
       *  possible return values for *result:
       *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
       *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
       *  - SCIP_DIDNOTRUN  : the heuristic was skipped
       *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
       *                      its frequency
       */

      virtual SCIP_RETCODE scip_exec(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur,               /**< the primal heuristic itself */
         SCIP_HEURTIMING    heurtiming,         /**< current point in the node solving loop */
	 unsigned int       nodeinfeasible,     /** indicates if the node is infeasible */
         SCIP_RESULT*       result              /**< pointer to store the result of the heuristic call */
         ); 
   };


#endif
