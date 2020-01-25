#include<iostream>
#include"objscip/objscip.h"
#include"ProblemePhoto.h"
#include"POP.h"
#include"POPHeur.h"
#include<algorithm>
#include"scip_exception.hpp"



using namespace std;

POPHeur::POPHeur(SCIP* scip, POP * IInstance)
  : ObjHeur(scip, "POPHeur", "Primal POP Heuristic", 'H', 1, 10, 0, -1, SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPLOOP, FALSE)
    //int  	priority, int  	freq, int  	freqofs, int  	maxdepth, unsigned int  	timingmask, SCIP_Bool  	usessubscip 

      {
	Instance=IInstance;
      }

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_RETCODE POPHeur::scip_free(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_RETCODE POPHeur::scip_init(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{

   return SCIP_OKAY;
}
   
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_RETCODE POPHeur::scip_exit(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{

   return SCIP_OKAY;
}
   
/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_RETCODE POPHeur::scip_initsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_RETCODE POPHeur::scip_exitsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
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

SCIP_RETCODE POPHeur::scip_exec(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur,               /**< the primal heuristic itself */
   SCIP_HEURTIMING    heurtiming,         /**< current point in the node solving loop */
   unsigned int       nodeinfeasible,     /** indicates if the node is infeasible */
   SCIP_RESULT*       result              /**< pointer to store the result of the heuristic call */
   )
{
	//cout<<"Lancement de l'heuristique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl<<endl;
	vector< vector< vector<double> > > fsolx;
	vector< vector< vector<int> > > intsolx;
     	double intsolz=0;
     	unsigned int maxj=0;
     	int maxt=0;
     	ProblemePhoto* _probData = Instance->K;
     	int n = _probData->nombreJobs();
     	SCIP_SOL* sol;

/**********************************Introduction de la meilleure solution connue***********************************************/
	if(Instance->reut == 'o' && Instance->lancementHeuristiquePrimale==0 && Instance->obj<=2){
		Instance->lancementHeuristiquePrimale=1;
		cout<<"premier appel de l'heuristique primale, qu'on reserve a la solution obtenue par la methode approchee"<<endl;
		ostringstream chemin;
		chemin<<"../Instances/"<<_probData->getNomInstance()<<"/Meilleure_Solution_Connue_"<<Instance->obj<<".txt";
		ifstream init(chemin.str().c_str(),ios::in);
		char buff[512];
		if(init){
			SCIP_SOL* newsol;
			SCIP_Bool success;
			SCIP_CALL( SCIPcreateSol (scip, &newsol, heur) );
			for(int i=0;i<n;i++){
				for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
					for(int t=0;t<=Instance->T-_probData->getExec()[i][j];t++){
         					SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varu[i][j][t], 0.0) );
					}
				}
			}
			init.getline(buff,512);
			init.getline(buff,512);
			for(int j=1;j<=_probData->nombreMachines();j++){
				int lot=0,temps=0;
				init.getline(buff,512,':');
				init.get();
				init.getline(buff,512,' ');
				while(buff[0] == 'T'){ //On a un lot pour la machine actuelle
					init>>lot;
					init.getline(buff,512,'s');
					init.get();
					init>>temps;
					for(unsigned int jj=0;jj<_probData->getCorresp()[lot-1].size();jj++)
						if(_probData->getCorresp()[lot-1][jj]==j)
							SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varu[lot-1][jj][temps], 1.0) );
					init.getline(buff,512,'|');
					init.get();
					init.getline(buff,512,' ');
				}
			}
			init.close();
			SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &success) );
			if( success )
				*result = SCIP_FOUNDSOL;
			else
				*result=SCIP_DIDNOTFIND;
			
			SCIP_CALL_EXC( SCIPfreeSol(scip, &newsol) );
			return SCIP_OKAY;
		}
	}
/******************************************************************************************************************/

	// only call heuristic, if an optimal LP solution is at hand 
     	if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
		return SCIP_OKAY;
	SCIP_CALL_EXC( SCIPcreateSol(scip, &sol, heur) );  // set up the working solution in your heuristic data
     	SCIP_CALL_EXC( SCIPlinkLPSol(scip, sol) );   // copy the current LP solution to the working solution
     	fsolx.resize(n);
     	for(int i=0;i<n;i++){
		fsolx[i].resize(_probData->getCorresp()[i].size());
		for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
			fsolx[i][j].resize(1+Instance->T-_probData->getExec()[i][j]);
		}
     	}
	intsolx.resize(n);
	for(int i=0;i<n;i++){
		intsolx[i].resize(_probData->getCorresp()[i].size());
		for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
			intsolx[i][j].resize(1+Instance->T-_probData->getExec()[i][j]);
		}
	}
     	//On recupere la solution fractionnaire dans fsolx et fsolz
     	////////////////////////////////////Gestion des variables fixees///////////////////////////////////////////
     	vector<int> tabTache(0); //tableau donnant l'ordre de parcours des taches
     	vector<int> membreTabTache(n,0);
     	for(int i=0;i<n;i++){
		for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
			for(int t=0;t<=Instance->T-_probData->getExec()[i][j];t++){
				fsolx[i][j][t]=SCIPgetSolVal(scip, sol, Instance->v_varu[i][j][t]);
				if(SCIPvarGetUbLocal(Instance->v_varu[i][j][t]) < SCIPepsilon(scip)){
					intsolx[i][j][t] = 0;
					if(membreTabTache[i]==0){
						membreTabTache[i]=1;
						tabTache.push_back(i+1);
					}
				}
				else if(SCIPvarGetLbLocal(Instance->v_varu[i][j][t]) > 0.9999){
					intsolx[i][j][t] = 1;
					if(membreTabTache[i]==0){
						membreTabTache[i]=1;
						tabTache.push_back(i+1);
					}
				}
				else
					intsolx[i][j][t]=-1; 
				//if(i==0) cout<<j<<" , "<<t<<" : "<<fsolx[i][j][t]<<endl;
			}
		}
     	}

     	for(int i=0;i<n;i++){
		if(membreTabTache[i]==0){
			membreTabTache[i]=1;
			tabTache.push_back(i+1);
		}
     	}
/**/	
	//Ici, on a deja pris la solution opt en relaxation dans fsolx et initialise a -1 les variables de la solution entiere 
	//(non fixe <=> -1)
	vector< vector<int> > invcorresp(_probData->nombreMachines()); //Donne les indices des taches pour lesquels chaque machine est qualifiee
	for(int i=0;i<n;i++){
		for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++)
			invcorresp[_probData->getCorresp()[i][j]-1].push_back(i+1);	
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pour toute tache i, on va choisir le couple (j,t) le plus adapte selon la relaxation continue
	//L'ordre de parcours des i sera donc celui donne par tabTache
	for(int ind=0;ind<n;ind++){
		int i = tabTache[ind]-1;
		unsigned int j = 0;
		int t=0;
		bool onreste=true;
		while(j < _probData->getCorresp()[i].size() && onreste==true){
			t=0;
			while((t<=Instance->T-_probData->getExec()[i][j]) && ((intsolx[i][j][t]==0) || (t==0 && _probData->getInitMask()[_probData->getPhi()[i]-1]!=_probData->getCorresp()[i][j])))
				t++;
			if(t<=Instance->T-_probData->getExec()[i][j])
				onreste=false;
			j++;
		}
		j--;
		maxj=j;
		maxt=t;
		//cout<<"couple initial (j,t) trouve : ("<<maxj<<" , "<<maxt<<")"<<endl;
		for(unsigned int j0=j;j0<_probData->getCorresp()[i].size();j0++){
			for(int t0=0;t0<=Instance->T-_probData->getExec()[i][j0];t0++){
				if((t0!=0 || _probData->getInitMask()[_probData->getPhi()[i]-1]==_probData->getCorresp()[i][j0]) && (intsolx[i][j0][t0]!=0) && (fsolx[i][j0][t0] - fsolx[i][maxj][maxt] > SCIPepsilon(scip)) && (SCIPvarGetUbLocal(Instance->v_varu[i][j0][t0]) > SCIPepsilon(scip))){
					maxj=j0;
					maxt=t0;
				}
			}
		}
		if(onreste==true){
			SCIP_CALL( SCIPfreeSol(scip, &sol) );
			*result=SCIP_DIDNOTFIND;
			return SCIP_OKAY;	
		}
		//Ici on a le maxj et le maxt : on met a 1 et on fixe a 0 ce qu'il faut pour maintenir la realisabilite	
		for(unsigned int jj=0;jj<_probData->getCorresp()[i].size();jj++){
			for(int tt=0;tt<=Instance->T-_probData->getExec()[i][jj];tt++){
				if(jj==maxj && tt==maxt)
					intsolx[i][jj][tt] = 1;
				else
					intsolx[i][jj][tt] = 0;
			}
		}
		int indicemachine = _probData->getCorresp()[i][maxj];
		for(unsigned int i0=1;i0<=invcorresp[indicemachine-1].size();i0++){
			int indiceLotCourant = invcorresp[indicemachine-1][i0-1];
			int indicejPouri0 = _probData->indiceMachineQualif(indiceLotCourant,indicemachine);
			if(indiceLotCourant != i+1 && indicejPouri0!=-1){//pour tout lot != i pouvant aller sur la machine
				int borne = ((int)(Instance->T-_probData->getExec()[indiceLotCourant-1][indicejPouri0])<=(int)(maxt+_probData->getExec()[i][maxj]-1+_probData->getSetup()[_probData->getFamilles()[indicemachine-1]-1][_probData->getBatch()[i][maxj]-1][_probData->getBatch()[indiceLotCourant-1][indicejPouri0]-1]))?Instance->T-_probData->getExec()[indiceLotCourant-1][indicejPouri0]:maxt+_probData->getExec()[i][maxj]-1+_probData->getSetup()[_probData->getFamilles()[indicemachine-1]-1][_probData->getBatch()[i][maxj]-1][_probData->getBatch()[indiceLotCourant-1][indicejPouri0]-1];
				int borneinf=(maxt-_probData->getExec()[indiceLotCourant-1][indicejPouri0]+1-_probData->getSetup()[_probData->getFamilles()[indicemachine-1]-1][_probData->getBatch()[indiceLotCourant-1][indicejPouri0]-1][_probData->getBatch()[i][maxj]-1]<0)?0:maxt-_probData->getExec()[indiceLotCourant-1][indicejPouri0]+1-_probData->getSetup()[_probData->getFamilles()[indicemachine-1]-1][_probData->getBatch()[indiceLotCourant-1][indicejPouri0]-1][_probData->getBatch()[i][maxj]-1];
				for(int t0=borneinf;t0<=borne;t0++){
					if(intsolx[indiceLotCourant-1][indicejPouri0][t0]==1){
						SCIP_CALL( SCIPfreeSol(scip, &sol) );
						*result=SCIP_DIDNOTFIND;
						return SCIP_OKAY;
					}
					else
						intsolx[indiceLotCourant-1][indicejPouri0][t0] = 0;
				}
			}
		}
		
		for(int i0=1;i0<=n;i0++){
			if((_probData->getPhi()[i0-1]==_probData->getPhi()[i])&&(i0!=i+1)){
				for(unsigned int j0=1;j0<=_probData->getCorresp()[i0-1].size();j0++){
					if(_probData->getCorresp()[i][maxj]!=_probData->getCorresp()[i0-1][j0-1]){
						int borne= ((int)(Instance->T-_probData->getExec()[i0-1][j0-1])<=(int)(maxt+_probData->getExec()[i][maxj]))?Instance->T-_probData->getExec()[i0-1][j0-1]:maxt+_probData->getExec()[i][maxj];
						int borneinf=(maxt-_probData->getExec()[i0-1][j0-1]<0)?0:maxt-_probData->getExec()[i0-1][j0-1];
						for(int t0=borneinf;t0<=borne;t0++){
								if(intsolx[i0-1][j0-1][t0]==1){
									SCIP_CALL( SCIPfreeSol(scip, &sol) );
									*result=SCIP_DIDNOTFIND;
									return SCIP_OKAY;
								}
								else
							  		intsolx[i0-1][j0-1][t0] = 0;
						}
					}
				}
			}	
		}
	}
		//Fin de la boucle principale : on a normalement fixe toutes les variables

	//cout<<"Fin de la boucle principale"<<endl;
	//On va donner la valeur de Z induite par les Uijt
	if(Instance->obj==4){
		int critere1 = 0;
		int critere2 = 0;
		//calcul des criteres 1 et 2
		for(int i=0;i<n;i++){
			for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
				for(int t=0;t<=Instance->T-_probData->getExec()[i][j];t++){
					if(t<=_probData->getH())
						critere1 += intsolx[i][j][t]*((t<=_probData->getH() - _probData->getExec()[i][j])?_probData->getW()[i]:_probData->getW()[i] * (_probData->getH()-t)/_probData->getExec()[i][j]);
					critere2 += intsolx[i][j][t]*(_probData->getC()[i]*(t+_probData->getExec()[i][j]));
				}
			}
		}
		
		double ecart1 = Instance->l1 * (Instance->r1 - critere1);
		double ecart2 = Instance->l2 * (critere2 - Instance->r2);
		intsolz = (ecart1>ecart2)?ecart1:ecart2;
		intsolz = intsolz + Instance->epsilon * (ecart1 + ecart2);  //Norme de Tchebychev pondérée augmentée
		if(intsolz > SCIPvarGetUbLocal(Instance->v_varz[0])){
			*result=SCIP_DIDNOTFIND;
			return SCIP_OKAY;		
		}
	}
/**/
      SCIP_SOL* newsol;
      SCIP_Bool success;
      // now create a solution out of the ...........and try to add it to SCIP
      SCIP_CALL( SCIPcreateSol (scip, &newsol, heur) );      
 //On met les valeurs de la solution entiere
	for(int i=0;i<n;i++){
		for(unsigned int j=0;j<_probData->getCorresp()[i].size();j++){
			for(int t=0;t<=Instance->T-_probData->getExec()[i][j];t++){
				if(intsolx[i][j][t]==1 /**&& (1.0-SCIPgetSolVal(scip, sol, Instance->v_varu[i][j][t])<SCIPepsilon(scip))**/){
         				SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varu[i][j][t], 1.0) );
				}
				else
         				SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varu[i][j][t], 0.0) );
       			}
		}
	}

	if(Instance->obj==4)
		SCIP_CALL( SCIPsetSolVal(scip, newsol, Instance->v_varz[0], intsolz) );

      //cout<<"Heuristic solution found with "<<cpt<<" nodes"<<endl;
            //SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, FALSE, FALSE, FALSE, &success) );
      SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &success) );

      //SCIP_Bool  	printreason, SCIP_Bool  	checkbounds, SCIP_Bool  	checkintegrality, SCIP_Bool  	checklprows, SCIP_Bool *  	stored 

      if( success ){
	//cout<<"trouvé mieux"<<endl;
        *result = SCIP_FOUNDSOL;  
      }
      else{
	//cout<<"pas mieux"<<endl;
	*result=SCIP_DIDNOTFIND;
      }

      SCIP_CALL_EXC( SCIPfreeSol(scip, &newsol) );
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}



