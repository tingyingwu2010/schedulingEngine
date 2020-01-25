#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<vector>

#define MASKCAPACITY 11
#define MAXWAFER 25

using namespace std;

// File opening function
string openFile(ostringstream *pathMS, string name, ofstream *f1){
	string Name;
	Name.append((*pathMS).str());
	Name.append(name);
	(*f1).open(Name.c_str());
	cout<<endl<<"Creating and opening file "<<Name<<endl;
	return Name;
}


// Function for writing line of ChangeOverTime.txt file (stream f1) giving the sequence-dependent setup time from j to k for machine type i. Here, limitConfigTime is the setup time limit
void writeSetupTimes(ofstream *f1, unsigned int nbType, unsigned int nbBatchConfig, unsigned int limitConfigTime){
	*f1<<"BATCH_CONFIG_TYPE;FROM_BATCH_CONFIG;TO_BATCH_CONFIG;CHANGEOVER_TIME"<<endl;
	for(int i=1;i<=nbType;i++){
		for(int j=1;j<=nbBatchConfig;j++){
			for(int k=1;k<=nbBatchConfig;k++)
				*f1<<"Batch_Config_Type"<<i<<";B"<<j<<";B"<<k<<";"<<((j==k)?0:(rand()%limitConfigTime+1))<<endl;
		
		}
	}
}


// Function for writing machines data. Here, to ensure that no setup time family (batch-config) is unused, we first assign one to each machine at first
void writeMachines(ofstream *f1, unsigned int nbMachines, unsigned int nbBatchConfigType){
	*f1<<"TOOL_NAME;LOAD_PORTS_COUNT;RETICLES_IN_LIBRARY_COUNT;BATCH_CONFIG_TYPE"<<endl;
	for(int i=1;i<=nbMachines;i++){
		*f1<<"MACHINE"<<i<<";4;"<<(MASKCAPACITY)<<";"<<((i<=nbBatchConfigType)?i:(rand()%nbBatchConfigType + 1))<<endl;
	}
}

//Function for writing the transport times of auxiliary resources. We ensure they are equal to zero between similar locations. And a ficticious location (different from machines) is defined
void writeTransportTime(ofstream *f1, unsigned int nbMachines, unsigned int transportLimit){
	*f1<<"FROM_TOOL;TO_TOOL;RETICLE_TRANSPORT_TIME"<<endl;
	for(int i=0;i<=nbMachines;i++){
		for(int j=0;j<=nbMachines;j++){
			*f1<<((i==0)?"DEPOT/STOCK":"MACHINE")<<i<<";"<<((j==0)?"DEPOT/STOCK":"MACHINE")<<j<<";"<<((i==j)?0:(rand()%transportLimit+1))<<endl;
		}
	}
}

//Function to write jobs data. A given percentage of jobs will have a number of wafers lower than 25 and maxPriority parameter gives the maximum authorized weight for jobs (as they are defined in the sum of weighted completion times minimization)
void writeJob(ofstream *f1, unsigned int nbJobs, unsigned int nbMasks, double percentage, int maxPriority){	*f1<<"JOB_ID;STEP_NAME;Planned_Arrival_Time3;Carrier2;LOT_QUANTITY;PIT2;PIT_CLASS2;Start_Operation_Time3;MES_PRIORITY;Weight2;Step_Start_Time3;State_Start_Time3;CONSTRAIN_START_TIME;CONSTRAIN_DURATION;LOT_TYPE;Location2;Status2;BRICK_NAME;GENERIC_OPERATION_NAME;PRODUCT_NAME;TECHNO_NAME;RESOLVED_RECIPE_NAME;TRACK_RECIPE_NAME;MASK_ID;CAPABILITY;USAGE;Report_Run_Time3"<<endl;
	for(int i=1;i<=nbJobs;i++){
		int probability = rand()%100+1;
		*f1<<"JOB"<<i<<";0;;;"<<((probability<100*(1.0-percentage))?MAXWAFER:rand()%MAXWAFER+1)<<";;;;;"<<(rand()%maxPriority+1)<<";;;;;;;Wait;;;;;;;"<<"MASK"<<((i<=nbMasks)?i:rand()%nbMasks+1)<<";;;"<<endl;	
	}
}

// Function that writes processability (which machines can process each job) and processing times. The processLimit parameter is the limit to the processing durations
void writeProcess(ofstream *f1, unsigned int nbJobs, unsigned int nbQualif, unsigned int nbMachines, unsigned int nbBatchConfig, unsigned int processLimit){
	*f1<<"JOB_ID;STEP_NAME;PPID;BATCH_CONFIG;TOOL_NAME;PROCESSING_TIME;CTIME"<<endl;
	for(int i=1;i<=nbJobs;i++){
		vector<bool> taken(nbMachines,false);
		vector<int> machinesQualif;
		int k;
		if(i<=nbMachines){
			machinesQualif.push_back(i);
			k = rand()%nbQualif;
			taken[i-1] = true;
		}
		else{
			k = rand()%nbQualif+1;		
		}
		for(int j=1;j<=k;j++){
			int l;
			do{
				l = rand()%nbMachines+1;
			}while(taken[l-1]);
			taken[l-1] = true;
			machinesQualif.push_back(l);
		}
		for(unsigned int j=1;j<=machinesQualif.size();j++){
			*f1<<"JOB"<<i<<";0;;"<<"B"<<(rand()%nbBatchConfig+1)<<";"<<"MACHINE"<<machinesQualif[j-1]<<";"<<(rand()%processLimit)+1<<";"<<(rand()%processLimit)+1<<endl;		
		}
	}
}

// Fonction that writes data about masks, like their initial location. Here, we ensure that each machine holds no more than the given mask capacity. The stock is a specific location, with infinite capacity
void writeMask(ofstream *f1, unsigned int nbMasks, unsigned int nbMachines){
	*f1<<"Mask_ID;State_Name;Equipment_Name;Availability_Date;Report_Run_Time"<<endl;
	vector<int> m(nbMachines,0);
	for(int i=1;i<=nbMasks;i++){
		int l = rand()%(nbMachines+1);
		while((l!=0) && (m[l]>MASKCAPACITY))
			l = rand()%(nbMachines+1);
		m[l]++;
		*f1<<"MASK"<<i<<";Production;"<<((l==0)?"DEPOT/STOCK":"MACHINE")<<l<<";;"<<endl;	
	}
}

// Function that writes ToolStatuses.txt. Not used
void writeToolStatuses(ofstream *f1){
	*f1<<"TOOL_NAME;LOAD_PORT;START_DATE3;END_DATE4;REPORT_RUN_TIME3"<<endl;
}

// Function that writes TTL.txt. Not used
void writeTTL(ofstream *f1){
	*f1<<"Lot_ID;Step_Name;Tool_Name;Rank;Report_Run_Time;Lot start date"<<endl;
}

// Function that writes Target.txt. Not used
void writeTarget(ofstream *f1){
	*f1<<"TechnoGroup;LotTypeGroup;BRICK_GENE;OPERATION_GENE;MovesObj;MovesDone;Begin_Date;End_Date;"<<endl;
}


// Main function of instance creation: Folders ModelData ans StatusData are created. ModelData contains ChangeOverTime.txt, LithoTools.txt and TransportTime.txt, StatusData contains Lot.txt, Processability.txt, Reticle.txt, ToolStatuses.txt, TTL.txt and Target.txt. An additional data file is created next to the two folders. It contains the initial setup family for each machine, which will allow to compute the initial setup time for each machine, at time 0, and the horizon H. They will be created in ../Instances/PHOTO_n1_n2_n3_n4_n5_n6_n7_n8_n9_n10 where: n1 is the special number for the instance (useful for random generator), n2 is the jobs count, n3 is the machines count, n4 is the masks count, n5 the number of setup family types, n6 the number of setup families (batch-configs), n7 the processing times limit, n8 the transport times limit, n9 the setup times limit, n10 the maximum number of qualified machines for a job. We also have two parameters, the time horizon and the maximum priority for the jobs
int main(int argc, char ** argv){
	if(argc!=13){
		cout<<"usage : genePhoto NInstance Nbjobs Nbmachines NbMasks NbType NbBatchConf ProcessTimeLimit TransportTimeLimit SetupTimeLimit NbMachineQualif Horizon MaxPriority"<<endl;
		exit(EXIT_FAILURE);
	}
	unsigned int H = atoi(argv[11]);
	ostringstream path;
	path<<("../Instances/PHOTO_");
	path<<atoi(argv[1])<<"_"<<atoi(argv[2])<<"_"<<atoi(argv[3])<<"_"<<atoi(argv[4])<<"_"<<atoi(argv[5])<<"_"<<atoi(argv[6])<<"_"<<atoi(argv[7]);
	path<<"_"<<atoi(argv[8])<<"_"<<atoi(argv[9])<<"_"<<atoi(argv[10]);
	ostringstream pathModel;
	pathModel<<path.str()<<"/ModelData";
	ostringstream pathStatus;
	pathStatus<<path.str()<<"/StatusData";
	
	// Creating instance directory
	ostringstream commandDirectory1;
	commandDirectory1<<"mkdir "<<path.str();
	cout<<commandDirectory1.str().c_str()<<endl;
	system(commandDirectory1.str().c_str());
	cout<<"Creating directory "<<path.str()<<endl;
	
	// Creating both sub-directories
	ostringstream commandSubDirectory1;
	commandSubDirectory1<<"mkdir "<<pathModel.str();
	cout<<commandSubDirectory1.str().c_str()<<endl;
	system(commandSubDirectory1.str().c_str());
	ostringstream commandSubDirectory2;
	commandSubDirectory2<<"mkdir "<<pathStatus.str();
	cout<<commandSubDirectory2.str().c_str()<<endl;
	system(commandSubDirectory2.str().c_str());

	srand(atoi(argv[1]));	
	
	ofstream f1;

	// Creating ModelData files

	unsigned int nbBatchConfigType = atoi(argv[5]);
	unsigned int nbBatchConfig = atoi(argv[6]);
	unsigned int setupTimeLimit = atoi(argv[9]);
	
	string strChgOT = openFile(&pathModel, "/ChangeOverTime.txt", &f1);	;
	
	writeSetupTimes(&f1, nbBatchConfigType, nbBatchConfig, setupTimeLimit);

	// Closing file
	f1.close();	
	cout<<"Closing file "<<strChgOT<<endl;

	string machinesFile = openFile(&pathModel, "/LithoTools.txt", &f1);	;
	
	unsigned int nbMachines = atoi(argv[3]);
	
	writeMachines(&f1, nbMachines, nbBatchConfigType);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<machinesFile<<endl;

	string transportTimeFile = openFile(&pathModel,"/TransportTime.txt",&f1);	

	unsigned int transportTimeLimit = atoi(argv[8]);

	writeTransportTime(&f1, nbMachines, transportTimeLimit);
	
	// Closing file
	f1.close();	
	cout<<"Closing file "<<transportTimeFile<<endl;
	
	// End of ModelData files creation

	// Creating StatusData files
	
	string jobFile = openFile(&pathStatus, "/Lot.txt", &f1);	

	unsigned int nbJobs = atoi(argv[2]);
	unsigned int nbMasks = atoi(argv[4]);
	int maxPriority = atoi(argv[12]);

	// We give it a 10% probability that a job has less than maximum number of wafers (only useful for wafers objective function)
	writeJob(&f1, nbJobs, nbMasks, 0.1, maxPriority);
	
	// Closing file
	f1.close();	
	cout<<"Closing file "<<jobFile<<endl;

	string processFile = openFile(&pathStatus, "/Processability.txt", &f1);
	
	unsigned int nbQualif = atoi(argv[10]);
	unsigned int processTimeLimit = atoi(argv[7]);

	writeProcess(&f1, nbJobs, nbQualif, nbMachines, nbBatchConfig, processTimeLimit);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<processFile<<endl;
	
	string maskFile = openFile(&pathStatus, "/Reticle.txt", &f1);

	writeMask(&f1, nbMasks, nbMachines);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<maskFile<<endl;

	string targetFile = openFile(&pathStatus, "/Target.txt", &f1);
	
	writeTarget(&f1);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<targetFile<<endl;
	
	string ttlFile = openFile(&pathStatus, "/TTL.txt", &f1);
	
	writeTTL(&f1);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<ttlFile<<endl;
	
	string toolStatusesFile = openFile(&pathStatus, "/ToolStatuses.txt", &f1);
	
	writeToolStatuses(&f1);
	
	// Closing file
	f1.close();
	cout<<"Closing file "<<toolStatusesFile<<endl;
	
	/******** File containing additional data, namely: 
			1) Dummy setup family (batch-config), at initial state, having the same setup time towards all other families,
			2) Time horizon H of the schedule (only used for the number of wafers maximization objective)
			3) Machines, jobs, masks and setup families counts.
	************************************************************/

	string compFile = openFile(&path ,"/AdditionalData.txt", &f1);
	
	f1<<"DEFAULT_CHANGEOVER_TIME "<<"0"<<endl<<"H "<<H<<endl;
	f1<<nbMachines<<" "<<nbJobs<<" "<<nbMasks<<" "<<nbBatchConfig<<" "<<nbBatchConfigType<<endl;

	// Closing file
	f1.close();
	cout<<"Closing file "<<compFile<<endl;

	/*********************************************************************************/

	cout<<endl<<endl<<"Instance data: "<<endl<<endl;
	cout<<"Jobs count: "<<nbJobs<<endl;
	cout<<"Machines count: "<<nbMachines<<endl;
	cout<<"Masks count: "<<nbMasks<<endl;
	cout<<"Setup family types count: "<<nbBatchConfigType<<endl;
	cout<<"Setup families count: "<<nbBatchConfig<<endl;
	cout<<"Maximal number of qualified machines per job: "<<nbQualif<<endl;	
		
	return 0;
}


