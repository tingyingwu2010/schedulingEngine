#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<climits>

using namespace std;

void Instance::setBound(int obj, double val){
	if(obj == 1)
		bestBound1 = val;
	if(obj == 2)
		bestBound2 = val;
}

int Instance::nombrePlaquettes(){
	int nb=0;
	for(int i=0;i<n;i++)
		nb+=w[i];
	return nb;
}

int Instance::member(vector<int> tab, int j){
	int ind = -1;
	for(unsigned int i=1;i<=tab.size();i++)
		if(tab[i-1]==j)
			ind = i-1;
	return ind;
}

void Instance::afficheTemps(){
	cout<<"Running time "<<endl;
	for(int i=1;i<=n;i++){
		for(int j=1;j<=m;j++){
			int memb=-1;
			for(unsigned int k=1;k<=corresp[i-1].size();k++)
				if(corresp[i-1][k-1]==j)
					memb=k-1;
			if(memb!=-1)
				cout<<" "<<exec[i-1][memb]<<" ";
			else
				cout<<"   ";
		}
	cout<<endl;
	}
}

void Instance::modifyHorizon(int nH){
	H=nH;
}

Instance::Instance(string name){
	instanceName.replace(instanceName.begin(), instanceName.end(), name.substr(name.find_last_of("/")+1));
	// Path where to look to read input files
	ostringstream path;
	path<<name<<"/";
	ostringstream pathModel;
	pathModel<<path.str()<<"ModelData/";
	ostringstream pathStatus;
	pathStatus<<path.str()<<"StatusData/";

	// Garnering additional data
	path<<"AdditionalData.txt";
	ifstream fAdditionalData(path.str().c_str());
	if(!fAdditionalData){
		cout<<"Error: File "<<path.str()<<" not found\n";
		exit(EXIT_FAILURE);
	}
	char buff[512];
	fAdditionalData.getline(buff,512);
	fAdditionalData.get();
	fAdditionalData>>H;
	fAdditionalData>>m;
	fAdditionalData>>n;
	fAdditionalData>>l;
	fAdditionalData>>nB;
	fAdditionalData>>f;
	fAdditionalData.close();
	
	// Garnering best known solution
	ostringstream best1;
	best1<<name<<"/"<<"Best_Known_Sol_1.txt";
	ifstream fBest1(best1.str().c_str(), ios::in);
	if(best1)
		fBest1>>bestKnownSol1;
	else
		bestKnownSol1 = -1;
	fBest1.close();

	ostringstream best2;
	best2<<name<<"/"<<"Best_Known_Sol_2.txt";
	ifstream fBest2(best2.str().c_str(),ios::in);
	if(best2)
		fBest2>>bestKnownSol2;
	else
		bestKnownSol2 = -1;
	fBest2.close();

	ostringstream best3;
	best3<<name<<"/"<<"Best_Known_Sol_3.txt";
	ifstream fBest3(best3.str().c_str(), ios::in);
	if(best3)
		fBest3>>bestKnownSol3;
	else
		bestKnownSol3 = -1;
	fBest3.close();	
	
	// Garnering best bounds

	ostringstream b1;
	b1<<name<<"/"<<"Best_Upper_Bound_1.txt";
	ifstream best11(b1.str().c_str(), ios::in);
	if(best11)
		best11>>bestBound1;
	else
		bestBound1 = -1;
	best11.close();

	ostringstream b2;
	b2<<name<<"/"<<"Best_Lower_Bound_2.txt";
	ifstream best12(b2.str().c_str(), ios::in);
	if(best12)
		best12>>bestBound2;
	else
		bestBound2 = -1;
	best12.close();

	// Garnering number of wafers per job + job priorities + required masks
	w.resize(n);
	c.resize(n);
	phi.resize(n);
	int k,k1,k2,k3;
	ostringstream pathStatusLot;
	pathStatusLot<<pathStatus.str()<<"Lot.txt";
	cout<<"Reading file "<<pathStatusLot.str()<<endl;
	ifstream fJob(pathStatusLot.str().c_str(),ios::in);
	fJob.getline(buff,512);
	int j = 0;
	for(int i=0;i<n;i++){
		fJob.getline(buff,512,';');
		fJob.get();fJob.get();fJob.get();fJob.get();
		fJob>>k1;
		w[i]=k1;
		fJob.get();fJob.get();fJob.get();fJob.get();fJob.get();
		fJob>>k1;
		fJob.getline(buff,256,'K'); // Get to the next field, which is mask ID
		c[i]=k1;
		fJob>>k1;
		phi[i]=k1;
		if(fJob.getline(buff,256))
			j++;
	}
	n=j;
	fJob.close();
	
	// Garnering processing times and qualified machines. First qualified machines
	corresp.resize(n);

	ostringstream pathStatusEx;
	pathStatusEx<<pathStatus.str()<<"Processability.txt";
	cout<<"Reading file "<<pathStatusEx.str()<<endl;
	ifstream fexec;
	fexec.open(pathStatusEx.str().c_str());
	fexec.getline(buff,512);
	while(fexec.good()){
		fexec.getline(buff,256,'B');
		fexec>>k;fexec.get();
		if(k<=n){
			fexec.getline(buff,256,'B');
			fexec>>k1; // Number of setup family (batch-config)
			fexec.getline(buff,256,'E');
			fexec>>k2; // Number of machine
			fexec.getline(buff,512);
			if(member(corresp[k-1],k2)==-1)
				corresp[k-1].push_back(k2);
		}
		else
			fexec.getline(buff,512);
	}
	fexec.close();
	exec.resize(n);
	batch.resize(n);
	
	// We read it a second time to get the processing times after we know the qualified machines
	fexec.open(pathStatusEx.str().c_str());
	fexec.getline(buff,512);
	while(fexec.good()){
		fexec.getline(buff,256,'B');
		fexec>>k;fexec.get();
		if(k<=n){
			fexec.getline(buff,256,'B');
			fexec>>k1; // Number of batch-config
			fexec.getline(buff,256,'E');
			fexec>>k2; // Number of machine
			fexec.get();
			fexec>>k3; // Processing time	
			batch[k-1].push_back(k1);
			exec[k-1].push_back(k3);
		}
		fexec.getline(buff,256);
	}
	fexec.close();

	// Garnering masks/auxiliary resources data

	initMask.resize(l);
	ostringstream pathStatusMask;
	pathStatusMask<<pathStatus.str()<<"Reticle.txt";
	cout<<"Reading file "<<pathStatusMask.str()<<endl;
	ifstream fMask(pathStatusMask.str().c_str(), ios::in);
	fMask.getline(buff,512);
	for(int i=0;i<l;i++){
		fMask.getline(buff,512,';');
		fMask.getline(buff,512,';');
		char carac;
		fMask>>carac;
		if(carac == 'M'){
			fMask.getline(buff,512,'K');
			fMask>>initMask[i];
		}
		else
			initMask[i]=0;
		fMask.getline(buff,512);
	}
	fMask.close();

	// Reading setup families and their types

	families.resize(m);
	ostringstream pathModelFamily;
	pathModelFamily<<pathModel.str()<<"LithoTools.txt";
	cout<<"Reading file "<<pathModelFamily.str()<<endl;
	ifstream fMof(pathModelFamily.str().c_str(), ios::in);
	fMof.getline(buff,512);	
	for(int i=0;i<m;i++){
		fMof.getline(buff,512,';');
		fMof.getline(buff,512,';');
		fMof.getline(buff,512,';');
		fMof>>families[i];
	}

	// Setup times reading
	
	ostringstream pathModelCOT;
	pathModelCOT<<pathModel.str()<<"ChangeOverTime.txt";
	cout<<"Reading file "<<pathModelCOT.str()<<endl;
	ifstream fCot(pathModelCOT.str().c_str(),ios::in);
	fCot.getline(buff,512);
	setup.resize(f);
	for(int f0=0;f0<f;f0++){
		setup[f0].resize(nB);
		for(int i=1;i<=nB;i++)
			setup[f0][i-1].resize(nB);
	}
	while(fCot.good()){
		fCot.getline(buff,512,'e');
		fCot>>k; // Family
		fCot.get();
		fCot.get();
		fCot>>k1; // Setup family 1
		fCot.get();
		fCot.get();
		fCot>>k2; // Setup family 2
		fCot.get();
		fCot>>setup[k-1][k1-1][k2-1]; 
	}

	// Colour initialization: HSV mode
	colors.resize(l);
	float sumtmp = 0;											
	for(int i=0;i<l;i++){
		colors[i].resize(3);
		int alea = rand(); // TODO: We want to avoid too similar colors
		sumtmp += alea;
		colors[i][0] = (float) (sumtmp/INT_MAX);
		colors[i][1] = (float) (rand())/RAND_MAX; // For saturation
		colors[i][2] = 0.25 + (float) (rand())/RAND_MAX; // Avoid too dark colors
	}
}

double Instance::getBound1(){
	return bestBound1;
}

double Instance::getBound2(){
	return bestBound2;
}

bool Instance::estQualifiee(int i, int j){
	bool estq = false;
	for(unsigned int k=0;k<corresp[i-1].size();k++){
		if(corresp[i-1][k]==j)
			estq=true;
	}
	return estq;
}

bool Instance::estTrainQualifiee(vector<int> train, int j){
	bool estq = true;
	for(unsigned int i=0;i<train.size();i++){
		estq = estq && estQualifiee(train[i],j);
	}
	return estq;
}

unsigned int Instance::indexMachineQualif(int job, int machine){ 
	unsigned int k = 0;
	while(k<corresp[job-1].size() && corresp[job-1][k]!=machine)
		k++;
	return k;
}

bool Instance::estBatchCompatible(vector<int> train ,int batchMachine, int j){
	bool comp = true;
	for(unsigned int i=1;i<train.size();i++){
		if(batch[train[i]-1][indexMachineQualif(train[i],j)]!=batch[train[i-1]-1][indexMachineQualif(train[i-1],j)]){
			comp=false;
		}
	}
	return (comp && ((batchMachine==0) || (batch[train[0]-1][indexMachineQualif(train[0],j)]==batchMachine)));
}

int Instance::nombreMachines(){
	return m;
}

int Instance::nombreMasques(){
	return l;
}

vector< vector<int> > Instance::getCorresp(){
	return corresp;
}

int Instance::jobsCount(){
	return n;
}

double Instance::getCost(int obj){
	if(obj == 1)
		return obj1;
	if(obj == 2)
		return obj2;
	else
		return obj3;
}

vector< vector<int> > Instance::getExec(){
	return exec;
}

vector< std::vector<int> > Instance::getBatch(){
	return batch;
}

vector< std::vector< std::vector<int> > > Instance::getSetup(){
	return setup;
}

vector<int> Instance::getPhi(){
	return phi;
}
	
vector<int> Instance::getInitMask(){
	return initMask;
}

vector<double> Instance::getW(){
	return w;
}

vector<double> Instance::getC(){
	return c;
}

vector<int> Instance::getFamilies(){
	return families;
}

int Instance::getH(){
	return H;
}

string Instance::getInstanceName(){
	return instanceName;
}

double Instance::getMSol1(){
	return bestKnownSol1;
}

double Instance::getMSol2(){
	return bestKnownSol2;
}

double Instance::getMSol3(){
	return bestKnownSol3;
}

void Instance::writeBestKnownSolution(int fctObj, double cc, vector< vector<vector<int> > > schedule){
	ostringstream path;
	path<<"../Instances/"<<instanceName<<"/Best_Known_Solution_"<<fctObj<<".txt";
	ofstream fm;
	fm.open(path.str().c_str());
	fm<<cc<<endl<<endl;
	for(int j=1;j<=m;j++){
		fm<<"Machine "<<j<<": ";
		for(unsigned int kk=1;kk<=schedule[j-1].size();kk++){	
			fm<<"Job "<<schedule[j-1][kk-1][0]<<" at time "<<schedule[j-1][kk-1][1]<<" | ";
		}
		fm<<endl;
	}
	fm.close();
	switch(fctObj){
		case 1: bestKnownSol1=cc;
			break;
		case 2: bestKnownSol2=cc;
			break;
		case 3: bestKnownSol3=cc;
			break;
	}
}

void Instance::printScheduleToConsole(vector<vector<vector<int> > > schedule){
	cout<<"Obtained schedule: "<<endl;	
	for(int j=1;j<=m;j++){
		cout<<"Machine "<<j<<" : ";
		for(unsigned int kk=1;kk<=schedule[j-1].size();kk++){	
			cout<<"Job "<<schedule[j-1][kk-1][0]<<" at time "<<schedule[j-1][kk-1][1]<<" | ";
		}
		cout<<endl;
	}
}

int Instance::visualization(vector<vector<vector<int> > > schedule, int fctObj, int meth){
	string s1 = (fctObj == 1)?"W":((fctObj == 2)?"C":((fctObj == 3)?"M":((fctObj == 4)?"Tch":"Tch3")));
	string s2 = (meth == 1)?"PLNE":((meth == 2)?"Meta":((meth == 3)?"App":((meth == 4)?"PLMO":"MetaTri")));
	string s3 = (fctObj == 1)?"the maximum number of wafers before horizon":((fctObj == 2)?"the minimum weighted sum of completion times":((fctObj == 3)?"the minimum number of mask moves":"the minimum Chebychev norm around the reference point"));
	string s4 = (meth == 1)?"ILP method":((meth == 2)?"metaheuristic":((meth == 3)?"approximation algorithm":((meth == 4)?"the ILPMO":"the multi-objective metaheuristic")));
	// Writing the .dot input of GraphViz
	ostringstream pathDot;
	pathDot<<"./DotFiles/"<<instanceName<<"_"<<s1<<"_"<<s2<<".dot";
	ofstream fd(pathDot.str().c_str());
	fd<<"digraph G{"<<endl;
	fd<<"label=\"Solution of the scheduling problem with "<<n<<" jobs, "<<m<<" machines and "<<l<<" masks for "<<s3<<", given by "<<s4<<".\""<<endl;
	for(int j=0;j<m;j++){ // For each machine
		fd<<"M"<<j+1<<"[pos=\"-1,"<<m-j<<"!\",shape=\"none\",size=2]"<<endl;
		for(unsigned int i=0;i<schedule[j].size();i++){
			int indj = schedule[j][i][0];
			int exc = exec[indj-1][indexMachineQualif(indj,j+1)];
			fd<<"J"<<indj<<"[pos=\""<<(float)(schedule[j][i][1]+exc/2.0)<<","<<m-j<<"!\",width=\""<<exc<<"\",height=\"1\",color=\""<<colors[phi[indj-1]-1][0]<<" "<<colors[phi[indj-1]-1][1]<<" "<<colors[phi[indj-1]-1][2]<<"\",shape=box,style=\"rounded,filled\",fillcolor=\""<<colors[phi[indj-1]-1][0]<<" "<<colors[phi[indj-1]-1][1]<<" "<<colors[phi[indj-1]-1][2] <<"\"]"<<endl;
			fd<<schedule[j][i][1]<<"[pos=\""<<schedule[j][i][1]<<",-0.5!\",shape=none]"<<endl;
			fd<<schedule[j][i][1]+exc<<"[pos=\""<<schedule[j][i][1]+exc<<",-0.5!\",shape=none]"<<endl;		
		}
	}
	fd<<"0[pos=\"0,-0.5!\",shape=none]"<<endl;
	
	if(fctObj==1 || fctObj==4 || fctObj==5)
		fd<<H<<"[pos=\""<<H<<",-0.5!\",shape=none,fontcolor=red]"<<endl; // Display of horizon H

	// Display the legend
	fd<<"rank = sink;"<<endl<<"Legend [pos=\"-2,-2!\",shape=none, margin=0, label=< <TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"4\">"<<endl<<"<TR>"<<endl<<"<TD COLSPAN=\"2\"><B>Legend</B>"<<endl<<"</TD>"<<endl<<"</TR>"<<endl;
	for(int i=0;i<l;i++){
		 fd<<"<TR>"<<endl<<"<TD>A"<<i+1<<"</TD>"<<endl<<"<TD BGCOLOR=\""<<colors[phi[i]-1][0]<<" "<<colors[phi[i]-1][1]<<" "<<colors[phi[i]-1][2]<<"\"></TD>"<<endl<<"</TR>"<<endl;	
	}
	fd<<endl<<" </TABLE>"<<endl<<">];"<<endl;

	// End of display	
	fd<<"}"<<endl;
	fd.close();
	// Building the schedule
	ostringstream command;
	command<<"neato ./DotFiles/"<<instanceName<<"_"<<s1<<"_"<<s2<<".dot -Tpdf -o ./Visualization/"<<instanceName<<"_"<<s1<<"_"<<s2<<".pdf";
	int x1 = system(command.str().c_str());
	ostringstream command2;
	command2<<"evince ./Visualization/"<<instanceName<<"_"<<s1<<"_"<<s2<<".pdf";
	x1 = system(command2.str().c_str());
	if(fctObj != 4){
		if((fctObj == 1 && (obj1 > bestKnownSol1 || bestKnownSol1 == -1 )) || (fctObj == 2 && (obj2 < bestKnownSol2 || bestKnownSol2 == -1)) || (fctObj == 3 && (obj3 < bestKnownSol3 || bestKnownSol3 == -1))){
			ostringstream command3;
			command3<<"neato ./DotFiles/"<<instanceName<<"_"<<s1<<"_"<<s2<<".dot -Tpdf -o ../Instances/"<<instanceName<<"/"<<instanceName<<"_"<<s1<<"_"<<s2<<".pdf";
			x1 = system(command3.str().c_str());
		}
	}
	if(fctObj == 4){
		if(obj1 > bestKnownSol1 || (obj2 < bestKnownSol2 || bestKnownSol2 == -1)){
			ostringstream command3;
			command3<<"neato ./DotFiles/"<<instanceName<<"_"<<s1<<"_"<<s2<<".dot -Tpdf -o ../Instances/"<<instanceName<<"/"<<instanceName<<"_"<<s1<<"_"<<s2<<".pdf";
			x1 = system(command3.str().c_str());
		}
	}
	if(fctObj == 5){
		if(obj1 > bestKnownSol1 || (obj2 < bestKnownSol2 || bestKnownSol2 == -1) || (obj3 < bestKnownSol3 || bestKnownSol3 == -1)){
			ostringstream command3;
			command3<<"neato ./DotFiles/"<<instanceName<<"_"<<s1<<"_"<<s2<<".dot -Tpdf -o ../Instances/"<<instanceName<<"/"<<instanceName<<"_"<<s1<<"_"<<s2<<".pdf";
			x1 = system(command3.str().c_str());
		}
	}
	return x1;
}


