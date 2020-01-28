#include"Instance.h"
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<queue>

using namespace std;

void Instance::printVector(vector<unsigned int> v){
	for(unsigned int i=0;i<v.size();i++)
		cout<<v[i]<<"\t";
	cout<<endl;
}

Individual Instance::crossoverMemetic(Individual parent1, Individual parent2){
	vector<int> pattern(n);
	vector<int> buf;
	vector<int> onesIndices;
	Individual child(parent1.fctObj,n,corresp);
	for(int i=0;i<n;i++){
		pattern[i]=rand()%2;
		if(pattern[i] == 1)
			onesIndices.push_back(i);
	}	
	int ii = 0;
	int j = 0;
	while(ii<n){
		if(pattern[ii] == 0){
			child.machines[ii]=parent1.machines[ii];
			child.jobs[ii]=parent1.jobs[ii];
		}
		else{
			child.machines[ii]=-1;
			buf.push_back(parent1.jobs[ii]);
		}
		ii++;
	}
	j = 0;
	for(int k=0;k<n;k++){
		int mem = member(buf,parent2.jobs[k]);
		if(mem != -1){
			child.machines[onesIndices[j]] = parent2.machines[k];
			child.jobs[onesIndices[j]] = parent2.jobs[k];
			buf.erase(buf.begin() + mem);
			j++;
		}
	}
	return child;
}

double Instance::newFitness(vector<vector<unsigned int> > &child, vector<unsigned int> &completionTimes, int indice, double fitness){	
	vector<int> tAvail(m,0); // Earliest availability date per machine (t=0 at start)
	vector<int> D(m,0); // Last job per machine (used to compute setup times), 0 if none
	for(int j=0;j<m;j++){
		int k(indice-1);
		while(k >= 0){
			if(child[0][k]-1 == (unsigned int)(j)){
				tAvail[j] = completionTimes[k];	
				D[j] = child[1][k];		
				k = 0;			
			}
			k--;
		}
	}
	vector<unsigned int> S(l); // Current machine per mask
	vector<int> tMask(l,0); // Earliest availability date per mask
	for(int i=0;i<l;i++){
		S[i] = initMask[i];
		int k(indice-1);
		while(k >= 0){
			if(phi[child[1][k]-1]-1 == i){
				S[i] = child[0][k];
				tMask[i] = completionTimes[k];				
				k = 0;			
			}
			k--;
		}
	}
	// Data prepared. Start at given index
	for(int j=indice;j<n;j++){ // Browse the individual
		if(S[phi[child[1][j]-1]-1]!=child[0][j]){
			tMask[phi[child[1][j]-1]-1]++; // If moving the mask is needed
			S[phi[child[1][j]-1]-1] = child[0][j];
		}
		int imQ(0);
		int imQ1 = indexMachineQualif(child[1][j],child[0][j]);
		if(D[child[0][j]-1]!=0)
			imQ = indexMachineQualif(D[child[0][j]-1],child[0][j]);
		int Beta=(D[child[0][j]-1] == 0)?0:setup[families[child[0][j]-1]-1][batch[D[child[0][j]-1]-1][imQ]-1][batch[child[1][j]-1][imQ1]-1];
		// Formula
		int end=max(tMask[phi[child[1][j]-1]-1],tAvail[child[0][j]-1]+Beta)+exec[child[1][j]-1][imQ1];
		//On actualise toutes les infos, sauf si on est à la end, auquel cas ça ne sert à rien
		if(j-n+1!=0){
			D[child[0][j]-1] = child[1][j];
			tAvail[child[0][j]-1] = end;
			tMask[phi[child[1][j]-1]-1] = end;
		}
		// Update fitness
		fitness = fitness+c[child[1][j]-1]*(end-completionTimes[j]);
		// End date is kept
		completionTimes[j] = end;
	}
	return fitness;
}

vector<vector<vector<int> > > Instance::HeuristicMemetic(Individual individual){
	vector<vector<vector<int> > > finalSchedule(m);
	vector<int> tAvail(m,0); // Earliest avalability date for each machine (t=0 at start)
	vector<int> S(l); // Current machine/location per mask
	vector<int> tMask(l,0); // Earliest availability date per mask
	vector<int> D(m,0); // Last job per machine (used to compute setup time) : 0 if none
	for(int i=1;i<=l;i++)
		S[i-1]=initMask[i-1];
	for(int j=0;j<n;j++){ // Browse the individual
		int currentJob=individual.jobs[j];
		int currentMachine=individual.machines[j];

		vector<int> elt(2);
		elt[0]=currentJob;

		if(S[phi[currentJob-1]-1]!=currentMachine){
			tMask[phi[currentJob-1]-1]++;
			S[phi[currentJob-1]-1] = currentMachine;
		}
		int imQ;
		if(D[currentMachine-1]!=0)
			imQ=indexMachineQualif(D[currentMachine-1],currentMachine);
		int Beta = (D[currentMachine-1]==0)?0:setup[families[currentMachine-1]-1][batch[D[currentMachine-1]-1][imQ]-1][batch[currentJob-1][indexMachineQualif(currentJob,currentMachine)]-1];
		// Formula
		int start =(tMask[phi[currentJob-1]-1] > tAvail[currentMachine-1]+Beta)?tMask[phi[currentJob-1]-1]:tAvail[currentMachine-1]+Beta;
		int end = start + exec[currentJob-1][indexMachineQualif(currentJob,currentMachine)];
		elt[1]=start;
		finalSchedule[currentMachine-1].push_back(elt);	
		// Update
		D[currentMachine-1] = currentJob;
		tAvail[currentMachine-1] = end;
		tMask[phi[currentJob-1]-1] = end;
	}
	return finalSchedule;
}


double Instance::evaluateSchedule(vector< vector< vector<int> > > schedule, int fct, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){ 
	//fct : 1->wafer, 2->wiCi, 3->mask moves, 4->Chebychev tricriteria
	if(fct == 1){ // Number of wafers
		double sum = 0;
		for(int j=0;j<m;j++){
			for(unsigned int i=0;i<schedule[j].size();i++){
				int job = schedule[j][i][0];
				int ti = schedule[j][i][1];
				int h = exec[job-1][indexMachineQualif(job,j+1)];
				double theta = (H<ti+h)?H:ti+h;
				theta = (theta-ti)/h;
				sum = sum + ((H<ti)?0:w[job-1]*theta);		
			}
		}
		obj1 = sum;
		return sum;
	}
	if(fct == 2){ // Weighted sum of completion times
		double sum=0;
		for(int j=0;j<m;j++){
			for(unsigned int i=0;i<schedule[j].size();i++){
				int job = schedule[j][i][0];
				int Ci = schedule[j][i][1] + exec[job-1][indexMachineQualif(job,j+1)];
				sum = sum + c[job-1] * Ci;		
			}
		}
		obj2 = sum;
		return sum;
	}
	else if(fct == 3){ // Number of mask moves
		vector< vector<int> > maskM(l);
		vector< vector<int> > lots(n);
		double deplmask = 0;
		int nbl = 0;
		for(int i=0;i<l;i++)
			maskM[i].push_back(initMask[i]);
		int k = 0;
		for(int j=0;j<m;j++){
			nbl+=schedule[j].size();
			for(unsigned int i=0;i<schedule[j].size();i++){
				lots[k].push_back(schedule[j][i][0]);
				lots[k].push_back(schedule[j][i][1]);
				lots[k].push_back(j+1);
				k++;
			}
		}
		// Sorting in increasing order of starting times
		for(int i=0;i<n;i++){
			int min = i;
			for(int j=i+1;j<n;j++)
				if(lots[j][1] < lots[min][1])
					min = j;		
			lots[i].swap(lots[min]);
		}
		// It is sorted. We can now find the number of mask moves
		for(int i=0;i<n;i++){
			if(maskM[phi[lots[i][0]-1]-1][maskM[phi[lots[i][0]-1]-1].size()-1]!=lots[i][2]){
				maskM[phi[lots[i][0]-1]-1].push_back(lots[i][2]);
				deplmask=deplmask+1;
			}
		}
		obj3 = deplmask;
		return deplmask;
	}
	else{ // Chebychev tricriteria
		double c1(evaluateSchedule(schedule,1,r1,r2,r3,l1,l2,l3,epsilon));
		double c2(evaluateSchedule(schedule,2,r1,r2,r3,l1,l2,l3,epsilon));
		double c3(evaluateSchedule(schedule,3,r1,r2,r3,l1,l2,l3,epsilon));
		c1 = l1 * (r1 - c1);
		c2 = l2 * (c2 - r2);
		c3 = l3 * (c3 - r3);
		return max(c1,max(c2,c3)) + epsilon * (c1+c2+c3);
	}
}

double Instance::fitnessMemetic(Individual individu, int fctObj, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	if(fctObj==4){
		double c1(fitnessMemetic(individu,1,r1,r2,r3,l1,l2,l3,epsilon));
		double c2(fitnessMemetic(individu,2,r1,r2,r3,l1,l2,l3,epsilon));
		double c3(fitnessMemetic(individu,3,r1,r2,r3,l1,l2,l3,epsilon));
		c1 = l1 * (r1 - c1);
		c2 = l2 * (c2 - r2);
		c3 = l3 * (c3 - r3);
		double comp(epsilon * (c1 + c2 + c3));
		if(c1 > c2){
			if(c1 > c3)
				return c1 + comp;
			else
				return c3 + comp;
		}
		if(c2 > c3){
			return c2 + comp;
		}
		else{
			return c3 + comp;
		}
	}
	double fitn(0.0);
	vector<int> tAvail(m,0); // Availability date per machine (t=0 at start)
	vector<int> S(l); // Current machine per mask
	vector<int> tMask(l,0); // Earliest availability date per mask
	vector<int> D(m,0); // Last job per machine: O if none
	for(int i=1;i<=l;i++)
		S[i-1] = initMask[i-1];
	for(int j=0;j<n;j++){
		int currentJob = individu.jobs[j];
		int currentMachine = individu.machines[j];
		if(S[phi[currentJob-1]-1] != currentMachine){
			tMask[phi[currentJob-1]-1]++;
			S[phi[currentJob-1]-1] = currentMachine;
			if(fctObj == 3)
				fitn++;
		}
		int imQ(0);
		int imQ1(indexMachineQualif(currentJob,currentMachine));
		if(D[currentMachine-1]!=0)
			imQ = indexMachineQualif(D[currentMachine-1],currentMachine);
		int Beta=(D[currentMachine-1] == 0)?0:setup[families[currentMachine-1]-1][batch[D[currentMachine-1]-1][imQ]-1][batch[currentJob-1][imQ1]-1];
		// Formula
		int start = (tMask[phi[currentJob-1]-1]>tAvail[currentMachine-1]+Beta)?tMask[phi[currentJob-1]-1]:tAvail[currentMachine-1]+Beta;
		int end = start + exec[currentJob-1][imQ1];		
		// Update
		D[currentMachine-1] = currentJob;
		tAvail[currentMachine-1] = end;
		tMask[phi[currentJob-1]-1] = end;
		if(fctObj == 1)
			if(H > start){
				if(H >= end)
				fitn+=w[currentJob-1];
			else
				fitn+=(float)(w[currentJob-1]*(H-start)/exec[currentJob-1][imQ1]);
		}
		if(fctObj == 2)
			fitn+=c[currentJob-1]*end;
	}
	return fitn;
}

double Instance::localSearchMemetic(Individual &child, int nbIter, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	// Variables to assess usefulness of meaningful neighbour
	int fctObj(child.fctObj);
	double nbAttempts(0.0);
	double delta(0.0);
	double fit(0.0);
	child.fitness = fitnessMemetic(child,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
	fit = child.fitness;
	// Best solution in local search
	Individual best(child);
	for(int i=0;i<2*nbIter;i++){
		int p = rand()%100;
		if(p < p1){ // V1 neighbourhood: Machine re-assignment
			int changeIndex=rand()%n;
			int oldMachineIndex = child.machines[changeIndex];
			int newMachineIndex = rand()%(corresp[child.jobs[changeIndex]-1].size()); // New machine must be qualified
			child.machines[changeIndex]=corresp[child.jobs[changeIndex]-1][newMachineIndex];

			double currentFitness(0.0);
			currentFitness = fitnessMemetic(child,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta = currentFitness - fit;
			double pr = (double)(rand() / ((double)(RAND_MAX)));
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ( (fctObj == 1 && pr > exp(delta / T)) || (fctObj != 1 && pr > exp(-delta / T)))){ // No improvement
				child.machines[changeIndex] = oldMachineIndex;
			}
			else{ // Improvement
				fit = currentFitness;
				child.fitness = currentFitness;
			}
		}
		else if(p < p1 + p2){ // V2 neighbourhood: Re-sequencing
			nbAttempts++;
			int moveIndex = rand()%n;
			int destinationIndex = rand() % (n+1); // We can move to n+1

			// Move now
			child.machines.insert(child.machines.begin()+destinationIndex,child.machines[moveIndex]);
			child.jobs.insert(child.jobs.begin()+destinationIndex,child.jobs[moveIndex]);
			if(destinationIndex <= moveIndex)
				moveIndex++;
			else
				destinationIndex--;
			child.machines.erase(child.machines.begin()+moveIndex);
			child.jobs.erase(child.jobs.begin()+moveIndex);

			// Test if the change is improved
			double currentFitness(0.0);
			currentFitness = fitnessMemetic(child,fctObj,r1,r2,r3,l1,l2,l3,epsilon);

			delta=currentFitness-fit;
			double pr = (double)(rand()/((double)(RAND_MAX)));
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ((fctObj == 1 && pr > exp(delta/T)) || (fctObj != 1 && pr > exp(-delta / T))) ){ //No improvement, back to original solution
				child.machines.insert(child.machines.begin()+moveIndex,child.machines[destinationIndex]);
				child.jobs.insert(child.jobs.begin()+moveIndex,child.jobs[destinationIndex]);
				if(destinationIndex > moveIndex)
					destinationIndex++;
				child.machines.erase(child.machines.begin() + destinationIndex);
				child.jobs.erase(child.jobs.begin() + destinationIndex);
				child.fitness = fit;
			}
			else{ // Improvement
				fit = currentFitness;
				child.fitness = currentFitness;
			}
		}
		else{ //Neighbourhood V3: Swap
			int moveIndex = rand() % n;
			int destinationIndex = rand() % n;
			// Now swap
			swap(child.machines[moveIndex],child.machines[destinationIndex]);
			swap(child.jobs[moveIndex],child.jobs[destinationIndex]);
			// Test if there is an improvement
			double currentFitness(0.0);
			currentFitness = fitnessMemetic(child,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta = currentFitness - fit;
			double pr = (double)(rand() / ((double)(RAND_MAX)));
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ((fctObj == 1 && pr > exp(delta / T)) || (fctObj != 1 && pr > exp(-delta / T)))){ // No improvement. Back to original solution
				swap(child.machines[moveIndex],child.machines[destinationIndex]);
				swap(child.jobs[moveIndex],child.jobs[destinationIndex]);
				child.fitness = fit;
			}
			else{ // Improvement
				fit = currentFitness;
				child.fitness = currentFitness;
			}
		}
		if(i % step == 0)
			T*=a; // Update temperature
		// Update best solution of local search
		if(child > best){
			best = child;
		}
	}
	child = best;
	return best.fitness;
}

double Instance::simulatedAnnealing(Individual &i1, int nbIter, int conv, int stagnation, double resetRatio, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, double epsilon){
	int fctObj(i1.fctObj);
	double Treset(T/10.0);
	i1.fitness = fitnessMemetic(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
	double fit(i1.fitness);
	double delta(0.0);
	int nbWithoutImprovement(0);
	int bestWithoutImprovement(0);
	// Best found solution
	Individual best(i1);
	for(int i=0;i<nbIter;i++){
		int p = rand() % 100;
		if(p < p1){ // V1
			int changeIndex = rand()%n;
			int oldMachineIndex = i1.machines[changeIndex];
			int newMachineIndex = rand() % (corresp[i1.jobs[changeIndex]-1].size());
			i1.machines[changeIndex]=corresp[i1.jobs[changeIndex]-1][newMachineIndex];
			double currentFitness(0.0);
			currentFitness = fitnessMemetic(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta = currentFitness - fit;
			double pr = (double)(rand() / ((double)(RAND_MAX)));
			if((fctObj == 1 && delta <= 0) || (fctObj != 1 && delta >= 0))
				nbWithoutImprovement++;
			else
				nbWithoutImprovement = 0;
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ((fctObj == 1 && pr > exp(delta / T)) || (fctObj != 1 && pr > exp(-delta / T)))){
				i1.machines[changeIndex] = oldMachineIndex;
				i1.fitness = fit;
			}
			else{
				fit = currentFitness;
				i1.fitness = currentFitness;
			}
		}
		else if(p < p1 + p2){ // V2
			int moveIndex = rand() % n;
			int destinationIndex = rand() % (n+1);
			i1.machines.insert(i1.machines.begin() + destinationIndex,i1.machines[moveIndex]);
			i1.jobs.insert(i1.jobs.begin() + destinationIndex,i1.jobs[moveIndex]);
			if(destinationIndex <= moveIndex)
				moveIndex++;
			else
				destinationIndex--;
			i1.machines.erase(i1.machines.begin() + moveIndex);
			i1.jobs.erase(i1.jobs.begin() + moveIndex);
			//Test if the change is improving
			double currentFitness(0.0);
			currentFitness = fitnessMemetic(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta = currentFitness - fit;
			double pr = (double)(rand() / ((double)(RAND_MAX)));
			if((fctObj == 1 && delta <= 0) || (fctObj != 1 && delta >= 0))
				nbWithoutImprovement++;
			else
				nbWithoutImprovement = 0;
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ((fctObj == 1 && pr > exp(delta / T)) || (fctObj != 1 && pr > exp(-delta / T)))){
				i1.machines.insert(i1.machines.begin() + moveIndex,i1.machines[destinationIndex]);
				i1.jobs.insert(i1.jobs.begin() + moveIndex,i1.jobs[destinationIndex]);
				if(destinationIndex > moveIndex)
					destinationIndex++;
				i1.machines.erase(i1.machines.begin() + destinationIndex);
				i1.jobs.erase(i1.jobs.begin() + destinationIndex);
				i1.fitness = fit;
			}
			else{
				fit = currentFitness;
				i1.fitness = currentFitness;
			}
		}
		else{ // V3
			int moveIndex = rand() % n;
			int destinationIndex = rand() % n;
			//Do the swap
			swap(i1.machines[moveIndex],i1.machines[destinationIndex]);
			swap(i1.jobs[moveIndex],i1.jobs[destinationIndex]);
			double currentFitness(0.0);
			currentFitness = fitnessMemetic(i1,fctObj,r1,r2,r3,l1,l2,l3,epsilon);
			delta = currentFitness - fit;
			double pr = (double)(rand() / ((double)(RAND_MAX)));
			if((fctObj == 1 && delta <= 0) || (fctObj != 1 && delta >= 0))
				nbWithoutImprovement++;
			else
				nbWithoutImprovement = 0;
			if(((fctObj == 1 && delta < 0) || (fctObj != 1 && delta > 0)) && ((fctObj == 1 && pr > exp(delta / T))||(fctObj != 1 && pr > exp(-delta / T)))){
				swap(i1.machines[moveIndex],i1.machines[destinationIndex]);
				swap(i1.jobs[moveIndex],i1.jobs[destinationIndex]);
				i1.fitness = fit;
			}
			else{
				fit = currentFitness;
				i1.fitness = currentFitness;
			}
		}
		if(i % step == 0)
			T=T*a;
		if(nbWithoutImprovement >= conv){
			T = Treset;
			nbWithoutImprovement = 0;
		}
		// Update best
		if(i1 > best){
			best = i1;
			bestWithoutImprovement = 0;
		}
		else
			bestWithoutImprovement++;
		// When best has not changed recently, too much decrease of initial temperature
		if(bestWithoutImprovement >= stagnation){ // Either we decrease Treset, or we decrease a
			Treset*=resetRatio;
			//a*=resetRatio;
			bestWithoutImprovement = 0;
		}
	}
	i1 = best;
	return best.fitness;
}

vector<vector<vector<int> > > Instance::POPMemetic(int Npop, int nbIter, int conv, int stagnation, double resetRatio, int fctObj, int Ktour, double Rep, double Ts, int ReinIter, int nbIterSA, int p1, int p2, double T, int step, double a, double r1, double r2, int r3, double l1, double l2, double l3, unsigned int duration, double epsilon){
	int Desc = Rep * Npop;
	int Surv = Desc * Ts;	
	int nbIterWithoutImprovement = 0;	
	cout<<"Number of iterations without improvement: "<<ReinIter<<endl<<endl;
	int withoutImprovement = 0;

	int iterLoc = nbIter;

	priority_queue<Individual> population;
	vector<Individual> pop;
	for(int i=0;i<Npop;i++){
		Individual i1(fctObj,n,corresp);
		//localSearchMemetic(i1,iterLoc,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3);
		simulatedAnnealing(i1,iterLoc,conv,stagnation,resetRatio,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3,epsilon);
		pop.push_back(i1);
	}
	for(int i=0;i<Npop;i++)
		population.push(pop[i]);

	Individual bestSol(population.top());

	int nIt(0);
	// currentBest is the best solution since last reset
	double currentBest(bestSol.fitness);
	double bestCost(bestSol.fitness);
	clock_t stopDuration = clock() + (duration * CLOCKS_PER_SEC); // Limit duration

	while(nbIterWithoutImprovement < nbIterSA && clock() < stopDuration){
		// Display procedure
		if(nIt % 20 == 0)
			cout<<"Best solution cost since: "<<endl<<"start\t(last reset)"<<endl;
		priority_queue<Individual> children;
		priority_queue<Individual> temp;
		vector<Individual> chd;

		// Selection, then crossover for the Npop*Rep individuals of next generation
		for(int j=0;j<Desc;j++){
			// Selection by K-tournament
			int selec1 = rand() % Npop;
			int ind(0);
			for(int i=0;i<Ktour-1;i++){
				ind = rand() % Npop;
				if(pop[ind] > pop[selec1])
					selec1=ind;
			}
			// selec1 is the index of the first individual selected by K-Tournament
			int selec2 = rand() % Npop;
			ind = 0;
			for(int i=0;i<Ktour;i++){
				ind = rand() % Npop;
				if(pop[ind] > pop[selec2])
					selec2 = ind;
			}
			// selec2 is the index of the second individual selected by K-Tournament
			// Crossover to form a child
			chd.push_back(crossoverMemetic(pop[selec1],pop[selec2]));
			//localSearchMemetic(chd[j],iterLoc,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3);
			simulatedAnnealing(chd[j],iterLoc,conv,stagnation,resetRatio,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3,epsilon);
			children.push(chd[j]);
		}
		if(Desc != 0){
			for(int j=0;j<Surv;j++){ // Take the Surv best children
				temp.push(children.top());
				children.pop();
			}
			for(int j=0;j<Npop-Surv;j++){ // Take the Npop-Surv best parents
				temp.push(population.top());
				population.pop();
			}
			while(!population.empty()) // Empty the population
				population.pop();		
			for(int j=0;j<Npop;j++){ // Gather everyone in the population
				population.push(temp.top());
				pop[j] = temp.top();
				temp.pop();		
			}
		}

		// Best global solution management
		if(fctObj == 1 && population.top().fitness > bestCost){
			bestCost = population.top().fitness;
			bestSol = population.top();
			nbIterWithoutImprovement = 0;
		}
		else if(fctObj != 1 && population.top().fitness < bestCost){
			bestCost = population.top().fitness;
			bestSol = population.top();
			nbIterWithoutImprovement = 0;
		}
		else{
			nbIterWithoutImprovement++;
		}
		// currentBest is the best solution since last reset
		if(fctObj == 1 && population.top().fitness > currentBest){
			currentBest = population.top().fitness;
			withoutImprovement = 0;
		}
		else if(fctObj != 1 && population.top().fitness < currentBest){
			currentBest = population.top().fitness;
			withoutImprovement = 0;
		}
		else{
			withoutImprovement++;
		}
		// If converged too early, reset
		if(withoutImprovement == ReinIter){
			cout<<" We reset the population: "<<endl<<endl<<endl;
			while(!population.empty()) // Empty the population
				population.pop();
			pop[0] = bestSol; // Keep the best solution or not?
			population.push(pop[0]);
			for(int ii=1;ii<Npop;ii++){
				Individual i1(fctObj,n,corresp);
				//localSearchMemetic(i1,iterLoc,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3);
				simulatedAnnealing(i1,iterLoc,conv,stagnation,resetRatio,p1,p2,T,step,a,r1,r2,r3,l1,l2,l3,epsilon);
				pop[ii] = i1;
			}
			for(int ii=1;ii<Npop;ii++)
				population.push(pop[ii]);
			// Update currentBest
			currentBest = population.top().fitness;
			withoutImprovement = 0;
		}
		nIt++;
		cout<<bestCost<<"\t\t("<<currentBest<<")"<<endl;
	}
	// Build a solution from the best found solution
	return HeuristicMemetic(bestSol);
}


