#include<iostream>
#include<cstdlib>
#include<vector>
#include"Individual.h"

using namespace std;

// Write an individual
ostream& operator<<(ostream &out, Individual const& i1){
	for(unsigned int i=0;i<i1.jobs.size();i++)
		out<<"("<<i1.jobs[i]<<","<<i1.machines[i]<<") ";
	out<<endl<<"Cout : "<<i1.fitness<<endl;
	return out;
}

// Comparaison
bool operator<(Individual const& i1, Individual const& i2){
	return (i1.fctObj == 1 && i1.fitness < i2.fitness) || (i1.fctObj != 1 && i1.fitness > i2.fitness);
}

// Comparison
bool operator>(Individual const& i1, Individual const& i2){
	return (i1.fctObj == 1 && i1.fitness > i2.fitness) || (i1.fctObj != 1 && i1.fitness < i2.fitness);
}

Individual::Individual(int n){
	machines.resize(n);
	jobs.resize(n);
}

Individual::Individual(unsigned int obj, int n, vector<vector<int> > corresp) : fitness(0), fctObj(obj) {
	machines.resize(n);
	jobs.resize(n);
	vector<int> nb(n);
	for(int i=0;i<n;i++)
		nb[i]=i+1;
	//On tire au hasard la permutation des jobs puis pour chaque job, on tire sa machine au hasard parmi son ensemble de machines qualifiées
	for(int i=0;i<n;i++){
		int h=rand()%nb.size();
		jobs[i]=nb[h];
		machines[i]=corresp[nb[h]-1][rand()%(corresp[nb[h]-1].size())]; //Un indice de machine au hasard dans corresp[nb[h]-1]
		nb.erase(nb.begin()+h);
	}
}




