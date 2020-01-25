#ifndef _INDIVIDU_
#define _INDIVIDU_

#include<vector>
#include<iostream>

class Individu{
	public : 
		double fitness;
		unsigned int fctObj;	
		std::vector<unsigned int> jobs;
		std::vector<unsigned int> machines;

	public :
		Individu(int n); //Constructor 
		Individu(unsigned int obj, int n, std::vector<std::vector<int> > corresp);
		
};

bool operator<(Individu const& i1, Individu const& i2);
bool operator>(Individu const& i1, Individu const& i2);
std::ostream& operator<<(std::ostream &out, Individu const& i1);

#endif
